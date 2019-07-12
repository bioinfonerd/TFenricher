
## Removes wells that have fewer than threshold counts (total, across all genes) in them
## Returns filtered Y data frame
filterWells <- function( Y, X, thresh = 1e5 )
{
    wKeep <- X %>% summarize_at( vars(-HUGO), sum ) %>% gather( Well, TotalCounts ) %>%
        filter( TotalCounts >= thresh ) %>% pull(Well)
    filter( Y, Well %in% wKeep )
}

## Custom parser for TF->Target relationships
parseTF <- function(fn)
{
    read_lines( fn ) %>% str_split( "\\t" ) %>%
        set_names( map_chr(., first) ) %>% map( ~.x[-1:-2] ) %>%
        keep( ~length(.)>=5 )
}

## Wrangles DGE data into common expression and metadata matrices
loadDGE <- function()
{
    ## DGE1
    X1 <- syn_csv( "syn18145755" ) %>% rename_at( vars(-HUGO), ~str_c("DGE1_",.) )
    Y1 <- syn_csv( "syn15673460" ) %>% mutate_at( "Well", ~str_c("DGE1_",.) ) %>%
        mutate_at( "Drug", str_to_lower ) ##%>% ##filterWells( X1, 3e4 ) %>%
##        filter( Concentration == 10 )

    ## DGE2
    X2 <- syn_csv( "syn17115624" ) %>% rename_at( vars(-HUGO), ~str_c("DGE2_",.) )
    Y2 <- syn_csv( "syn17115625" ) %>% mutate_at( "Well", ~str_c("DGE2_",.) ) %>%
        mutate_at( "Drug", str_to_lower ) ##%>% ##filterWells( X2, 8.5e4 ) %>%
##        filter( Concentration == 10 )

    ## Combine everything into a single set of matrices
    Y <- bind_rows(Y1, Y2)
    X <- inner_join(X1, X2, by="HUGO") %>% select( HUGO, Y$Well )

    list( X=X, Y=Y )
}

## Computes differential gene expression for binders vs. non-binders using edgeR
myEdgeR <- function( target, DGE, TAS )
{
    ## Identify binders and non-binders
    dPos <- TAS %>% filter( symbol == target, tas %in% c(1,2,3) ) %>% pull(name) %>% unique()
    dNeg <- TAS %>% filter( symbol == target, tas == 10 ) %>% pull(name) %>% unique()

    ## Pull corresponding wells
    wPos <- DGE$Y %>% filter( Drug %in% dPos ) %>% select(Well) %>% mutate( Binder = "Yes" )
    wNeg <- DGE$Y %>% filter( Drug %in% dNeg ) %>% select(Well) %>% mutate( Binder = "No" )

    ## Compose inputs for edgeR
    eRY <- bind_rows(wPos, wNeg) %>% mutate_at("Binder", factor, levels=c("No","Yes")) %>%
        as.data.frame %>% column_to_rownames("Well")
    eRX <- DGE$X %>% select(HUGO, rownames(eRY)) %>%
        as.data.frame %>% column_to_rownames("HUGO")

    ## Create the design matrix and estimate dispersion
    dl <- edgeR::DGEList( counts = eRX, samples = eRY )
    dl <- edgeR::calcNormFactors( dl )
    mmx <- model.matrix( ~Binder, data = dl$samples )
    dl <- edgeR::estimateDisp( dl, mmx )

    ## Compute differential gene expression
    gf <- edgeR::glmFit( dl, mmx ) %>% edgeR::glmLRT( coef = 2 )
    edgeR::topTags( gf, nrow(eRX) ) %>% as.data.frame %>% rownames_to_column( "Gene" )
}

## Enrichment analysis to infer TF activity based on target expression
myGSEA <- function( DFX, TF )
{
    fgsea <- function(v) { suppressWarnings(fgsea::fgsea(TF, v, 10000)) }
    with( DFX, set_names(logFC, Gene) ) %>% fgsea() %>% as_tibble() %>% arrange(pval)
##    with( DFX, set_names(log2FoldChange, HUGO) ) %>%
##        fgsea() %>% as_tibble() %>% arrange(pval)
}

main <- function()
{
    ## Download and wrangle DGE data
    DGE <- loadDGE()

    ## Download ISG list (positive control)
    ISG <- syn("syn11629935") %>% scan( what=character() )

    ## TAS scores
    TAS <- syn_csv("syn18268627") %>% mutate_at("name", str_to_lower)

    ## Transcription Factor -> Target relationships
    fnTF <- str_c("https://raw.githubusercontent.com/bioinfonerd/",
                  "Transcription-Factor-Databases/master/Ttrust_v2/human_TF_GSEA/all.gmt")
    TF <- parseTF(fnTF) %>% c( list(ISG=ISG) )

    ## JAK2
    dfxJAK2 <- myEdgeR( "JAK2", DGE, TAS )
    tfJAK2 <- dfxJAK2 %>% myGSEA( TF )
}
