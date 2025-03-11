### Generate volcano plots for rnaseq ### 
### Input : output from deseq2 `results` function
  ### excepts the following columns : 
  ###   1. log2foldchange : log2 foldchanges 
  ###   2. qvalue : adjusted p-value
  ###   3. gene_id : ensembl id 
  ###   4. gene_name : gene symbol
  ###   5. geneID : string concatenating both id and column (Serves as a unique identifier to count the # DEGs)
### Output : 
### Volcano plot with & without gene labels of interest 
### Data used to generate the volcano plot 

scriptn = "001_RNASeq_VolcanoPlots.R"

.libPaths( c( .libPaths() , "/home/yn9w/R/x86_64-pc-linux-gnu-library/3.6" ) )
options( stringsAsFactors = F )
set.seed(47) 

requiredPackages = c("tidyverse","DESeq2","rtracklayer","optparse","ggrepel")

if (length(setdiff(requiredPackages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(requiredPackages, rownames(installed.packages())))  
} else { 
  lapply( requiredPackages , require, character.only = TRUE) 
}

option_list <- list( 
                    make_option( c("-i", "--inputFile") ,  type = "character" , action = "store" , help = "Enter the deseq2File" ) , 
                    make_option( c("-f" , "--log2FoldChange" ) , type = "character" , action = "store" , help = "Enter the log2FoldChange cutoff") , 
                    make_option( c("-q" , "--qvalue" ) , type = "character" , action = "store" , help = "Enter the qvalue cutoff") , 
                    make_option( c("-l" , "--outLabel" ) , type = "character" , action = "store" , help = "Enter the name of output file") , 
                    make_option( c("-o" , "--outdir" ) , type = "character" , action = "store" , help = "Enter the name of the output directory") , 
                    make_option( c("-g" , "--geneList" ) , type = "character" , default = NULL , action = "store" , help = "Enter a gene list to label") ,
                    make_option( c("-c" , "--colors" ) , type = "character" , default = NULL , action = "store" , help = "Enter the colors for up and down genes respectively") , 
                    make_option( c("-n" , "--topn" ) , type = "character" , default = 20 , action = "store" , help = "Enter the topn genes to show on volcano plot") )  

opt = parse_args( OptionParser( option_list = option_list,
                  description = "This script plots the mutation frequencies in patientclasses specified") )

deseq2File = opt$inputFile
deseq2Result = read_tsv( deseq2File )
log2FC_thresh = as.numeric(opt$log2FoldChange)
qvalue_thresh = as.numeric(opt$qvalue)
outLabel = opt$outLabel
outDIR = opt$outdir
geneList = NULL
if( ! is.null(opt$geneList) ){ 
geneList = unlist( strsplit( opt$geneList , "," ) )
}
plot_colors = setNames( unlist(strsplit( opt$colors , "," )) , c("UpRegulated","DownRegulated") )
topn = opt$topn

argLine = paste0( unlist( lapply( seq_along( opt[ ! names(opt) %in% 'help']) , function(idx,optList) { paste( names(optList)[[idx]] , optList[[idx]] , sep = " " ) } , opt ) ) , collapse = " " )
headerLine = paste( "#" , scriptn , " " , argLine , sep = "" )
writeLines( headerLine , con = file.path( outDIR , paste0( outLabel , "_VolcanoPlot_Params_logFC_" , log2FC_thresh , "_qvalue_" , qvalue_thresh , ".txt" ) ) )

## Generating plotDat with up/ down regulated annotaitons based on the log2FC and qvalue threshold
## genes with NA padj values are filtered out 
## genes with padj greater than qvalue threshold are annotated as NonDEG
plotDat = deseq2Result %>% 
              dplyr::filter( ! is.na(log2FoldChange) ) %>% 
              mutate( group = case_when( log2FoldChange > log2FC_thresh & padj < qvalue_thresh ~ "UpRegulated" , 
                                         log2FoldChange < -log2FC_thresh & padj < qvalue_thresh ~ "DownRegulated" , 
                                         padj > qvalue_thresh ~ "NonDEG") , 
                      negLog10padj = -log10(padj) ) 
## 
deg_count = plotDat %>% 
                filter( group %in% c("UpRegulated","DownRegulated") ) %>% 
                dplyr::select( group , geneID ) %>% 
                group_by( group ) %>% 
                tally() %>% 
                mutate( label = paste0( gsub("Regulated","",group) , "(n=" , n ,")") )


pdf( file.path( outDIR , paste0( outLabel , "_VolcanoPlot_logFC_" , log2FC_thresh , "_qvalue_" , qvalue_thresh , ".pdf" ) ) )

## Plot params 
hline_intercept = -log10(qvalue_thresh)
vline_intercepts = c( -log2FC_thresh , log2FC_thresh ) 
maxFC_in_dat = ceiling( max( abs(plotDat$log2FoldChange) ) )
maxPvalue = ceiling( max( abs(plotDat$negLog10padj) , na.rm = T ) )
downRegulated_text_coord = -log2FC_thresh - 2
upRegulated_text_coord = log2FC_thresh + 2 
## We want to name the top n up and down regulated genes in the plot. To do so , 
## we are identifying top and bottom genes and adding a column to label them
## if the user provides a genelist of interest to name , we will override the top n genes
TopUP = plotDat %>% filter( group == "UpRegulated" ) %>% arrange( desc(log2FoldChange) ) %>% head( n = topn ) %>% .$gene_id
TopDOWN = plotDat %>% filter( group == "DownRegulated" ) %>% arrange( log2FoldChange ) %>% head( n = topn ) %>% .$gene_id

if( !is.null(geneList) & length(geneList) > 0 ) {
  plotDat = plotDat %>% mutate( pointLabel = if_else( is.element( gene_name , geneList ) , gene_name , "" ) )
} else {
  plotDat = plotDat %>% mutate( pointLabel = if_else( is.element( gene_id , c(TopUP,TopDOWN) ) , gene_name , "" ) )
}


#volcanoPlot = ggplot( plotDat , aes( x = log2FoldChange , y = negLog10padj  ) ) + 
#                geom_point( aes( fill = factor(group) , color=factor(group) , size = abs(log2FoldChange) ) , alpha = 0.5  )
#volcanoPlot = volcanoPlot + scale_color_manual( values = c(plot_colors,"NonDEG"="grey") ) 
volcanoPlot = ggplot( plotData , aes( x = log2FoldChange , y = negLog10padj ) ) + geom_point( shape = 19 , aes(color=factor(group)) , alpha = 0.7 )
volcanoPlot = volcanoPlot + scale_color_manual( values = c(plot_colors,"NonDEG"="grey") )  
#volcanoPlot = volcanoPlot + scale_size_continuous( range = c(0.2,3) )
volcanoPlot = volcanoPlot + theme_bw() + xlim( -maxFC_in_dat , maxFC_in_dat )
volcanoPlot = volcanoPlot + theme( axis.text = element_text( size = 15 , face = "bold" ) , 
                                   axis.title = element_text( size = 15 , face = "bold" ) ,
                                   panel.grid.major = element_blank() ,   
                                   panel.grid.minor = element_blank() , 
                                   legend.position = "none" )
volcanoPlot = volcanoPlot + geom_hline( yintercept = 1.30103 , linetype = "dashed" , color = "darkgrey" )
volcanoPlot = volcanoPlot + geom_vline( xintercept = vline_intercepts , linetype = "dashed" , color = "darkgrey" )
volcanoPlot = volcanoPlot + labs( y = "-log10( Adjusted P-value )" ) 
#volcanoPlot = volcanoPlot + geom_text_repel( aes(label = pointLabel) , 
#                                            size = 2  , max.overlaps = Inf , force=1, point.padding=unit(1,'lines'),
#                                              direction='y', nudge_x=0.1,segment.size=0.2)

print( volcanoPlot )

# vp_tmp = volcanoPlot + coord_flip()
# print(vp_tmp)

volcanoPlot = volcanoPlot + geom_text_repel( aes(label = pointLabel) , nudge_x = .15,
                  box.padding = 0.5,
                  nudge_y = 1,
                  segment.curvature = -0.1,
                  segment.ncp = 3,  , color = "red"  , segment.color = "black" , 
                  segment.angle = 20 , size = 4 , max.overlaps = Inf )

volcanoPlot = volcanoPlot + 
                  annotate("text", x = downRegulated_text_coord , y = maxPvalue , 
                           label = deg_count %>% filter( group == "DownRegulated" ) %>% .$label , 
                           label.padding=unit(4, "lines") ) + 
                  annotate("text", x = upRegulated_text_coord , y = maxPvalue , 
                           label =  deg_count %>% filter( group == "UpRegulated" ) %>% .$label  )

print( volcanoPlot )

dev.off()

write.table( plotDat , file = file.path( outDIR , paste0( outLabel , "_VolcanoPlot_logFC_" , log2FC_thresh , "_qvalue_" , qvalue_thresh ,  ".txt" ) ) , quote = F , sep = "\t" , row.names = F , col.names = T )
# pdf( file.path( outDIR , paste0( outLabel , "_VolcanoPlot_logFC_" , log2FC_thresh , "_qvalue_" , qvalue_thresh , "_" , Sys.Date() , ".pdf" ) ) )
# print( volcanoPlot )

# barplot_dat = plotDat %>% filter( gene_name %in% deg_symbols ) %>% 
# 						group_by( status ) %>% summarize( n = n() ) %>% mutate( freq = n / sum(n) ) %>% 
# 						mutate( status = factor( status , levels = c("UpRegulated","DownRegulated") ) )
# barplot = ggplot( barplot_dat , aes( x = status , y = freq , fill = status ) ) + geom_bar( stat = "identity" , width = 0.5 ) +
# 						scale_fill_manual( values = c("UpRegulated"="darkmagenta","DownRegulated"="darkgreen") ) + 
# 						theme_classic() + theme( legend.position = "none" ) + labs( x = "" )
# print( barplot )

# barplot_dat = barplot_dat %>% arrange( status )
# barplot_dat$pieLabels = paste0( barplot_dat$status , "( " , round( barplot_dat$freq*100 , 2 ) , " % )" )
# pie( barplot_dat$freq , labels = barplot_dat$pieLabels , col = alpha( c("UpRegulated"="darkmagenta","DownRegulated"="darkgreen"),0.6) , border = c("UpRegulated"="darkmagenta","DownRegulated"="darkgreen") )
# volcanoPlot = volcanoPlot + scale_color_manual( values = c("UpRegulated"="darkmagenta","DownRegulated"="darkgreen","NonDEG"="grey") ) 


