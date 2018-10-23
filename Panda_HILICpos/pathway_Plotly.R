library('tidyverse')
library('plotly')

pathway_Plotly <- function(data, pathway_name){
  complete_pathway <- separate_rows(data, pathway, sep = "\\$")
  
  pathway_list <- complete_pathway %>% 
    filter(pathway == pathway_name)
  
  g <- ggplot(pathway_list, aes(y = id, x = foldchange)) + 
    geom_point(aes(col = foldchange,size = -rank, text = name)) +
    scale_colour_gradient2(low = "green1", high = "deeppink1", mid = "grey",
                          space = "Lab", na.value = "grey50") +
    geom_vline(xintercept = 0, size = 1, col = "yellow3") +
    theme(legend.position="none") + 
    labs(title=pathway_list$pathway[1], y = NULL, x = "Log(Fold Change)",
         caption = pathway_list$name)
  
  ggplotly(g, tooltip = c("text"))
}

