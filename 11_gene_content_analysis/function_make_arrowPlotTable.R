
library(tidyverse)

make_arrowPlotTable = function(T_gene) {
  
  T_gene %>%
    split(f = 1:nrow(.)) %>%
    lapply(FUN = function(row) {
      if (row$direction == "+") {
        x1 = row$start
        x2 = row$stop - 700
        x3 = row$stop
        if ((x3 - x1) < 700) {
          x = c(x1, x3, x1)
          y = c(-1.5, 0, 1.5)
        } else {
          x = c(x1, x2, x2, x3, x2, x2, x1)
          y = c(-1, -1, -1.5, 0, 1.5, 1, 1)
        }
      } else {
        x1 = row$start
        x2 = row$start + 700
        x3 = row$stop
        if ((x3 - x1) < 700) {
          x = c(x1, x3, x3)
          y = c(0, -1.5, 1.5)
        } else {
          x = c(x1, x2, x2, x3, x3, x2, x2)
          y = c(0, -1.5, -1, -1, 1, 1, 1.5)
        }
      }
      data.frame(x = x, y = y, gene = row$gene, stringsAsFactors = F) %>%
        left_join(row) %>%
        return()
    }) %>% bind_rows() %>%
    return()
  
}