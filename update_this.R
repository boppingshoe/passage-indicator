
# commends to run the Rmd file that contains counts data and passage indicator updaters
# schedule to run this file during migration season to monitor passage

setwd("G:\\STAFF\\Bobby\\adult passage\\passage-indicator\\counts_update")
library(rmarkdown)
rmarkdown::render("counts_update.Rmd")






