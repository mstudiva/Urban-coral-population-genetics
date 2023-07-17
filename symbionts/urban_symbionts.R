#### PACKAGES ####

library(tidyverse) # install.packages('tidyverse')
# library(devtools)  # install.packages("devtools")
# devtools::install_github("jrcunning/steponeR")
library(steponeR)
library(reshape2)
library(ggplot2)


#### DATA IMPORT ####

# Get the list of files with qPCR data in the directory 'data'
Plates <- list.files(path="data", pattern=".csv", full.names=T)
Plates # Plates 1-8 and 11-16 have Sample 12/NTC split across plates, so results .csv's were modified manually
# Plate 7 was rerun entirely, so the corresponding csv has been removed and the metadata csv was modified accordingly

# Run the steponeR program following your plate target labels
Urban_Out <- steponeR(files=Plates, target.ratios=c("B.A", "C.A", "D.A", "A.B", "C.B", "D.B", "A.C", "B.C", "D.C", "A.D", "B.D", "C.D"), 
                      fluor.norm=list(A=0, B=0, C=0, D=0),
                      copy.number=list(A=1, B=1, C=1, D=1), # look up copy numbers
                      ploidy=list(A=1, B=1, C=1, D=1),
                      extract=list(A=0.813, B=0.813, C= 0.813, D=0.813))

# Target ratio results
Urban<-Urban_Out$result
