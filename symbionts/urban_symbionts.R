#### PACKAGES ####

library(tidyverse) # install.packages('tidyverse')
# library(devtools)  # install.packages("devtools")
# devtools::install_github("jrcunning/steponeR")
library(steponeR)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(gridExtra) # install.packages("gridExtra")
library(vegan) # install.packages("vegan")
library(ape) # install.packages("ape")
library(pairwiseAdonis) # install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
set.seed(666) # Setting seed allows randomized processes (PERMANOVA) to be repeated later


#### QAQC ####

# DATA IMPORT
# Get the list of files with qPCR data in the directory 'data'
Plates <- list.files(path="data", pattern=".csv", full.names=T)
Plates # Plates 1-8 and 11-16 have Sample 12/NTC split across plates, so results .csv's were modified manually
# Plate 7 was rerun entirely, so the corresponding csv has been removed and the metadata csv was modified accordingly

# Run the steponeR program following your plate target labels
Urban_Out <- steponeR(files=Plates, target.ratios=c("A.B","A.C","A.D","B.A","B.C","B.D","C.A","C.B","C.D","D.A","D.B","D.C"), 
                      fluor.norm=list(A=0, B=0, C=0, D=0),
                      copy.number=list(A=9, B=1, C=23, D=3), 
                      ploidy=list(A=1, B=1, C=1, D=1),
                      extract=list(A=0.813, B=0.813, C= 0.813, D=0.813))

# Target ratio results
Urban<-Urban_Out$result

# METADATA (not needed here since it's already included in the symbiont abundance data)
# Metadata<-read.csv("metadata.csv")
# Metadata$Species <- as.factor(Metadata$Species)
# Metadata$Region <- as.factor(Metadata$Region)
# Metadata$Site <- as.factor(Metadata$Site)
# Metadata$CollectionDate <- as.factor(Metadata$CollectionDate)
# Metadata$Season <- as.factor(Metadata$Season)
# Metadata$Rain <- as.factor(Metadata$Rain)

# Create unique sample ID+FileName to merge with symbiont abundance data  
Urban$Sample.Plate<-paste(Urban$Sample.Name, Urban$File.Name, sep = "_" )

# Importing the symbiont abundance data from Rich
Urban_Prop <- read.csv("symProps_3Apr2024.csv", head=T)
Urban_Prop$Species <- as.factor(Urban_Prop$Species)
Urban_Prop$Region <- as.factor(Urban_Prop$Region)
Urban_Prop$CollectionDate <- as.factor(Urban_Prop$CollectionDate)
Urban_Prop$Season <- as.factor(Urban_Prop$Season)
Urban_Prop$Rain <- as.factor(Urban_Prop$Rain)

# Joining metadata with qPCR data
Urban_Join<-left_join(Urban, Urban_Prop, by="Sample.Plate")


# DATA CLEANING 1 
# Identifies all samples where only one symbiont technical replicate amplifies
One_A<- Urban_Join[which(Urban_Join$A.reps==1),]
One_B<- Urban_Join[which(Urban_Join$B.reps==1),]
One_C<- Urban_Join[which(Urban_Join$C.reps==1),]
One_D<- Urban_Join[which(Urban_Join$D.reps==1),]

# If only one technical replicate amplifies, set all corresponding symbiont ratios to NA
Urban_Join$B.A[which(Urban_Join$A.reps==1 | Urban_Join$B.reps==1)] <- NA
Urban_Join$C.A[which(Urban_Join$A.reps==1 | Urban_Join$C.reps==1)] <- NA
Urban_Join$D.A[which(Urban_Join$A.reps==1 | Urban_Join$D.reps==1)] <- NA
Urban_Join$A.B[which(Urban_Join$B.reps==1 | Urban_Join$A.reps==1)] <- NA
Urban_Join$C.B[which(Urban_Join$B.reps==1 | Urban_Join$C.reps==1)] <- NA
Urban_Join$D.B[which(Urban_Join$B.reps==1 | Urban_Join$D.reps==1)] <- NA
Urban_Join$A.C[which(Urban_Join$C.reps==1 | Urban_Join$A.reps==1)] <- NA
Urban_Join$B.C[which(Urban_Join$C.reps==1 | Urban_Join$B.reps==1)] <- NA
Urban_Join$D.C[which(Urban_Join$C.reps==1 | Urban_Join$D.reps==1)] <- NA
Urban_Join$A.D[which(Urban_Join$D.reps==1 | Urban_Join$A.reps==1)] <- NA
Urban_Join$B.D[which(Urban_Join$D.reps==1 | Urban_Join$B.reps==1)] <- NA
Urban_Join$C.D[which(Urban_Join$D.reps==1 | Urban_Join$C.reps==1)] <- NA

# Now set all NAs to 0
Urban_Join$B.A[is.na(Urban_Join$B.A)] <- 0
Urban_Join$C.A[is.na(Urban_Join$C.A)] <- 0
Urban_Join$D.A[is.na(Urban_Join$D.A)] <- 0
Urban_Join$A.B[is.na(Urban_Join$A.B)] <- 0
Urban_Join$C.B[is.na(Urban_Join$C.B)] <- 0
Urban_Join$D.B[is.na(Urban_Join$D.B)] <- 0
Urban_Join$A.C[is.na(Urban_Join$A.C)] <- 0
Urban_Join$B.C[is.na(Urban_Join$B.C)] <- 0
Urban_Join$D.C[is.na(Urban_Join$D.C)] <- 0
Urban_Join$A.D[is.na(Urban_Join$A.D)] <- 0
Urban_Join$B.D[is.na(Urban_Join$B.D)] <- 0
Urban_Join$C.D[is.na(Urban_Join$C.D)] <- 0


# FILTERING BY SPECIES 
Urban_Join %>%
  filter(Species == "Mcav") -> Mcav
Urban_Join %>%
  filter(Species == "Ofav") -> Ofav
Urban_Join %>%
  filter(Species == "Cnat") -> Cnat
Urban_Join %>%
  filter(Species == "Pseu") -> Pseu
Urban_Join %>%
  filter(Species == "Dlab") -> Dlab

# Determining mean symbiont abundance by species and region for downstream QAQC
Urban_Join %>%
  filter(Species == "Mcav") %>%
  group_by(Region) %>%
  summarize(across(propA:propD, mean, na.rm=TRUE)) -> Mcav_dom
Mcav_dom # reef: AC, urban: AD
Urban_Join %>%
  filter(Species == "Ofav") %>%
  group_by(Region) %>%
  summarize(across(propA:propD, mean, na.rm=TRUE)) -> Ofav_dom
Ofav_dom # reef: ABD, urban: D
Urban_Join %>%
  filter(Species == "Cnat") %>%
  group_by(Region) %>%
  summarize(across(propA:propD, mean, na.rm=TRUE)) -> Cnat_dom
Cnat_dom # reef: BD, urban: BD
Urban_Join %>%
  filter(Species == "Pseu") %>%
  group_by(Region) %>%
  summarize(across(propA:propD, mean, na.rm=TRUE)) -> Pseu_dom
Pseu_dom # reef: BD, urban: BD
Urban_Join %>%
  filter(Species == "Dlab") %>%
  group_by(Region) %>%
  summarize(across(propA:propD, mean, na.rm=TRUE)) -> Dlab_dom
Dlab_dom # reef: ABD, urban: D

# Further filtering species by site based on dominant symbionts
Mcav %>%
  filter(grepl('reef', Region)) -> Mcav_reef
Mcav %>%
  filter(grepl('urban', Region)) -> Mcav_urban
Ofav %>%
  filter(grepl('reef', Region)) -> Ofav_reef
Ofav %>%
  filter(grepl('urban', Region)) -> Ofav_urban
# Not needed since dominant symbionts are the same between regions
# Cnat %>%
#   filter(grepl('reef', Region)) -> Cnat_reef
# Cnat %>%
#     filter(grepl('urban', Region)) -> Cnat_urban
# Pseu %>% 
#   filter(grepl('reef', Region)) -> Pseu_reef
# Pseu %>%
#   filter(grepl('urban', Region)) -> Pseu_urban
Dlab %>%
  filter(grepl('reef', Region)) -> Dlab_reef
Dlab %>%
  filter(grepl('urban', Region)) -> Dlab_urban


# DATA CLEANING 2 
# Makes a list of samples with only one technical replicate of the dominant symbiont type
Reps_Mcav_reef <- Mcav_reef[which(Mcav_reef$A.reps<2 | Mcav_reef$C.reps<2), ] 
Reps_Mcav_reef %>%
  mutate(violation = "reps", .before = Sample.Name) -> Reps_Mcav_reef
Reps_Mcav_urban <- Mcav_urban[which(Mcav_urban$A.reps<2 | Mcav_urban$D.reps<2), ] 
Reps_Mcav_urban %>%
  mutate(violation = "reps", .before = Sample.Name) -> Reps_Mcav_urban

Reps_Ofav_reef <- Ofav_reef[which(Ofav_reef$A.reps<2 | Ofav_reef$B.reps<2 | Ofav_reef$D.reps<2), ] 
Reps_Ofav_reef %>%
  mutate(violation = "reps", .before = Sample.Name) -> Reps_Ofav_reef
Reps_Ofav_urban <- Ofav_urban[which(Ofav_urban$D.reps<2), ] 
Reps_Ofav_urban %>%
  mutate(violation = "reps", .before = Sample.Name) -> Reps_Ofav_urban

Reps_Cnat <- Cnat[which(Cnat$B.reps<2 | Cnat$D.reps<2), ] 
Reps_Cnat %>%
  mutate(violation = "reps", .before = Sample.Name) -> Reps_Cnat

Reps_Pseu <- Pseu[which(Pseu$B.reps<2 | Pseu$D.reps<2), ] 
Reps_Pseu %>%
  mutate(violation = "reps", .before = Sample.Name) -> Reps_Pseu

Reps_Dlab_reef <- Dlab_reef[which(Dlab_reef$A.reps<2 | Dlab_reef$B.reps<2 | Dlab_reef$D.reps<2), ] 
Reps_Dlab_reef %>%
  mutate(violation = "reps", .before = Sample.Name) -> Reps_Dlab_reef
Reps_Dlab_urban <- Dlab_urban[which(Dlab_urban$D.reps<2), ] 
Reps_Dlab_urban %>%
  mutate(violation = "reps", .before = Sample.Name) -> Reps_Dlab_urban

# Makes a list of samples where technical replicates of dominant symbiont types had standard deviation >1.5
StDe1.5_Mcav_reef <- Mcav_reef[which(Mcav_reef$A.CT.sd>1.5 | Mcav_reef$C.CT.sd>1.5), ] 
StDe1.5_Mcav_reef %>%
  mutate(violation = "stdv", .before = Sample.Name) -> StDe1.5_Mcav_reef
StDe1.5_Mcav_urban <- Mcav_urban[which(Mcav_urban$A.CT.sd>1.5 | Mcav_urban$D.CT.sd>1.5), ] 
StDe1.5_Mcav_urban %>%
  mutate(violation = "stdv", .before = Sample.Name) -> StDe1.5_Mcav_urban

StDe1.5_Ofav_reef <- Ofav_reef[which(Ofav_reef$A.CT.sd>1.5 | Ofav_reef$B.CT.sd>1.5 | Ofav_reef$D.CT.sd>1.5), ] 
StDe1.5_Ofav_reef %>%
  mutate(violation = "stdv", .before = Sample.Name) -> StDe1.5_Ofav_reef
StDe1.5_Ofav_urban <- Ofav_urban[which(Ofav_urban$D.CT.sd>1.5), ] 
StDe1.5_Ofav_urban %>%
  mutate(violation = "stdv", .before = Sample.Name) -> StDe1.5_Ofav_urban

StDe1.5_Cnat <- Cnat[which(Cnat$B.CT.sd>1.5 | Cnat$D.CT.sd>1.5), ] 
StDe1.5_Cnat %>%
  mutate(violation = "stdv", .before = Sample.Name) -> StDe1.5_Cnat

StDe1.5_Pseu <- Pseu[which(Pseu$B.CT.sd>1.5 | Pseu$D.CT.sd>1.5), ] 
StDe1.5_Pseu %>%
  mutate(violation = "stdv", .before = Sample.Name) -> StDe1.5_Pseu

StDe1.5_Dlab_reef <- Dlab_reef[which(Dlab_reef$A.CT.sd>1.5 | Dlab_reef$B.CT.sd>1.5 | Dlab_reef$D.CT.sd>1.5), ] 
StDe1.5_Dlab_reef %>%
  mutate(violation = "stdv", .before = Sample.Name) -> StDe1.5_Dlab_reef
StDe1.5_Dlab_urban <- Dlab_urban[which(Dlab_urban$D.CT.sd>1.5), ] 
StDe1.5_Dlab_urban %>%
  mutate(violation = "stdv", .before = Sample.Name) -> StDe1.5_Dlab_urban

# Makes a list of samples where the Ct mean of dominant symbiont types was >25 (late amplification)
Late_Mcav_reef<-Mcav_reef[which(Mcav_reef$A.CT.mean>35 & Mcav_reef$C.CT.mean>35), ] 
Late_Mcav_reef %>%
  mutate(violation = "late", .before = Sample.Name) -> Late_Mcav_reef
Late_Mcav_urban<-Mcav_urban[which(Mcav_urban$A.CT.mean>35 & Mcav_urban$D.CT.mean>35), ] 
Late_Mcav_urban %>%
  mutate(violation = "late", .before = Sample.Name) -> Late_Mcav_urban

Late_Ofav_reef<-Ofav_reef[which(Ofav_reef$A.CT.mean>35 & Ofav_reef$B.CT.mean>35 & Ofav_reef$D.CT.mean>35), ] 
Late_Ofav_reef %>%
  mutate(violation = "late", .before = Sample.Name) -> Late_Ofav_reef
Late_Ofav_urban<-Ofav_urban[which(Ofav_urban$D.CT.mean>35), ] 
Late_Ofav_urban %>%
  mutate(violation = "late", .before = Sample.Name) -> Late_Ofav_urban

Late_Cnat<-Cnat[which(Cnat$B.CT.mean>35 & Cnat$D.CT.mean>35), ] 
Late_Cnat %>%
  mutate(violation = "late", .before = Sample.Name) -> Late_Cnat

Late_Pseu<-Pseu[which(Pseu$B.CT.mean>35 & Pseu$D.CT.mean>35), ] 
Late_Pseu %>%
  mutate(violation = "late", .before = Sample.Name) -> Late_Pseu

Late_Dlab_reef<-Dlab_reef[which(Dlab_reef$A.CT.mean>35 & Dlab_reef$B.CT.mean>35 & Dlab_reef$D.CT.mean>35), ] 
Late_Dlab_reef %>%
  mutate(violation = "late", .before = Sample.Name) -> Late_Dlab_reef
Late_Dlab_urban<-Dlab_urban[which(Dlab_urban$D.CT.mean>35), ] 
Late_Dlab_urban %>%
  mutate(violation = "late", .before = Sample.Name) -> Late_Dlab_urban


# RE-RUNS
# Combines all lists above by species and finds distinct samples
ToReRun_Mcav<-rbind(Reps_Mcav_reef, Reps_Mcav_urban, StDe1.5_Mcav_reef, StDe1.5_Mcav_urban, Late_Mcav_reef, Late_Mcav_urban)
ToReRun_Mcav<-ToReRun_Mcav %>% distinct_at(vars(-violation), .keep_all=T)
ToReRun_Ofav<-rbind(Reps_Ofav_reef, Reps_Ofav_urban, StDe1.5_Ofav_reef, StDe1.5_Ofav_urban, Late_Ofav_reef, Late_Ofav_urban)
ToReRun_Ofav<-ToReRun_Ofav %>% distinct_at(vars(-violation), .keep_all=T)
ToReRun_Cnat<-rbind(Reps_Cnat, StDe1.5_Cnat, Late_Cnat)
ToReRun_Cnat<-ToReRun_Cnat %>% distinct_at(vars(-violation), .keep_all=T)
ToReRun_Pseu<-rbind(Reps_Pseu, StDe1.5_Pseu, Late_Pseu)
ToReRun_Pseu<-ToReRun_Pseu %>% distinct_at(vars(-violation), .keep_all=T)
ToReRun_Dlab<-rbind(Reps_Dlab_reef, Reps_Dlab_urban, StDe1.5_Dlab_reef, StDe1.5_Dlab_urban, Late_Dlab_reef, Late_Dlab_urban)
ToReRun_Dlab<-ToReRun_Dlab %>% distinct_at(vars(-violation), .keep_all=T)

# Write all files to csv, look at each sample that fails, and add a column called 'redo' with a 'y' or 'n' value for each sample
# Commented out so it's not accidentally overwritten
# write.csv(ToReRun_Mcav, file = "ToReRun_Mcav.csv")
# write.csv(ToReRun_Ofav, file = "ToReRun_Ofav.csv")
# write.csv(ToReRun_Cnat, file = "ToReRun_Cnat.csv")
# write.csv(ToReRun_Pseu, file = "ToReRun_Pseu.csv")
# write.csv(ToReRun_Dlab, file = "ToReRun_Dlab.csv")

# Re-importing to make a master list of reruns
ToReRun_Mcav <- read.csv(file = "ToReRun_Mcav.csv", head = T)
ToReRun_Ofav <- read.csv(file = "ToReRun_Ofav.csv", head = T)
ToReRun_Cnat <- read.csv(file = "ToReRun_Cnat.csv", head = T)
ToReRun_Pseu <- read.csv(file = "ToReRun_Pseu.csv", head = T)
ToReRun_Dlab <- read.csv(file = "ToReRun_Dlab.csv", head = T)

# Combining only the samples with 'y' under the column 'failed'
ToReRun <- rbind(ToReRun_Mcav,ToReRun_Ofav,ToReRun_Cnat,ToReRun_Pseu,ToReRun_Dlab)
ToReRun %>%
  filter(failed == 'y') -> ToReRun_Failed

# Adding a new column that counts row numbers, creating a list of samples that have multiple qPCR runs, counting how many times, then keeping only the last attempt
Urban_Prop %>%
  mutate(sort = 1:n()) %>%
  group_by(ID) %>%
  filter(duplicated(ID)|n()>1) %>%
  mutate(count = n()) %>%
  filter(sort == max(sort)) %>%
  distinct(ID, .keep_all = TRUE) -> ReRuns_Done
write.csv(ReRuns_Done, file = "ReRuns_Done.csv")

# Now matching this with the list of failed samples to determine if the last attempt also failed
ReRuns_Done %>%
  inner_join(select(ToReRun_Failed, ID, failed), by = "ID") %>%
  distinct(ID, .keep_all = TRUE) -> ReRuns_Failed
write.csv(ReRuns_Failed, file = "ReRuns_Failed.csv")

ReRuns_Done %>%
  left_join(select(ToReRun_Failed, ID, failed), by = "ID") %>%
  filter(is.na(failed)) -> ReRuns_Worked
# Sanity check: The number of rows in ReRuns_Failed and ReRuns_Worked should sum to ReRuns_Done
ReRuns_Worked %>%
  filter(count==2) -> ReRuns_2Worked # now filtering out any samples run thrice
write.csv(ReRuns_2Worked, file = "ReRuns_Worked.csv")

# Did any samples that were run 3 times actually work the 2nd time?
ReRuns_Failed %>%
  filter(count==3) -> ReRuns_3Attempts

# Now compare it to the other attempts
# Commented out so it's not accidentally overwritten
# ToReRun %>%
#   inner_join(select(ReRuns_3Attempts, ID, sort), by = "ID") %>%
#   arrange(ID) -> ReRuns_Check
# write.csv(ReRuns_Check, file = "ReRuns_Check.csv")
# Add a column 'use' and one run (if any) that worked per sample

# Now read it back in, and modify it to look like the others below
ReRuns_Check <- read.csv("ReRuns_Check.csv", head = T)
ReRuns_Check %>%
  select(ID, Sample.Plate, 32:44,use) %>%
  filter(use == "y") -> Urban_Prop_3Good

# One-and-done samples, add a column called 'use' with all rows having 'y'
Urban_Prop %>%
  anti_join(ToReRun, by = "Sample.Plate") %>%
  mutate(use = "y") -> Urban_Prop_Good

# Rerun samples that worked the 2nd time
ReRuns_2Worked %>%
  select(1:15) %>%
  mutate(use = "y") -> Urban_Prop_2Good

# Merge all the 'good' datasets
Urban_Prop_Merged <- rbind(Urban_Prop_Good, Urban_Prop_2Good, Urban_Prop_3Good)

# Creating a filtered list of remaining reruns
# ToReRun %>%
#   anti_join(Urban_ReRuns, by = "ID") -> ReRuns_Remain

# But what if we already plated some of these samples? Let's check
# AlreadyPlated <- read.csv("AlreadyPlated.csv", head = T)
# AlreadyPlated %>%
#   group_by(ID) -> AlreadyPlated
# ReRuns_Remain %>%
#   anti_join(AlreadyPlated, by = "ID") -> ReRuns_Remain_ForReal
# write.csv(ReRuns_Remain_ForReal, file = "ReRuns_Remain.csv")


#### STATISTICS ####

# Finds duplicated sample IDs with good data
# Commented out so it's not accidentally overwritten
# Urban_Prop_Merged %>%
  # group_by(ID) %>%
  # mutate(count = n()) %>%
#   filter(count > 1) %>%
#   distinct(Sample.Plate, .keep_all = T) %>%
#   arrange(ID) -> Urban_Prop_Duplicates
# write.csv(Urban_Prop_Duplicates, file = "Duplicates_Worked.csv")

# Read back in and filter out any runs that you modified to 'use' = 'n'
Urban_Prop_Duplicates <- read.csv("Duplicates_Worked.csv", head = T)
Urban_Prop_Duplicates %>%
  filter(use == "y") -> Urban_Prop_Duplicates_Use

# Filter out the rest of the data to only include samples with one run
Urban_Prop_Merged %>%
  group_by(ID) %>%
  mutate(count = n()) %>%
  filter(count == 1) -> Urban_Prop_Use

# Now combine the two for a set of unique sample rows
Urban_Prop_Unique <- rbind(Urban_Prop_Use, Urban_Prop_Duplicates_Use)
# Sanity check: do any samples drop when you filter by unique sample IDs?
Urban_Prop_Unique %>%
  distinct(ID, .keep_all = TRUE) -> Urban_Prop_Unique # 452 unique samples
# 2nd Sanity check: How many unique samples should there be?
Urban_Prop %>%
  distinct(ID, .keep_all = TRUE) -> Urban_Prop_Count # 835 unique samples
ReRuns_Failed %>%
  distinct(ID, .keep_all = TRUE) -> ReRuns_Failed_Count # 69 unique samples
# Huh, so where are the ~300 missing samples?

# Ungroups ID since you don't need it as a grouping variable
Urban_Prop_Unique %>%
  ungroup() -> Urban_Prop_Unique

# Filtering out any samples that did not amplify
Urban_Prop_Unique %>%
  filter(totalSym != 0) -> Urban_Prop_Final

# Filtering by species
Urban_Prop_Final %>%
  filter(Species == "Mcav") -> Mcav_Prop
Urban_Prop_Final %>%
  filter(Species == "Ofav") -> Ofav_Prop
Urban_Prop_Final %>%
  filter(Species == "Cnat") -> Cnat_Prop
Urban_Prop_Final %>%
  filter(Species == "Pseu") -> Pseu_Prop
Urban_Prop_Final %>%
  filter(Species == "Dlab") -> Dlab_Prop

# Identifying A-dominated samples for QAQC with TaqMan
Urban_Prop_Final %>% 
  filter(dom == "A" | grepl('A', coDom)) -> domA

Urban_Prop_Final %>% 
  filter(grepl('A', background)) %>%
  group_by(Species) %>%
  sample_n(1) %>%
  filter(Species != "Cnat") -> backA

testA <- rbind(domA, backA)
testA %>%
  arrange(Species, Region, Site, ID) -> testA

write.csv(testA, file = "qaqcA.csv")

# Randomly selecting samples by sites for QAQC with TaqMan
Urban_Prop_Final %>% 
  filter(dom != "A" | grepl('A', coDom))
  group_by(Species, Region) %>%
  sample_n(3) %>%
  ungroup() %>%
  sample_n(28) -> testABCD

testABCD %>%
  arrange(Species, Region, Site, ID) -> testABCD

write.csv(testABCD, file = "qaqcABCD.csv")


#### MCAV PERMANOVA/PCoA ####

# making a matrix of just symbiont proportion data for each species
Mcav_Matrix <- Mcav_Prop[c(10:13)]
Mcav_Site <- make.unique(Mcav_Prop$Site, sep = "_") # making unique row names for dissimilarity matrix
rownames(Mcav_Matrix) <- rownames(Mcav_Site) # setting row names
Mcav_Matrix <- sqrt(Mcav_Matrix) # sqrt transformation for zero-inflated data

# dissimilarity matrix in vegan
Mcav_Diss <- vegdist(Mcav_Matrix, "bray") # using Bray-Curtis dissimilarity
Mcav_Pcoa <- pcoa(Mcav_Diss) # running the PCoA
Mcav_Vectors <- Mcav_Pcoa$vectors # vectors of datapoints for plot

# setting site as a factor for signficance testing
Mcav_Prop$Site <- factor(Mcav_Prop$Site, levels=c("Emerald Reef","Rainbow Reef","Star Island","Belle Isle", "MacArthur North"))
fill.color<-c("#01665e","#5ab4ac","#f6e8c3","#d8b365","#8c510a") # custom color palette

# PCoA plot
pdf(file="Miami Mcav PCoA.pdf", width=8, height=6)
plot(Mcav_Vectors[,1], Mcav_Vectors[,2],col=fill.color[as.numeric(as.factor(Mcav_Prop$Site))], pch=16, xlab=(paste("Coordinate 1 (", round((Mcav_Pcoa$values$Relative_eig[[1]] * 100), 2), "%)", sep = "")), ylab=(paste("Coordinate 2 (", round((Mcav_Pcoa$values$Relative_eig[[2]] * 100), 2), "%)", sep = "")), main="Montastraea cavernosa")
legend("bottomright", legend=c("Emerald Reef","Rainbow Reef","Star Island","Belle Isle", "MacArthur North"), fill = fill.color, bty="n")
dev.off()

# testing for homogeneity of variance among species
Mcav_Disp <- betadisper(Mcav_Diss, group=Mcav_Prop$Site)
permutest(Mcav_Disp, bias.adjust = TRUE, perm = 9999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999
# 
# Response: Distances
# Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
# Groups     4 0.39763 0.099407 2.2481   9999 0.0775 .
# Residuals 48 2.12251 0.044219                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# variance is homogeneous among sites, proceed with PERMANOVA

# PERMANOVA
Mcav_Perm <- adonis2(Mcav_Diss ~ Site + Rain + Depth, Mcav_Prop, permutations = 9999, parallel = getOption("mc.cores"))
Mcav_Perm
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = Mcav_Diss ~ Site + Rain + Depth, data = Mcav_Prop, permutations = 9999, parallel = getOption("mc.cores"))
# Df SumOfSqs      R2      F Pr(>F)    
# Site      4   1.6539 0.34138 6.0516 0.0004 ***
#   Rain      1   0.0070 0.00145 0.1025 0.8401    
# Depth     1   0.0409 0.00844 0.5985 0.5455    
# Residual 46   3.1430 0.64873                  
# Total    52   4.8448 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Only site factor significant

# pairwise PERMANOVA
Mcav_Pair <- pairwise.adonis2(Mcav_Diss ~ Site, data = Mcav_Prop)
Mcav_Pair 

# creating dataframe of pairwise results for plotting below
Mcav_Pair_Out <- bind_rows(Mcav_Pair$`Emerald Reef_vs_Rainbow Reef`, Mcav_Pair$`Emerald Reef_vs_Belle Isle`, Mcav_Pair$`Emerald Reef_vs_MacArthur North`, Mcav_Pair$`Emerald Reef_vs_Star Island`, Mcav_Pair$`Rainbow Reef_vs_Belle Isle`, Mcav_Pair$`Rainbow Reef_vs_MacArthur North`, Mcav_Pair$`Rainbow Reef_vs_Star Island`, Mcav_Pair$`Belle Isle_vs_MacArthur North`, Mcav_Pair$`Belle Isle_vs_Star Island`, Mcav_Pair$`MacArthur North_vs_Star Island`,  .id = "Comparison")
Mcav_Pair_Out


#### OFAV PERMANOVA/PCoA ####

# making a matrix of just symbiont proportion data for each species
Ofav_Matrix <- Ofav_Prop[c(10:13)]
Ofav_Site <- make.unique(Ofav_Prop$Site, sep = "_") # making unique row names for dissimilarity matrix
rownames(Ofav_Matrix) <- rownames(Ofav_Site) # setting row names
Ofav_Matrix <- sqrt(Ofav_Matrix) # sqrt transformation for zero-inflated data

# dissimilarity matrix in vegan
Ofav_Diss <- vegdist(Ofav_Matrix, "bray") # using Bray-Curtis dissimilarity
Ofav_Pcoa <- pcoa(Ofav_Diss) # running the PCoA
Ofav_Vectors <- Ofav_Pcoa$vectors # vectors of datapoints for plot

# setting site as a factor for signficance testing
Ofav_Prop$Site <- factor(Ofav_Prop$Site, levels=c("Emerald Reef","Rainbow Reef","Star Island","MacArthur North"))
fill.color2<-c("#01665e","#5ab4ac","#f6e8c3","#8c510a") # custom color palette

# PCoA plot
pdf(file="Miami Ofav PCoA.pdf", width=8, height=6)
plot(Ofav_Vectors[,1], Ofav_Vectors[,2],col=fill.color2[as.numeric(as.factor(Ofav_Prop$Site))], pch=16, xlab=(paste("Coordinate 1 (", round((Ofav_Pcoa$values$Relative_eig[[1]] * 100), 2), "%)", sep = "")), ylab=(paste("Coordinate 2 (", round((Ofav_Pcoa$values$Relative_eig[[2]] * 100), 2), "%)", sep = "")), main="Orbicella faveolata")
legend("topright", legend=c("Emerald Reef","Rainbow Reef","Star Island", "MacArthur North"), fill = fill.color2, bty="n")
dev.off()

# testing for homogeneity of variance among species
Ofav_Disp <- betadisper(Ofav_Diss, group=Ofav_Prop$Site)
permutest(Ofav_Disp, bias.adjust = TRUE, perm = 9999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups     3 0.99483 0.33161 10.802   9999  1e-04 ***
#   Residuals 65 1.99540 0.03070                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# variance is heterogeneous among sites, so need to determine whether highest variance is associated with site with the lowest sample size

Ofav_Disp_Site <- TukeyHSD(Ofav_Disp)
Ofav_Disp_Site
boxplot(Ofav_Disp) # plots showing  variance among sites
Ofav_Prop %>% 
  group_by(Site) %>%
  tally() # Emerald and Rainbow have the highest variance but also the highest sample size, so good to proceed with PERMANOVA

# PERMANOVA
Ofav_Perm <- adonis2(Ofav_Diss ~ Site + Rain + Depth, Ofav_Prop, permutations = 9999, parallel = getOption("mc.cores"))
Ofav_Perm
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = Ofav_Diss ~ Site + Rain + Depth, data = Ofav_Prop, permutations = 9999, parallel = getOption("mc.cores"))
# Df SumOfSqs      R2       F Pr(>F)    
# Site      3   4.3806 0.39236 16.5960 0.0001 ***
#   Rain      1   1.1625 0.10412 13.2120 0.0003 ***
#   Depth     1   0.0787 0.00705  0.8945 0.4010    
# Residual 63   5.5431 0.49648                   
# Total    68  11.1649 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Site and Rain factors signficant

# pairwise PERMANOVA
Ofav_Pair <- pairwise.adonis2(Ofav_Diss ~ Site, data = Ofav_Prop)
Ofav_Pair 

# creating dataframe of pairwise results for plotting below
Ofav_Pair_Out <- bind_rows(Ofav_Pair$`Emerald Reef_vs_Rainbow Reef`, Ofav_Pair$`Emerald Reef_vs_MacArthur North`, Ofav_Pair$`Emerald Reef_vs_Star Island`, Ofav_Pair$`Rainbow Reef_vs_MacArthur North`, Ofav_Pair$`Rainbow Reef_vs_Star Island`, Ofav_Pair$`MacArthur North_vs_Star Island`, .id = "Comparison")
Ofav_Pair_Out


#### CNAT PERMANOVA/PCoA ####

# Need to first combine Emerald and Rainbow as reef since sample size is low
Cnat_Prop %>% 
  mutate(across('Site', str_replace, 'Emerald Reef|Rainbow Reef', 'reef')) -> Cnat_Prop

# making a matrix of just symbiont proportion data for each species
Cnat_Matrix <- Cnat_Prop[c(10:13)]
Cnat_Site <- make.unique(Cnat_Prop$Site, sep = "_") # making unique row names for dissimilarity matrix
rownames(Cnat_Matrix) <- rownames(Cnat_Site) # setting row names
Cnat_Matrix <- sqrt(Cnat_Matrix) # sqrt transformation for zero-inflated data

# dissimilarity matrix in vegan
Cnat_Diss <- vegdist(Cnat_Matrix, "bray") # using Bray-Curtis dissimilarity
Cnat_Pcoa <- pcoa(Cnat_Diss) # running the PCoA
Cnat_Vectors <- Cnat_Pcoa$vectors # vectors of datapoints for plot

# setting site as a factor for signficance testing
Cnat_Prop$Site <- factor(Cnat_Prop$Site, levels=c("reef","Star Island","Belle Isle", "MacArthur North"))
fill.color3<-c("#003c30","#f6e8c3","#d8b365","#8c510a") # custom color palette

# PCoA plot
pdf(file="Miami Cnat PCoA.pdf", width=8, height=6)
plot(Cnat_Vectors[,1], Cnat_Vectors[,2],col=fill.color3[as.numeric(as.factor(Cnat_Prop$Site))], pch=16, xlab=(paste("Coordinate 1 (", round((Cnat_Pcoa$values$Relative_eig[[1]] * 100), 2), "%)", sep = "")), ylab=(paste("Coordinate 2 (", round((Cnat_Pcoa$values$Relative_eig[[2]] * 100), 2), "%)", sep = "")), main="Colpophyllia natans")
legend("bottomleft", legend=c("reef","Star Island","Belle Isle", "MacArthur North"), fill = fill.color3, bty="n")
dev.off()

# testing for homogeneity of variance among species
Cnat_Disp <- betadisper(Cnat_Diss, group=Cnat_Prop$Site)
permutest(Cnat_Disp, bias.adjust = TRUE, perm = 9999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999
# 
# Response: Distances
# Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups     3 0.83445 0.278150 13.826   9999  2e-04 ***
#   Residuals 64 1.28757 0.020118                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# variance is heterogeneous among sites, so need to determine whether highest variance is associated with site with the lowest sample size

Cnat_Disp_Site <- TukeyHSD(Cnat_Disp)
Cnat_Disp_Site
boxplot(Cnat_Disp) # plots showing  variance among sites
Cnat_Prop %>% 
  group_by(Site) %>%
  tally() # Star Island has the highest variance but also a high sample size, so good to proceed with PERMANOVA

# PERMANOVA
Cnat_Perm <- adonis2(Cnat_Diss ~ Site + Rain + Depth, Cnat_Prop, permutations = 9999, parallel = getOption("mc.cores"))
Cnat_Perm
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = Cnat_Diss ~ Site + Rain + Depth, data = Cnat_Prop, permutations = 9999, parallel = getOption("mc.cores"))
# Df SumOfSqs       R2       F Pr(>F)    
# Site      3   1.4054  0.29643  8.7382 0.0004 ***
#   Rain      1   0.0186  0.00392  0.3470 0.6318    
# Depth     1  -0.0068 -0.00143 -0.1261 0.9823    
# Residual 62   3.3238  0.70108                   
# Total    67   4.7411  1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Only site factor significant

# pairwise PERMANOVA
Cnat_Pair <- pairwise.adonis2(Cnat_Diss ~ Site, data = Cnat_Prop)
Cnat_Pair 

# creating dataframe of pairwise results for plotting below
Cnat_Pair_Out <- bind_rows(Cnat_Pair$`reef_vs_Belle Isle`, Cnat_Pair$`reef_vs_MacArthur North`, Cnat_Pair$`reef_vs_Star Island`, Cnat_Pair$`Belle Isle_vs_MacArthur North`, Cnat_Pair$`Belle Isle_vs_Star Island`, Cnat_Pair$`MacArthur North_vs_Star Island`, .id = "Comparison")
Cnat_Pair_Out


#### PSEU PERMANOVA/PCoA ####

# making a matrix of just symbiont proportion data for each species
Pseu_Matrix <- Pseu_Prop[c(10:13)]
Pseu_Site <- make.unique(Pseu_Prop$Site, sep = "_") # making unique row names for dissimilarity matrix
rownames(Pseu_Matrix) <- rownames(Pseu_Site) # setting row names
Pseu_Matrix <- sqrt(Pseu_Matrix) # sqrt transformation for zero-inflated data

# dissimilarity matrix in vegan
Pseu_Diss <- vegdist(Pseu_Matrix, "bray") # using Bray-Curtis dissimilarity
Pseu_Pcoa <- pcoa(Pseu_Diss) # running the PCoA
Pseu_Vectors <- Pseu_Pcoa$vectors # vectors of datapoints for plot

# setting site as a factor for signficance testing
Pseu_Prop$Site <- factor(Pseu_Prop$Site, levels=c("Emerald Reef","Rainbow Reef","Star Island","Belle Isle", "MacArthur North"))

# PCoA plot
pdf(file="Miami Pseu PCoA.pdf", width=8, height=6)
plot(Pseu_Vectors[,1], Pseu_Vectors[,2],col=fill.color[as.numeric(as.factor(Pseu_Prop$Site))], pch=16, xlab=(paste("Coordinate 1 (", round((Pseu_Pcoa$values$Relative_eig[[1]] * 100), 2), "%)", sep = "")), ylab=(paste("Coordinate 2 (", round((Pseu_Pcoa$values$Relative_eig[[2]] * 100), 2), "%)", sep = "")), main="Pseudodiploria spp.")
legend("topright", legend=c("Emerald Reef","Rainbow Reef","Star Island","Belle Isle", "MacArthur North"), fill = fill.color, bty="n")
dev.off()

# testing for homogeneity of variance among species
Pseu_Disp <- betadisper(Pseu_Diss, group=Pseu_Prop$Site)
permutest(Pseu_Disp, bias.adjust = TRUE, perm = 9999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups      4  1.7167 0.42918 6.1812   9999  5e-04 ***
#   Residuals 151 10.4843 0.06943                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# variance is heterogeneous among sites, so need to determine whether highest variance is associated with site with the lowest sample size

Pseu_Disp_Site <- TukeyHSD(Pseu_Disp)
Pseu_Disp_Site
boxplot(Pseu_Disp) # plots showing  variance among sites
Pseu_Prop %>% 
  group_by(Site) %>%
  tally() # Sample sizes roughly equal, so good to proceed with PERMANOVA

# PERMANOVA
Pseu_Perm <- adonis2(Pseu_Diss ~ Site + Rain + Depth, Pseu_Prop, permutations = 9999, parallel = getOption("mc.cores"))
Pseu_Perm
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = Pseu_Diss ~ Site + Rain + Depth, data = Pseu_Prop, permutations = 9999, parallel = getOption("mc.cores"))
# Df SumOfSqs      R2       F Pr(>F)    
# Site       4    7.236 0.22652 11.5825 0.0001 ***
#   Rain       1    1.002 0.03136  6.4145 0.0081 ** 
#   Depth      1    0.435 0.01361  2.7845 0.0864 .  
# Residual 149   23.271 0.72850                   
# Total    155   31.943 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Site and Rain factors signficant

# pairwise PERMANOVA
Pseu_Pair <- pairwise.adonis2(Pseu_Diss ~ Site, data = Pseu_Prop)
Pseu_Pair 

# creating dataframe of pairwise results for plotting below
Pseu_Pair_Out <- bind_rows(Pseu_Pair$`Emerald Reef_vs_Rainbow Reef`, Pseu_Pair$`Emerald Reef_vs_Belle Isle`, Pseu_Pair$`Emerald Reef_vs_MacArthur North`, Pseu_Pair$`Emerald Reef_vs_Star Island`, Pseu_Pair$`Rainbow Reef_vs_Belle Isle`, Pseu_Pair$`Rainbow Reef_vs_MacArthur North`, Pseu_Pair$`Rainbow Reef_vs_Star Island`, Pseu_Pair$`Belle Isle_vs_MacArthur North`, Pseu_Pair$`Belle Isle_vs_Star Island`, Pseu_Pair$`MacArthur North_vs_Star Island`,  .id = "Comparison")
Pseu_Pair_Out


#### DLAB PERMANOVA/PCoA ####

# making a matrix of just symbiont proportion data for each species
Dlab_Matrix <- Dlab_Prop[c(10:13)]
Dlab_Site <- make.unique(Dlab_Prop$Site, sep = "_") # making unique row names for dissimilarity matrix
rownames(Dlab_Matrix) <- rownames(Dlab_Site) # setting row names
Dlab_Matrix <- sqrt(Dlab_Matrix) # sqrt transformation for zero-inflated data

# dissimilarity matrix in vegan
Dlab_Diss <- vegdist(Dlab_Matrix, "bray") # using Bray-Curtis dissimilarity
Dlab_Pcoa <- pcoa(Dlab_Diss) # running the PCoA
Dlab_Vectors <- Dlab_Pcoa$vectors # vectors of datapoints for plot

# setting site as a factor for signficance testing
Dlab_Prop$Site <- factor(Dlab_Prop$Site, levels=c("Emerald Reef","Rainbow Reef","Belle Isle"))
fill.color4<-c("#01665e","#5ab4ac","#d8b365") # custom color palette

# PCoA plot
pdf(file="Miami Dlab PCoA.pdf", width=8, height=6)
plot(Dlab_Vectors[,1], Dlab_Vectors[,2],col=fill.color4[as.numeric(as.factor(Dlab_Prop$Site))], pch=16, xlab=(paste("Coordinate 1 (", round((Dlab_Pcoa$values$Relative_eig[[1]] * 100), 2), "%)", sep = "")), ylab=(paste("Coordinate 2 (", round((Dlab_Pcoa$values$Relative_eig[[2]] * 100), 2), "%)", sep = "")), main="Diploria labyrinthiformis")
legend("topleft", legend=c("Emerald Reef","Rainbow Reef","Belle Isle"), fill = fill.color4, bty="n")
dev.off()

# testing for homogeneity of variance among species
Dlab_Disp <- betadisper(Dlab_Diss, group=Dlab_Prop$Site)
permutest(Dlab_Disp, bias.adjust = TRUE, perm = 9999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999
# 
# Response: Distances
# Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
# Groups     2 0.41539 0.207695 2.9523   9999  0.079 .
# Residuals 19 1.33666 0.070351                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# variance is homogeneous among sites, proceed with PERMANOVA

# PERMANOVA
Dlab_Perm <- adonis2(Dlab_Diss ~ Site + Rain + Depth, Dlab_Prop, permutations = 9999, parallel = getOption("mc.cores"))
Dlab_Perm
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = Dlab_Diss ~ Site + Rain + Depth, data = Dlab_Prop, permutations = 9999, parallel = getOption("mc.cores"))
# Df SumOfSqs      R2      F Pr(>F)
# Site      2   0.3864 0.08478 0.8933 0.4860
# Depth     1   0.2784 0.06109 1.2874 0.2942
# Residual 18   3.8925 0.85413              
# Total    21   4.5573 1.00000   
# No factors signficant


#### RELATIVE ABUNDANCE ####

# creating a custom color palette
getPalette = colorRampPalette(brewer.pal(6, "Accent"))
values = colorRampPalette(brewer.pal(6, "Accent"))(6)


#### MCAV ABUNDANCE ####

# transposing and reformatting dataframes to make abundance a single column
Mcav_Perc <- dplyr::select(Mcav_Prop, 5, 10:13)
Mcav_Perc <- reshape2::melt(Mcav_Perc, id = "Site")
Mcav_Perc$Site=factor(Mcav_Perc$Site, levels=c("Emerald Reef","Rainbow Reef","Star Island","Belle Isle", "MacArthur North")) 
Mcav_Perc <- dcast(Mcav_Perc, Site~variable, mean)
Mcav_Perc <- melt(Mcav_Perc, id = "Site")
Mcav_Perc <- dplyr::rename(Mcav_Perc, symtype = variable, abundance = value)

# creating dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
Mcav_Letters <- data.frame(cbind(Mcav_Pair_Out$Comparison,Mcav_Pair_Out$'Pr(>F)'))
Mcav_Letters <- na.omit(Mcav_Letters)
Mcav_Letters <- dplyr::rename(Mcav_Letters, comparison = X1, p.adj = X2)
Mcav_Letters$comparison = c("Emerald Reef-Rainbow Reef", "Emerald Reef-Belle Isle", "Emerald Reef-MacArthur North", "Emerald Reef-Star Island", "Rainbow Reef-Belle Isle", "Rainbow Reef-MacArthur North", "Rainbow Reef-Star Island", "Belle Isle-MacArthur North", "Belle Isle-Star Island", "MacArthur North-Star Island")
Mcav_Letters$p.adj <- as.numeric(paste(Mcav_Letters$p.adj))
Mcav_Letters

# creates compact letter display of significant pairwise differences for figure
Mcav_Cld <- cldList(p.adj ~ comparison, data = Mcav_Letters, threshold = 0.05, remove.zero = FALSE, remove.space = FALSE)
Mcav_Cld=Mcav_Cld[order(Mcav_Cld$Group),] 
Mcav_Cld$Site <- Mcav_Cld$Group
Mcav_Cld <- Mcav_Cld %>% mutate(symtype="Symbiodinium")
Mcav_Cld

# percent stacked barplot
Mcav_Plot <- ggplot(Mcav_Perc, aes(fill=symtype, y=abundance, x=Site)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Site",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Montastraea cavernosa") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  geom_text(data = Mcav_Cld, aes(y=-0.05, label=Letter)) +
  geom_text(x = 1.8, y=1.035, label = "Site: F4,52 = 6.1, R2 = 0.341, p < 0.001") +
  theme(plot.title = element_text(hjust=0.5))
Mcav_Plot 

# takes the legend and saves it as a separate object (grob)
get_legend <-function(Mcav_Plot){
  tmp <- ggplot_gtable(ggplot_build(Mcav_Plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend=get_legend(Mcav_Plot)

# replotting without legend and x axis label
Mcav_Plot <- ggplot(Mcav_Perc, aes(fill=symtype, y=abundance, x=Site)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Site",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Montastraea cavernosa") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  geom_text(data = Mcav_Cld, aes(y=-0.05, label=Letter)) +
  geom_text(x = 1.8, y=1.035, label = "Site: F4,52 = 6.1, R2 = 0.341, p < 0.001") +
  theme(plot.title = element_text(hjust=0.5), legend.position = "none")
Mcav_Plot 

ggsave("Miami Mcav symbionts.pdf", plot= Mcav_Plot, width=6, height=4, units="in", dpi=300)


#### OFAV ABUNDANCE ####

# transposing and reformatting dataframes to make abundance a single column
Ofav_Perc <- dplyr::select(Ofav_Prop, 5, 10:13)
Ofav_Perc <- reshape2::melt(Ofav_Perc, id = "Site")
Ofav_Perc$Site=factor(Ofav_Perc$Site, levels=c("Emerald Reef","Rainbow Reef","Star Island","Belle Isle", "MacArthur North")) 
Ofav_Perc <- dcast(Ofav_Perc, Site~variable, mean)
Ofav_Perc <- melt(Ofav_Perc, id = "Site")
Ofav_Perc <- dplyr::rename(Ofav_Perc, symtype = variable, abundance = value)

# creating dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
Ofav_Letters <- data.frame(cbind(Ofav_Pair_Out$Comparison,Ofav_Pair_Out$'Pr(>F)'))
Ofav_Letters <- na.omit(Ofav_Letters)
Ofav_Letters <- dplyr::rename(Ofav_Letters, comparison = X1, p.adj = X2)
Ofav_Letters$comparison = c("Emerald Reef-Rainbow Reef", "Emerald Reef-MacArthur North", "Emerald Reef-Star Island", "Rainbow Reef-MacArthur North", "Rainbow Reef-Star Island", "MacArthur North-Star Island")
Ofav_Letters$p.adj <- as.numeric(paste(Ofav_Letters$p.adj))
Ofav_Letters

# creates compact letter display of significant pairwise differences for figure
Ofav_Cld <- cldList(p.adj ~ comparison, data = Ofav_Letters, threshold = 0.05, remove.zero = FALSE, remove.space = FALSE)
Ofav_Cld=Ofav_Cld[order(Ofav_Cld$Group),] 
Ofav_Cld$Site <- Ofav_Cld$Group
Ofav_Cld <- Ofav_Cld %>% mutate(symtype="Symbiodinium")
Ofav_Cld

# percent stacked barplot
Ofav_Plot <- ggplot(Ofav_Perc, aes(fill=symtype, y=abundance, x=Site)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Site",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Orbicella faveolata") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  geom_text(data = Ofav_Cld, aes(y=-0.05, label=Letter)) +
  geom_text(x = 1.8, y=1.035, label = "Site: F3,68 = 16.6, R2 = 0.392, p < 0.001") +
  theme(plot.title = element_text(hjust=0.5), legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())
Ofav_Plot 

ggsave("Miami Ofav symbionts.pdf", plot= Ofav_Plot, width=6, height=4, units="in", dpi=300)


#### CNAT ABUNDANCE ####

# transposing and reformatting dataframes to make abundance a single column
Cnat_Perc <- dplyr::select(Cnat_Prop, 5, 10:13)
Cnat_Perc <- reshape2::melt(Cnat_Perc, id = "Site")
Cnat_Perc$Site=factor(Cnat_Perc$Site, levels=c("reef","Star Island","Belle Isle", "MacArthur North")) 
Cnat_Perc <- dcast(Cnat_Perc, Site~variable, mean)
Cnat_Perc <- melt(Cnat_Perc, id = "Site")
Cnat_Perc <- dplyr::rename(Cnat_Perc, symtype = variable, abundance = value)

# creating dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
Cnat_Letters <- data.frame(cbind(Cnat_Pair_Out$Comparison,Cnat_Pair_Out$'Pr(>F)'))
Cnat_Letters <- na.omit(Cnat_Letters)
Cnat_Letters <- dplyr::rename(Cnat_Letters, comparison = X1, p.adj = X2)
Cnat_Letters$comparison = c("reef-Belle Isle", "reef-MacArthur North", "reef-Star Island", "Belle Isle-MacArthur North", "Belle Isle-Star Island", "MacArthur North-Star Island")
Cnat_Letters$p.adj <- as.numeric(paste(Cnat_Letters$p.adj))
Cnat_Letters

# creates compact letter display of significant pairwise differences for figure
Cnat_Cld <- cldList(p.adj ~ comparison, data = Cnat_Letters, threshold = 0.05, remove.zero = FALSE, remove.space = FALSE)
Cnat_Cld=Cnat_Cld[order(Cnat_Cld$Group),] 
Cnat_Cld$Site <- Cnat_Cld$Group
Cnat_Cld <- Cnat_Cld %>% mutate(symtype="Symbiodinium")
Cnat_Cld

# percent stacked barplot
Cnat_Plot <- ggplot(Cnat_Perc, aes(fill=symtype, y=abundance, x=Site)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Site",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Colpophyllia natans") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  geom_text(data = Cnat_Cld, aes(y=-0.05, label=Letter)) +
  geom_text(x = 1.8, y=1.035, label = "Site: F3,67 = 8.7, R2 = 0.296, p < 0.001") +
  theme(plot.title = element_text(hjust=0.5), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank())
Cnat_Plot 

ggsave("Miami Cnat symbionts.pdf", plot= Cnat_Plot, width=6, height=4, units="in", dpi=300)


#### PSEU ABUNDANCE ####

# transposing and reformatting dataframes to make abundance a single column
Pseu_Perc <- dplyr::select(Pseu_Prop, 5, 10:13)
Pseu_Perc <- reshape2::melt(Pseu_Perc, id = "Site")
Pseu_Perc$Site=factor(Pseu_Perc$Site, levels=c("Emerald Reef","Rainbow Reef","Star Island","Belle Isle", "MacArthur North")) 
Pseu_Perc <- dcast(Pseu_Perc, Site~variable, mean)
Pseu_Perc <- melt(Pseu_Perc, id = "Site")
Pseu_Perc <- dplyr::rename(Pseu_Perc, symtype = variable, abundance = value)

# creating dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
Pseu_Letters <- data.frame(cbind(Pseu_Pair_Out$Comparison,Pseu_Pair_Out$'Pr(>F)'))
Pseu_Letters <- na.omit(Pseu_Letters)
Pseu_Letters <- dplyr::rename(Pseu_Letters, comparison = X1, p.adj = X2)
Pseu_Letters$comparison = c("Emerald Reef-Rainbow Reef", "Emerald Reef-Belle Isle", "Emerald Reef-MacArthur North", "Emerald Reef-Star Island", "Rainbow Reef-Belle Isle", "Rainbow Reef-MacArthur North", "Rainbow Reef-Star Island", "Belle Isle-MacArthur North", "Belle Isle-Star Island", "MacArthur North-Star Island")
Pseu_Letters$p.adj <- as.numeric(paste(Pseu_Letters$p.adj))
Pseu_Letters

# creates compact letter display of significant pairwise differences for figure
Pseu_Cld <- cldList(p.adj ~ comparison, data = Pseu_Letters, threshold = 0.05, remove.zero = FALSE, remove.space = FALSE)
Pseu_Cld=Pseu_Cld[order(Pseu_Cld$Group),] 
Pseu_Cld$Site <- Pseu_Cld$Group
Pseu_Cld <- Pseu_Cld %>% mutate(symtype="Symbiodinium")
Pseu_Cld

# percent stacked barplot
Pseu_Plot <- ggplot(Pseu_Perc, aes(fill=symtype, y=abundance, x=Site)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Site",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Pseudodiploria spp.") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  geom_text(data = Pseu_Cld, aes(y=-0.05, label=Letter)) +
  geom_text(x = 1.8, y=1.035, label = "Site: F4,155 = 11.6, R2 = 0.227, p < 0.001") +
  theme(plot.title = element_text(hjust=0.5), legend.position = "none", axis.title.x = element_blank())
Pseu_Plot 

ggsave("Miami Pseu symbionts.pdf", plot= Pseu_Plot, width=6, height=4, units="in", dpi=300)


#### DLAB ABUNDANCE ####

# transposing and reformatting dataframes to make abundance a single column
Dlab_Perc <- dplyr::select(Dlab_Prop, 5, 10:13)
Dlab_Perc <- reshape2::melt(Dlab_Perc, id = "Site")
Dlab_Perc$Site=factor(Dlab_Perc$Site, levels=c("Emerald Reef","Rainbow Reef","Belle Isle")) 
Dlab_Perc <- dcast(Dlab_Perc, Site~variable, mean)
Dlab_Perc <- melt(Dlab_Perc, id = "Site")
Dlab_Perc <- dplyr::rename(Dlab_Perc, symtype = variable, abundance = value)

# percent stacked barplot
Dlab_Plot <- ggplot(Dlab_Perc, aes(fill=symtype, y=abundance, x=Site)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Site",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Diploria labyrinthiformis") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  geom_text(x = 1.8, y=1.035, label = "Site: F2,21 = 0.9, R2 = 0.085, p = ns") +
  theme(plot.title = element_text(hjust=0.5), legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())
Dlab_Plot 

ggsave("Miami Dlab symbionts.pdf", plot= Dlab_Plot, width=6, height=4, units="in", dpi=300)

#### ABUNDANCE MULTIPLOT ####

Perc_Multiplot <- grid.arrange(Pseu_Plot, Cnat_Plot, Dlab_Plot, Mcav_Plot, Ofav_Plot, legend, ncol=3, nrow=2, widths=c(6,4.5,3.5), heights=c(3.75,4))

ggsave("Miami multiplot symbionts.pdf", plot= Perc_Multiplot, width=14, height=8, units="in", dpi=300)
