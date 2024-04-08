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
library(rcompanion)
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

# Create unique sample ID+FileName to merge with symbiont abundance data  
Urban$Sample.Plate<-paste(Urban$Sample.Name, Urban$File.Name, sep = "_" )

# Importing the symbiont abundance data from Rich
Urban_Prop <- read.csv("symProps_7Apr2024.csv", head=T)
Urban_Prop$Species <- as.factor(Urban_Prop$Species)
Urban_Prop$Site <- as.factor(Urban_Prop$Site)
Urban_Prop$CollectionDate <- as.factor(Urban_Prop$CollectionDate)
Urban_Prop$Season <- as.factor(Urban_Prop$Season)
Urban_Prop$Rain <- as.factor(Urban_Prop$Rain)

# Joining metadata with qPCR data
Urban_Join<-left_join(Urban, Urban_Prop, by="Sample.Plate")
Urban_Join %>%
  filter(!is.na(ID)) -> Urban_Join
# Did we lose any samples due to metadata mismatch?
Urban_Join %>%
  filter(is.na(ID)) -> Urban_Join_Check # 0 samples

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
  summarize(across(Symbiodinium:Durusdinium, mean, na.rm=TRUE)) -> Mcav_dom
Mcav_dom # reef: AC, urban: AD
Urban_Join %>%
  filter(Species == "Ofav") %>%
  group_by(Region) %>%
  summarize(across(Symbiodinium:Durusdinium, mean, na.rm=TRUE)) -> Ofav_dom
Ofav_dom # reef: ABD, urban: D
Urban_Join %>%
  filter(Species == "Cnat") %>%
  group_by(Region) %>%
  summarize(across(Symbiodinium:Durusdinium, mean, na.rm=TRUE)) -> Cnat_dom
Cnat_dom # reef: BD, urban: BD
Urban_Join %>%
  filter(Species == "Pseu") %>%
  group_by(Region) %>%
  summarize(across(Symbiodinium:Durusdinium, mean, na.rm=TRUE)) -> Pseu_dom
Pseu_dom # reef: BD, urban: BD
Urban_Join %>%
  filter(Species == "Dlab") %>%
  group_by(Region) %>%
  summarize(across(Symbiodinium:Durusdinium, mean, na.rm=TRUE)) -> Dlab_dom
Dlab_dom # reef: ABD, urban: D

# Further filtering species by region based on dominant symbionts
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
ToReRun %>%
  inner_join(select(ReRuns_3Attempts, ID, sort), by = "ID") %>%
  arrange(ID) -> ReRuns_Check
# Commented out so it's not accidentally overwritten
# write.csv(ReRuns_Check, file = "ReRuns_Check.csv")
# Add a column 'use' and one run (if any) that worked per sample

# Now read it back in, and modify it to look like the others below
ReRuns_Check <- read.csv("ReRuns_Check.csv", head = T)
ReRuns_Check %>%
  select(ID, Sample.Plate, 32:45, use) %>%
  filter(use == "y") -> Urban_Prop_3Good

# One-and-done samples, add a column called 'use' with all rows having 'y'
Urban_Prop %>%
  anti_join(ToReRun_Failed, by = "Sample.Plate") %>%
  mutate(use = "y") -> Urban_Prop_Good

# Rerun samples that worked the 2nd time
ReRuns_2Worked %>%
  select(1:16) %>%
  mutate(use = "y") -> Urban_Prop_2Good

# Merge all the 'good' datasets
Urban_Prop_Merged <- rbind(Urban_Prop_Good, Urban_Prop_2Good, Urban_Prop_3Good)

# Old code from the samples shipped to Lorelei
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

# Identifying A-dominated samples for QAQC with TaqMan
# Urban_Prop_Final %>% 
#   filter(dom == "A" | grepl('A', coDom)) -> domA
# 
# Urban_Prop_Final %>% 
#   filter(grepl('A', background)) %>%
#   group_by(Species) %>%
#   sample_n(1) %>%
#   filter(Species != "Cnat") -> backA
# 
# testA <- rbind(domA, backA)
# testA %>%
#   arrange(Species, Region, Region, ID) -> testA
# 
# write.csv(testA, file = "qaqcA.csv")

# Randomly selecting samples by regions for QAQC with TaqMan
# Urban_Prop_Final %>% 
#   filter(dom != "A" | grepl('A', coDom))
#   group_by(Species, Region) %>%
#   sample_n(3) %>%
#   ungroup() %>%
#   sample_n(28) -> testABCD
# 
# testABCD %>%
#   arrange(Species, Region, Region, ID) -> testABCD
# 
# write.csv(testABCD, file = "qaqcABCD.csv")

# Finds duplicated sample IDs with good data
Urban_Prop_Merged %>%
group_by(ID) %>%
mutate(count = n()) %>%
  filter(count > 1) %>%
  distinct(Sample.Plate, .keep_all = T) %>%
  arrange(ID) -> Urban_Prop_Duplicates
# Commented out so it's not accidentally overwritten
# write.csv(Urban_Prop_Duplicates, file = "Duplicates_Worked.csv")
# Modify column 'use' to select samples to drop

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
  group_by(ID) %>% 
  filter(n()>1) # 0 samples
Urban_Prop_Unique %>%
  distinct(ID, .keep_all = TRUE) -> Urban_Prop_Unique # 802 unique samples
# 2nd Sanity check: How many unique samples should there be?
Urban_Prop_Merged %>%
  distinct(ID, .keep_all = TRUE) -> Urban_Prop_Merged_Count # 802 unique samples
Urban_Prop_Merged_Count %>%
  anti_join(Urban_Prop_Unique, by = "ID") -> Urban_Prop_Check # 0 samples
Urban_Prop %>%
  distinct(ID, .keep_all = TRUE) -> Urban_Prop_Count # 836 unique samples
nrow(ReRuns_Failed) # 104 unique samples
# 3rd sanity check: How many samples are in the original dataset, but not in the filtered dataset?
Urban_Prop %>%
  distinct(ID, .keep_all = TRUE) %>%
  # filter(totalSym != 0) %>%
  anti_join(Urban_Prop_Unique, by = "ID") -> Urban_Prop_Missing # 34 unique samples

# Ungroups ID since you don't need it as a grouping variable
Urban_Prop_Unique %>%
  ungroup() -> Urban_Prop_Unique

# Filtering out any samples that did not amplify
Urban_Prop_Unique %>%
  filter(totalSym != 0) -> Urban_Prop_Final
write.csv(Urban_Prop_Final, file = "symProps_final.csv")


#### STATISTICS ####

Urban_Prop_Final <- read.csv("symProps_final.csv", head=T)
Urban_Prop_Final$Species <- as.factor(Urban_Prop_Final$Species)
Urban_Prop_Final$Site <- as.factor(Urban_Prop_Final$Site)
Urban_Prop_Final$Season <- as.factor(Urban_Prop_Final$Season)
Urban_Prop_Final$Rain <- as.factor(Urban_Prop_Final$Rain)

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

# What is the sample size distribution among sites?
Mcav_Prop %>% 
  group_by(Region,Site) %>%
  tally() 
Ofav_Prop %>% 
  group_by(Region,Site) %>%
  tally() 
Cnat_Prop %>% 
  group_by(Region,Site) %>%
  tally() 
Pseu_Prop %>% 
  group_by(Region,Site) %>%
  tally() 
Dlab_Prop %>% 
  group_by(Region,Site) %>%
  tally() 
# A lot of sites have low sample sizes, so combining by region for statistical testing
# Will also create abundance figures by site

# Need to ungroup Region for PCoAs below
Mcav_Prop %>%
  ungroup() -> Mcav_Prop
Ofav_Prop %>%
  ungroup() -> Ofav_Prop
Cnat_Prop %>%
  ungroup() -> Cnat_Prop
Pseu_Prop %>%
  ungroup() -> Pseu_Prop
Dlab_Prop %>%
  ungroup() -> Dlab_Prop


#### MCAV PERMANOVA/PCoA ####

# making a matrix of just symbiont proportion data for each species
Mcav_Matrix <- Mcav_Prop[c(12:14)]
Mcav_Region <- make.unique(Mcav_Prop$Region, sep = "_") # making unique row names for dissimilarity matrix
rownames(Mcav_Matrix) <- Mcav_Region # setting row names
Mcav_Matrix <- sqrt(Mcav_Matrix) # sqrt transformation for zero-inflated data

# dissimilarity matrix in vegan
Mcav_Diss <- vegdist(Mcav_Matrix, "bray") # using Bray-Curtis dissimilarity
Mcav_Pcoa <- pcoa(Mcav_Diss) # running the PCoA
Mcav_Vectors <- Mcav_Pcoa$vectors # vectors of datapoints for plot

# setting region as a factor for significance testing
Mcav_Prop$Region <- factor(Mcav_Prop$Region, levels=c("Palm Beach urban", "Broward reef", "North Miami reef", "North Miami urban", "Miami reef", "Miami urban"))
fill.color<-c("#fddbc7", "#e0e0e0", "#999999", "#ef8a62", "#4d4d4d", "#b2182b") # custom color palette

# PCoA plot
pdf(file="Mcav PCoA.pdf", width=8, height=6)
plot(Mcav_Vectors[,1], Mcav_Vectors[,2],col=fill.color[as.numeric(as.factor(Mcav_Prop$Region))], pch=16, xlab=(paste("Coordinate 1 (", round((Mcav_Pcoa$values$Rel_corr_eig[[1]] * 100), 2), "%)", sep = "")), ylab=(paste("Coordinate 2 (", round((Mcav_Pcoa$values$Rel_corr_eig[[2]] * 100), 2), "%)", sep = "")), main="Montastraea cavernosa")
legend("bottomright", legend=c("Palm Beach urban", "Broward reef", "North Miami reef", "North Miami urban", "Miami reef", "Miami urban"), fill = fill.color, bty="n")
dev.off()

# testing for homogeneity of variance among species
Mcav_Disp <- betadisper(Mcav_Diss, group=Mcav_Prop$Region)
permutest(Mcav_Disp, bias.adjust = TRUE, perm = 9999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999
# 
# Response: Distances
# Df Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups      5 1.4271 0.285411 5.9691   9999  5e-04 ***
#   Residuals 125 5.9768 0.047815                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# variance is heterogeneous among regions, so need to determine whether highest variance is associated with region with the lowest sample size
Mcav_Disp_Region <- TukeyHSD(Mcav_Disp)
Mcav_Disp_Region
boxplot(Mcav_Disp) # plots showing  variance among regions
Mcav_Prop %>% 
  group_by(Region) %>%
  tally() # North Miami reef has the highest variance but also a decent sample size, so good to proceed with PERMANOVA

# PERMANOVA
Mcav_Perm <- adonis2(Mcav_Diss ~ Region + Rain + Depth, Mcav_Prop, permutations = 9999, parallel = getOption("mc.cores"))
Mcav_Perm
write.csv(Mcav_Perm, file = "Mcav PERMANOVA.csv")
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = Mcav_Diss ~ Region + Rain + Depth, data = Mcav_Prop, permutations = 9999, parallel = getOption("mc.cores"))
# Df SumOfSqs      R2       F Pr(>F)    
# Region     5  19.6861 0.71342 64.1905 0.0001 ***
#   Rain       1   0.2902 0.01052  4.7319 0.0279 *  
#   Depth      1   0.0734 0.00266  1.1973 0.2768    
# Residual 123   7.5444 0.27341                   
# Total    130  27.5942 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Region and Rain factors significant

# pairwise PERMANOVA
Mcav_Pair <- pairwise.adonis2(Mcav_Diss ~ Region, data = Mcav_Prop)
Mcav_Pair 

# creating dataframe of pairwise results for plotting below
Mcav_Pair_Out <- bind_rows(Mcav_Pair$`Miami reef_vs_Miami urban`, Mcav_Pair$`Miami reef_vs_Broward reef`, Mcav_Pair$`Miami reef_vs_North Miami reef`, Mcav_Pair$`Miami reef_vs_Palm Beach urban`, Mcav_Pair$`Miami reef_vs_North Miami urban`, Mcav_Pair$`Miami urban_vs_Broward reef`, Mcav_Pair$`Miami urban_vs_North Miami reef`, Mcav_Pair$`Miami urban_vs_Palm Beach urban`, Mcav_Pair$`Miami urban_vs_North Miami urban`, Mcav_Pair$`Broward reef_vs_North Miami reef`, Mcav_Pair$`Broward reef_vs_Palm Beach urban`, Mcav_Pair$`Broward reef_vs_North Miami urban`, Mcav_Pair$`North Miami reef_vs_Palm Beach urban`, Mcav_Pair$`North Miami reef_vs_North Miami urban`, Mcav_Pair$`Palm Beach urban_vs_North Miami urban`,  .id = "Comparison")
Mcav_Pair_Out
write.csv(Mcav_Pair_Out, file = "Mcav PERMANOVA pairwise.csv")


#### OFAV PERMANOVA/PCoA ####

# since Ofav doesn't have a lot of samples in some regions, need to create a new grouping variable for testing
Ofav_Prop %>% mutate(RegionNew = case_when(
  Region == "Palm Beach urban"  ~ "urban",
  Region == "Broward reef"  ~ "Broward reef",
  Region == "North Miami reef"  ~ "Miami reef",
  Region == "North Miami urban"  ~ "urban",
  Region == "Miami reef"  ~ "Miami reef",
  Region == "Miami urban"  ~ "urban")) -> Ofav_Prop

# making a matrix of just symbiont proportion data for each species
Ofav_Matrix <- Ofav_Prop[c(12:14)]
Ofav_Region <- make.unique(Ofav_Prop$RegionNew, sep = "_") # making unique row names for dissimilarity matrix
rownames(Ofav_Matrix) <- Ofav_Region # setting row names
Ofav_Matrix <- sqrt(Ofav_Matrix) # sqrt transformation for zero-inflated data

# dissimilarity matrix in vegan
Ofav_Diss <- vegdist(Ofav_Matrix, "bray") # using Bray-Curtis dissimilarity
Ofav_Pcoa <- pcoa(Ofav_Diss) # running the PCoA
Ofav_Vectors <- Ofav_Pcoa$vectors # vectors of datapoints for plot

# setting region as a factor for signficance testing
Ofav_Prop$RegionNew <- factor(Ofav_Prop$RegionNew, levels=c("Broward reef", "Miami reef", "urban"))
fill.color2<-c("#999999", "#4d4d4d", "#b2182b") # custom color palette

# PCoA plot
pdf(file="Ofav PCoA.pdf", width=8, height=6)
plot(Ofav_Vectors[,1], Ofav_Vectors[,2],col=fill.color2[as.numeric(as.factor(Ofav_Prop$RegionNew))], pch=16, xlab=(paste("Coordinate 1 (", round((Ofav_Pcoa$values$Rel_corr_eig[[1]] * 100), 2), "%)", sep = "")), ylab=(paste("Coordinate 2 (", round((Ofav_Pcoa$values$Rel_corr_eig[[2]] * 100), 2), "%)", sep = "")), main="Orbicella faveolata")
legend("topright", legend=c("Broward reef", "Miami reef", "urban"), fill = fill.color2, bty="n")
dev.off()

# testing for homogeneity of variance among species
Ofav_Disp <- betadisper(Ofav_Diss, group=Ofav_Prop$RegionNew)
permutest(Ofav_Disp, bias.adjust = TRUE, perm = 9999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups      2 1.6944 0.84720 16.151   9999  1e-04 ***
#   Residuals 120 6.2946 0.05246                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# variance is heterogeneous among regions, so need to determine whether highest variance is associated with region with the lowest sample size
Ofav_Disp_Region <- TukeyHSD(Ofav_Disp)
Ofav_Disp_Region
boxplot(Ofav_Disp) # plots showing  variance among regions
Ofav_Prop %>% 
  group_by(RegionNew) %>%
  tally() # Miami reef has the highest variance but also a high sample size, so good to proceed with PERMANOVA

# PERMANOVA
Ofav_Perm <- adonis2(Ofav_Diss ~ RegionNew + Rain + Depth, Ofav_Prop, permutations = 9999, parallel = getOption("mc.cores"))
Ofav_Perm
write.csv(Ofav_Perm, file = "Ofav PERMANOVA.csv")
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = Ofav_Diss ~ RegionNew + Rain + Depth, data = Ofav_Prop, permutations = 9999, parallel = getOption("mc.cores"))
# Df SumOfSqs      R2       F Pr(>F)    
# RegionNew   2  10.5178 0.48821 56.7010 0.0001 ***
#   Rain        1   0.0166 0.00077  0.1788 0.6940    
# Depth       1   0.0648 0.00301  0.6988 0.4118    
# Residual  118  10.9442 0.50801                   
# Total     122  21.5434 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# RegionNew factor signficant

# pairwise PERMANOVA
Ofav_Pair <- pairwise.adonis2(Ofav_Diss ~ RegionNew, data = Ofav_Prop)
Ofav_Pair 

# creating dataframe of pairwise results for plotting below
Ofav_Pair_Out <- bind_rows(Ofav_Pair$`Miami reef_vs_urban`, Ofav_Pair$`Miami reef_vs_Broward reef`, Ofav_Pair$`urban_vs_Broward reef`, .id = "Comparison")
Ofav_Pair_Out
write.csv(Ofav_Pair_Out, file = "Ofav PERMANOVA pairwise.csv")


#### CNAT PERMANOVA/PCoA ####

# since Cnat doesn't have a lot of samples in some regions, need to create a new grouping variable for testing
Cnat_Prop %>% mutate(RegionNew = case_when(
  Region == "Palm Beach urban"  ~ "Palm Beach urban",
  Region == "Broward reef"  ~ "reef",
  Region == "North Miami reef"  ~ "reef",
  Region == "North Miami urban"  ~ "Miami Dade urban",
  Region == "Miami reef"  ~ "reef",
  Region == "Miami urban"  ~ "Miami Dade urban")) -> Cnat_Prop

# making a matrix of just symbiont proportion data for each species
Cnat_Matrix <- Cnat_Prop[c(12:14)]
Cnat_Region <- make.unique(Cnat_Prop$RegionNew, sep = "_") # making unique row names for dissimilarity matrix
rownames(Cnat_Matrix) <- Cnat_Region # setting row names
Cnat_Matrix <- sqrt(Cnat_Matrix) # sqrt transformation for zero-inflated data

# dissimilarity matrix in vegan
Cnat_Diss <- vegdist(Cnat_Matrix, "bray") # using Bray-Curtis dissimilarity
Cnat_Pcoa <- pcoa(Cnat_Diss) # running the PCoA
Cnat_Vectors <- Cnat_Pcoa$vectors # vectors of datapoints for plot

# setting region as a factor for signficance testing
Cnat_Prop$RegionNew <- factor(Cnat_Prop$RegionNew, levels=c("reef", "Palm Beach urban", "Miami Dade urban"))
fill.color3<-c("#4d4d4d","#fddbc7", "#b2182b") # custom color palette

# PCoA plot
pdf(file="Cnat PCoA.pdf", width=8, height=6)
plot(Cnat_Vectors[,1], Cnat_Vectors[,2],col=fill.color3[as.numeric(as.factor(Cnat_Prop$RegionNew))], pch=16, xlab=(paste("Coordinate 1 (", round((Cnat_Pcoa$values$Rel_corr_eig[[1]] * 100), 2), "%)", sep = "")), ylab=(paste("Coordinate 2 (", round((Cnat_Pcoa$values$Rel_corr_eig[[2]] * 100), 2), "%)", sep = "")), main="Colpophyllia natans")
legend("topright", legend=c("reef", "Palm Beach urban", "Miami Dade urban"), fill = fill.color3, bty="n")
dev.off()

# testing for homogeneity of variance among species
Cnat_Disp <- betadisper(Cnat_Diss, group=Cnat_Prop$RegionNew)
permutest(Cnat_Disp, bias.adjust = TRUE, perm = 9999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)  
# Groups      2  0.6553 0.32767 3.9445   9999 0.0222 *
#   Residuals 135 11.2145 0.08307                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# variance is heterogeneous among regions, so need to determine whether highest variance is associated with region with the lowest sample size
Cnat_Disp_Region <- TukeyHSD(Cnat_Disp)
Cnat_Disp_Region
boxplot(Cnat_Disp) # plots showing  variance among regions
Cnat_Prop %>% 
  group_by(RegionNew) %>%
  tally() # reef has the highest variance but does not have the smallest sample size, so good to proceed with PERMANOVA

# PERMANOVA
Cnat_Perm <- adonis2(Cnat_Diss ~ RegionNew + Rain + Depth, Cnat_Prop, permutations = 9999, parallel = getOption("mc.cores"))
Cnat_Perm
write.csv(Cnat_Perm, file = "Cnat PERMANOVA.csv")
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = Cnat_Diss ~ RegionNew + Rain + Depth, data = Cnat_Prop, permutations = 9999, parallel = getOption("mc.cores"))
# Df SumOfSqs       R2       F Pr(>F)    
# RegionNew   2   4.0949  0.26996 31.1861 0.0001 ***
#   Rain        1   2.3633  0.15581 35.9978 0.0001 ***
#   Depth       1  -0.0217 -0.00143 -0.3309 0.9994    
# Residual  133   8.7318  0.57566                   
# Total     137  15.1683  1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# RegionNew and Rain factors significant

# pairwise PERMANOVA
Cnat_Pair <- pairwise.adonis2(Cnat_Diss ~ RegionNew, data = Cnat_Prop)
Cnat_Pair 

# creating dataframe of pairwise results for plotting below
Cnat_Pair_Out <- bind_rows(Cnat_Pair$`reef_vs_Miami Dade urban`, Cnat_Pair$`reef_vs_Palm Beach urban`, Cnat_Pair$`Miami Dade urban_vs_Palm Beach urban`, .id = "Comparison")
Cnat_Pair_Out
write.csv(Cnat_Pair_Out, file = "Cnat PERMANOVA pairwise.csv")


#### PSEU PERMANOVA/PCoA ####

# making a matrix of just symbiont proportion data for each species
Pseu_Matrix <- Pseu_Prop[c(12:14)]
Pseu_Region <- make.unique(Pseu_Prop$Region, sep = "_") # making unique row names for dissimilarity matrix
rownames(Pseu_Matrix) <- Pseu_Region # setting row names
Pseu_Matrix <- sqrt(Pseu_Matrix) # sqrt transformation for zero-inflated data

# dissimilarity matrix in vegan
Pseu_Diss <- vegdist(Pseu_Matrix, "bray") # using Bray-Curtis dissimilarity
Pseu_Pcoa <- pcoa(Pseu_Diss) # running the PCoA
Pseu_Vectors <- Pseu_Pcoa$vectors # vectors of datapoints for plot

# setting region as a factor for signficance testing
Pseu_Prop$Region <- factor(Pseu_Prop$Region, levels=c("Palm Beach urban", "Broward reef", "North Miami reef", "North Miami urban", "Miami reef", "Miami urban"))

# PCoA plot
pdf(file="Pseu PCoA.pdf", width=8, height=6)
plot(Pseu_Vectors[,1], Pseu_Vectors[,2],col=fill.color[as.numeric(as.factor(Pseu_Prop$Region))], pch=16, xlab=(paste("Coordinate 1 (", round((Pseu_Pcoa$values$Rel_corr_eig[[1]] * 100), 2), "%)", sep = "")), ylab=(paste("Coordinate 2 (", round((Pseu_Pcoa$values$Rel_corr_eig[[2]] * 100), 2), "%)", sep = "")), main="Pseudodiploria spp.")
legend("topleft", legend=c("Palm Beach urban", "Broward reef", "North Miami reef", "North Miami urban", "Miami reef", "Miami urban"), fill = fill.color, bty="n")
dev.off()

# testing for homogeneity of variance among species
Pseu_Disp <- betadisper(Pseu_Diss, group=Pseu_Prop$Region)
permutest(Pseu_Disp, bias.adjust = TRUE, perm = 9999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
# Groups      5  1.775 0.35492 2.5382   9999 0.0291 *
#   Residuals 317 44.326 0.13983                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# variance is heterogeneous among regions, so need to determine whether highest variance is associated with region with the lowest sample size
Pseu_Disp_Region <- TukeyHSD(Pseu_Disp)
Pseu_Disp_Region
boxplot(Pseu_Disp) # plots showing  variance among regions
Pseu_Prop %>% 
  group_by(Region) %>%
  tally() # Highest variance is not necessarily associated with lowest sample sizes, so good to proceed with PERMANOVA

# PERMANOVA
Pseu_Perm <- adonis2(Pseu_Diss ~ Region + Rain + Depth, Pseu_Prop, permutations = 9999, parallel = getOption("mc.cores"))
Pseu_Perm
write.csv(Pseu_Perm, file = "Pseu PERMANOVA.csv")
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = Pseu_Diss ~ Region + Rain + Depth, data = Pseu_Prop, permutations = 9999, parallel = getOption("mc.cores"))
# Df SumOfSqs      R2       F Pr(>F)    
# Region     5    9.518 0.13790 10.2326 0.0001 ***
#   Rain       1    0.105 0.00153  0.5663 0.4627    
# Depth      1    0.797 0.01155  4.2868 0.0407 *  
#   Residual 315   58.598 0.84902                   
# Total    322   69.018 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Region and Depth factors significant

# pairwise PERMANOVA
Pseu_Pair <- pairwise.adonis2(Pseu_Diss ~ Region, data = Pseu_Prop)
Pseu_Pair 

# creating dataframe of pairwise results for plotting below
Pseu_Pair_Out <- bind_rows(Pseu_Pair$`Miami reef_vs_Miami urban`, Pseu_Pair$`Miami reef_vs_Broward reef`, Pseu_Pair$`Miami reef_vs_Palm Beach urban`, Pseu_Pair$`Miami reef_vs_North Miami reef`, Pseu_Pair$`Miami reef_vs_North Miami urban`, Pseu_Pair$`Miami urban_vs_Broward reef`, Pseu_Pair$`Miami urban_vs_Palm Beach urban`, Pseu_Pair$`Miami urban_vs_North Miami reef`, Pseu_Pair$`Miami urban_vs_North Miami urban`, Pseu_Pair$`Broward reef_vs_Palm Beach urban`, Pseu_Pair$`Broward reef_vs_North Miami reef`, Pseu_Pair$`Broward reef_vs_North Miami urban`, Pseu_Pair$`Palm Beach urban_vs_North Miami reef`, Pseu_Pair$`Palm Beach urban_vs_North Miami urban`, Pseu_Pair$`North Miami reef_vs_North Miami urban`,  .id = "Comparison")
Pseu_Pair_Out
write.csv(Pseu_Pair_Out, file = "Pseu PERMANOVA pairwise.csv")


#### DLAB PERMANOVA/PCoA ####

# since Dlab doesn't have a lot of samples in some regions, need to create a new grouping variable for testing
Dlab_Prop %>% mutate(RegionNew = case_when(
  Region == "Palm Beach urban"  ~ "Palm Beach urban",
  Region == "Broward reef"  ~ "reef",
  Region == "North Miami reef"  ~ "reef",
  Region == "Miami reef"  ~ "reef",
  Region == "Miami urban"  ~ "Miami urban")) -> Dlab_Prop

# making a matrix of just symbiont proportion data for each species
Dlab_Matrix <- Dlab_Prop[c(12:14)]
Dlab_Region <- make.unique(Dlab_Prop$RegionNew, sep = "_") # making unique row names for dissimilarity matrix
rownames(Dlab_Matrix) <- Dlab_Region # setting row names
Dlab_Matrix <- sqrt(Dlab_Matrix) # sqrt transformation for zero-inflated data

# dissimilarity matrix in vegan
Dlab_Diss <- vegdist(Dlab_Matrix, "bray") # using Bray-Curtis dissimilarity
Dlab_Pcoa <- pcoa(Dlab_Diss) # running the PCoA
Dlab_Vectors <- Dlab_Pcoa$vectors # vectors of datapoints for plot

# setting region as a factor for signficance testing
Dlab_Prop$RegionNew <- factor(Dlab_Prop$RegionNew, levels=c("reef","Palm Beach urban","Miami urban"))

# PCoA plot
pdf(file="Dlab PCoA.pdf", width=8, height=6)
plot(Dlab_Vectors[,1], Dlab_Vectors[,2],col=fill.color3[as.numeric(as.factor(Dlab_Prop$RegionNew))], pch=16, xlab=(paste("Coordinate 1 (", round((Dlab_Pcoa$values$Rel_corr_eig[[1]] * 100), 2), "%)", sep = "")), ylab=(paste("Coordinate 2 (", round((Dlab_Pcoa$values$Rel_corr_eig[[2]] * 100), 2), "%)", sep = "")), main="Diploria labyrinthiformis")
legend("topleft", legend=c("reef","Palm Beach urban","Miami urban"), fill = fill.color3, bty="n")
dev.off()

# testing for homogeneity of variance among species
Dlab_Disp <- betadisper(Dlab_Diss, group=Dlab_Prop$RegionNew)
permutest(Dlab_Disp, bias.adjust = TRUE, perm = 9999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups     2 1.6812 0.84059 18.901   9999  1e-04 ***
#   Residuals 80 3.5578 0.04447                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# variance is heterogeneous among regions, so need to determine whether highest variance is associated with region with the lowest sample size
Dlab_Disp_Region <- TukeyHSD(Dlab_Disp)
Dlab_Disp_Region
boxplot(Dlab_Disp) # plots showing  variance among regions
Dlab_Prop %>% 
  group_by(RegionNew) %>%
  tally() # reef has the highest variance but not the lowest sample size, so good to proceed with PERMANOVA

# PERMANOVA
Dlab_Perm <- adonis2(Dlab_Diss ~ RegionNew + Rain + Depth, Dlab_Prop, permutations = 9999, parallel = getOption("mc.cores"))
Dlab_Perm
write.csv(Dlab_Perm, file = "Dlab PERMANOVA.csv")
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = Dlab_Diss ~ RegionNew + Rain + Depth, data = Dlab_Prop, permutations = 9999, parallel = getOption("mc.cores"))
# Df SumOfSqs      R2       F Pr(>F)    
# RegionNew  2   1.3994 0.22875 12.0764 0.0002 ***
#   Rain       1   0.1751 0.02863  3.0228 0.0672 .  
# Depth      1   0.0237 0.00387  0.4085 0.6034    
# Residual  78   4.5193 0.73875                   
# Total     82   6.1176 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Region factor significant, Rain marginal

# pairwise PERMANOVA
Dlab_Pair <- pairwise.adonis2(Dlab_Diss ~ RegionNew, data = Dlab_Prop)
Dlab_Pair 

# creating dataframe of pairwise results for plotting below
Dlab_Pair_Out <- bind_rows(Dlab_Pair$`reef_vs_Miami urban`, Dlab_Pair$`reef_vs_Palm Beach urban`, Dlab_Pair$`Miami urban_vs_Palm Beach urban`, .id = "Comparison")
Dlab_Pair_Out
write.csv(Dlab_Pair_Out, file = "Dlab PERMANOVA pairwise.csv")


#### RELATIVE ABUNDANCE ####

# creating a custom color palette
getPalette = colorRampPalette(brewer.pal(6, "Accent"))
values = colorRampPalette(brewer.pal(6, "Accent"))(6)


#### MCAV ABUNDANCE ####

# transposing and reformatting dataframes to make abundance a single column
Mcav_Perc <- dplyr::select(Mcav_Prop, 5, 11:14)
Mcav_Perc <- reshape2::melt(Mcav_Perc, id = "Region")
Mcav_Perc$Region=factor(Mcav_Perc$Region, levels=c("Palm Beach urban", "Broward reef", "North Miami reef", "North Miami urban", "Miami reef", "Miami urban")) 
Mcav_Perc <- dcast(Mcav_Perc, Region~variable, mean)
Mcav_Perc <- melt(Mcav_Perc, id = "Region")
Mcav_Perc <- dplyr::rename(Mcav_Perc, symtype = variable, abundance = value)

# creating dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
Mcav_Letters <- data.frame(cbind(Mcav_Pair_Out$Comparison,Mcav_Pair_Out$'Pr(>F)'))
Mcav_Letters <- na.omit(Mcav_Letters)
Mcav_Letters %>%
  add_row(X1 = "15", X2 = "1") -> Mcav_Letters
Mcav_Letters <- dplyr::rename(Mcav_Letters, comparison = X1, p.adj = X2)
Mcav_Letters$comparison = c("Miami reef-Miami urban", "Miami reef-Broward reef", "Miami reef-North Miami reef", "Miami reef-Palm Beach urban", "Miami reef-North Miami urban", "Miami urban-Broward reef", "Miami urban-North Miami reef", "Miami urban-Palm Beach urban", "Miami urban-North Miami urban", "Broward reef-North Miami reef", "Broward reef-Palm Beach urban", "Broward reef-North Miami urban", "North Miami reef-Palm Beach urban", "North Miami reef-North Miami urban", "Palm Beach urban-North Miami urban")
Mcav_Letters$p.adj <- as.numeric(paste(Mcav_Letters$p.adj))
Mcav_Letters

# creates compact letter display of significant pairwise differences for figure
Mcav_Cld <- cldList(p.adj ~ comparison, data = Mcav_Letters, threshold = 0.05, remove.zero = FALSE, remove.space = FALSE)
Mcav_Cld=Mcav_Cld[order(Mcav_Cld$Group),] 
Mcav_Cld$Region <- Mcav_Cld$Group
Mcav_Cld <- Mcav_Cld %>% mutate(symtype="Durusdinium")
Mcav_Cld

# percent stacked barplot
Mcav_Plot <- ggplot(Mcav_Perc, aes(fill=symtype, y=abundance, x=Region)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Region",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Montastraea cavernosa") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  # ylim(-0.05,1.1) +
  geom_text(data = Mcav_Cld, aes(y=-0.05, label=Letter)) +
  geom_text(x = 1.85, y=1.035, label = "Region: F5,130 = 29.1, R2 = 0.515, p < 0.001") +
  # geom_text(x = 5.25, y=1.035, label = "Rain: F1,130 = 12.9, R2 = 0.046, p < 0.001") +
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
Mcav_Plot <- ggplot(Mcav_Perc, aes(fill=symtype, y=abundance, x=Region)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Region",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Montastraea cavernosa") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  # ylim(-0.05,1.1) +
  geom_text(data = Mcav_Cld, aes(y=-0.05, label=Letter)) +
  geom_text(x = 1.85, y=1.035, label = "Region: F5,130 = 29.1, R2 = 0.515, p < 0.001") +
  # geom_text(x = 5.25, y=1.035, label = "Rain: F1,130 = 12.9, R2 = 0.046, p < 0.001") +
  theme(plot.title = element_text(hjust=0.5), legend.position = "none")
Mcav_Plot 

ggsave("Mcav symbionts.pdf", plot= Mcav_Plot, width=8, height=4, units="in", dpi=300)

# Symbiont abundance by site
Mcav_Perc_Site <- dplyr::select(Mcav_Prop, 6, 11:14)
Mcav_Perc_Site <- reshape2::melt(Mcav_Perc_Site, id = "Site")
Mcav_Perc_Site$Site=factor(Mcav_Perc_Site$Site, levels=c("Peanut Island","T328","BC1","FTL4","South Canyon","Pillars","Graceland","Haulover Inlet","FIU Biscayne Bay","Arch Creek","Emerald Reef","Rainbow Reef","Bill Baggs","Seaquarium","Fisher Island","Coral City Camera","Miami Beach Marina","Star Island","MacArthur North","Belle Isle","FEC Slip Downtown"))
Mcav_Perc_Site <- dcast(Mcav_Perc_Site, Site~variable, mean)
Mcav_Perc_Site <- melt(Mcav_Perc_Site, id = "Site")
Mcav_Perc_Site <- dplyr::rename(Mcav_Perc_Site, symtype = variable, abundance = value)

# Adding region for facet_wrap
Mcav_Perc_Site %>%
  ungroup() %>%
  mutate(Region = case_when(
  Site == "Peanut Island"  ~ "Palm Beach urban",
  Site == "T328"  ~ "Broward reef",
  Site == "BC1"  ~ "Broward reef",
  Site == "FTL4"  ~ "Broward reef",
  Site == "South Canyon"  ~ "North Miami reef",
  Site == "Pillars"  ~ "North Miami reef",
  Site == "Graceland"  ~ "North Miami reef",
  Site == "Haulover Inlet"  ~ "North Miami urban",
  Site == "FIU Biscayne Bay"  ~ "North Miami urban",
  Site == "Arch Creek"  ~ "North Miami urban",
  Site == "Emerald Reef"  ~ "Miami reef",
  Site == "Rainbow Reef"  ~ "Miami reef",
  Site == "Bill Baggs"  ~ "Miami urban",
  Site == "Seaquarium"  ~ "Miami urban",
  Site == "Fisher Island"  ~ "Miami urban",
  Site == "Coral City Camera"  ~ "Miami urban",
  Site == "Miami Beach Marina"  ~ "Miami urban",
  Site == "Star Island"  ~ "Miami urban",
  Site == "MacArthur North"  ~ "Miami urban",
  Site == "Belle Isle"  ~ "Miami urban",
  Site == "FEC Slip Downtown"  ~ "Miami urban")) -> Mcav_Perc_Site
Mcav_Perc_Site$Region=factor(Mcav_Perc_Site$Region, levels=c("Palm Beach urban","Broward reef","North Miami reef","North Miami urban","Miami reef","Miami urban"))

# percent stacked barplot by site
Mcav_Plot_Site <- ggplot(Mcav_Perc_Site, aes(fill=symtype, y=abundance, x=Site)) + 
  geom_bar(position="fill", stat="identity") +
    labs(x = "Site",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Montastraea cavernosa") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5)) +
  facet_grid(~Region, scales = "free", space='free') 
Mcav_Plot_Site 

ggsave("Mcav symbionts by site.pdf", plot= Mcav_Plot_Site, width=20, height=4, units="in", dpi=300)


#### OFAV ABUNDANCE ####

# transposing and reformatting dataframes to make abundance a single column
Ofav_Perc <- dplyr::select(Ofav_Prop, 21, 11:14)
Ofav_Perc <- reshape2::melt(Ofav_Perc, id = "RegionNew")
Ofav_Perc$RegionNew=factor(Ofav_Perc$RegionNew, levels=c("Broward reef", "Miami reef", "urban")) 
Ofav_Perc <- dcast(Ofav_Perc, RegionNew~variable, mean)
Ofav_Perc <- melt(Ofav_Perc, id = "RegionNew")
Ofav_Perc <- dplyr::rename(Ofav_Perc, symtype = variable, abundance = value)

# creating dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
Ofav_Letters <- data.frame(cbind(Ofav_Pair_Out$Comparison,Ofav_Pair_Out$'Pr(>F)'))
Ofav_Letters <- na.omit(Ofav_Letters)
Ofav_Letters <- dplyr::rename(Ofav_Letters, comparison = X1, p.adj = X2)
Ofav_Letters$comparison = c("Miami reef-urban", "Miami reef-Broward reef", "urban-Broward reef")
Ofav_Letters$p.adj <- as.numeric(paste(Ofav_Letters$p.adj))
Ofav_Letters

# creates compact letter display of significant pairwise differences for figure
Ofav_Cld <- cldList(p.adj ~ comparison, data = Ofav_Letters, threshold = 0.05, remove.zero = FALSE, remove.space = FALSE)
Ofav_Cld=Ofav_Cld[order(Ofav_Cld$Group),] 
Ofav_Cld$RegionNew <- Ofav_Cld$Group
Ofav_Cld <- Ofav_Cld %>% mutate(symtype="Durusdinium")
Ofav_Cld

# percent stacked barplot
Ofav_Plot <- ggplot(Ofav_Perc, aes(fill=symtype, y=abundance, x=RegionNew)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Region",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Orbicella faveolata") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  geom_text(data = Ofav_Cld, aes(y=-0.05, label=Letter)) +
  geom_text(x = 1.85, y=1.035, label = "Region: F2,122 = 60.0, R2 = 0.492, p < 0.001") +
  # geom_text(x = 2, y=1.035, label = "Rain: F1,122 = 4.7, R2 = 0.019, p = 0.026") +
  theme(plot.title = element_text(hjust=0.5), legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())
Ofav_Plot 

ggsave("Ofav symbionts.pdf", plot= Ofav_Plot, width=4, height=4, units="in", dpi=300)

# Symbiont abundance by site
Ofav_Perc_Site <- dplyr::select(Ofav_Prop, 6, 11:14)
Ofav_Perc_Site <- reshape2::melt(Ofav_Perc_Site, id = "Site")
Ofav_Perc_Site$Site=factor(Ofav_Perc_Site$Site, levels=c("Peanut Island","T328","BC1","FTL4","South Canyon","Pillars","Graceland","Haulover Inlet","FIU Biscayne Bay","Arch Creek","Emerald Reef","Rainbow Reef","Bill Baggs","Seaquarium","Fisher Island","Coral City Camera","Miami Beach Marina","Star Island","MacArthur North","Belle Isle","FEC Slip Downtown"))
Ofav_Perc_Site <- dcast(Ofav_Perc_Site, Site~variable, mean)
Ofav_Perc_Site <- melt(Ofav_Perc_Site, id = "Site")
Ofav_Perc_Site <- dplyr::rename(Ofav_Perc_Site, symtype = variable, abundance = value)

# Adding region for facet_wrap
Ofav_Perc_Site %>%
  ungroup() %>%
  mutate(Region = case_when(
    Site == "Peanut Island"  ~ "Palm Beach urban",
    Site == "T328"  ~ "Broward reef",
    Site == "BC1"  ~ "Broward reef",
    Site == "FTL4"  ~ "Broward reef",
    Site == "South Canyon"  ~ "North Miami reef",
    Site == "Pillars"  ~ "North Miami reef",
    Site == "Graceland"  ~ "North Miami reef",
    Site == "Haulover Inlet"  ~ "North Miami urban",
    Site == "FIU Biscayne Bay"  ~ "North Miami urban",
    Site == "Arch Creek"  ~ "North Miami urban",
    Site == "Emerald Reef"  ~ "Miami reef",
    Site == "Rainbow Reef"  ~ "Miami reef",
    Site == "Bill Baggs"  ~ "Miami urban",
    Site == "Seaquarium"  ~ "Miami urban",
    Site == "Fisher Island"  ~ "Miami urban",
    Site == "Coral City Camera"  ~ "Miami urban",
    Site == "Miami Beach Marina"  ~ "Miami urban",
    Site == "Star Island"  ~ "Miami urban",
    Site == "MacArthur North"  ~ "Miami urban",
    Site == "Belle Isle"  ~ "Miami urban",
    Site == "FEC Slip Downtown"  ~ "Miami urban")) -> Ofav_Perc_Site
Ofav_Perc_Site$Region=factor(Ofav_Perc_Site$Region, levels=c("Palm Beach urban","Broward reef","North Miami reef","North Miami urban","Miami reef","Miami urban"))

# percent stacked barplot by site
Ofav_Plot_Site <- ggplot(Ofav_Perc_Site, aes(fill=symtype, y=abundance, x=Site)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Site",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Orbicella faveolata") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5)) +
  facet_grid(~Region, scales = "free", space='free') 
Ofav_Plot_Site 

ggsave("Ofav symbionts by site.pdf", plot= Ofav_Plot_Site, width=20, height=4, units="in", dpi=300)


#### CNAT ABUNDANCE ####

# transposing and reformatting dataframes to make abundance a single column
Cnat_Perc <- dplyr::select(Cnat_Prop, 21, 11:14)
Cnat_Perc <- reshape2::melt(Cnat_Perc, id = "RegionNew")
Cnat_Perc$RegionNew=factor(Cnat_Perc$RegionNew, levels=c("reef", "Palm Beach urban", "Miami Dade urban")) 
Cnat_Perc <- dcast(Cnat_Perc, RegionNew~variable, mean)
Cnat_Perc <- melt(Cnat_Perc, id = "RegionNew")
Cnat_Perc <- dplyr::rename(Cnat_Perc, symtype = variable, abundance = value)

# creating dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
Cnat_Letters <- data.frame(cbind(Cnat_Pair_Out$Comparison,Cnat_Pair_Out$'Pr(>F)'))
Cnat_Letters <- na.omit(Cnat_Letters)
Cnat_Letters <- dplyr::rename(Cnat_Letters, comparison = X1, p.adj = X2)
Cnat_Letters$comparison = c("reef-Miami Dade urban", "reef-Palm Beach urban", "Miami Dade urban-Palm Beach urban")
Cnat_Letters$p.adj <- as.numeric(paste(Cnat_Letters$p.adj))
Cnat_Letters

# creates compact letter display of significant pairwise differences for figure
Cnat_Cld <- cldList(p.adj ~ comparison, data = Cnat_Letters, threshold = 0.05, remove.zero = FALSE, remove.space = FALSE)
Cnat_Cld=Cnat_Cld[order(Cnat_Cld$Group),] 
Cnat_Cld$RegionNew <- Cnat_Cld$Group
Cnat_Cld <- Cnat_Cld %>% mutate(symtype="Durusdinium")
Cnat_Cld

# percent stacked barplot
Cnat_Plot <- ggplot(Cnat_Perc, aes(fill=symtype, y=abundance, x=RegionNew)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Region",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Colpophyllia natans") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  geom_text(data = Cnat_Cld, aes(y=-0.05, label=Letter)) +
  geom_text(x = 1.85, y=1.035, label = "Region: F2,137 = 31.2, R2 = 0.270, p < 0.001") +
  # geom_text(x = 2, y=1.035, label = "Rain: F1,137 = 36.0, R2 = 0.156, p < 0.001") +
  theme(plot.title = element_text(hjust=0.5), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank())
Cnat_Plot 

ggsave("Cnat symbionts.pdf", plot= Cnat_Plot, width=4, height=4, units="in", dpi=300)

# Symbiont abundance by site
Cnat_Perc_Site <- dplyr::select(Cnat_Prop, 6, 11:14)
Cnat_Perc_Site <- reshape2::melt(Cnat_Perc_Site, id = "Site")
Cnat_Perc_Site$Site=factor(Cnat_Perc_Site$Site, levels=c("Peanut Island","T328","BC1","FTL4","South Canyon","Pillars","Graceland","Haulover Inlet","FIU Biscayne Bay","Arch Creek","Emerald Reef","Rainbow Reef","Bill Baggs","Seaquarium","Fisher Island","Coral City Camera","Miami Beach Marina","Star Island","MacArthur North","Belle Isle","FEC Slip Downtown"))
Cnat_Perc_Site <- dcast(Cnat_Perc_Site, Site~variable, mean)
Cnat_Perc_Site <- melt(Cnat_Perc_Site, id = "Site")
Cnat_Perc_Site <- dplyr::rename(Cnat_Perc_Site, symtype = variable, abundance = value)

# Adding region for facet_wrap
Cnat_Perc_Site %>%
  ungroup() %>%
  mutate(Region = case_when(
    Site == "Peanut Island"  ~ "Palm Beach urban",
    Site == "T328"  ~ "Broward reef",
    Site == "BC1"  ~ "Broward reef",
    Site == "FTL4"  ~ "Broward reef",
    Site == "South Canyon"  ~ "North Miami reef",
    Site == "Pillars"  ~ "North Miami reef",
    Site == "Graceland"  ~ "North Miami reef",
    Site == "Haulover Inlet"  ~ "North Miami urban",
    Site == "FIU Biscayne Bay"  ~ "North Miami urban",
    Site == "Arch Creek"  ~ "North Miami urban",
    Site == "Emerald Reef"  ~ "Miami reef",
    Site == "Rainbow Reef"  ~ "Miami reef",
    Site == "Bill Baggs"  ~ "Miami urban",
    Site == "Seaquarium"  ~ "Miami urban",
    Site == "Fisher Island"  ~ "Miami urban",
    Site == "Coral City Camera"  ~ "Miami urban",
    Site == "Miami Beach Marina"  ~ "Miami urban",
    Site == "Star Island"  ~ "Miami urban",
    Site == "MacArthur North"  ~ "Miami urban",
    Site == "Belle Isle"  ~ "Miami urban",
    Site == "FEC Slip Downtown"  ~ "Miami urban")) -> Cnat_Perc_Site
Cnat_Perc_Site$Region=factor(Cnat_Perc_Site$Region, levels=c("Palm Beach urban","Broward reef","North Miami reef","North Miami urban","Miami reef","Miami urban"))

# percent stacked barplot by site
Cnat_Plot_Site <- ggplot(Cnat_Perc_Site, aes(fill=symtype, y=abundance, x=Site)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Site",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Colpophyllia natans") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5)) +
  facet_grid(~Region, scales = "free", space='free') 
Cnat_Plot_Site 

ggsave("Cnat symbionts by site.pdf", plot= Cnat_Plot_Site, width=20, height=4, units="in", dpi=300)


#### PSEU ABUNDANCE ####

# transposing and reformatting dataframes to make abundance a single column
Pseu_Perc <- dplyr::select(Pseu_Prop, 5, 11:14)
Pseu_Perc <- reshape2::melt(Pseu_Perc, id = "Region")
Pseu_Perc$Region=factor(Pseu_Perc$Region, levels=c("Palm Beach urban", "Broward reef", "North Miami reef", "North Miami urban", "Miami reef", "Miami urban")) 
Pseu_Perc <- dcast(Pseu_Perc, Region~variable, mean)
Pseu_Perc <- melt(Pseu_Perc, id = "Region")
Pseu_Perc <- dplyr::rename(Pseu_Perc, symtype = variable, abundance = value)

# creating dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
Pseu_Letters <- data.frame(cbind(Pseu_Pair_Out$Comparison,Pseu_Pair_Out$'Pr(>F)'))
Pseu_Letters <- na.omit(Pseu_Letters)
Pseu_Letters <- dplyr::rename(Pseu_Letters, comparison = X1, p.adj = X2)
Pseu_Letters$comparison = c("Miami reef-Miami urban", "Miami reef-Broward reef", "Miami reef-Palm Beach urban", "Miami reef-North Miami reef", "Miami reef-North Miami urban", "Miami urban-Broward reef", "Miami urban-Palm Beach urban", "Miami urban-North Miami reef", "Miami urban-North Miami urban", "Broward reef-Palm Beach urban", "Broward reef-North Miami reef", "Broward reef-North Miami urban", "Palm Beach urban-North Miami reef", "Palm Beach urban-North Miami urban", "North Miami reef-North Miami urban")
Pseu_Letters$p.adj <- as.numeric(paste(Pseu_Letters$p.adj))
Pseu_Letters

# creates compact letter display of significant pairwise differences for figure
Pseu_Cld <- cldList(p.adj ~ comparison, data = Pseu_Letters, threshold = 0.05, remove.zero = FALSE, remove.space = FALSE)
Pseu_Cld=Pseu_Cld[order(Pseu_Cld$Group),] 
Pseu_Cld$Region <- Pseu_Cld$Group
Pseu_Cld <- Pseu_Cld %>% mutate(symtype="Durusdinium")
Pseu_Cld

# percent stacked barplot
Pseu_Plot <- ggplot(Pseu_Perc, aes(fill=symtype, y=abundance, x=Region)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Region",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Pseudodiploria spp.") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  geom_text(data = Pseu_Cld, aes(y=-0.05, label=Letter)) +
  geom_text(x = 1.85, y=1.035, label = "Region: F5,322 = 10.2, R2 = 0.140, p < 0.001") +
  # geom_text(x = 2, y=1.035, label = "Depth: F1,322 = 4.3, R2 = 0.012, p = 0.041") +
  theme(plot.title = element_text(hjust=0.5), legend.position = "none", axis.title.x = element_blank())
Pseu_Plot 

ggsave("Pseu symbionts.pdf", plot= Pseu_Plot, width=8, height=4, units="in", dpi=300)

# Symbiont abundance by site
Pseu_Perc_Site <- dplyr::select(Pseu_Prop, 6, 11:14)
Pseu_Perc_Site <- reshape2::melt(Pseu_Perc_Site, id = "Site")
Pseu_Perc_Site$Site=factor(Pseu_Perc_Site$Site, levels=c("Peanut Island","T328","BC1","FTL4","South Canyon","Pillars","Graceland","Haulover Inlet","FIU Biscayne Bay","Arch Creek","Emerald Reef","Rainbow Reef","Bill Baggs","Seaquarium","Fisher Island","Coral City Camera","Miami Beach Marina","Star Island","MacArthur North","Belle Isle","FEC Slip Downtown"))
Pseu_Perc_Site <- dcast(Pseu_Perc_Site, Site~variable, mean)
Pseu_Perc_Site <- melt(Pseu_Perc_Site, id = "Site")
Pseu_Perc_Site <- dplyr::rename(Pseu_Perc_Site, symtype = variable, abundance = value)

# Adding region for facet_wrap
Pseu_Perc_Site %>%
  ungroup() %>%
  mutate(Region = case_when(
    Site == "Peanut Island"  ~ "Palm Beach urban",
    Site == "T328"  ~ "Broward reef",
    Site == "BC1"  ~ "Broward reef",
    Site == "FTL4"  ~ "Broward reef",
    Site == "South Canyon"  ~ "North Miami reef",
    Site == "Pillars"  ~ "North Miami reef",
    Site == "Graceland"  ~ "North Miami reef",
    Site == "Haulover Inlet"  ~ "North Miami urban",
    Site == "FIU Biscayne Bay"  ~ "North Miami urban",
    Site == "Arch Creek"  ~ "North Miami urban",
    Site == "Emerald Reef"  ~ "Miami reef",
    Site == "Rainbow Reef"  ~ "Miami reef",
    Site == "Bill Baggs"  ~ "Miami urban",
    Site == "Seaquarium"  ~ "Miami urban",
    Site == "Fisher Island"  ~ "Miami urban",
    Site == "Coral City Camera"  ~ "Miami urban",
    Site == "Miami Beach Marina"  ~ "Miami urban",
    Site == "Star Island"  ~ "Miami urban",
    Site == "MacArthur North"  ~ "Miami urban",
    Site == "Belle Isle"  ~ "Miami urban",
    Site == "FEC Slip Downtown"  ~ "Miami urban")) -> Pseu_Perc_Site
Pseu_Perc_Site$Region=factor(Pseu_Perc_Site$Region, levels=c("Palm Beach urban","Broward reef","North Miami reef","North Miami urban","Miami reef","Miami urban"))

# percent stacked barplot by site
Pseu_Plot_Site <- ggplot(Pseu_Perc_Site, aes(fill=symtype, y=abundance, x=Site)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Site",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Pseudodiploria spp.") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5)) +
  facet_grid(~Region, scales = "free", space='free') 
Pseu_Plot_Site 

ggsave("Pseu symbionts by site.pdf", plot= Pseu_Plot_Site, width=27, height=4, units="in", dpi=300)


#### DLAB ABUNDANCE ####

# transposing and reformatting dataframes to make abundance a single column
Dlab_Perc <- dplyr::select(Dlab_Prop, 21, 11:14)
Dlab_Perc <- reshape2::melt(Dlab_Perc, id = "RegionNew")
Dlab_Perc$RegionNew=factor(Dlab_Perc$RegionNew, levels=c("reef","Palm Beach urban","Miami urban")) 
Dlab_Perc <- dcast(Dlab_Perc, RegionNew~variable, mean)
Dlab_Perc <- melt(Dlab_Perc, id = "RegionNew")
Dlab_Perc <- dplyr::rename(Dlab_Perc, symtype = variable, abundance = value)

# creating dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
Dlab_Letters <- data.frame(cbind(Dlab_Pair_Out$Comparison,Dlab_Pair_Out$'Pr(>F)'))
Dlab_Letters <- na.omit(Dlab_Letters)
Dlab_Letters <- dplyr::rename(Dlab_Letters, comparison = X1, p.adj = X2)
Dlab_Letters$comparison = c("reef-Miami urban", "reef-Palm Beach urban", "Miami urban-Palm Beach urban")
Dlab_Letters$p.adj <- as.numeric(paste(Dlab_Letters$p.adj))
Dlab_Letters

# creates compact letter display of significant pairwise differences for figure
Dlab_Cld <- cldList(p.adj ~ comparison, data = Dlab_Letters, threshold = 0.05, remove.zero = FALSE, remove.space = FALSE)
Dlab_Cld=Dlab_Cld[order(Dlab_Cld$Group),] 
Dlab_Cld$RegionNew <- Dlab_Cld$Group
Dlab_Cld <- Dlab_Cld %>% mutate(symtype="Durusdinium")
Dlab_Cld

# percent stacked barplot
Dlab_Plot <- ggplot(Dlab_Perc, aes(fill=symtype, y=abundance, x=RegionNew)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Region",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Diploria labyrinthiformis") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  geom_text(data = Dlab_Cld, aes(y=-0.05, label=Letter)) +
  geom_text(x = 1.85, y=1.035, label = "Region: F2,82 = 12.1, R2 = 0.229, p < 0.001") +
  # geom_text(x = 2, y=1.035, label = "Rain: F1,82 = 3.0, R2 = 0.029, p = 0.067") +
  theme(plot.title = element_text(hjust=0.5), legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank()) 
Dlab_Plot 

ggsave("Dlab symbionts.pdf", plot= Dlab_Plot, width=4, height=4, units="in", dpi=300)

# Symbiont abundance by site
Dlab_Perc_Site <- dplyr::select(Dlab_Prop, 6, 11:14)
Dlab_Perc_Site <- reshape2::melt(Dlab_Perc_Site, id = "Site")
Dlab_Perc_Site$Site=factor(Dlab_Perc_Site$Site, levels=c("Peanut Island","T328","BC1","FTL4","South Canyon","Pillars","Graceland","Haulover Inlet","FIU Biscayne Bay","Arch Creek","Emerald Reef","Rainbow Reef","Bill Baggs","Seaquarium","Fisher Island","Coral City Camera","Miami Beach Marina","Star Island","MacArthur North","Belle Isle","FEC Slip Downtown"))
Dlab_Perc_Site <- dcast(Dlab_Perc_Site, Site~variable, mean)
Dlab_Perc_Site <- melt(Dlab_Perc_Site, id = "Site")
Dlab_Perc_Site <- dplyr::rename(Dlab_Perc_Site, symtype = variable, abundance = value)

# Adding region for facet_wrap
Dlab_Perc_Site %>%
  ungroup() %>%
  mutate(Region = case_when(
    Site == "Peanut Island"  ~ "Palm Beach urban",
    Site == "T328"  ~ "Broward reef",
    Site == "BC1"  ~ "Broward reef",
    Site == "FTL4"  ~ "Broward reef",
    Site == "South Canyon"  ~ "North Miami reef",
    Site == "Pillars"  ~ "North Miami reef",
    Site == "Graceland"  ~ "North Miami reef",
    Site == "Haulover Inlet"  ~ "North Miami urban",
    Site == "FIU Biscayne Bay"  ~ "North Miami urban",
    Site == "Arch Creek"  ~ "North Miami urban",
    Site == "Emerald Reef"  ~ "Miami reef",
    Site == "Rainbow Reef"  ~ "Miami reef",
    Site == "Bill Baggs"  ~ "Miami urban",
    Site == "Seaquarium"  ~ "Miami urban",
    Site == "Fisher Island"  ~ "Miami urban",
    Site == "Coral City Camera"  ~ "Miami urban",
    Site == "Miami Beach Marina"  ~ "Miami urban",
    Site == "Star Island"  ~ "Miami urban",
    Site == "MacArthur North"  ~ "Miami urban",
    Site == "Belle Isle"  ~ "Miami urban",
    Site == "FEC Slip Downtown"  ~ "Miami urban")) -> Dlab_Perc_Site
Dlab_Perc_Site$Region=factor(Dlab_Perc_Site$Region, levels=c("Palm Beach urban","Broward reef","North Miami reef","North Miami urban","Miami reef","Miami urban"))

# percent stacked barplot by site
Dlab_Plot_Site <- ggplot(Dlab_Perc_Site, aes(fill=symtype, y=abundance, x=Site)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Site",
       y = "Relative Abundance",
       fill = 'Symbiodiniaceae Genus',
       title = "Diploria labyrinthiformis") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(6)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5)) +
  facet_grid(~Region, scales = "free", space='free') 
Dlab_Plot_Site 

ggsave("Dlab symbionts by site.pdf", plot= Dlab_Plot_Site, width=16, height=4, units="in", dpi=300)


#### ABUNDANCE MULTIPLOTS ####

Perc_Multiplot <- grid.arrange(Pseu_Plot, Cnat_Plot, Dlab_Plot, Mcav_Plot, Ofav_Plot, legend, ncol=3, nrow=2, widths=c(6,3,3), heights=c(3.75,4))

ggsave("Multiplot symbionts.pdf", plot= Perc_Multiplot, width=16, height=8, units="in", dpi=300)
