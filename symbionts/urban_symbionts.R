#### PACKAGES ####

library(tidyverse) # install.packages('tidyverse')
# library(devtools)  # install.packages("devtools")
# devtools::install_github("jrcunning/steponeR")
library(steponeR)
library(reshape2)
library(ggplot2)


#### DATA IMPORT 1 ####

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


#### DATA CLEANING 1 ####

# Identifies all samples where only one symbiont technical replicate amplifies
One_A<- Urban[which(Urban$A.reps==1),]
One_B<- Urban[which(Urban$B.reps==1),]
One_C<- Urban[which(Urban$C.reps==1),]
One_D<- Urban[which(Urban$D.reps==1),]

# If only one technical replicate amplifies, set all corresponding symbiont ratios to NA
Urban$B.A[which(Urban$A.reps==1 | Urban$B.reps==1)] <- NA
Urban$C.A[which(Urban$A.reps==1 | Urban$C.reps==1)] <- NA
Urban$D.A[which(Urban$A.reps==1 | Urban$D.reps==1)] <- NA
Urban$A.B[which(Urban$B.reps==1 | Urban$A.reps==1)] <- NA
Urban$C.B[which(Urban$B.reps==1 | Urban$C.reps==1)] <- NA
Urban$D.B[which(Urban$B.reps==1 | Urban$D.reps==1)] <- NA
Urban$A.C[which(Urban$C.reps==1 | Urban$A.reps==1)] <- NA
Urban$B.C[which(Urban$C.reps==1 | Urban$B.reps==1)] <- NA
Urban$D.C[which(Urban$C.reps==1 | Urban$D.reps==1)] <- NA
Urban$A.D[which(Urban$D.reps==1 | Urban$A.reps==1)] <- NA
Urban$B.D[which(Urban$D.reps==1 | Urban$B.reps==1)] <- NA
Urban$C.D[which(Urban$D.reps==1 | Urban$C.reps==1)] <- NA

# Now set all NAs to 0
Urban$B.A[is.na(Urban$B.A)] <- 0
Urban$C.A[is.na(Urban$C.A)] <- 0
Urban$D.A[is.na(Urban$D.A)] <- 0
Urban$A.B[is.na(Urban$A.B)] <- 0
Urban$C.B[is.na(Urban$C.B)] <- 0
Urban$D.B[is.na(Urban$D.B)] <- 0
Urban$A.C[is.na(Urban$A.C)] <- 0
Urban$B.C[is.na(Urban$B.C)] <- 0
Urban$D.C[is.na(Urban$D.C)] <- 0
Urban$A.D[is.na(Urban$A.D)] <- 0
Urban$B.D[is.na(Urban$B.D)] <- 0
Urban$C.D[is.na(Urban$C.D)] <- 0


#### METADATA ####

Metadata<-read.csv("metadata.csv")
Metadata$Species <- as.factor(Metadata$Species)
Metadata$Region <- as.factor(Metadata$Region)
Metadata$Site <- as.factor(Metadata$Site)
Metadata$CollectionDate <- as.factor(Metadata$CollectionDate)
Metadata$Season <- as.factor(Metadata$Season)
Metadata$Rain <- as.factor(Metadata$Rain)

# Create unique sample ID+FileName to merge with metadata  
Urban$Sample.Plate<-paste(Urban$Sample.Name, Urban$File.Name, sep = "_" )

# Joining metadata with qPCR data
Urban_Join<-left_join(Urban, Metadata, by="Sample.Plate")


#### FILTERING BY SPECIES ####

# Filtering by species, then determining Ct means for each symbiont by site
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

# The lowest mean Ct per site indicates the dominant symbiont genus
Urban_Join %>%
  filter(Species == "Mcav") %>%
  group_by(Site) %>%
  summarize(across(A.CT.mean:D.CT.mean, mean, na.rm=TRUE)) -> Mcav_dom
Urban_Join %>%
  filter(Species == "Ofav") %>%
  group_by(Site) %>%
  summarize(across(A.CT.mean:D.CT.mean, mean, na.rm=TRUE)) -> Ofav_dom
Urban_Join %>%
  filter(Species == "Cnat") %>%
  group_by(Site) %>%
  summarize(across(A.CT.mean:D.CT.mean, mean, na.rm=TRUE)) -> Cnat_dom
Urban_Join %>%
  filter(Species == "Pseu") %>%
  group_by(Site) %>%
  summarize(across(A.CT.mean:D.CT.mean, mean, na.rm=TRUE)) -> Pseu_dom
Urban_Join %>%
  filter(Species == "Dlab") %>%
  group_by(Site) %>%
  summarize(across(A.CT.mean:D.CT.mean, mean, na.rm=TRUE)) -> Dlab_dom

# Further filtering species by site based on dominant symbionts
Mcav %>%
  filter(Region == "Miami offshore") -> Mcav_offshore
Mcav %>%
  filter(Region == "Miami urban") -> Mcav_urban
Ofav %>%
  filter(Region == "Miami offshore") -> Ofav_offshore
Ofav %>%
  filter(Region == "Miami urban") -> Ofav_urban
# Cnat %>% # Not needed since Cnat is dominated by D across all sites
#   filter(Region == "Miami offshore") -> Cnat_offshore
# Cnat %>%
#     filter(Region == "Miami urban") -> Cnat_urban
Pseu %>%
  filter(Region == "Miami offshore") -> Pseu_offshore
Pseu %>%
  filter(Region == "Miami urban") -> Pseu_urban
# Dlab %>% # Not needed since Dlab is a small dataset
#   filter(Region == "Miami offshore") -> Dlab_offshore
# Dlab %>%
#   filter(Region == "Miami urban") -> Dlab_urban


#### DATA CLEANING 2 ####

# Makes a list of samples with only one technical replicate of the dominant symbiont type
ReRun_Cnat <- Cnat[which(Cnat$D.reps==1), ] # 1 sample violates, but not necessary to rerun
ReRun_Ofav_urban <- Ofav_urban[which(Ofav_urban$D.reps==1), ] # no samples violate
ReRun_Ofav_offshore <- Ofav_offshore[which(Ofav_offshore$A.reps==1 | Ofav_offshore$B.reps==1 | Ofav_offshore$D.reps==1), ] # 10 samples violate, but not all necessary to rerun
ReRun_Mcav_offshore <- Mcav_offshore[which(Mcav_offshore$C.reps==1), ] # no samples violate
ReRun_Mcav_urban <- Mcav_urban[which(Mcav_urban$D.reps==1 | Mcav_urban$A.reps==1), ] # 2 samples violate, but not necessary to rerun
ReRun_Pseu_offshore <- Pseu_offshore[which(Pseu_offshore$B.reps==1), ] # rerun Plate 20 Sample 10
ReRun_Pseu_urban <- Pseu_urban[which(Pseu_urban$D.reps==1), ] # 2 samples violate, but not necessary to rerun
ReRun_Dlab <- Dlab[which(Dlab$A.reps==1 | Dlab$C.reps==1 |Dlab$D.reps==1), ] # 7 sample violates, only 1 necessary to rerun

# Makes a list of samples where technical replicates of dominant symbiont types had standard deviation >1.5
StDe1.5_Cnat <- Cnat[which(Cnat$D.CT.sd>1.5), ] # no samples violate
StDe1.5_Ofav_urban <- Ofav_urban[which(Ofav_urban$D.CT.sd>1.5), ] # no samples violate
StDe1.5_Ofav_offshore <- Ofav_offshore[which(Ofav_offshore$A.CT.sd>1.5 | Ofav_offshore$B.CT.sd>1.5 | Ofav_offshore$B.CT.sd>1.5), ] # 13 samples violate, but not all necessary to rerun
StDe1.5_Mcav_offshore <- Mcav_offshore[which(Mcav_offshore$C.CT.sd>1.5), ] # no samples violate
StDe1.5_Mcav_urban <- Mcav_urban[which(Mcav_urban$D.CT.sd>1.5 | Mcav_urban$A.CT.sd>1.5), ] # 13 samples violate, but not all necessary to rerun
StDe1.5_Pseu_offshore <- Pseu_offshore[which(Pseu_offshore$B.CT.sd>1.5), ] # no samples violate
StDe1.5_Pseu_urban <- Pseu_urban[which(Pseu_urban$D.CT.sd>1.5), ] # no samples violate
StDe1.5_Dlab <- Dlab[which(Dlab$D.CT.sd>1.5 | Dlab$A.CT.sd>1.5 | Dlab$C.CT.sd>1.5), ] # 5 samples violate, not necessary to rerun

# Makes a list of samples where the Ct mean of dominant symbiont types was >25 (late amplification)
Late_Cnat<-Cnat[which(Cnat$D.CT.mean>25), ] # 3 samples violate, no need to rerun
Late_Ofav_urban<-Ofav_urban[which(Ofav_urban$D.CT.mean>25), ] # 1 sample violates, no need to rerun
Late_Ofav_offshore<-Ofav_offshore[which(Ofav_offshore$A.CT.mean>25 & Ofav_offshore$B.CT.mean>25 & Ofav_offshore$D.CT.mean>25), ] # 4 samples violate, no need to rerun
Late_Mcav_urban<-Mcav_urban[which(Mcav_urban$D.CT.mean>25 & Mcav_urban$A.CT.mean>25), ] # 1 sample violates, no need to rerun
Late_Mcav_offshore<-Mcav_offshore[which(Mcav_offshore$C.CT.mean>25), ] # 0 samples violate
Late_Pseu_urban<-Pseu_urban[which(Pseu_urban$D.CT.mean>25), ] # 7 sample violates, no need to rerun
Late_Pseu_offshore<-Pseu_offshore[which(Pseu_offshore$B.CT.mean>25), ] # 13 samples violate, no need to rerun
Late_Dlab<-Dlab[which(Dlab$D.CT.mean>25 & Dlab$A.CT.mean>25 & Dlab$C.CT.mean>25), ] # 2 samples violate, not necessary to rerun

# Combines all lists above by species and finds distinct samples
ToReRun_Cnat<-rbind(ReRun_Cnat, StDe1.5_Cnat, Late_Cnat)
ToReRun_Cnat<-ToReRun_Cnat %>% distinct()
ToReRun_Ofav<-rbind(ReRun_Ofav_urban, ReRun_Ofav_offshore, StDe1.5_Ofav_urban, StDe1.5_Ofav_offshore, Late_Ofav_urban, Late_Ofav_offshore)
ToReRun_Ofav<-ToReRun_Ofav %>% distinct()
ToReRun_Mcav<-rbind(ReRun_Mcav_urban, ReRun_Mcav_offshore, StDe1.5_Mcav_urban, StDe1.5_Mcav_offshore, Late_Mcav_urban, Late_Mcav_offshore)
ToReRun_Mcav<-ToReRun_Mcav %>% distinct()
ToReRun_Pseu<-rbind(ReRun_Pseu_urban, ReRun_Pseu_offshore, StDe1.5_Pseu_urban, StDe1.5_Pseu_offshore, Late_Pseu_urban, Late_Pseu_offshore)
ToReRun_Pseu<-ToReRun_Pseu %>% distinct()
ToReRun_Dlab<-rbind(ReRun_Dlab, StDe1.5_Dlab, Late_Dlab)
ToReRun_Dlab<-ToReRun_Dlab %>% distinct()


#### REMOVE DUPLICATES ####

# Finds and removes the first instance of duplicated sample IDs (samples that were rerun)
Urban_Join %>%
  group_by(ID) %>%
  filter(duplicated(ID)|n()==1) -> Urban_DeDup

# Does the same with metadata
Metadata %>%
  group_by(ID) %>%
  filter(duplicated(ID)|n()==1) -> Metadata_DeDup


#### REFILTERING BY SPECIES ####

# Refiltering by species using the deduplicated dataset
Urban_DeDup %>%
  filter(Species == "Mcav") -> Mcav_DeDup
Urban_DeDup %>%
  filter(Species == "Ofav") -> Ofav_DeDup
Urban_DeDup %>%
  filter(Species == "Cnat") -> Cnat_DeDup
Urban_DeDup %>%
  filter(Species == "Pseu") -> Pseu_DeDup
Urban_DeDup %>%
  filter(Species == "Dlab") -> Dlab_DeDup

# Further filtering species by site based on dominant symbionts
Mcav_DeDup %>%
  filter(Region == "Miami offshore") -> Mcav_DeDup_offshore
Mcav_DeDup %>%
  filter(Region == "Miami urban") -> Mcav_urban
Ofav_DeDup %>%
  filter(Region == "Miami offshore") -> Ofav_DeDup_offshore
Ofav_DeDup %>%
  filter(Region == "Miami urban") -> Ofav_DeDup_urban
# Cnat_DeDup %>% # Not needed since Cnat is dominated by D across all sites
#   filter(Region == "Miami offshore") -> Cnat_DeDup_offshore
# Cnat_DeDup %>%
#     filter(Region == "Miami urban") -> Cnat_DeDup_urban
Pseu_DeDup %>%
  filter(Region == "Miami offshore") -> Pseu_DeDup_offshore
Pseu_DeDup %>%
  filter(Region == "Miami urban") -> Pseu_DeDup_urban
# Dlab_DeDup %>% # Not needed since Dlab is a small dataset
#   filter(Region == "Miami offshore") -> Dlab_DeDup_offshore
# Dlab_DeDup %>%
#   filter(Region == "Miami urban") -> Dlab_DeDup_urban


#### SYMBIONT PROPORTIONS ####

# Determines the column with the lowest Ct value (most abundant symbiont) per sample
# the 'replace' string replaces NAs with Inf for the sake of the ranking, but does not change the raw data
Urban_DeDup$DomSym <- max.col(-replace(Urban_DeDup[3:6], is.na(Urban_DeDup[3:6]), Inf))
# Now replacing the column numbers with symbiont type in a new column
Urban_DeDup %>%
  mutate(DomSymType = recode(DomSym, '1' = 'A', '2' = 'B', '3' = 'C', '4' = 'D')) -> Urban_DeDup



Urban_DeDup$DomSym <- max.col(is.na(-Urban_DeDup[3:6]))

Urban_DeDup$TotalnoD <- (Urban_DeDup$A.D + Urban_DeDup$B.D+ Urban_DeDup$C.D)
Urban_DeDup$TotalD <- (1-Urban_DeDup$TotalnoD)

# Log 10
    Urban$logA.SH <- log10(Urban$A.D)
    Urban$logB.SH <- log10(Urban$B.D)
    Urban$logC.SH <- log10(Urban$C.D)
    Urban$logSH<-log10(Urban$TotalSH)

    Urban$logA.SH[which(Urban$A.D==0)] <- NA
    Urban$logB.SH[which(Urban$B.D==0)] <- NA
    Urban$logC.SH[which(Urban$C.D==0)] <- NA
    

    # Clade Proportion
    # D Proportion
    Urban_DeDup$D.Prp <- (Urban_DeDup$TotalD/Urban_DeDup$TotalnoD)
    # C Proportion
    #Ofav$C.Prp<-(Ofav$C.SH/Ofav$tot.SH)
    
    hist(Urban_DeDup$D.Prp)
    Dproportion<-Ofav[(Ofav$D.Prp<1 & Ofav$D.Prp>0.01),]
    hist(Dproportion$D.Prp)