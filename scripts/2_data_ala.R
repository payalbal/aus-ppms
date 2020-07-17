## Download and clean ALA data

## Set working environment
rm(list = ls())
gc()
system("ps")
system("pkill -f R")

setwd("./aus-ppms/")
data_path <- "./data" #"/Volumes/discovery_data/aus-ppms_data/"
ala_folder <- "AVES"
ala_path <- file.path(data_path, "ala_data", ala_folder)

## Get download and cleaning scripts from github repo
## Ref: https://github.com/cwarecsiro/gdmEngine/tree/master/gdmEngine
system(paste0("curl https://raw.githubusercontent.com/cwarecsiro/gdmEngine/master/gdmEngine/R/download_ala.R -o ", "./scripts/download_ala.R"))
system(paste0("curl https://raw.githubusercontent.com/cwarecsiro/gdmEngine/master/gdmEngine/R/download_taxalist.R -o ", "./scripts/download_ala_bytaxa.R"))
## Change line in download_ala_bytaxa.R
x <- readLines("./scripts/download_ala_bytaxa.R")
x[179] <- "        print(paste0('Searching the ALA for records of ', spp))"
cat(x, file="./scripts/download_ala_bytaxa.R", sep="\n")
system(paste0("curl https://raw.githubusercontent.com/cwarecsiro/gdmEngine/master/gdmEngine/R/merge_downloads.R -o ", "./scripts/merge_ala_downloads.R"))
system(paste0("curl https://raw.githubusercontent.com/cwarecsiro/gdmEngine/master/gdmEngine/R/filter_ALA_data.R -o ", "./scripts/filter_ala_data.R"))

## Download species list
## 1. Go to https://www.ala.org.au/
## 2. Click on total occurrences
## 3. Click download
## 4 Select 'Species checklist' and specify reason for download
## 5. Download list 

## List names of native species 
data_sp <- read.csv(file.path(data_path, "all_sp_records-2020-03-19.csv"))
data_sp <- data_sp[data_sp$Invasive %in% "",]
data_sp <- data_sp[,c("Species.Name", "Kingdom", 
                      "Phylum", "Class", "Order", 
                      "Family", "Genus", "Conservation")]

## Remove duplicates
if(!(dim(data_sp)[1] == length(unique(data_sp$Species.Name)))){
  data_sp <- data_sp[!duplicated(data_sp$Species.Name),]
}

## [Optional] Subset data_sp
sort(unique(data_sp$Class))
data_sp <- data_sp[data_sp$Class == "AVES",]

## Check that all names are unique
dim(data_sp)[1] == length(unique(data_sp$Species.Name))

## Download ALA data using species list
library(httr)
library(assertthat)
# library(parallel)
source("./scripts/download_ala.R")
source("./scripts/download_ala_bytaxa.R")

splist <- unique(data_sp$Species.Name)

download_taxalist(specieslist = splist, #trial: splist[111:112], 
                  dst = ala_path, 
                  parallel = FALSE, 
                  background = FALSE)
## parallel = TRUE gives error because package ‘gdmEngine’ is not available (for R version 3.6.2)

## Merge downloaded data
## ***NOTE: For gsdms, think of merging subsets of data (e.g. by taxa/class/# species) to run filter function on and merge final filtered data.
library(data.table)
library(dplyr)
source("./scripts/merge_ala_downloads.R")
ala.data <- merge_downloads(src = file.path(ala_path, "raw_files/"), output.folder = ala_path,
                            output.name = paste0("merged_data_", Sys.Date()),
                            keep_unzip = FALSE,
                            parallel = FALSE, 
                            verbose = TRUE)

## Filter merged data
library(sp)
library(raster)
source("./scripts/filter_ala_data.R")
ala.data <- read.csv(file.path(data_path, "merged_data_2020-04-01.csv"))
aus.mask <- readRDS(file.path(data_path,"mask_aus.rds")) # after covariate data prep

filtered.data <- filter_ALA_data(ALA.download.data = ala.data$data,             
                                 output.folder = ala_path,       
                                 output.name = "filtered_data_",  
                                 domain.mask = aus.mask,                   
                                 earliest.year = 1950,
                                 spatial.uncertainty.m = 2000,
                                 select.fields = NULL,
                                 verbose=TRUE)

## Subset filtered data to species with >=20 records
dat <- as.data.table(filtered.data)
dat.counts <- dat[,.N, by = scientificName]
sp.names <- dat.counts[N >= 20]$scientificName
dat <- dat[scientificName %in% sp.names]

saveRDS(dat, file.path(ala_path, "ala_data.rds"))
write.csv(dat, file.path(ala_path, "ala_data.csv"), row.names=FALSE) ## larger file.  



