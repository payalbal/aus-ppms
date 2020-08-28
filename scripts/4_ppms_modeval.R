## Set up work environment ####
rm(list = ls())
gc()
Sys.setenv(TZ='Australia/Melbourne')

# system("ps")
# system("pkill -f R")

x <- c('sp', 'raster', 'data.table', 'spatstat.data', 'nlme', 'rpart', 'spatstat', 'ppmlasso', 
       'parallel', 'doParallel', 'doMC', 'rgdal', 'plyr', 'dismo', 'future.apply', 'future')
lapply(x, require, character.only = TRUE)
source("/tempdata/workdir/aus-ppms/scripts/0_functions.R")

## *** NOTE: all files within to be recreated after BB's landuse maps (currently using regSSP_birds files)
ppm_path <- "/tempdata/workdir/aus-ppms/data" # "/Volumes/discovery_data/aus-ppms_data"
rdata_path <- file.path(ppm_path, "RData")
dir.create("/tempdata/workdir/aus-ppms/output")
output_path <- paste0("/tempdata/workdir/aus-ppms/output/output_", gsub("-", "", format(Sys.time(), "%F_%H%M")))
datestamp <- gsub("-", "", format(Sys.time(), "%F")) ## used later to extract results
dir.create(output_path)


## Global parameters ####
bkpts <- 10000 # number of background points
n.fits <- 20


## Master log file ####
job_start <- Sys.time()
masterlog <- paste0(output_path, "/ppm_run_bkp_",bkpts, "_fits_", n.fits, ".txt")
writeLines(c(""), masterlog)
cat(paste0(">> Job start = ", job_start, "\n"), file = masterlog, append = TRUE)
cat(paste0("---------------------------------------------"), file = masterlog, append = TRUE)


## I. Data preperation ####
## Load biodiversity data ####
data_start <- Sys.time()
cat(paste0("\n\nData prep start = ", data_start, "\n"), file = masterlog, append = TRUE)

occdat <- readRDS(file.path(ppm_path, "ala_data/ala_data.rds"))
occdat <- occdat[,.(decimalLongitude, decimalLatitude, scientificName)]
names(occdat) <- c("X", "Y", "species")


## Test for correlations in covariates ####
covariates <- readRDS(file.path(rdata_path, "covariates_aus.rds"))
all_covs <- names(covariates)
covs_exclude <- c()
# covs_exclude <- all_covs[grepl("di", all_covs)]
# covs_bias <- c("dibu", "diro", "pa", "popd", "roughness")
if(!is.null(covs_exclude)){
  covariates <- covariates[[-which(names(covariates)%in%covs_exclude)]]
}else{
  covariates <- covariates
}
cov_df <- getValues(covariates)
cov_df <- na.omit(cov_df)
preds <- colnames(correlations(cov_df))
saveRDS(preds, file = file.path(rdata_path, "preds_aus.rds"))
rm(cov_df, covariates)
# preds <- readRDS(file.path(rdata_path, "preds_aus.rds"))


## Load covariate data ####
## *** NOTE: To update based on scearios used for BB's lu maps
## Write covariates into temp folder structure
temp <-  file.path(ppm_path, "temp")
ssps <- paste0("ssp", 1:3)
rcps <- c("45", "60", "85")
quartiles <- "q2" #c("q2", "q1", "q3")
scens <- sort(apply(expand.grid(quartiles, ssps), 1, paste0, collapse="_"))
scens_rcps <- sort(apply(expand.grid(quartiles, rcps), 1, paste0, collapse="_"))

type <- c("covariates", "predlu", "bio")
dir.create(temp)
for(k in 1:length(type)){
  if(type[k]%in%type[c(2,3)]){
    for(i in 1:length(scens)){
      if(type[k] == "predlu"){
        pattern <- paste0(type[k], "_", scens[i])
      }
      if(type[k] == "bio"){
        pattern <- paste0(type[k], scens_rcps[i])
      }
      covs <- readRDS(list.files(rdata_path, pattern = pattern, full.names = T))
      covs <- covs[[which(names(covs)%in%c(preds, "landuse"))]]
      writeToDisk(covs, folder = file.path(temp, pattern))
    }
  }else{
    covs <- readRDS(list.files(rdata_path,
                               pattern = "covariates",
                               full.names = T))
    covs <- covs[[which(names(covs)%in%c(preds, "landuse"))]]
    writeToDisk(covs, folder = file.path(temp, "models"))
  }
}
rm(covs, type, pattern, i, k)

## Covariate stack for model building
covariates_model <- stack(list.files(file.path(temp, "models"), full.names = T))
static_var <- covariates_model[[-which(grepl(paste0(c("landuse", "bio"),
                                                    collapse = "|"), names(covariates_model)))]]
dynamic_var <- covariates_model[[which(grepl("bio", names(covariates_model)))]]
landuse_var <- covariates_model[[which(grepl("landuse", names(covariates_model)))]]
covs_mod <- stack(static_var, dynamic_var, landuse_var)
names(covs_mod)

## Covariate stack for predictions (as per scenario)
# ...for(j in 1:length(scens)){...
j=3 
dynamic_var <- stack(list.files(file.path(temp,
                                          paste0("bio", scens_rcps[j])), full.names = TRUE))
landuse_var <- stack(stack(), list.files(file.path(temp,
                                                   paste0("predlu_", scens[j])), full.names = TRUE))
covs_pred <- stack(static_var, dynamic_var, landuse_var)

## Check that names are the same
all(names(covs_pred) == names(covs_mod))

rm(covariates_model, static_var, dynamic_var, landuse_var)

## Estimate samplig bias ####
## Here, using all records from ALA. Other optiosn can be explored later
## Author: Skiptin Wooley
## Create empty raster
reg_mask <- readRDS(file.path(rdata_path, "mask_aus.rds"))
r <- reg_mask
r[] <- 0

## Count the number of records per cell
tab <- table(cellFromXY(r, occdat[,1:2]))
r[as.numeric(names(tab))] <- tab
r <- mask(r,reg_mask)

## Make zeros a very small number otherwise issues with log(0).
r[r[]==0] <- 1e-6
arear <- raster::area(r)

## Calculate the probability that a cell has been sampled while accounting for area differences in lat/lon
off.bias <- (-log(1-exp(-r*arear)) - log( arear))
names(off.bias) <- "off.bias"

## Add offset to covariate stacks for model and prediction
covs_mod <- stack(covs_mod, off.bias)
covs_pred <- stack(covs_pred, off.bias)


## Sync NA across all input rasters & save aligned_mask and covariate stacks ####
## Find min non-NA set values across mask and covariates and sync NAs
##  see 0_functions.R
## *** NOTE: mask is too slow; replace with jgarber's glda_calc.R function...
reg_mask <- align.maskNA(covs_mod, reg_mask)
reg_mask <- align.maskNA(covs_pred, reg_mask)
covs_mod <- mask(covs_mod, reg_mask)
covs_pred <- mask(covs_pred, reg_mask)
saveRDS(reg_mask, file = file.path(rdata_path, "aligned_mask_aus.rds"))
saveRDS(readAll(covs_mod), file = file.path(rdata_path, "covs_mod_aus.rds"))
saveRDS(readAll(covs_pred), file = file.path(rdata_path, "covs_pred_aus.rds"))

covs_mod <- readRDS(file.path(rdata_path, "covs_mod_aus.rds"))
covs_pred <- readRDS(file.path(rdata_path, "covs_pred_aus.rds"))

## Log run time
data_runtime <- Sys.time()-data_start
cat(paste0("Data prep stop = ", Sys.time(), "\n"), file = masterlog, append = TRUE)
cat(paste0(">> Data prep runtime = ", data_runtime, " ", 
           attr(data_runtime, "units"), "\n"), file = masterlog, append = TRUE)


## II. Define model parameters ####
model_prep_start <- Sys.time()
cat(paste0("\n\nModel prep start = ", model_prep_start, "\n"), file = masterlog, append = TRUE)

## Generate quadrature (background) points
## FOR LATER: https://github.com/skiptoniam/qrbp
reg_mask0 <- readRDS(file.path(rdata_path, "aligned_mask_aus.rds"))
reg_mask0[which(is.na(reg_mask0[]))] <- 0
rpts <- rasterToPoints(reg_mask0, spatial=TRUE)
backxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])     
backxyz <- backxyz[,-1]
backxyz <- na.omit(cbind(backxyz, as.matrix(covs_mod)))
backxyzK <- backxyz[sample(nrow(backxyz), bkpts), ]

# ## Checks
# summary(covs_mod)
# summary(backxyzK)
# plot(reg_mask0, legend = FALSE)
# plot(rasterFromXYZ(backxyzK[,1:3]), col = "black", add = TRUE, legend=FALSE)

## Prediction points
predxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])     
predxyz <- predxyz[,-1]
predxyz <- na.omit(cbind(predxyz, as.matrix(covs_pred)))

## Specify covariates with interactions
interaction_terms <- c("X", "Y")

## Specify factor variable
factor_terms <- c("landuse", "pa")

## Define ppm variables to be used in the model 
## NOTE: X, Y & landuse are called in the ppm formula directly, hence removed here
ppm_terms <- names(backxyzK)[!grepl("off.bias", names(backxyzK))]
ppm_terms <- ppm_terms[!(ppm_terms %in% interaction_terms)]
ppm_terms <- ppm_terms[!(ppm_terms %in% factor_terms)]

# ## Fix max.lambda = 100
# lambdaseq <- round(sort(exp(seq(-10, log(max.lambda + 1e-05), length.out = n.fits)), 
#                         decreasing = TRUE),5)

## Estimate weights - by Skipton Wooley
## See ppmlasso::ppmdat 
ppmlasso_weights <- function (sp.xy, quad.xy, coord = c("X", "Y")){
  sp.col = c(which(names(sp.xy) == coord[1]), which(names(sp.xy) ==
                                                      coord[2]))
  quad.col = c(which(names(quad.xy) == coord[1]), which(names(quad.xy)
                                                        == coord[2]))
  X.inc = sort(unique(quad.xy[, quad.col[1]]))[2] -
    sort(unique(quad.xy[, quad.col[1]]))[1]
  Y.inc = sort(unique(quad.xy[, quad.col[2]]))[2] -
    sort(unique(quad.xy[, quad.col[2]]))[1]
  quad.0X = min(quad.xy[, quad.col[1]]) - floor(min(quad.xy[,
                                                            quad.col[1]])/X.inc) * X.inc
  quad.0Y = min(quad.xy[, quad.col[2]]) - floor(min(quad.xy[,
                                                            quad.col[2]])/Y.inc) * Y.inc
  X = c(sp.xy[, quad.col[1]], quad.xy[, quad.col[1]])
  Y = c(sp.xy[, quad.col[2]], quad.xy[, quad.col[2]])
  round.X = round((X - quad.0X)/X.inc) * X.inc
  round.Y = round((Y - quad.0Y)/Y.inc) * Y.inc
  round.id = paste(round.X, round.Y)
  round.table = table(round.id)
  wt = X.inc * Y.inc/as.numeric(round.table[match(round.id,
                                                  names(round.table))])
}

## Previously...
# spdat <- gbif
# species_names <- levels(factor(spdat$species))
# spwts <- list()
# for(i in seq_along(species_names)){
#       print(i)
#       spxy <- spdat[spdat$species %in% species_names[i], c(4,3)]
#       names(spxy) <- c("X", "Y")
#       cellNo <- cellFromXY(ar,spxy)
#       cellNoCounts <- table(cellNo)
#       tmp_cell_area <- extract(ar, spxy, fun = mean, na.rm = TRUE)*1000
#       tmp_dat <- data.frame(area=tmp_cell_area,cell_number=cellNo)
#       tmp_wts <- ((tmp_dat$area*1000)/cellNoCounts[as.character(cellNo)])/totarea
#       spwts[[i]] <- tmp_wts
# }

# ## Estimate weights for background points - NOT REQUIRED WITH ppmlasso_weights()
# ## calculate the area of the study area and then 
# ##  work out the weights based on the total area divided by number of points
# ar <- raster::area(reg_mask0)
# ar <- mask(ar,reg_mask)
# totarea <- cellStats(ar,'sum')*1000 ## in meters^2
# area_offset <- extract(ar, backxyzK[,c('X','Y')], small = TRUE, fun = mean, na.rm = TRUE)*1000 ## in meters
# bkgrd_wts <- c(totarea/area_offset)

# ## Estimate offset for unequal area due to projection [FOR GLOBAL ANALYSIS ONLY]
# temp.ala<-SpatialPoints(occdat[,c(1, 2)], proj4string=crs(reg_mask))
# 
# ## Craete empty raster to catch data
# r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
#             crs=crs(reg_mask), resolution=res(reg_mask))
# r[]<-0
# 
# ## Calculate area to adjust for 
# ## QUESTION: Difference between lambda_offset and area_offset?
# arear <- raster::area(r)
# lambda <- rasterize(temp.ala, arear, fun = function(x,...){length(x)})
# l1 <- crop(lambda, extent(-180, 0, -90, 90))
# l2 <- crop(lambda, extent(0, 180, -90, 90))   
# extent(l1) <- c(180, 360, -90, 90)
# lm <- merge(l1, l2)
# lm[is.na(lm)]<-0
# lm1 <- lm
# lambda_mask <- mask(lm,reg_mask)
# lambda_offset <- 1-exp(-1 * lambda_mask)
# lambda_offset[lambda_offset == 0] <- 1e-6
# area_rast_mask <- mask(raster::area(lm), reg_mask)
# ## If offset as just area: spp ~ blah + foo + offset(log(area_rast))
# ## IF offset as area + bias: spp ~ blah + foo + offset(log(off.area_rast)-log(off.bias))
# 
# writeRaster(lambda_offset, 
#             filename = file.path(output_path, "effort_offset_lambda_0360.tif"), 
#             format = "GTiff",overwrite = TRUE)
# saveRDS(lambda_offset, file.path(output_path, "effort_offset_lambda_0360.rds"))
# 
# writeRaster(area_rast_mask,
#             filename = file.path(output_path, "area_offset_0360.tif"), 
#             format = "GTiff", overwrite = TRUE)
# saveRDS(area_rast_mask,filename = file.path(output_path, "area_offset_0360.rds"))


## Initialise log file for ppm computation times comparison ####
timelog <- paste0(output_path, "/ppm_timelog_",Sys.Date(), ".txt")
writeLines(c(""), timelog)
cat(paste(c("i", "species_name", "pr_pts", "back_pts",	"ppmfit_0folds_min",	
            "ppmfit_1fold_avg_min",	"ppmfit_5folds_min",	"ppmeval_train_sec",	
            "ppmeval_test_sec", "distribution", "\n"), collapse = ","), 
    file = timelog, append = T)


## Define model function ####
fit_ppms_apply <- function(i, spdat, bkdat, interaction_terms, 
                           ppm_terms, species_names, mask_path, 
                           n.fits=50, min.obs = 60, 
                           modeval = TRUE, 
                           output_list = FALSE) {
  
  ## Initialise log file for species
  cat("Fitting a ppm to", species_names[i],"\nThis is the", i,"^th model of",length(species_names),"\n")
  specieslog <- paste0(output_path, "/ppm_log_",gsub(" ","_",species_names[i]),"_",gsub("-", "", Sys.Date()), ".txt")
  writeLines(c(""), specieslog)
  cat(paste("Fitting a ppm to", species_names[i],"\nThis is the", i,"^th model of",length(species_names),"\n"),
      file = specieslog, append = T)
  cat(paste("-------------------------------------------\n"),
      file = specieslog, append = T)
  ppm_dat_start <- Sys.time()
  
  ## Define species specific data & remove NA (if any)
  spxy <- spdat[spdat$species %in% species_names[i], c(1,2)]
  names(spxy) <- c("X", "Y")
  
  # # ## Check if points fall outside the mask
  # plot(reg_mask)
  # points(spxy)
  
  ## Move species occurrence points falling off the mask to nearest 'land' cells
  ## First find points which fall in NA areas on/off the raster
  reg_mask <- readRDS(mask_path)
  vals <- extract(reg_mask, spxy)
  outside_mask <- is.na(vals)
  if(sum(outside_mask) > 0){
    outside_pts <- spxy[outside_mask, ]
    ## find the nearest land within 5 decimal degrees of these
    land <- nearestLand(outside_pts, reg_mask, 1000000)
    ## replace points falling in NA with new points on nearest land
    spxy[outside_mask, ] <- land
    ## count how many were moved
    sum(!is.na(land[, 1]))
  }
  # ## WARNING: We lose data again i.e. number of unique locations is reduced. This can be problematic for ppms...
  # nrow(unique(outside_pts))
  # nrow(unique(land))
  
  
  ## Extract covariates for presence points
  r <- dropLayer(covs_mod, "landuse")
  spxy_nolu <- extract(r, spxy,
                       fun = mean, na.rm = TRUE)
  # ## NOTE: This step takes too logn to run ...
  # spxyz_nolu <- extract(r, spxy, buffer = 1000000, small = TRUE, 
  #                       fun = mean, na.rm = TRUE)
  
  ## For landuse: Take the non-NA value at shortest distance from point
  r = covs_mod$landuse ## landuse raster
  spxy_lu <- extract(r, spxy, method='simple', na.rm=TRUE, factor = TRUE) 
  ## NOTE: extract(covs_mod[[1]], spxy, method = 'simple') gives NAs
  ## WARNING: Even with the na.rm argumennt, we may still give NAs.
  ## Alternate approach (below) but it takes too long to run; no good, abort!
  # spxy_lu <- apply(X = spxy, MARGIN = 1, FUN = function(X) r@data@values[which.min(replace(distanceFromPoints(r, X), is.na(r), NA))])
  
  spxyz <- cbind(spxy, spxy_lu, spxy_nolu)
  names(spxyz)[grep("lu", names(spxyz))] <- "landuse"
  spxyz <- na.omit(spxyz)
  
  ## Check for NAs in data and return error if any
  if (all(is.na(spxyz))) {
    stop("NA values not allowed in model data.")
  }
  
  ## Check if number of observationns for a species is > 20, return NULL
  nk <- nrow(spxyz)
  
  if (nk < 20) {
    return(NULL)
  } else {
    
    ## Build model
    ## Add 'Pres' column and calculate weights for PPM data
    ppmxyz <- rbind(cbind(spxyz, Pres = 1),
                    cbind(bkdat, Pres = 0))
    wts <- ppmlasso_weights(sp.xy = spxyz, quad.xy = bkdat)
    ppmxyz <- cbind(ppmxyz, wt = wts)
    
    ## Log run time for data preparationn
    ppm_dat_end <- Sys.time() - ppm_dat_start
    cat(paste("\nRun time for data preparation: ", ppm_dat_end, attr(ppm_dat_end, "units"), "\n"), 
        file = specieslog, append = T)
    
    ## Specify PPM formula based on number of observatios
    ## If number of observations < 20, do not fit a model.
    ## If number of observations >= 20, fit a model accordig to the followig rule:
    ##  1. If number of observations <= min.obs (default as 60), 
    ##    use only the spatial covariates i.e. lat and long and offset. This gives 6 
    ##    terms: X, Y, X^2, Y^2, X:Y (linear, quadratic with interaction), offset
    ##  2. If number of observations is between min.obs and min.obs+10*, add landuse covariate.
    ##    *10 as per the one-in-ten rule (see point 3).Factor variables are fit as 1 varible only, 
    ##    irrespective of the number of levels within 
    ##  3. If number of observatons > min.obs, 
    ##    use the one-in-ten rule (https://en.wikipedia.org/wiki/One_in_ten_rule)
    ##    where, an additional covariate is added for every 10 additonal obeservations.
    ##    Because we fit poly with degree 2, adding a covariate will add two 
    ##    terms x and x^2 to the model. Therefore, to operationalize this rule, 
    ##    we add 1 raw covariate for every 20 additional observations. 
    ## NOTE: if area offset considered, use "+ offset(log(off.area)-log(off.bias))"
    
    ## *** ISSUE: There is inherent bias in how covariates are added i.e.
    ##  which covariates are included in the model as they are selected 
    ##  based on their order in the dataframe; pa being added last.
    ## *** ISSUE: factor(pa) only included when sufficient obs to fit all variables...
    
    if(nk <= min.obs){
      ppmform <- formula(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                ", degree = 2, raw = FALSE) + offset(log(off.bias))",collapse =""))
    } else {
      if(nk > min.obs & nk <= min.obs + 10){
        ppmform <- formula(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                  ", degree = 2, raw = FALSE) + factor(landuse) + offset(log(off.bias))",collapse ="")) 
      } else {
        extra_covar <- ceiling((nk - min.obs)/20) ## upto nk = 360 obs
        if(extra_covar <= length(ppm_terms)) { 
          ppmform <- formula(paste0(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                           ", degree = 2, raw = FALSE) + factor(landuse)"), 
                                    paste0(" + poly(", ppm_terms[1:extra_covar],
                                           ", degree = 2, raw = FALSE)", 
                                           collapse = ""), " + offset(log(off.bias))" ,
                                    collapse =""))
        } else {
          if(extra_covar > length(ppm_terms)) { ## Fit all covariates + pa
            extra_covar <- length(ppm_terms)
            ppmform <- formula(paste0(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                             ", degree = 2, raw = FALSE) + factor(landuse)"), 
                                      paste0(" + poly(", ppm_terms[1:extra_covar],
                                             ", degree = 2, raw = FALSE)",
                                             collapse = ""), " + factor(pa) + offset(log(off.bias))" ,
                                      collapse =""))
          }
        }
      }
    }
    
    ## Fit ppm & save output
    cat(paste("\nFitting ppm model for species ",i , ": ", spp[i], " @ ", Sys.time(), "\n"), 
        file = specieslog, append = T)
    cat(paste("   # presence points for species (original) = ", 
              dim(spxy)[1], "\n"), file = specieslog, append = T)
    cat(paste("   # presence points for species (extracted) = ", 
              dim(spxyz)[1], "\n\n"), file = specieslog, append = T)
    cat(paste("   # total data points for species (presence + background) = ", 
              dim(ppmxyz)[1], "\n\n"), file = specieslog, append = T)
    
    if(!modeval){ ## if model_eval == FALSE
      ## Fit model
      ppm_mod_start <- Sys.time()
      
      mod <- tryCatch(ppmlasso(formula = ppmform, data = ppmxyz, n.fits = n.fits, 
                               criterion = "bic", standardise = FALSE),
                      error = function(e){ 
                        cat(paste("\nModel ",i," for species ", spp[i], " has errors. \n"),
                            file = specieslog, 
                            append = T)
                        return(NA) 
                      })
      
      if(is.na(mod)){
        cat(paste("\nERROR: Output for model ",i , " for species ", spp[i], " is NA. \n"), 
            file = specieslog, append = T)
      }
      
      # ## Print warnings to screen (for sequence runs only)
      # cat('Warnings for ', species_names[i],':\n')
      # warnings()
      
      # ## Capture messages and errors in file (for sequence runs only)
      # sink(specieslog, type = "message", append = TRUE, split = FALSE)
      # try(warnings())
      # ## reset message sink and close the file connection
      # sink(type="message")
      # close(specieslog)
      
      ## Log run time for model run
      ppm_mod_end <- Sys.time() - ppm_mod_start
      cat(paste(">> Run time for ppm model fitting: ", ppm_mod_end, 
                attr(ppm_mod_end, "units"), "\n"), 
          file = specieslog, append = T)
      
      ## Log time for comparison
      time_vec <- c(dim(spxy)[1], dim(bkdat)[1], ppm_mod_end, rep(NA, 5))
      cat(paste(c(i, trimws(gsub(" ", "_", spp[i])), time_vec, "\n"), collapse = ","),
          file = timelog, append = T)
      
      if (output_list == TRUE) {
        return(mod)
      } else {
        saveRDS(mod, paste0(output_path, "/ppm_out_",gsub(" ","_",species_names[i]),"_",gsub("-", "", Sys.Date()), ".rds"))
      }
      
    } else { ## if model_eval == TRUE
      
      ppm_mod_start <- Sys.time()
      
      ## Run a basic k-fold for a ppmlasso
      ## ... Brendan might want to check this
      nk <- k <- 5 # number of folds
      kfold_train <- list()
      kfold_test <- list()
      rm(folds)
      
      ## Sample without replacement 
      createFolds <- function(x,k){
        n <- nrow(x)
        x$folds <- rep(1:k,length.out = n)[sample(n,n)]
        x
      }
      
      folds <- plyr::ddply(ppmxyz,.(Pres),createFolds,k = nk)
      for(ii in 1:nk){
        kfold_test[[ii]]<-folds[folds$folds==ii,] ## dim(folds)[1]/k rows in kfolds_test[[ii]]
        kfold_train[[ii]]<-folds[folds$folds!=ii,] ## all rows not in kfolds_test[[ii]]
        # print(colSums(folds[folds$folds==ii,],na.rm = T))
      }
      
      ppmCV <- list()
      ppmCV_dats <- list()
      ppmCV_end_kmodels <- c() ## to log time per run
      
      ## Fit model
      ## fit k models with a random (without replacement) basic K-fold
      for (ii in seq_len(k)){
        ppmCV_dats[[ii]] <- kfold_train[[ii]]
        ppmCV[[ii]] <- tryCatch(ppmlasso(formula = ppmform, 
                                         data = ppmCV_dats[[ii]], 
                                         n.fits = n.fits, 
                                         criterion = "bic", 
                                         standardise = FALSE),
                                error = function(e){ 
                                  cat(paste("\n Model",ii,"/",k," for species ", i , ": ", spp[i], " has errors. \n"),
                                      file = specieslog, 
                                      append = T)
                                  return(NA) 
                                })
        
        ppmCV_end <- Sys.time()
        
        if(is.na(ppmCV[[ii]])){
          cat(paste("\nERROR: Output for model ",ii,"/",k," for species ", i , ": ", spp[i], " is NA. \n"), 
              file = specieslog, append = T)
          stop("\nERROR: Output for model ",ii,"/",k," for species ", i , ": ", spp[i], " is NA \n")
        } else {
          cat(paste("\nModel",ii,"/",k, " done @ ", ppmCV_end, "\n"),
              file = specieslog, append = T)
        }
        ppmCV_end_kmodels[[i]] <- ppmCV_end
      }
      
      ## Log run time for k models run
      ppm_mod_end <- Sys.time() - ppm_mod_start
      cat(paste(">> Run time for all", k," models for species ",
                i , ": ", ppm_mod_end, 
                attr(ppm_mod_end, "units"), "\n"), 
          file = specieslog, append = T)
      
      ## Model evaluation with training data
      ppm_traineval_start <- Sys.time()
      
      cat(paste("\n\nModel evaluation with training data: ", Sys.time(), "\n"), 
          file = specieslog, append = T)
      
      cell_area <- prod(res(reg_mask))
      train.preds <- lapply(1:k,function(x)predict(ppmCV[[x]],
                                                   newdata=ppmCV_dats[[x]])*cell_area)
      model.evaluations <- lapply(1:k,function(x)dismo::evaluate(train.preds[[x]][ppmCV_dats[[x]]$Pres==1],
                                                                 train.preds[[x]][ppmCV_dats[[x]]$Pres==0]))
      
      # mean AUC from k-folds the model
      train.meanAUC <- mean(sapply(model.evaluations,function(x)x@auc))
      
      # ## ROC curves
      # par(mfrow=c(3,2))
      # sapply(1:k,function(ii)plot(model.evaluations[[ii]],"ROC"))
      # 
      # ## TPR plots
      # par(mfrow=c(3,2))
      # sapply(1:k,function(ii)plot(model.evaluations[[ii]],"TPR"))
      # 
      # ## Density plots
      # par(mfrow=c(3,2))
      # sapply(1:k,function(ii)density(model.evaluations[[ii]]))
      
      ## Log run time for evaluation using training data
      ppm_traineval_end <- Sys.time() - ppm_traineval_start
      cat(paste(">> Run time for model evaluation with training data: ", 
                ppm_traineval_end, attr(ppm_traineval_end, "units"), "\n"), 
          file = specieslog, append = T)
      
      
      ## Model evaluations with test data
      ppm_testeval_start <- Sys.time()
      
      cat(paste("\n\nModel evaluation with test data: ", Sys.time(), "\n"), 
          file = specieslog, append = T)
      
      test.preds <- lapply(1:k,function(x)predict(ppmCV[[x]],
                                                  newdata=kfold_test[[x]])*cell_area)
      test.evaluations <- lapply(1:k,function(x)dismo::evaluate(test.preds[[x]][kfold_test[[x]]$Pres==1],
                                                                test.preds[[x]][kfold_test[[x]]$Pres==0]))
      
      
      ## Mean AUC from k-folds the model
      test.meanAUC <- mean(sapply(test.evaluations,function(x)x@auc))
      
      # ## ROC curves
      # par(mfrow=c(3,2))
      # sapply(1:k,function(ii)plot(test.evaluations[[ii]],"ROC"))
      # 
      # ## TPR plots
      # par(mfrow=c(3,2))
      # sapply(1:k,function(ii)plot(test.evaluations[[ii]],"TPR"))
      # 
      # ## Density plots
      # par(mfrow=c(3,2))
      # sapply(1:k,function(ii)density(model.evaluations[[ii]]))
      
      ## Log run time for evaluation using test data
      ppm_testeval_end <- Sys.time() - ppm_testeval_start
      cat(paste(">> Run time for model evaluation with test data: ", 
                ppm_testeval_end, attr(ppm_testeval_end, "units"), "\n"), 
          file = specieslog, append = T)
      
      ## Save output
      mod <- list(ppmCV = ppmCV, 
                  train_eval = model.evaluations, 
                  train_AUC = train.meanAUC, 
                  test_eval = test.evaluations,
                  test_AUC = test.meanAUC)
      
      ## Log time for comparison
      time_vec <- c(dim(spxy)[1], dim(bkdat)[1], NA, mean(ppmCV_end_kmodels), 
                    ppm_mod_end, ppm_traineval_end, ppm_testeval_end, NA)
      cat(paste(c(i, trimws(gsub(" ", "_", spp[i])), time_vec, "\n"), collapse = ","),
          file = timelog, append = T)
      
      gc()
      
      if (output_list == TRUE) {
        return(mod)
      } else {
        saveRDS(mod, paste0(output_path, "/ppm_out_",gsub(" ","_",species_names[i]),"_",gsub("-", "", Sys.Date()), ".rds"))
        out <- list(likelihood = sapply(mod$ppmCV, `[[`, 8), train_AUC = mod$train_AUC, test_AUC = mod$test_AUC)
        return(out)
      }
    }
  }
}

## Log run time
model_prep_runtime <- Sys.time()-model_prep_start
cat(paste0("Model prep stop = ", Sys.time(), "\n"), file = masterlog, append = TRUE)
cat(paste0(">> Model prep runtime = ", model_prep_runtime, " ",
           attr(model_prep_runtime, "units"), "\n"), file = masterlog, append = TRUE)

## III. Fit models & log info ####
ppm_start <- Sys.time()
cat(paste0("\n\nModel runs start = ", ppm_start, "\n"), file = masterlog, append = TRUE)

spdat <- as.data.frame(occdat)
bkdat <- backxyzK
modeval <- TRUE
spp <- unique(spdat$species)[1:10]
mc.cores <- 10
seq_along(spp)
ppm_models <- list()

plan(multiprocess, workers = mc.cores)
options(future.globals.maxSize = +Inf) ## CAUTION: Set this to a value, e.g. availablecores-1?/RAM-10?
ppm_models <- future_lapply(1:length(spp), fit_ppms_apply, spdat,
                            bkdat, interaction_terms, ppm_terms,
                            species_names = spp,
                            mask_path = file.path(rdata_path, "aligned_mask_aus.rds"),
                            n.fits = n.fits, min.obs = 60,
                            modeval = modeval, output_list = FALSE)
names(ppm_models) <- tolower(gsub(" ", "_", spp))
# ppm_models <- parallel::mclapply(1:length(spp), fit_ppms_apply, spdat,
#                                  bkdat, interaction_terms, ppm_terms,
#                                  species_names = spp,
#                                  mask_path = file.path(rdata_path, "aligned_mask_aus.rds"),
#                                  n.fits = n.fits, min.obs = 60, mc.cores = mc.cores,
#                                  modeval = modeval, output_list = FALSE)

# ppm_models <- lapply(1:length(spp), fit_ppms_apply, spdat,
#                      bkdat, interaction_terms, ppm_terms,
#                      species_names = spp,
#                      mask_path = file.path(rdata_path, "aligned_mask_aus.rds"),
#                      n.fits=10, min.obs = 60, modeval = TRUE)

# ppm_models <- fit_ppms_apply(1, spdat, bkdat, interaction_terms,
#                              ppm_terms, species_names = spp,
#                              mask_path = file.path(rdata_path, "aligned_mask_aus.rds"),
#                              n.fits=5, min.obs = 60, modeval = TRUE)


saveRDS(ppm_models, paste0(output_path, "/ppm_out_", gsub("-", "", Sys.Date()), ".rds"))

ppm_runtime <- Sys.time()-ppm_start
cat(paste0("\n\nfit_ppms_apply() run time for ", length(spp), " species: ", 
           ppm_runtime, " ", attr(ppm_runtime, "units"), "\n\n"))
cat(paste0("Model runs stop = ", Sys.time(), "\n"), file = masterlog, append = TRUE)
cat(paste0(">> Model runtime = ", ppm_runtime, " ",
           attr(ppm_runtime, "units"), "\n"), file = masterlog, append = TRUE)

cat(paste0("\n\n---------------------------------------------\n"), 
    file = masterlog, append = TRUE)
cat(paste0("Run details: \n", 
           "modeval = ", modeval, "\n",
           "# species = ", length(spp), "\n",
           "# background points = ", bkpts, "\n",
           "# n.fits = ", n.fits, "\n"),
    file = masterlog, append = TRUE)
cat(paste0("\n---------------------------------------------\n\n"), 
    file = masterlog, append = TRUE)

cat("Session details: \n", file = masterlog, append = TRUE)
cat(paste0("Date:\n"), file = masterlog, append = TRUE)
sink(file = masterlog, append = TRUE)
Sys.Date()
sink(NULL)
cat("\n\nSession Info:\n", file = masterlog, append = TRUE)
sink(file = masterlog, append = TRUE)
sessionInfo()
sink(NULL)


## Extract mean likelihood
meanlik <- sapply(ppm_models, `[[`, 1)
as.data.frame(colMeans(meanlik))

gc()


## Load ppm_models from rds files ####
## Species ppm files
ppm_files <- list.files(output_path, pattern = "ppm_out", full.names = TRUE)

## Extract species names from files
ppm_sp <- regmatches(ppm_files,regexec(paste0("ppm_out_(.*?)_", datestamp), ppm_files))
ppm_sp <- sapply(ppm_sp, `[[`, 2)

## Create list of results
ppm_out <- list()
for(i in 1:length(ppm_files)){
  ppm_out[[i]] <- readRDS(ppm_files[i])
}
names(ppm_out) <- ppm_sp

## Check model fits/predictions
names(ppm_out[[1]])
ppm_test_auc <- sapply(ppm_out, `[[`, 5)


## Check for species distributions
reg_mask <- readRDS(file.path(rdata_path, "aligned_mask_aus.rds"))
par(mfrow = c(2,5))
for (i in 1:length(spp)){
  ## Define species specific data & remove NA (if any)
  spxy <- spdat[spdat$species %in% spp[i], c(1,2)]
  names(spxy) <- c("X", "Y")
  
  ## Check if points fall outside the mask
  plot(reg_mask, axes = FALSE, legend = FALSE, box = FALSE, main = spp[i])
  points(spxy)
}

## Check for NULL models
length(ppm_models[!sapply(ppm_models,is.null)])
length(ppm_models[sapply(ppm_models,is.null)])
which(sapply(ppm_models,is.null))


## IV. Catch errors in models & save outputs ####
error_list <- list()
model_list <- list()
n <- 1
m <- 1
errorfile <- paste0(output_path, "/errorfile_", gsub("-", "", Sys.Date()), ".txt")

for(i in 1:length(ppm_models)){
  if(!class(ppm_models[[i]])[1] == "try-error") {
    model_list[[n]] <- ppm_models[[i]]
    n <- n+1
  }else{
    print(paste0("Model ",i, " for '", spp[i], "' has errors"))
    cat(paste(i, ",", spp[i], "\n"),
        file = errorfile, append = T)
    error_list[[m]] <- ppm_models[[i]]
    m <- m+1
  }
}

saveRDS(model_list, file = paste0(output_path, "/modlist_",  gsub("-", "", Sys.Date()), ".rds"))
saveRDS(error_list, file = paste0(output_path, "/errlist_",  gsub("-", "", Sys.Date()), ".rds"))



## V. Predict and save output ####
## Define prediction function ####
predict_ppms_apply <- function(i, models_list, newdata, bkdat, RCPs = c(26, 85)){
  
  cat("Predicting ", i,"^th model\n")
  
  if(class(models_list)== "try-error") { ## redundant because error models are removed
    return(NULL)
  } else {
    
    predmu <- list()
    
    ## Current predictions
    preddat <- newdata[which(names(newdata) %in% 
                               names(newdata)[-grep("26|85", names(newdata))])]
    predmu$current <- predict.ppmlasso(models_list[[i]], newdata = preddat)
    
    ## Future predictions
    for (rcp in RCPs) {
      if (rcp == 26) {
        preddat <- newdata[which(names(newdata) %in% 
                                   names(newdata)[-grep("85|bio", names(newdata))])]
        names(preddat) <- names(bkdat)[1:(length(names(bkdat))-1)]
        predmu$rcp26 <- predict.ppmlasso(models_list[[i]], newdata = preddat)
        # predmu <- 1-exp(-predmu) ## gives reative probabilities
        
      } else {
        preddat <- newdata[which(names(newdata) %in% 
                                   names(newdata)[-grep("26|bio", names(newdata))])]
        names(preddat) <- names(bkdat)[1:(length(names(bkdat))-1)]
        predmu$rcp85 <- predict.ppmlasso(models_list[[i]], newdata = preddat)
        # predmu <- 1-exp(-predmu) ## gives reative probabilities
      }
      rm(preddat)
    }
    return(predmu)
  }
}

# ## Bioregions layer _ TO USE FOR CLIPPING...
# bioreg <- readRDS(file.path(rdata_path, "bioregions_aus.rds")) 

newdata <- predxyz
bkdat <- backxyzK
RCPs <- c(26, 85)
now <- Sys.time()
prediction_list <- parallel::mclapply(seq_along(model_list), predict_ppms_apply,
                                      model_list, newdata, bkdat, RCPs, mc.cores = mc.cores)
pred_time <- Sys.time()-now
saveRDS(prediction_list, file = paste0(output_path, "/predlist_",  gsub("-", "", Sys.Date()), ".rds"))



## Locate errors and rerun analysis for species with errors ###
##  To be automated if error problem is not solved by species grouping
##  At the moment, it appears that error might be when species data is spatially restricted.
error_species <- read.table(output_path, "/errorfile_1_20190725.txt", header = FALSE, sep = ",")
colnames(error_species) <- c("index", "species")
error_index <- error_species$index

names(model_list) <- tolower(gsub(" ","_", levels(factor(spdat$species))[-error_index]))
names(prediction_list) <- tolower(gsub(" ","_", levels(factor(spdat$species))[-error_index]))

spdat <- gbif
bkdat <- backxyzK
bkwts <- bkgrd_wts
spp <- levels(factor(spdat$species))[error_index]
seq_along(spp)
mc.cores <- 1
error_models <- parallel::mclapply(1:length(spp), fit_ppms_apply, spdat, #spwts,
                                   bkdat, bkwts, interaction_terms, ppm_terms,
                                   species_names = spp, n.fits=100, min.obs = 50, mc.cores = mc.cores)
names(error_models) <- tolower(gsub(" ","_", levels(factor(spdat$species))[error_index]))

newdata <- predxyz
bkdat <- backxyzK
RCPs <- c(26, 85)
error_pred <- parallel::mclapply(seq_along(error_models), predict_ppms_apply,
                                 error_models, newdata, bkdat, RCPs, mc.cores = mc.cores)
names(error_pred) <- tolower(gsub(" ","_", levels(factor(spdat$species))[error_index]))

errorfile <- paste0(output_path, "/errorfile_2_", gsub("-", "", Sys.Date()), ".txt")
n <- length(model_list)+1
m <- length(errlist)+1
error_list <- list()
for (i in 1:length(error_models)){
  if(!class(error_models[[i]])[1] == "try-error") {
    model_list[[n]] <- error_models[[i]]
    n <- n+1
  }else{
    print(paste0("Model ",i, " for '", spp[i], "' has errors"))
    cat(paste(i, ",", spp[i], "\n"),
        file = errorfile, append = T)
    error_list[[m]] <- error_models[[i]]
    m <- m+1
  }
}
names(model_list)[((length(model_list)-length(error_models))+1):length(model_list)] <- names(error_models)

n <- length(prediction_list) + 1
for (i in 1:length(error_pred)) {
  prediction_list[n] <- error_pred[i]  
  n <- n + 1
}
names(prediction_list)[((length(prediction_list)-length(error_pred))+1):length(prediction_list)] <- names(error_pred)


## Catch remianing errors ####
n <- 1
m <- 1
catch_errors(seq_along(error_models), ppm_models = model_list, species_names = spp, errorfile = errorfile)



## Log job duration
job_end <- Sys.time()
job_duration = job_end - job_start
cat(paste0(">> Job end = ", job_end, "\n\n"), file = masterlog, append = TRUE)
cat(paste0("Job duration = ", job_duration, "\n\n"), file = masterlog, append = TRUE)
cat(paste0("Date:\n"), file = masterlog, append = TRUE)
sink(file = masterlog, append = TRUE)
Sys.Date()
sink(NULL)
cat("\n\nSession Info:\n", file = masterlog, append = TRUE)
sink(file = masterlog, append = TRUE)
sessionInfo()
sink(NULL)

rm(list=ls())
gc()



## ------ EXTRAS ----------------

## Testing for corerelations in vcovariates ####
## There is no hard and fast rule about how many covariates to fit to data, and it will change depending on the data and the amount of information in each observation and how they vary with the covariates, how the covariates are correlated ect... But to start I'd only fit sqrt(n_observations) covariates. So if you have 20 occurrences that's 4-5 covariates (including the polynomials) so that's only two variables! You might have to identify the most important variable and start from there.

## see slide 14 & 26 in: http://www.bo.astro.it/~school/school09/Presentations/Bertinoro09_Jasper_Wall_3.pdf

## let's do a PCA on the data to work out which are the most variable coefs.
Xoriginal=t(as.matrix(backxyzK))

# Center the data so that the mean of each row is 0
rowmns <- rowMeans(Xoriginal)
X <-  Xoriginal - matrix(rep(rowmns, dim(Xoriginal)[2]), nrow=dim(Xoriginal)[1])

# Calculate P
A <- X %*% t(X)
E <- eigen(A,TRUE)
P <- t(E$vectors)

dimnames(P) <- list(colnames(backxyzK),paste0("PCA",1:ncol(P)))
df <- as.data.frame(t(P[,1:5]))
df$row.names<-rownames(df)
library(reshape2)
library(ggplot2)
long.df<-reshape2::melt(df,id=c("row.names"))
pca_plot <- ggplot2::ggplot(long.df,aes(x=row.names,y=variable,fill=value))+
  geom_raster()+
  scale_fill_viridis_c()+
  theme_minimal()

## not surprising that elevation, slope and aspect are all correlated (choose one).
pca_plot


##ALTERNATIVE - Specify ppm formula based on number of observations for a species
nk <- ceiling(sqrt(nrow(spxy)))
## Specify ppm formula - independent with quadratic terms
## if low number of observations, (i.e. =< 25), just use the covariates interaction terms for the spatial variables
if(nk <= 5){ ## because linear and quad for X ad Y gives 5 terms
  ppmform <- formula(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                            ", degree = 2, raw = FALSE)",collapse =""))
  ## including land-use change for species with low number of observations
  # paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
  #        ", degree = 2, raw = FALSE) + poly(landuse, degree = 2, raw = FALSE)",collapse ="")
  
} else  {
  ## if number of observations, (i.e. > 25), fit independent with linear and quadratic terms startign with landuse
  extra_covar <- ceiling((nk - 5)/2) ## fit additional covariates based on 
  if(extra_covar > 10) extra_covar <- 10 ## if nrow(spxy) > 630 and nk > 25, only then all parameters are used!! ,...fix
  ppmform <- formula(paste0(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                   ", degree = 2, raw = FALSE)"), 
                            paste0(" + poly(", ppm_terms[1:extra_covar],
                                   ", degree = 2, raw = FALSE)", 
                                   collapse = ""), collapse =""))
}

## Explore relationship between a method to specify cut-off and 3obs
n.obs <- seq(20,1000, 10)
cutoff <- ceiling(sqrt(n.obs))
add.params <- ceiling((sqrt(n.obs) - 5)/2)
max.params <- length(names(covs_mod))
obs.cutoff <- n.obs[which(add.params > max.params)[1]]
add.params[which(add.params > max.params)] <- max.params
plot(n.obs,add.params, type ="l")
abline(v=obs.cutoff)
text(obs.cutoff - 20, 4, paste0("max obs for params cut off: ", obs.cutoff), srt=90)

