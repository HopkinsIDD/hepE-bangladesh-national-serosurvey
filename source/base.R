## utility and base functions for HEV serosurvey analyses
##
reload_source <- function(){

  ## load all packages
  p_load(dplyr,ggplot2,readxl,stringr,readr,sf,mgcv,INLA,raster,conflicted,survey,srvyr,janitor)

  source(here("source","base.R"))
}

## aligns Zila names with gadm name standards
fix_adm2names <- function(dat){

  dat %>% mutate(adm2name = recode(adm2name,
                                   "Barguna" = "Borgona",
                                   "Munshiganj" = "Munshigonj",
                                   "Mymensingh" = "Nasirabad",
                                   "Netrokona" = "Netrakona",
                                   "Narsingdi" = "Narshingdi",
                                   "Kushtia" = "Kustia",
                                   "Chapai Nawabganj" = "Nawabganj",
                                   "Gaibandha" = "Gaibanda",
                                   "Rangpur" = "Rongpur",
                                   "Habiganj" = "Hobiganj",
                                   "Maulvibazar" = "Moulvibazar"))
}

## transforms sf file to Bangladesh Transverse Mercator projection
## https://spatialreference.org/ref/epsg/gulshan-303-bangladesh-transverse-mercator/proj4/
## transforms sf file to Bangladesh Transverse Mercator projection
transform_to_btm <- function(my_obj){

  if ("sf" %in% class(my_obj)){
    rc <- st_transform(my_obj,
                       crs="+proj=tmerc +lat_0=0 +lon_0=90 +k=0.9996 +x_0=500000 +y_0=0 +a=6377276.345 +b=6356075.41314024 +towgs84=283.7,735.9,261.1,0,0,0,0 +units=m +no_defs")
  } else if(class(my_obj) == "RasterLayer"){
    rc <- raster::projectRaster(my_obj,
                                crs="+proj=tmerc +lat_0=0 +lon_0=90 +k=0.9996 +x_0=500000 +y_0=0 +a=6377276.345 +b=6356075.41314024 +towgs84=283.7,735.9,261.1,0,0,0,0 +units=m +no_defs")
  } else {
    stop("can only transform rasterlayers and sf objects at the moment")
  }

}

##' loads GADM data
##' could have used raster::getData('GADM',country="BGD",level=level)
##' but this seems more robust to any changes in GADM
##' @title
##' @return
##' @author Andrew A
load_bgd_shp <- function(level=2){

  if(level==2){
    rc <- readRDS(here("data","BGD_adm2.rds")) %>% sf::st_as_sf() %>% mutate(adm2name = NAME_2)
  } else if(level == 1){
    rc <- readRDS(here("data","BGD_adm1.rds")) %>% sf::st_as_sf() %>% mutate(adm1name = NAME_1)
  } else if(level == 0){
    rc <- readRDS(here("data","BGD_adm0.rds")) %>% sf::st_as_sf()
  }
  return(rc)
}

## makes age asex pyramid for surveyed population
make_age_pyramid <- function(dat){
  cdat <- read_csv(here("data","census_data_sex_age.csv"),skip=1) %>% filter(Age != "100+" & Age != "95-99") %>%
    mutate(prop_male = 100*`Male Population`/sum(`Both Sexes Population`),
           prop_female = 100*`Female Population`/sum(`Both Sexes Population`)
    )

  dat <- dat %>% mutate(age_cat = cut(age,seq(0,100,by=5),right = FALSE) %>% ordered %>% droplevels)

  ind_summary <- dat %>% summarize(mean_age = mean(age,na.rm=T),
                                   prop_u5 = mean(age<5,na.rm=T))

  age_sex_dat <- dat %>% filter(!is.na(age_cat)) %>% group_by(sex,age_cat) %>%
    summarize(count = n()) %>%
    ungroup %>%
    mutate(prop = 100* count / sum(count))


  cdat$age_cat <-  ordered(levels(age_sex_dat$age_cat),levels = levels(age_sex_dat$age_cat))

  age_sex_dat %>% ggplot(aes(x=age_cat)) +
    geom_histogram(data = age_sex_dat %>% filter(sex=="Male"),aes(x=age_cat,y=prop),
                   binwidth = 1,fill="steelblue",stat="identity") +
    geom_histogram(data = age_sex_dat %>% filter(sex=="Female") %>% filter(!is.na(age_cat)),aes(y=-prop),
                   binwidth = 1,fill="pink",stat="identity") +
    geom_point(data = cdat,aes(x=age_cat,y=prop_male)) +
    geom_point(data = cdat,aes(x=age_cat,y=-prop_female)) +
    coord_flip() +
    scale_y_continuous(breaks = -10:10,labels = abs(-10:10),limits=c(-7,7)) +
    theme_minimal() +
    xlab("age category") +
    ylab("percent of total population")

}



#' loads HEV lab data and merges with
#' survey data
#' also does a little relabeling etc
load_and_clean_hev_data <- function(){

  dat <- read_xlsx(here("data","HEV ELISA_national serosurvey_Result.xlsx"),na = c("OVERFLOW","OVRFLO"),
            col_types = c('text','text',rep('numeric',11),rep('text',2)) ) %>%
    setNames(c('date','sample_id','pc_1','pc_2','pc_avg','nc_1','nc_2','nc_3','nc_avg','nc_adj','od','cutoff','absorb_over_cutoff','seropos_yn','seropos')) %>%
    mutate(
        pos = ifelse(absorb_over_cutoff>1.1,1,ifelse(absorb_over_cutoff<0.9,0,NA))) %>%
  select(sample_id,od,cutoff,absorb_over_cutoff,pos) %>%
  filter(!sample_id %in% c("430101-B","740103-A")) %>% ## two aliqots of two samples run. same +/- results and similar ODs so choosing one of each
    mutate(sample_id = str_replace(sample_id,"-C",""))

  ## the only NAs in OD are for overflow so we will assign those as positive
  dat <- dat %>% mutate(pos = ifelse(is.na(od),1,pos))

  sample_meta_dat <- read_csv(here("data","CombinedDataY1Y2.csv")) %>%
    rename(sample_id = indid,community_id = CommunityID) %>%
    mutate(sample_id = as.character(sample_id),
           year = substr(indidyr,nchar(indidyr),nchar(indidyr))) %>%
    filter(year ==2) %>%
    mutate(age = ifelse(age == 999,NA,age),
           incomeSumm = incomeSumm %>% factor(levels=1:4),
           educationSumm = educationSumm %>% factor(levels=1:4),
           travel = travel %>% factor(levels=c(4,1,2,3)),
           sex = ifelse(sex == 1,"Female","Male"),
           age_cat_1 = cut(age,c(0,seq(5,45,by=5),110)))

  ## now deal with animal codes
  sample_meta_dat <- sample_meta_dat %>% rename(animal_chickens = animal1,
                                                animal_pigs = animal2,
                                                animal_asses = animal3,
                                                animal_cattle = animal4,
                                                animal_goats = animal5,
                                                animal_buffaloes = animal6,
                                                animal_horses = animal7,
                                                animal_rabbits = animal8,
                                                animal_sheep = animal9,
                                                animal_ducks = animal10,
                                                animal_geese = animal11,
                                                animal_pigeons = animal12,
                                                animal_other = animal13,
                                                animal_none = animal14)

  rc <- left_join(dat,sample_meta_dat) %>%
    mutate(adm2name = str_replace_all(Zila_name," Zila","")) %>% fix_adm2names()

  ## adding BTM coordinate data
  sampling_locs <- sf::st_as_sf(rc %>% filter(!is.na(lat)) %>% dplyr::select(lon,lat),coords = c("lon","lat"),crs="+proj=longlat +datum=WGS84 +no_defs")
  sampling_locs_btm <- sampling_locs %>% transform_to_btm %>% st_coordinates

  rc$x_btm <- NA
  rc$y_btm <- NA
  rc$x_btm[which(!is.na(rc$lat))] <- sampling_locs_btm[,"X"]
  rc$y_btm[which(!is.na(rc$lat))] <- sampling_locs_btm[,"Y"]

  return(rc)
}

##gets corrected seropositive
get_corrected_proportion_positive <- function(){
  return(NULL)
}

recode_risk_factors <- function(x) {

  factor_levels <-
    c('age_group_code1',
      'age_group_code3',
      'animal_chickensTRUE',
      'animal_pigsTRUE',
      'animal_cattleTRUE',
      'animal_goatsTRUE',
      'animal_rabbitsTRUE',
      'animal_ducksTRUE',
      'pigs_rabbitsTRUE',
      'non_pig_rabbit_animalsTRUE',
      'sexMale',
      'alt',
      'travel1',
      'travel2',
      'travel3',
      'elecTRUE',
      'ownerHHTRUE',
      'landTRUE',
      'urbanTRUE',
      'educationSumm2',
      'educationSumm3',
      'educationSumm1',
      'incomeSumm2',
      'incomeSumm3',
      'incomeSumm4',
      'HHSize',
      'logpop',
      'travel_time',
      'poverty',
      'distance_to_water')


  ## label lookup function for plotting and tables
  risk_factor_label_lookup <-c(
    `age_group_code1` = "0-4 years",
    `age_group_code3` = ">14 years",
     `sexMale` = "male",
    `animal_chickensTRUE` = "chickens",
    `animal_pigsTRUE` = "pigs",
    `animal_cattleTRUE` = "cattle",
    `animal_goatsTRUE` = "goats",
    `animal_rabbitsTRUE` = "rabbits",
    `animal_ducksTRUE` = "ducks",
    `pigs_rabbitsTRUE`="pigs or rabbits",
    `non_pig_rabbit_animalsTRUE` = "other animals",
    `alt10` = "altitude (x10m)",
    `alt` = "altitude (meters)",
    `altitude` = "altitude (meters)",
     `travel1` = "trav. last week",
     `travel2` = "trav. last month",
     `travel3` = "trav. last 6 months",
     `elecTRUE` = "electricity in house",
     `ownerHHTRUE` = "owns home",
     `landTRUE` = "owns land",
     `urbanTRUE` = "urban",
    `anyMosqContTRUE` = "mosquito control",
    `educationSumm1` = "post-secondary education (head hhl)",
    `educationSumm2` = "secondary school (head hhl)",
    `educationSumm3` = "primary school (head hhl)",
    `educationSumm4` = "no school (head hhl)",
    `incomeSumm1` = "<7,000TK",
    `incomeSumm2` = "7,000-9,999TK",
    `incomeSumm3` = "10,000-20,000TK",
    `incomeSumm4` = ">20,000TK",
    `logpop` = "population density (log)",
    `logtravel_time` = "travel time to nearest city (log)",
    `travel_time` = "travel time to nearest city (min)",
    `aedesAeg` = "Aedes aegypti mosquitos captured",
    `aedesAlbo` = "Aedes albopictus mosquitos captured",
    `poverty` = "poverty index",
    `distance_to_water` = "distance to major water body (per 10km)")


  ## factor(x,levels=factor_levels) %>%
  ##   recode(.,!!!risk_factor_label_lookup)

    factor(x) %>%
    recode(.,!!!risk_factor_label_lookup)

}


## makes a bootstrap dataset mimicing our sampling
make_boot_data <- function(dat){

  n_communities <- dat$community_id %>% unique %>% length

  ## resample 70 communities
  coms <- sample(dat$community_id,n_communities,replace = TRUE)

  ## in each community resample the same number of houeholds we actually sampled
  ## 10-12
  hhls <- sapply(coms,function(x) {
    hhl_list <- dat %>% filter(community_id == x) %>% dplyr::select(hhid) %>% unlist %>% unique
    sample(hhl_list,size = length(hhl_list),replace=TRUE)
  })

  ## for each household, resample people
  ind_dat <- sapply(hhls %>% unlist,function(hh){
    tmp_dat <- dat %>% filter(hhid == hh)
    tmp_dat %>% sample_n(size=nrow(tmp_dat),replace=TRUE)
  },simplify = FALSE) %>%
    bind_rows

  return(ind_dat)
}

## summary statistics for bootstraps
bs_stats <- function(bs_rep){

  sero_p <- mean(bs_rep$pos,na.rm=T)
  sero_p_by_age <- bs_rep %>% filter(!is.na(age_cat_1)) %>% group_by(age_cat_1) %>% summarize(sero_p_a = mean(pos,na.rm=T)) %>% dplyr::select(sero_p_a) %>% unlist
  sero_p_by_sex <- bs_rep %>% filter(!is.na(sex)) %>% group_by(sex) %>% summarize(sero_p_sex = mean(pos,na.rm=T)) %>% dplyr::select(sero_p_sex) %>% unlist
  sero_p_by_urban <- bs_rep %>% filter(!is.na(urban)) %>% group_by(urban) %>% summarize(sero_p_urban = mean(pos,na.rm=T)) %>% dplyr::select(sero_p_urban) %>% unlist


  sero_p_age_sex <- bs_rep %>%
    filter(!is.na(sex),!is.na(age_cat_1)) %>%
    group_by(age_cat_1,sex) %>%
    summarize(sero_p_sex_a = mean(pos,na.rm=T))



  labels <- sero_p_age_sex %>% mutate(lab = paste0(sex,"-",age_cat_1)) %>% ungroup %>% dplyr::select(lab) %>% unlist %>% unname

  age_sex_out <- sero_p_age_sex %>% ungroup %>% dplyr::select(sero_p_sex_a) %>% unlist %>% unname %>% setNames(labels)

  c(sero_p=sero_p,sero_p_by_age,sero_p_by_sex,age_sex_out,sero_p_by_urban)
}



#' Preps data for inla model
#' @param eqn 
#' @param est_dat 
#' @param loc_btm 
#' @param ran_eff 
#' @param pred_dat 
#'
#' @return
#' @export
#'
#' @examples
prep_inla_pred <- function(eqn="logpop + distance_to_water + poverty + travel_time + altitude",
                           est_dat,
                           loc_btm,
                           ran_eff,
                           pred_dat=my_grid){
  ## prepare the estimation data
  inla_prep_est <- prep_inla_est(eqn=eqn,
                                 dat=est_dat,
                                 loc_btm=loc_btm,
                                 ran_eff=ran_eff)
  ## find centroids for prediction data
  pred_cents <- pred_dat %>%
    st_centroid() %>%
    st_coordinates
  
  ## calculate a lattice projection from the estimated mesh to the
  ## predicted centroids
  pred_projector <- inla.mesh.projector(inla_prep_est$mesh,
                                        loc=pred_cents)
  ## get A matrix
  A_pred <- pred_projector$proj$A
  
  ## create prediction stack
  stack_pred <- inla.stack(
    data=list(y=NA,Ntrials=1),
    A=list(A_pred,1),
    effects=list(inla_prep_est$s_index,
                 data.frame(Intercept=rep(1,nrow(pred_dat)))),
    tag="pred")
  
  ## combine stacks
  est_and_pred_stack <- inla.stack(stack_pred,inla_prep_est$stack)
  
  ## get stack data
  est_and_pred_data <- inla.stack.data(est_and_pred_stack)
  
  return(list(est_prep=inla_prep_est,
              pred_projector=pred_projector,
              A_pred=A_pred,
              stack_pred=stack_pred,
              est_and_pred_stack=est_and_pred_stack,
              est_and_pred_data=est_and_pred_data))
}



#----------------------------------
# TAKEN from Ben Arnold: https://github.com/ben-arnold/enterics-seroepi/blob/master/R/SI-File6-Fig3-Fig4-agecurve-analyses.Rmd
# simulataneous CIs for GAMs
# estimated by resampling the
# Baysian posterior estimates of
# the variance-covariance matrix
# assuming that it is multivariate normal
# the function below also estimates
# the unconditional variance-covariance
# matrix, Vb=vcov(x,unconditional=TRUE),
# which allows for undertainty in the actual
# estimated mean as well
# (Marra & Wood 2012 Scandinavian Journal of Statistics,
#  Vol. 39: 53–74, 2012, doi: 10.1111/j.1467-9469.2011.00760.x )
# simultaneous CIs provide much better coverage than pointwise CIs
# see: http://www.fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
#----------------------------------
gamCI <- function(m,newdata,nreps=10000) {
  require(mgcv)
  require(dplyr)
  Vb <- vcov(m,unconditional = TRUE)
  pred <- predict(m, newdata, se.fit = TRUE)
  fit <- pred$fit
  se.fit <- pred$se.fit
  BUdiff <- MASS::mvrnorm(n=nreps, mu = rep(0, nrow(Vb)), Sigma = Vb)
  Cg <- predict(m, newdata, type = "lpmatrix")
  simDev <- Cg %*% t(BUdiff)
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  masd <- apply(absDev, 2L, max)
  crit <- quantile(masd, prob = 0.95, type = 8)
  pred <- data.frame(newdata,fit=pred$fit,se.fit=pred$se.fit)
  pred <- mutate(pred,
                 uprP = fit + (2 * se.fit),
                 lwrP = fit - (2 * se.fit),
                 uprS = fit + (crit * se.fit),
                 lwrS = fit - (crit * se.fit)
  )
  return(pred)
}

#' Preps for INLA
#'
#' @param eqn 
#' @param dat 
#' @param loc_btm 
#' @param ran_eff 
#'
#' @return
#' @export
#'
#' @examples
prep_inla_est <- function(eqn="logpop + sex + distance_to_water + poverty + travel_time + altitude",
                          dat,
                          loc_btm,
                          ran_eff){
  ## MESH
  ## Build mesh around household coordinates
  bnd <- inla.nonconvex.hull(loc_btm,convex=-0.1) ## inla mesh segment
  
  #meshbuilder() ## shiny app
  mesh <- inla.mesh.2d(boundary = bnd,
                       loc=loc_btm,
                       offset=c(-0.05, -0.05),
                       cutoff = 2000, ## if households are less than 2km apart, builds only a single vertex
                       max.edge=c(30000,50000))
  
  ## let's have a quick look at the mesh
  # plot(mesh,asp=1,main="com")
  # points(loc_btm[,"X"],loc_btm[,"Y"],pch=19,cex=0.5,col="orange")
  
  ## now make the SPDE model on the mesh
  spde <- inla.spde2.matern(mesh=mesh, alpha=2) # alpha is Fractional operator order
  s_index <- inla.spde.make.index(name="spatial.field",n.spde=spde$n.spde)
  
  ## observation/prediction matrix
  A_est_multi <- inla.spde.make.A(mesh=mesh,
                                  loc=loc_btm)
  
  ## create model matrix of data
  model_mtx <- model.matrix(as.formula(paste("y~",eqn)), dat)[,-1, drop=F]
  
  ## create stack for estimates
  effect_list <- list(s_index,
                      data.frame(Intercept=rep(1,nrow(model_mtx)),
                                 model_mtx))
  if(!missing(ran_eff)){
    for(i in 1:length(ran_eff)){
      effect_list[[length(effect_list)+1]] <-
        list(as.vector(unlist(dat[,ran_eff[i]])) %>% as.factor())
      names(effect_list[[length(effect_list)]]) <- ran_eff[[i]]
    }
  }
  A_list <- list(A_est_multi,1)
  if(!missing(ran_eff)){
    for(i in 1:length(ran_eff)){
      A_list[length(A_list)+1] <- 1
    }
  }
  est_stack = inla.stack(
    data=list(y=dat[,'y'] %>% unlist,
              Ntrials=dat[,'obs'] %>% unlist), ## response
    A=A_list, ## projector matrix
    effects=effect_list,
    tag="est")
  ## get stack data
  est_stack_data <- inla.stack.data(est_stack)
  
  return(list(bnd=bnd,
              mesh=mesh,
              spde=spde,
              s_index=s_index,
              A_est_multi=A_est_multi,
              model_mtx=model_mtx,
              stack=est_stack,
              stack_data=est_stack_data))
}




#' @param name_start character string, first letters in file name
#' @param path character string, path to folder of interest, end with "/"
#' @param exclude character string, patterns to exclude from the file names of interest
#'
#' @return character string, path to most recent file
#' @export
#'
#' @examples
find_recent_file <- function(name_start, path, exclude=NULL){
    if(substring(path, nchar(path))!="/")
        warning("Path does not end with a '/', problems may ensue.")
    ## view all files of that name at that path
    file_list <- list.files(path=path,
                            pattern=paste0(name_start, "*"))
    ## remove files with unwanted patterns
    if(!is.null(exclude)){
        for(i in 1:length(exclude))
            file_list <- file_list[!grepl(pattern = exclude[i], file_list)]
    }
    if(length(file_list)==0){
        warning('File not found')
        return(NA)
    }
    ## view file info
    file_info <- file.info(paste0(path, file_list))
    ## find most recent file
    most_recent_file <- paste0(path,
                               file_list[which.max(file_info$mtime)])

    cat(sprintf("Loaded file: \n %s last updated on \n %s \n",most_recent_file,file_info$mtime[which.max(file_info$mtime)]))

    return(most_recent_file)
}

## make pretty CI
make_pretty_ci <- function(x,multiplier=100){
  sprintf("%.2f%% (95%%CI %.2f-%.2f)",round(multiplier*x[1],2),round(multiplier*x[2],2),round(multiplier*x[3],2))
}
