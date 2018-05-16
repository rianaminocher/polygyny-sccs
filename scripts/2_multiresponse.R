
### Fit multi-response model to 100 imputed datasets ###

rm(list=ls())
setwd("/Volumes/riana_minocher/sccs-polygyny") # direct to the github repo

library(ape)
library(MuMIn)
library(MCMCglmm)
library(mice)
library(rethinking)
library(phytools)
library(phangorn)
library(wesanderson)
library(mapdata)


        ### load data ###

  load("data/Cross-cultural data.RData")
  updated_SCCS <- read.csv("data/SCCS-var1-2000.csv")
  updated_SCCS <- updated_SCCS[-187,] # remove extra "NA" row on bottom
  rownames(updated_SCCS) <- updated_SCCS$socname
  
  # reduce dataset to variables we use and impute missing values
  data.multi <- SCCS[ , c("v79", "v860", "v871", "v872", 
                          "v1260", "v740", "v64", "v53", 
                          "v186", "v158", "v2021", "v1649", 
                          "v1666", "v890", "v714")]
  data.multi$v186 <- updated_SCCS$v186 
  
  # any rows with no outcome data?
  data.multi[is.na(data.multi$v860) & is.na(data.multi$v871), ] # Ajie, Gond
  # --> drop from dataset
  data.multi <- data.multi[rownames(data.multi)!="Ajie" & rownames(data.multi)!="Gond", ]
  
       ### impute missing data ###

  set.seed(1)
  imputed_dfs_multi <- mice(data.multi, m = 100)

  # extract dfs
  list_nos <- seq(1:100)
  fun_extract_dfs <- function(x) {
    complete(imputed_dfs_multi, action = x) 
  }

list_dfs_multi <- lapply(list_nos, fun_extract_dfs)



        ### load phylogeny ###
  
  # get socnames and phylo
  socnamesphylo <- read.csv("data/socnamesphylo.csv")
  phylonames <- socnamesphylo$phylo
  
  # read tree
  phylo <- read.nexus("data/Time-calibrated SCCS supertree ultrametric.nex")
  is.ultrametric(phylo, option = 2)
  inv.phylo <- inverseA(phylo, nodes = "TIPS", scale = TRUE)
  (max(vcv(phylo)) - min(vcv(phylo)) / max(vcv(phylo))) # 140

  # compute the NNLS ultrametric tree
  nnls <- nnls.tree(cophenetic(phylo), phylo, rooted = TRUE)
  is.ultrametric(nnls)
  # compare distances
  tips <- phylo$tip.label
  plot(cophenetic(phylo)[tips, tips], cophenetic(nnls)[tips, tips],
       xlab="original distances",ylab="distances on NNLS tree")
  cor(as.vector(cophenetic(phylo)[tips,tips]),
      as.vector(cophenetic(nnls)[tips,tips]))
  # 1
  
  inv.phylo <- inverseA(nnls, nodes = "TIPS", scale = TRUE)	
  
  # check names
  setdiff(phylonames, phylo$tip.label)
  setdiff(phylo$tip.label, phylonames)
  phylonames<- as.character(phylonames)
  # change to match
  phylonames[phylonames=="Kazakh"]<- "Kazak"
  phylonames[phylonames=="MISSING_Tobelorese"]<- "Tobelorese"
  phylonames[phylonames=="Twana"]<- "MISSING_Twana"
  phylonames[phylonames=="'Paiute_North.'"]<- "Paiute_North"
  phylonames[phylonames=="Omaha"]<- "MISSING_Omaha"
  phylonames[phylonames=="Huron"]<- "MISSING_Huron"
  phylonames[phylonames=="Nambicuara"]<- "MISSING_Nambicuara"
  phylonames[phylonames=="Botocudo"]<- "MISSING_Botocudo"




        ### fit multi-response model ###

# % men (v871) and % women (v872) very highly correlated (r=0.97, P<0.0001) 
# --> use v871 (fewer NA's)
# % men (v871), cultural rules (v860) as response 

  
  ### all predictors, no phylo ###
  
  fun_full_model_nophy <- function(df) {

  # reorganise variables 
  
    # v890 (female contribution to agriculture) from the subsistence block
  df$v890 <- as.factor(df$v890)
  levels(df$v890) <- list("2"="0", "9.5"="1", "19.5"="2", "29.5"="3",
                          "39.5"="4", "49.5"="5", "59.5"="6", "69.5"="7",
                          "79.5"="8")
  df$v890 <- as.numeric(as.character(df$v890)) - mean(as.numeric(as.character(df$v890)))
  
    # v1649 (internal warfare) and v1666 (assault frequency) from aggressiveness block
  df$v1649 <- df$v1649 - mean(df$v1649)
  df$v1666 <- df$v1666 - mean(df$v1666)
  
    # v158 (stratification) and v2021 (wealth inequalities) from stratification block
  df$v158 <- as.factor(df$v158)
  levels(df$v158) <- list("1"= c("1","2"), 
                          "2"=  c("3","4"), 
                          "3" = "5")
  df$v2021 <- as.factor(df$v2021)
  levels(df$v2021) <- list("1"= c("1","1.5"), 
                           "2"=c("2","2.5"),
                           "3"="3")
  
    # v186 (temperature) from climate block
  df$v186 <- scale(df$v186)
  
    # v53 (father-inf proximity)
  df$v53 <- df$v53 - mean(df$v53)
  
    # v64 (population density)
  df$v64 <- as.factor(df$v64)
  levels(df$v64) <- list("0.2"="1", "0.4"="2", "3"="3", 
                         "13"="4", "63"="5", "300.5"="6", "500"="7")
  df$v64 <- as.numeric(as.character(df$v64)) - mean(as.numeric(as.character(df$v64)))
  
    # v740 (female arranged marriages) - change to 2 level
  df$v740 <- as.factor(df$v740)
  levels(df$v740) <- list("1" = c("1","2","3","4","5"), 
                          "2"="6")
  
    # v1260 (pathogen stress)
  df$v1260 <- scale(log(df$v1260))
  
  # - v714 sex ratio
  ## center on even sex ratio
  df$v714<- df$v714-2
  
  # polygyny = stress v1260 + femarr v740 + density v64 + proxi v53
  # + temp v186 + strat v158 + ineq v2021 + intwar 1649 + assfreq v1666 + fcontmean v890 + sexratio v714
  
  prior.multi <- list(B=list(mu=c(rep(0,28)),  V=diag(28)*5), R=list(V=diag(2),nu=0.02,fix=2)) 
  # weakly informative prior
  fullmodel.multi <- MCMCglmm(cbind(v871, v860) ~ trait -1 + trait:v1260 + trait:v740 + trait:v64 + trait:v53 + trait:v186 + trait:v158 + trait:v2021 + trait:v1649 + trait:v1666 + trait:v890 + trait:v714,
                          family=c("poisson","ordinal"), rcov=~us(trait):units, 
                          prior=prior.multi, data=df, nitt=5e5, burnin=1e4, thin=1e2)
  
  #plot(fullmodel.multi) # all good
  #summary(fullmodel.multi) # all eff. samp >500
  
  #fullmodel.multib <- MCMCglmm(cbind(v871, v860) ~ trait -1 + trait:v1260 + trait:v740 + trait:v64 + trait:v53 + trait:v186 + trait:v158 + trait:v2021 + trait:v1649 + trait:v1666 + trait:v890 + trait:v714,
  #                        family=c("poisson","ordinal"), rcov=~us(trait):units, 
  #                        prior=prior.multi, data=df, nitt=5e5, burnin=1e4, thin=1e2)
  #fullmodel.multic <- MCMCglmm(cbind(v871, v860) ~ trait -1 + trait:v1260 + trait:v740 + trait:v64 + trait:v53 + trait:v186 + trait:v158 + trait:v2021 + trait:v1649 + trait:v1666 + trait:v890 + trait:v714,
  #                       family=c("poisson","ordinal"), rcov=~us(trait):units, 
  #                        prior=prior.multi, data=df, nitt=5e5, burnin=1e4, thin=1e2)
  #gelman.diag(list(fullmodel.multi$Sol, fullmodel.multib$Sol, fullmodel.multic$Sol)) ## all 1 (CI up to 1.01)
  #gelman.diag(list(fullmodel.multi$VCV[,1:2], fullmodel.multib$VCV[,1:2], fullmodel.multic$VCV[,1:2]), transform=TRUE) # 1
  }  
   
  
  # fit nophy model to 100 dfs
  
  list_mods_nophy <- lapply(list_dfs_multi[1:100], fun_full_model_nophy)

  
  # summarize no phy models

  combined_mod_nophy <- list_mods_nophy[[1]]
  for(i in 2:100) { 
    combined_mod_nophy$Sol <- rbind(combined_mod_nophy$Sol, list_mods_nophy[[i]]$Sol)
    combined_mod_nophy$VCV <- rbind(combined_mod_nophy$VCV, list_mods_nophy[[i]]$VCV)
    combined_mod_nophy$CP <- rbind(combined_mod_nophy$CP, list_mods_nophy[[i]]$CP)
  }
  
  
  # summarize fixed effects
  combined_mod_nophy_fixed_means <- colMeans(combined_mod_nophy$Sol[ , 1:28])
  combined_mod_nophy_fixed_hpdi <- cbind(rep(NA, 28), rep(NA, 28)) 
  for (i in 1:28) {
    combined_mod_nophy_fixed_hpdi[i, 1] <- HPDI(combined_mod_nophy$Sol[ ,i], prob = 0.95)[1]
    combined_mod_nophy_fixed_hpdi[i,2] <- HPDI(combined_mod_nophy$Sol[,i], prob=0.95)[2]
  }
  
    combined_mod_nophy_fixed_hpdi
    round(cbind(combined_mod_nophy_fixed_means, combined_mod_nophy_fixed_hpdi),3)
    combined_mod_nophy_fixed_means

    
  # proportion of samples on the same side of 0 as the mean
  
  p.combined_mod_nophy <- data.frame(p = rep(NA, 28))
  rownames(p.combined_mod_nophy) <- names(combined_mod_nophy_fixed_means)
  for (i in 1:nrow(p.combined_mod_nophy)) {
    if (sum(combined_mod_nophy$Sol[ ,i] > 0) / length(combined_mod_nophy$Sol[ ,i]) > 0.5) 
      p.combined_mod_nophy[i, ] <- sum(combined_mod_nophy$Sol[ ,i] > 0) / length(combined_mod_nophy$Sol[ ,i])
    else p.combined_mod_nophy[i, ] <- 1 - sum(combined_mod_nophy$Sol[ ,i] > 0) / length(combined_mod_nophy$Sol[,i])
  }
  
  
  # variance components
  
  combined_mod_nophy_random_means <- colMeans(combined_mod_nophy$VCV[ ,1:4])
  combined_mod_nophy_random_hpdi <- cbind(rep(NA, 4), rep(NA, 4)) 
  for (i in 1:4) {
    combined_mod_nophy_random_hpdi[i, 1] <- HPDI(combined_mod_nophy$VCV[ ,i], prob = 0.95)[1]
    combined_mod_nophy_random_hpdi[i, 2] <- HPDI(combined_mod_nophy$VCV[ ,i], prob = 0.95)[2]
  }
  
  round(cbind(combined_mod_nophy_random_means, combined_mod_nophy_random_hpdi), 3)

  
  # cut points
  
  combined_mod_nophy_cp_means <- colMeans(combined_mod_nophy$CP[ ,1:3])
  combined_mod_nophy_cp_hpdi <- cbind(rep(NA, 3), rep(NA, 3)) 
  for (i in 1:3) {
    combined_mod_nophy_cp_hpdi[i, 1] <- HPDI(combined_mod_nophy$CP[,i], prob = 0.95)[1]
    combined_mod_nophy_cp_hpdi[i, 2] <- HPDI(combined_mod_nophy$CP[,i], prob = 0.95)[2]
  }
  
  round(cbind(combined_mod_nophy_cp_means, combined_mod_nophy_cp_hpdi), 3)

  


   ### all predictors, phylo ###
  
fun_full_model_phy <- function(df) {
  
  # add phylo
  
  df$phylo <- phylonames[phylonames!="Ajie" & phylonames!="Gond"]
  
  # reorganise variables
  
    # v890
  df$v890 <- as.factor(df$v890)
  levels(df$v890) <- list("2"="0","9.5"="1","19.5"="2","29.5"="3",
                          "39.5"="4", "49.5"="5", "59.5"="6", "69.5"="7",
                          "79.5"="8")
  df$v890 = as.numeric(as.character(df$v890))-mean(as.numeric(as.character(df$v890)))
  
    # v1649
  df$v1649<- df$v1649 - mean(df$v1649)
  df$v1666<- df$v1666 - mean(df$v1666)
  
    # v158 and v2021
  df$v158 <- as.factor(df$v158)
  levels(df$v158) <- list("1"=c("1","2"), 
                          "2"=c("3","4"), 
                          "3"="5")
  df$v2021 <- as.factor(df$v2021)
  levels(df$v2021) <- list("1"=c("1","1.5"), 
                           "2"=c("2","2.5"), 
                           "3"="3")
  
    # v186 
  df$v186 <- scale(df$v186)
  
    # v53
  df$v53 <- df$v53 - mean(df$v53)
  
    # v64
  df$v64 <- as.factor(df$v64)
  levels(df$v64) <- list("0.2"="1", "0.4"="2", "3"="3", 
                         "13"="4", "63"="5", "300.5"="6", "500"="7")
  df$v64 <- as.numeric(as.character(df$v64)) - mean(as.numeric(as.character(df$v64)))
  
    # v740
  df$v740 <- as.factor(df$v740)
  levels(df$v740) <- list("1"=c("1","2","3","4","5"), 
                          "2"="6")
  
    # v1260
  ### log transform, SCALE
  df$v1260 <- scale(log(df$v1260))
  
    # - v714 sex ratio
  ## center on even sex ratio
  df$v714<- df$v714-2
  
  # polygyny = stress v1260 + femarr v740 + density v64 + proxi v53
  # + temp v186 + strat v158 + ineq v2021 + intwar 1649 + assfreq v1666 + fcontmean v890 + sexratio v714
  
  prior.multi.phy <- list(B=list(mu=c(rep(0,28)),  V=diag(28)*5), R=list(V=diag(2),nu=2,fix=2), G=list(G1=list(V=diag(2),nu=2)))
  # weakly informative prior
  fullmodel.multi.phy <- MCMCglmm(cbind(v871, v860) ~ trait -1 + trait:v1260 + trait:v740 + trait:v64 + trait:v53 + trait:v186 + trait:v158 + trait:v2021 + trait:v1649 + trait:v1666 + trait:v890 + trait:v714,
                          random=~us(trait):phylo, ginverse=list(phylo=inv.phylo$Ainv), 
						  family=c("poisson","ordinal"), rcov=~us(trait):units, 
                          prior=prior.multi.phy, data=df, nitt=100000, burnin=10000, thin=10)
  ## model checks
  #plot(fullmodel.multi.phy) # 
  #summary(fullmodel.multi.phy) # 
  
  #fullmodel.multi.phyb <- MCMCglmm(cbind(v871, v860) ~ trait -1 + trait:v1260 + trait:v740 + trait:v64 + trait:v53 + trait:v186 + trait:v158 + trait:v2021 + trait:v1649 + trait:v1666 + trait:v890 + trait:v714,
  #                       random=~us(trait):phylo, ginverse=list(phylo=inv.phylo$Ainv), 
	#					  family=c("poisson","ordinal"), rcov=~us(trait):units, 
  #                       prior=prior.multi.phy, data=df, nitt=100000, burnin=10000, thin=10)
  #fullmodel.multi.phyc <- MCMCglmm(cbind(v871, v860) ~ trait -1 + trait:v1260 + trait:v740 + trait:v64 + trait:v53 + trait:v186 + trait:v158 + trait:v2021 + trait:v1649 + trait:v1666 + trait:v890 + trait:v714,
  #                        random=~us(trait):phylo, ginverse=list(phylo=inv.phylo$Ainv), 
	#					  family=c("poisson","ordinal"), rcov=~us(trait):units, 
  #                        prior=prior.multi.phy, data=df, nitt=100000, burnin=10000, thin=10)
  #gelman.diag(list(fullmodel.multi.phy$Sol, fullmodel.multi.phyb$Sol, fullmodel.multi.phyc$Sol)) ## all 1 (CI up to 1.01)
  #gelman.diag(list(fullmodel.multi.phy$VCV[,1:2], fullmodel.multi.phyb$VCV[,1:2], fullmodel.multi.phyc$VCV[,1:2]), transform=TRUE) ## all 1 (CI up to 1.01)
  }


# list_mods_phy <- lapply(list_dfs_multi[1:100], fun_full_model_phy) 
# might not work if at least one model has singular equations --> need to do one by one, and restart if chains get stuck

  list_mods_phy <- list_mods_nophy
  # replace no phy models one by one
  # test: summary(list_mods_phy[[1]]); list_mods_phy[[1]]<- fun_full_model_phy(list_dfs_multi[[1]]); summary(list_mods_phy[[1]])
  # save as it happens in case R crashes
  # load object if R crashed: load("list_mods_phy.robj")

  for(i in 1:100){
  list_mods_phy[[i]] <- fun_full_model_phy(list_dfs_multi[[i]])
  save(list_mods_phy, file="list_mods_phy.robj")
  }


  # summarize phy models
  
  
  combined_mod_phy <- list_mods_phy[[1]]
  for(i in 2:100) {
    combined_mod_phy$Sol <- rbind(combined_mod_phy$Sol, list_mods_phy[[i]]$Sol)
    combined_mod_phy$VCV <- rbind(combined_mod_phy$VCV, list_mods_phy[[i]]$VCV)
    combined_mod_phy$CP <- rbind(combined_mod_phy$CP, list_mods_phy[[i]]$CP)
  }
  
    
  # summarize fixed effects
    
  combined_mod_phy_fixed_means <- colMeans(combined_mod_phy$Sol[ , 1:28])
  combined_mod_phy_fixed_hpdi <- cbind(rep(NA, 28), rep(NA, 28)) 
  for (i in 1:28) {
    combined_mod_phy_fixed_hpdi[i, 1] <- HPDI(combined_mod_phy$Sol[ , i], prob = 0.95)[1]
    combined_mod_phy_fixed_hpdi[i, 2] <- HPDI(combined_mod_phy$Sol[ , i], prob = 0.95)[2]
  } 

  round(cbind(combined_mod_phy_fixed_means, combined_mod_phy_fixed_hpdi), 3)

  
  # proportion of samples on the same side of 0 as the mean
 
  p.combined_mod_phy <- data.frame(p = rep(NA, 28))
  rownames(p.combined_mod_phy) <- names(combined_mod_phy_fixed_means)
  for (i in 1:nrow(p.combined_mod_phy)) {
    if 
    (sum(combined_mod_phy$Sol[ , i] > 0) / length(combined_mod_phy$Sol[ , i]) > 0.5) 
      p.combined_mod_phy[i, ] <- sum(combined_mod_phy$Sol[ , i] > 0) / length(combined_mod_phy$Sol[ , i]) 
    else 
    p.combined_mod_phy[i, ] <- 1-sum(combined_mod_phy$Sol[ , i] > 0) / length(combined_mod_phy$Sol[ , i])
  }
  
  round(p.combined_mod_phy, 2)

  
  # variance components
  
  combined_mod_phy_random_means <- colMeans(combined_mod_phy$VCV[ ,1:8])
  combined_mod_phy_random_hpdi <- cbind(rep(NA, 8), rep(NA, 8)) 
  for (i in 1:8) {
    combined_mod_phy_random_hpdi[i, 1] <- HPDI(combined_mod_phy$VCV[ , i], prob = 0.95)[1]
    combined_mod_phy_random_hpdi[i, 2] <- HPDI(combined_mod_phy$VCV[ , i], prob = 0.95)[2]
  }
  
  round(cbind(combined_mod_phy_random_means, combined_mod_phy_random_hpdi), 3)

  
  # cut points
  
  combined_mod_phy_cp_means <- colMeans(combined_mod_phy$CP[ , 1:3])
  combined_mod_phy_cp_hpdi <- cbind(rep(NA, 3), rep(NA, 3)) 
  for (i in 1:3) {
    combined_mod_phy_cp_hpdi[i, 1] <- HPDI(combined_mod_phy$CP[ , i], prob = 0.95)[1]
    combined_mod_phy_cp_hpdi[i, 2] <- HPDI(combined_mod_phy$CP[ , i], prob = 0.95)[2]
  }
  
  round(cbind(combined_mod_phy_cp_means, combined_mod_phy_cp_hpdi), 3)

  
  
  
  ## Table means, HDPI, PP
  
  tab_fixed <- 
    data.frame(
    round(cbind(combined_mod_nophy_fixed_means, combined_mod_nophy_fixed_hpdi), 3),
    p.combined_mod_nophy,
    round(cbind(combined_mod_phy_fixed_means, combined_mod_phy_fixed_hpdi), 3),
    p.combined_mod_phy
  )
  colnames(tab_fixed) <- c("Mean", "95% HDPI", "", "PP", "Mean", "95% HDPI", "", "PP")
  tab_fixed <- tab_fixed[order(rownames(tab_fixed), decreasing = TRUE), ]
  write.csv(tab_fixed, "table2.csv")
  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          