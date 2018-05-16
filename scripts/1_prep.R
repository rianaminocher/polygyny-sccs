
### Prep & variable selection ###

rm(list=ls())
setwd("/Volumes/riana_minocher/sccs-polygyny/sccs-polygyny") # direct to the github repo

library(foreign)
library(car)

  ### load SCCS workspace ###

  load("data/Cross-cultural data.RData")

  ### load updated SCCS vars ###
  
  updated_SCCS <- read.csv("data/SCCS-var1-2000.csv")
  updated_SCCS <- updated_SCCS[-187,] # remove extra "NA" row on bottom
  rownames(updated_SCCS) <- updated_SCCS$socname
  
  ### load Low's marriage code and Ember's non-sororal polygyny, mortality & plow ###
  
  lowmarrcode <- read.csv("data/lowmarrcode.csv")
  malemort <- as.data.frame(read.spss("data/male_mort_sh.sav"))
  malemort <- malemort[!is.na(malemort$id), ] # remove extra rows
  identical(malemort$name, SCCS$socname) # no, but names identical and SCCS id is correct (1:186)
  rownames(malemort) <- rownames(SCCS) # can assign rownames to match SCCS
  
  ### variables previously used in polygyny analyses ###
 
  all_var_list <- c("v11", "v1260", "v158", "v1649", "v1666", 
                    "v179", "v180", "v186", "v200", "v2021", 
                    "v203", "v204", "v205", "v2176", "v53", 
                    "v586", "v625", "v64", "v7", "v714", "v740", "v79", 
                    "v833", "v857", "v858", "v860", "v871", 
                    "v885", "v890", "v9")
  
  all.equal(SCCS$socname, as.character(updated_SCCS$socname)) 
  
  # check if there are any diffs with updated codebook
 
  names(setdiff(updated_SCCS[ , colnames(updated_SCCS) %in% all_var_list],
                SCCS[ , colnames(SCCS) %in% all_var_list]))
  # 4 variables are different v186 temp, v860 cultural rules, v1649 intwar, v1666 assault freq
  
  ### combine 'updated' and original SCCS vars in df "polygyny_SCCS" ###
  
  polygyny_SCCS <- updated_SCCS[ , colnames(updated_SCCS) %in% all_var_list] # take updated SCCS vars
  
  # v2021, v586, v64, v79 are not in updated SCCS list
  polygyny_SCCS$v2021 <- SCCS[ , "v2021"]
  polygyny_SCCS$v586 <- SCCS[ , "v586"] 
  polygyny_SCCS$v64 <- SCCS[ , "v64"] 
  polygyny_SCCS$v79 <- SCCS[ , "v79"] 
  
  # use original SCCS v1649, v1666 - updated list has '88' vals where NA was
  polygyny_SCCS$v1649 <- SCCS[ , "v1649"] 
  polygyny_SCCS$v1666 <- SCCS[ , "v1666"]
  
  # use original v860, updated has 4 more NAs
  polygyny_SCCS$v860 <- SCCS[ , "v860"] 
  
  # Low and Ember vars not in SCCS
  polygyny_SCCS$lowcode <- lowmarrcode$BSLmarrCode
  polygyny_SCCS$nsp <- malemort$gen_ns_poly
  polygyny_SCCS$mort <- malemort$malemortnew
  polygyny_SCCS$plow <- malemort$pres_plow
  
  
  ### Consider variables in categories ###
  
        # Subsistence

  # remove hunt, fish, gather non EA, female contribution 1 (we use mean), subsistence mode dominant classification
  polygyny_SCCS <- polygyny_SCCS[, !colnames(polygyny_SCCS) %in% c("v7","v9","v11","v885", "v833")] 
  
  ## Time in subsistence
  levels(SCCSfact$v586)
  
  ## Mode of subsistence
  polygyny_SCCS$v858 <- as.factor(polygyny_SCCS$v858)
  levels(SCCSfact$v858) # 1:11
  levels(polygyny_SCCS$v858) <- list("forager"=c("1","2","3"), 
                                     "complex forager"=c("4","5"),
                                     "pastoralist"="6", 
                                     "horticulturalist"=c("7","8","9"),
                                     "intensive agriculture"=c("11","12"))
  # the SCCS variable skips num "10"
  
  ## Female contribution to agriculture
  polygyny_SCCS$v890 <- as.factor(polygyny_SCCS$v890)
  levels(SCCSfact$v890)
  levels(polygyny_SCCS$v890) <- list("2"="0","9.5"="1","19.5"="2","29.5"="3",
                                     "39.5"="4", "49.5"="5", "59.5"="6", "69.5"="7",
                                     "79.5"="8")
  polygyny_SCCS$v890 <- as.numeric(as.character(polygyny_SCCS$v890))
  
  ## Dependence on hunting reassign
  polygyny_SCCS$v204 <- as.factor(polygyny_SCCS$v204)
  levels(SCCSfact$v204)
  levels(polygyny_SCCS$v204) <- list("2.5"="0","10.5"="1","20.5"="2","30.5"="3",
                                     "40.5"="4", "50.5"="5", "60.5"="6", "70.5"="7",
                                     "80.5"="8","93"="9")
  polygyny_SCCS$v204 <- as.numeric(as.character(polygyny_SCCS$v204))
  polygyny_SCCS$v204 <- log(polygyny_SCCS$v204)
  
  ## Dependence on fishing reassign
  polygyny_SCCS$v205 <- as.factor(polygyny_SCCS$v205)
  levels(SCCSfact$v205)
  levels(polygyny_SCCS$v205) <- list("2.5"="0","10.5"="1","20.5"="2","30.5"="3",
                                     "40.5"="4", "50.5"="5", "60.5"="6", "70.5"="7",
                                     "80.5"="8","93"="9")
  polygyny_SCCS$v205 <- as.numeric(as.character(polygyny_SCCS$v205))
  polygyny_SCCS$v205 <- log(polygyny_SCCS$v205)
  
  ## Dependence on gathering reassign, log, scale
  polygyny_SCCS$v203 <- as.factor(polygyny_SCCS$v203)
  levels(SCCSfact$v203)
  levels(polygyny_SCCS$v203) <- list("2.5"="0","10.5"="1","20.5"="2","30.5"="3",
                                     "40.5"="4", "50.5"="5", "60.5"="6", "70.5"="7",
                                     "80.5"="8","93"="9")
  polygyny_SCCS$v203 <- as.numeric(as.character(polygyny_SCCS$v203))
  polygyny_SCCS$v203 <- log(polygyny_SCCS$v203)
  
  ## Presence of plow
  polygyny_SCCS$plow
  
  vif(lm(v871 ~ v586 + v858 + v890
         + v204 + v205 + v203 + plow, polygyny_SCCS))
  # removing 858 mode of subsistence
  vif(lm(polygyny_SCCS$v871 ~ polygyny_SCCS$v586 + polygyny_SCCS$v890
         + polygyny_SCCS$v204 + polygyny_SCCS$v205 + polygyny_SCCS$v203, polygyny_SCCS))
  # or remove fish, hunt
  vif(lm(v871 ~ v586 + v858 + v890
         + v205, polygyny_SCCS))
  
  prcomp( ~ polygyny_SCCS$v586 + polygyny_SCCS$v890
          + polygyny_SCCS$v204 + polygyny_SCCS$v205 + polygyny_SCCS$v203 + plow, polygyny_SCCS)
  summary(prcomp( ~ polygyny_SCCS$v586 + polygyny_SCCS$v890
                  + polygyny_SCCS$v204 + polygyny_SCCS$v205 + polygyny_SCCS$v203 + plow, polygyny_SCCS))
  plot(prcomp( ~ polygyny_SCCS$v586 + polygyny_SCCS$v890
               + polygyny_SCCS$v204 + polygyny_SCCS$v205 + polygyny_SCCS$v203 + plow, polygyny_SCCS))
  biplot(prcomp( ~ polygyny_SCCS$v586 + polygyny_SCCS$v890
                 + polygyny_SCCS$v204 + polygyny_SCCS$v205 + polygyny_SCCS$v203 + plow, polygyny_SCCS))
  
  
  # Aggressiveness
  
  ## Internal warfare
  
  ## Assault frequency 
  
  ## Value placed on aggressiveness
  levels(SCCSfact$v625)
  
  ## Male mortality
  range(polygyny_SCCS$mort, na.rm=T)
  
  vif(lm(v871 ~ v1649 + v1666 + v625 + mort, polygyny_SCCS))
  # remove mort
  vif(lm(v871 ~ v1649 + v1666 + v625, polygyny_SCCS))
  
  prcomp(~ v1649 + v1666 + v625, polygyny_SCCS)
  summary(prcomp(~ v1649 + v1666 + v625, polygyny_SCCS))
  plot(prcomp(~ v1649 + v1666 + v625, polygyny_SCCS))
  biplot(prcomp(~ v1649 + v1666 + v625, polygyny_SCCS))
  
  # Stratification
  
  ## Social stratification
  levels(SCCSfact$v158)
  polygyny_SCCS$v158 <- as.factor(polygyny_SCCS$v158)
  levels(polygyny_SCCS$v158) <- list("1"=c("1","2"), "2"=c("3","4"), "3"="5")
  polygyny_SCCS$v158 <- as.numeric(as.character(polygyny_SCCS$v158))
  
  ## Distribution of wealth
  levels(SCCSfact$v2021)
  polygyny_SCCS$v2021 <- as.factor(polygyny_SCCS$v2021)
  levels(polygyny_SCCS$v2021) <- list("1"=c("1","1.5"), "2"=c("2","2.5"), "3"="3")
  polygyny_SCCS$v2021 <- as.numeric(as.character(polygyny_SCCS$v2021))
  # 1 = general equality + minor inequality
  # 2 = some differences in wealth
  # 3 = considerable differences in wealth
  
  vif(lm(v871 ~ v158 + v2021, polygyny_SCCS))
  # strat vars are ok
  prcomp(~ v158 + v2021, polygyny_SCCS)
  summary(prcomp(~ v158 + v2021, polygyny_SCCS))
  plot(prcomp(~ v158 + v2021, polygyny_SCCS))
  biplot(prcomp(~ v158 + v2021, polygyny_SCCS))
  
  # Climate
  
  ## Latitude
  polygyny_SCCS$v179 <- log(polygyny_SCCS$v179)
  
  ## Hemisphere factor (binom)
  levels(SCCSfact$v180)
  polygyny_SCCS$v180 <- as.factor(polygyny_SCCS$v180)
  
  ## Temperature
  
  ## Region
  levels(SCCSfact$v200)
  polygyny_SCCS$v200 <- as.factor(polygyny_SCCS$v200) # convert to categorical, ref Africa
  
  ## Climate type 
  polygyny_SCCS$v857 <- as.factor(polygyny_SCCS$v857)
  levels(SCCSfact$v857)
  
  vif(lm(v871 ~ v186 + v200 + v857 + v179 + v180, polygyny_SCCS))
  # region and climate type are collinear
  vif(lm(v871 ~ v186 + v857 + v179 + v180, polygyny_SCCS))
  # climate type > 3
  vif(lm(v871 ~ v186 + v200 + v179 + v180, polygyny_SCCS))
  # region > 3
  
  prcomp(~ v186 + v179 + as.numeric(v180), polygyny_SCCS)
  summary(prcomp(~ v186 + v179 + as.numeric(v180), polygyny_SCCS))
  plot(prcomp(~ v186 + v179 + as.numeric(v180), polygyny_SCCS))
  biplot(prcomp(~ v186 + v179 + as.numeric(v180), polygyny_SCCS))
  
  # Father-inf proximity
  levels(SCCSfact$v53)
  
  # Population density
  polygyny_SCCS$v64 <- as.factor(polygyny_SCCS$v64)
  levels(SCCSfact$v64)
  levels(polygyny_SCCS$v64) <- list("0.2"="1", "0.4"="2", "3"="3", 
                                    "13"="4", "63"="5", "300.5"="6", "500"="7")
  polygyny_SCCS$v64 <- as.numeric(as.character(polygyny_SCCS$v64))
  # 1 = 0.2 people per sq mile
  # 2 = 0.4 people per sq mile
  # 3 = 3 people per sq mile
  # 4 = 13 people per sq mile
  # 5 = 63 people per sq mile
  # 6 = 300.5 people per sq mile
  # 7 = 500 people per sq mile
  
  # Female arranged marriages reassign, factor (binom)
  levels(SCCSfact$v740)
  polygyny_SCCS$v740 <- as.factor(polygyny_SCCS$v740)
  levels(polygyny_SCCS$v740) <- list("1"=c("1","2","3","4","5"), "2"="6")
  
  # Pathogen stress
  
    vif(lm(v871 ~ v1260 + v53 + v64
         + v740 + v714, polygyny_SCCS))
	# all <3
  
  # summarise variables used
  
  data_summary <- polygyny_SCCS[ , c("v79", "v860", "v871", 
                                     "v1260", "v740", "v64", "v53", 
                                     "v186", "v158", "v2021", "v1649", 
                                     "v1666", "v890", "v714")]
  str(data_summary)
  data_summary$v860 <- as.factor(data_summary$v860)
  tab_summarise_vars <- data.frame(
    colnames(data_summary),
    round(apply(data_summary, 2, sd, na.rm = TRUE),2),
    paste(apply(data_summary, 2, range, na.rm = TRUE)[1, ],
          apply(data_summary, 2, range, na.rm = TRUE)[2, ],
          sep = "-"),
    apply(data_summary, 2, function(x) length(which(is.na(x) == TRUE)))
  )
  colnames(tab_summarise_vars) <- c("Variable", "SD", "Range", "# missing")
  write.csv(tab_summarise_vars, "table1.csv")
  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             