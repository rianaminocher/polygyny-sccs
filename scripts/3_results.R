
### Extract results from fitted models ###

# load workspace 
load("data/Adrians SCCS workspace w sex ratio.RData")

library(phytools)
library(rethinking)
library(wesanderson)
library(mapdata)
        


        ### Figure Posterior distributions ###

        # compare degree of support for predictors 
        # plot posterior distribution of each parameter, for phy and nophy

  # plotting function
  plot_post_dist <- function(sol, var_names) {
    dens(combined_mod_nophy$Sol[ , sol], 
         show.zero = TRUE, 
         col = "white", 
         xlab = "", 
         ylab = "", 
         main = var_names)
    polygon(density(combined_mod_nophy$Sol[ , sol]), 
            col = adjustcolor(wes_palettes$Zissou[1], alpha.f = 0.5))
    dens(combined_mod_phy$Sol[ , sol], 
         show.zero = TRUE, 
         col = "white", 
         xlab = "", 
         ylab = "", 
         main = var_names, 
         add = TRUE)
    polygon(density(combined_mod_phy$Sol[ , sol]), 
            col = adjustcolor(wes_palettes$Zissou[4], alpha.f = 0.5))
    legend("topright", 
           legend = c(round(p.combined_mod_nophy[sol, ], 2), round(p.combined_mod_phy[sol, ], 2)),
           border = "white", 
           bty = "n", 
           text.col = c(wes_palettes$Zissou[1], wes_palettes$Zissou[4]))
  }

  var_names <- c("Pathogens", 
                 "Arranged Marriages", 
                 "Pop. Density", 
                 "Father Role", 
                 "Temperature", 
                 "Stratification = 2", 
                 "Stratification = 3",
                 "Inequality = 2", 
                 "Inequality = 3", 
                 "Warfare", 
                 "Assault", 
                 "Female Agriculture", 
                 "Sex Ratio")
  
  # build plots
  
    ### v871
  sol <- seq(3, 27, 2) # solution indices for v871
  
  tiff("Figure 3.tif", compression = "lzw", 
       height = 6.0, width = 12.0, units = "cm", res = 600, pointsize = 5)
  par(mfrow = c(3, 5), mar = c(2, 2, 4, 2))
  
  mapply(plot_post_dist, sol = sol, var_names = var_names)
 
  par(xpd=NA)
  legend(x=5, y=1, fill=c(adjustcolor(wes_palettes$Zissou[1], alpha.f=0.5),adjustcolor(wes_palettes$Zissou[4], alpha.f=0.5)), border="black", legend=c("No phylogeny", "Phylogeny"), bty="n", cex=1.5)
 
  dev.off()

    ### v860
  sol <- seq(4, 28, 2) # solution indices for v860
  
  tiff("Figure S1.tif", compression = "lzw", 
       height = 6.0, width = 12.0, units = "cm", res = 600, pointsize = 5)
  par(mfrow = c(3, 5), mar = c(2, 2, 4, 2))
  
  mapply(plot_post_dist, sol = sol, var_names = var_names)
  
  par(xpd=NA)
  legend(x=5, y=1, fill=c(adjustcolor(wes_palettes$Zissou[1], alpha.f=0.5),adjustcolor(wes_palettes$Zissou[4], alpha.f=0.5)), border="black", legend=c("No phylogeny", "Phylogeny"), bty="n", cex=1.5)
 
  dev.off()

  
  
        
        ### Figure R2 distributions ###
  
        # proportion of variance explained by fixed effects, phylo; covariance for phylo of different outcomes
        # for fixed effects, have to do for each of the 100 models, because the input data are different
  
  # calculate R2 function
  r2_v871 <- function(mod)
  {
    X<- as.matrix(mod$X)[1:184,]
    sf <- 0
    for (i in seq(3, 27, 2)) {
      fixed <- coefTable(mod)[i,1] * X[,i]
      sf <- sf + fixed
    }
    VarF <- var(sf) # variance in fitted values
    VarG <- mean(mod$VCV[,'traitv871:traitv871.phylo']) # variance in random effects
    VarR <- mean(mod$VCV[,'traitv871:traitv871.units']) # residual variance
    VarD <- log(1 + 1/exp(mean(sf))) # distribution-specific variance for poisson, see Johnson 2014
    R2m <- VarF/ (VarF + VarG + VarR + VarD) ## Marginal R2
    R2c <- (VarF + VarG) / (VarF + VarG + VarR + VarD) ## Conditional R2
    R2g <- VarG/ (VarF + VarG + VarR + VarD) ## Proportion of total variance explained by phylo
    output <- list(R2m, R2c, R2g)
    output
  }
  results_r2_v871 <- sapply(list_mods_phy, r2_v871)
  rownames(results_r2_v871) <- c("R2m", "R2c", "R2g")

  r2_v860 <- function(mod)
  {
    X<- as.matrix(mod$X)[185:368,]
    sf <- 0
    for (i in seq(4, 28, 2)) {
      fixed <- coefTable(mod)[i,1] * X[,i]
      sf <- sf + fixed
    }
    VarF <- var(sf) # variance in fitted values
    VarG <- mean(mod$VCV[,'traitv860:traitv860.phylo']) # variance in random effects
    VarR <- 1 # residual variance (fixed to 1 in this case)
    VarD <- 1 # distribution-specific variance for probit (ordinal), see Hadfield course notes
    R2m <- VarF/ (VarF + VarG + VarR + VarD) ## Marginal R2
    R2c <- (VarF + VarG) / (VarF + VarG + VarR + VarD) ## Conditional R2
    R2g <- VarG/ (VarF + VarG + VarR + VarD) ## Proportion of total variance explained by phylo
    output <- list(R2m, R2c, R2g)
    output
  }
  results_r2_v860 <- sapply(list_mods_phy, r2_v860)
  rownames(results_r2_v860) <- c("R2m", "R2c", "R2g")


      # plot density distribution of R2m, R2g
  
  # plotting function
  plot_dens_dist_R2 <- function(results, dens_mains) {
  dens(as.numeric(results[1, ]), 
       show.zero = FALSE, 
       main = dens_mains, 
       col = "white", 
       xlab = "Proportion of variance explained", 
       ylab = "Density across imputation", 
       xlim = c(0,1))
  polygon(density(as.numeric(results[1, ])),
          col = adjustcolor(wes_palettes$Zissou[1], alpha.f = 0.5))
  polygon(density(as.numeric(results[3, ])), 
          col = adjustcolor(wes_palettes$Zissou[4], alpha.f = 0.5))
  polygon(density(as.numeric(results[2, ])), 
          col = adjustcolor(wes_palettes$Zissou[5], alpha.f = 0.5))
  
  legend(0,11, adj = c(0.15, 0.5), 
         legend = "Predictors", 
         bty = "n", 
         fill = adjustcolor(wes_palettes$Zissou[1], alpha.f = 0.5))
  legend(0.3,11, 
         adj = c(0.15, 0.5), 
         legend = "Phylogeny", 
         bty = "n", 
         fill = adjustcolor(wes_palettes$Zissou[4], alpha.f = 0.5))
  legend(0.6,11, adj = c(0.1,0.5), 
         legend = "Pred + Phylo", 
         bty = "n", 
         fill = adjustcolor(wes_palettes$Zissou[5], alpha.f = 0.5))
  }
  
  results <- list(results_r2_v871, results_r2_v860)
  dens_mains <- c("% Married men polygynous (v871)",
                "Cultural rules constraining polygyny (v860)") 
  
  # build plot
  
  tiff("Figure 4.tif", compression = "lzw", 
       height = 6.0, width = 6.0, units = "cm", res = 600, pointsize = 5)
  par(mfrow = c(2, 1), mar = c(4, 3, 4, 2))
  
  mapply(plot_dens_dist_R2, results = results, dens_mains = dens_mains)

  dev.off()


      # for phylogenetic signal, for v871 have to calculate VarD foe each df, then combine lambda samples
  lambda.v871 <- 1:900000
  for(i in 1:100){
    X<- as.matrix(list_mods_phy[[i]]$X)[1:184,]
    sf <- 0
    for (j in c(3,5,7,9,11,13,15,17,19,21,23,25,27)) {
      fixed <- coefTable(list_mods_phy[[i]])[j,1] * X[,j]
      sf <- sf + fixed
  }
  lambda.v871[(1+(i-1)*9000):(i*9000)] <- list_mods_phy[[i]]$VCV[,'traitv871:traitv871.phylo'] / 
    (list_mods_phy[[i]]$VCV[,'traitv871:traitv871.phylo'] 
     + list_mods_phy[[i]]$VCV[,'traitv871:traitv871.units'] 
     + log(1 + 1/exp(mean(sf)))) 
  }
  
  hist(lambda.v871)
  l.v871.mean <- mean(lambda.v871) # 0.45
  l.v871.HPDI<- HPDI(lambda.v871, prob=0.95) # 0.18 0.73

  # for v860 just calculate ICC from combined model
  lambda.v860 <- combined_mod_phy$VCV[,'traitv860:traitv860.phylo']/
    (combined_mod_phy$VCV[,'traitv860:traitv860.phylo']+combined_mod_phy$VCV[,'traitv860:traitv860.units']+1)
  hist(lambda.v860)
  l.v860.mean <- mean(lambda.v860) # 0.56
  l.v860.HPDI<- HPDI(lambda.v860, prob=0.95) # 0.27 0.84



  ## calculate correlation for phylogenetic covariance
  cor(combined_mod_phy$VCV) # overall 0.9272267
  cor.v871v860.phylo <- combined_mod_phy$VCV[,2] / sqrt(combined_mod_phy$VCV[,1]*combined_mod_phy$VCV[,4])
  # by sample -> cor12 = cov12 / sqrt(var1*var2)
  median(cor.v871v860.phylo) # 0.85
  mean(cor.v871v860.phylo) # 0.81
  HPDI(cor.v871v860.phylo, prob=0.95) # 0.55 0.98


  
  
          ### Figure Tree of SCCS populations ###

  # drop tips not in data/polygyny
  data.multi$phylo <- phylonames[phylonames!="Ajie" & phylonames!="Gond"]
  v871 <- na.omit(data.multi[ , c("v871", "phylo")])
  
  phylo_v871 <- drop.tip(phylo, setdiff(phylo$tip.label, v871$phylo))
  
  # create named vector of variable to be mapped onto tree
  polygyny <- v871$v871
  names(polygyny) <- v871$phylo
  
  # create a named vector to match the tiplabels
  obj <- contMap(phylo_v871, polygyny, plot = FALSE, type = "fan", method = "fastAnc")
  
  # change color palette --> see http://blog.phytools.org/2014/05/changing-color-ramp-in-contmap-or.html
  n <- length(obj$cols)
  obj$cols[1:n] <- colorRampPalette(c("#f0f9e8", "#ccebc5", "#7bccc4", "#43a2ca", "#0868ac"), space="Lab")(n)
  
  # build plot
  
  tiff("Figure 2.tif", compression = "lzw", 
       height = 6.0, width = 6.0, units = "cm", res = 600, pointsize = 5)
  
  plot(obj, type = "fan", cex = 0.25, mar = c(0, 0, 0, 0),
       legend = FALSE, lwd = 2, outline = FALSE, fsize = c(0.6, 1), sig = 2)
  add.color.bar(100, obj$cols, title = "% Polygyny",
                lims = obj$lims, digits = 2, prompt = FALSE, x = 1,
                y = -30, lwd = 2, fsize = 1, subtitle = "length = 100ky")
  dev.off()
  
  
  
        ### Figure Map of SCCS populations ###
  
  ## Loading food sharing data and sccs coordinates
  d <- read.csv("polygyny_data.csv")
  d1 <- read.csv("dplace1.csv")
  d1 <- subset(d1, d1$Source == "Standard cross-cultural sample")
  d1$id <- as.numeric(substring(as.character(d1$Society.id), 5))
  d1 <- d1[,c(8,9,25)]
  
  d <- d[!is.na(d$men), ]
  colnames(d)[1] <- "id"
  d <- merge(d, d1, by.y = "id")
  
  # custom color range
  cols <- colorRampPalette(c("#f0f9e8", "#ccebc5", "#7bccc4", "#43a2ca", "#0868ac"))(100)
  
  tiff("Figure 1.tif", compression = "lzw", 
       height = 7.0, width = 12.0, units = "cm", res=600, pointsize=5)
  map(database = 'worldHires', fill = TRUE, col = "white", bg = "#969696", border = NA) # plots worldmap
  points(d$Revised.longitude, d$Revised.latitude, cex = 1.5, pch = 16, col = cols[d$men])
  scale = (length(cols) - 1) / 100
  for (i in 1:(length(cols) - 1)) {
    x = 2*(i - 1)/scale
    rect(x - 100, -90, (x + 2 / scale) - 100, -75, col = cols[i], border = NA)
  }
  arrows(x0 = -100, x1 = 100, y0 = -87, y1 = -87, length = 0, col = "white")
  mtext("0", 1, adj = 0.225, col = "white")
  mtext("50", 1, adj = 0.5, col = "white")
  mtext("100", 1, adj = 0.75, col = "white")
  mtext("% Married Men Polygynous (v871)", 1, line = 1, adj = 0.5, col = "white", cex = 1.25)
  dev.off()

  

        ### Figure prediction plots ###
        
  # predictive plots for the predictors with the clearest effects on both polygyny measures

        # v871; v1260 pathogen stress & v1666 assault frequency
  
  set.seed(13)
  tiff("Figure 5.tif", compression = "lzw", 
       height = 6.0, width = 12.0, units = "cm", res = 600, pointsize = 5)
  par(mfrow = c(1, 2))
  
  plot(v871 ~ scale(log(v1260)), 
       data.multi, 
       xaxt = "n", 
       col = "white",
       xlab = "Pathogen Stress (v1260)", 
       ylab = "% Men Married Polygynously (v871)",
       ylim = c(0, 45))
  
  v1260 <- scale(seq(log(min(data.multi$v1260, na.rm = TRUE)), 
                     log(max(data.multi$v1260, na.rm = TRUE)),
                     length = 100))
  
  mu.link.phy <- function(v1260) as.numeric(combined_mod_phy$Sol[,"traitv871"]) +
    as.numeric(combined_mod_phy$Sol[ , "traitv871:v1260"])*v1260
  mu.pathogen.phy <- sapply(v1260, mu.link.phy) # this is same as link() in rethinking package
  mu.pathogen.phy.mean <- apply(mu.pathogen.phy, 2, mean)
  for (i in sample(1:nrow(mu.pathogen.phy), 100)) {
    lines(exp(mu.pathogen.phy[i, ]) ~ v1260, 
          col = col.alpha(wes_palettes$Zissou[1], 0.2))
  }
  lines(exp(mu.pathogen.phy.mean) ~ v1260, lwd = 2)
  axis(1, at = c(-1.65697167, -0.55123585, 0.70575543, 1.74886040), labels=c("7", "10", "15", "21"))
  mtext("A", adj = -0.15, line = 2, cex = 2)
  
  plot(v871 ~ v1666,
       data.multi, 
       col = "white", 
       xlab = "Assault Frequency (v1666)", 
       ylab = "% Men Married Polygynously (v871)", 
       ylim = c(0, 45))
  
  v1666 <- seq(min(data.multi$v1666, na.rm = TRUE), 
               max(data.multi$v1666, na.rm = TRUE), 
               length = 100) 
  
  mu.link.phy <- function(v1666) as.numeric(combined_mod_phy$Sol[ , "traitv871"]) +
    as.numeric(combined_mod_phy$Sol[ , "traitv871:v1666"])*v1666
  mu.assault.phy <- sapply(v1666, mu.link.phy) 
  mu.assault.phy.mean <- apply(mu.assault.phy, 2, mean)
  for (i in sample(1:nrow(mu.assault.phy), 100)) {
    lines(exp(mu.assault.phy[i, ]) ~ v1666, 
          col = col.alpha(wes_palettes$Zissou[1], 0.2))
  }
  lines(exp(mu.assault.phy.mean) ~ v1666, lwd = 2)
  mtext("B", adj = -0.15, line = 2, cex = 2)
  
  dev.off()
  
  # summarize range of predictions:
  exp(min(mu.pathogen.phy.mean)); exp(max(mu.pathogen.phy.mean))
  # 1.630345 - 6.36252
  
  # --> going from min to max pathogen stress, % polygynously married men predicted to increase from 1.7% to 7.0%
  exp(max(mu.pathogen.phy.mean)) / exp(min(mu.pathogen.phy.mean))
  # --> 3.902559 fold increase in polygynously married men
 
  # summarize range of predictions:
  exp(min(mu.assault.phy.mean)); exp(max(mu.assault.phy.mean))
  # 3.565793 - 8.049515
  
  # --> going from min to max assault, % polygynously married men predicted to increase from 3.8% to 8.9%
  exp(max(mu.assault.phy.mean))/exp(min(mu.assault.phy.mean))
  # --> ~2 (2.257426) fold increase in polygynously married men
  
  
  
        # v860; v1260 pathogen stress
  
  ## predicted probabilities for ordinal (see https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q1/019668.html):
  # Imagine a latent variable (l) that conforms to the  standard linear model: l = Xb+Zu+e (forget about Zu, that's the random effects)
  # The probabilities of falling into each of the five categories are:
  # pnorm(-l)
  # pnorm(cp[1]-l)-pnorm(-l)
  # pnorm(cp[2]-l)-pnorm(cp[1]-l)
  # pnorm(cp[3]-l)-pnorm(cp[2]-l)
  # 1-pnorm(cp[3]-l)
  # and get predicted probabilities, using the results from MCMCglmm, and the "hunch" approach Jarrod mentioned in the thread above:
  # if the residual variance is set to one and there are no random effects
  # pnorm(-Xb, 0, sqrt(2)),
  # pnorm(cp[1] - Xb,0, sqrt(2)) - pnorm(-Xb,0, sqrt(2))
  # pnorm(cp[2] - Xb,0, sqrt(2)) - pnorm(cp[1] - Xb,0, sqrt(2))
  # 1- pnorm(cp[2] - Xb,0, sqrt(2))

  mu.link.phy <- function(v1260) as.numeric(combined_mod_phy$Sol[ , "traitv860"]) +
    as.numeric(combined_mod_phy$Sol[ , "traitv860:v1260"])*v1260

  v1260 <- scale(seq(log(min(data.multi$v1260, na.rm = TRUE)), 
                     log(max(data.multi$v1260, na.rm = TRUE)), 
                     length = 100))
  
  mu.pathogen.phy <- sapply(v1260, mu.link.phy) 
  mu.pathogen.phy.mean <- apply(mu.pathogen.phy, 2, mean)

  p1.pathogen.phy <- pnorm(-mu.pathogen.phy.mean, 0, sqrt(2)) # probability of level 1
  p2.pathogen.phy <- pnorm(combined_mod_phy_cp_means[1] - mu.pathogen.phy.mean,0, sqrt(2)) - pnorm(-mu.pathogen.phy.mean,0, sqrt(2)) # probability of level 2
  p3.pathogen.phy <- pnorm(combined_mod_phy_cp_means[2] - mu.pathogen.phy.mean,0, sqrt(2)) - pnorm(combined_mod_phy_cp_means[1] - mu.pathogen.phy.mean,0, sqrt(2)) # probability of level 3
  p4.pathogen.phy <- pnorm(combined_mod_phy_cp_means[3] - mu.pathogen.phy.mean,0, sqrt(2)) - pnorm(combined_mod_phy_cp_means[2] - mu.pathogen.phy.mean,0, sqrt(2)) # probability of level 4
  p5.pathogen.phy <- 1 - pnorm(combined_mod_phy_cp_means[3] - mu.pathogen.phy.mean,0, sqrt(2)) # probability of level 5
  
  # build plot
  ## use same color scale as in tree and map: "#f0f9e8", "#ccebc5", "#7bccc4", "#43a2ca", "#0868ac"
  
  tiff("Figure 6.tif", compression = "lzw", 
       height = 6.0, width = 9.0, units = "cm", res = 600, pointsize = 5)
  par(mfrow = c(2,3))
  
  plot(v860 ~ scale(log(v1260)), 
       data.multi, 
       main = "Monogamy prescribed", 
       xaxt = "n", 
       col = "white", 
       ylim = c(0, 0.7), 
       xlab = "Pathogen Stress (v1260)",
       ylab = "Probability of cultural rules (v860)")
  axis(1, 
       at = c(-1.65697167, -0.55123585, 0.70575543, 1.74886040), 
       labels = c("7", "10", "15", "21"))
    for (i in sample(1:nrow(mu.pathogen.phy), 100)){
    lines(pnorm(-mu.pathogen.phy[i, ], 0, sqrt(2)) ~ v1260, col = col.alpha("#f0f9e8", 0.2), lwd = 0.5)
  }
  lines(p1.pathogen.phy ~ v1260)
  
  plot(v860 ~ scale(log(v1260)), 
       data.multi, 
       main = "Monogamy preferred", 
       xaxt = "n", 
       col = "white",
       ylim = c(0, 0.7), 
       xlab = "Pathogen Stress (v1260)",
       ylab = "Probability of cultural rules (v860)")
  axis(1, 
       at = c(-1.65697167, -0.55123585, 0.70575543, 1.74886040), 
       labels = c("7", "10", "15", "21"))
    for (i in sample(1:nrow(mu.pathogen.phy), 100)){
    lines(pnorm(combined_mod_phy$CP[i, 1] -mu.pathogen.phy[i, ], 0, sqrt(2)) - pnorm(-mu.pathogen.phy[i, ], 0, sqrt(2)) ~ v1260, col = col.alpha("#ccebc5", 0.2), lwd = 0.5)
  }
  lines(p2.pathogen.phy ~ v1260)
  
  plot(v860 ~ scale(log(v1260)), 
       data.multi, 
       main = "Polygyny for leaders", 
       xaxt = "n",
       col = "white", 
       ylim = c(0,0.7), 
       xlab = "Pathogen Stress (v1260)", 
       ylab = "Probability of cultural rules (v860)")
  axis(1, 
       at = c(-1.65697167, -0.55123585, 0.70575543, 1.74886040), 
       labels = c("7", "10", "15", "21"))
    for (i in sample(1:nrow(mu.pathogen.phy), 100)){
    lines(pnorm(combined_mod_phy$CP[i, 2] - mu.pathogen.phy[i, ], 0, sqrt(2)) - pnorm(combined_mod_phy$CP[i, 1] - mu.pathogen.phy[i, ], 0, sqrt(2)) ~ v1260, col = col.alpha("#7bccc4", 0.2), lwd = 0.5)
  }
  lines(p3.pathogen.phy ~ v1260)
  
  plot(v860 ~ scale(log(v1260)),
       data.multi, 
       main = "Polygyny for upper class",
       xaxt = "n",
       col = "white",
       ylim = c(0,0.7),
       xlab = "Pathogen Stress (v1260)", 
       ylab = "Probability of cultural rules (v860)")
  axis(1, 
       at = c(-1.65697167, -0.55123585, 0.70575543, 1.74886040), 
       labels = c("7", "10", "15", "21"))
  for (i in sample(1:nrow(mu.pathogen.phy), 100)){
    lines(pnorm(combined_mod_phy$CP[i, 3] - mu.pathogen.phy[i, ], 0, sqrt(2)) - pnorm(combined_mod_phy$CP[i, 2] - mu.pathogen.phy[i, ], 0, sqrt(2)) ~ v1260, col = col.alpha("#43a2ca", 0.2), lwd = 0.5)
  }
    lines(p4.pathogen.phy ~ v1260)
  
  plot(v860 ~ scale(log(v1260)),
       data.multi,
       main = "Polygyny prevalent", 
       xaxt = "n", 
       col = "white",
       ylim = c(0, 0.7),
       xlab = "Pathogen Stress (v1260)",
       ylab = "Probability of cultural rules (v860)")
  axis(1, 
       at = c(-1.65697167, -0.55123585, 0.70575543, 1.74886040), 
       labels = c("7", "10", "15", "21"))
  for (i in sample(1:nrow(mu.pathogen.phy), 100)){
    lines(1 - pnorm(combined_mod_phy$CP[i, 3] - mu.pathogen.phy[i, ], 0, sqrt(2)) ~ v1260, col = col.alpha("#0868ac", 0.2), lwd = 0.5)
	}
	lines(p5.pathogen.phy ~ v1260)
  dev.off()
  

        # v860; v1666 assault frequency

  mu.link.phy <- function(v1666) as.numeric(combined_mod_phy$Sol[ ,"traitv860"]) +
    as.numeric(combined_mod_phy$Sol[ ,"traitv860:v1666"])*v1666
  
  v1666 <- seq(min(data.multi$v1666, na.rm = TRUE),
               max(data.multi$v1666, na.rm = TRUE),
               length = 100) 
  
  mu.assault.phy <- sapply(v1666, mu.link.phy) # this is same as link() in rethinking package
  mu.assault.phy.mean <- apply(mu.assault.phy, 2, mean)

  p1.assault.phy <- pnorm(-mu.assault.phy.mean, 0, sqrt(2)) # probability of level 1
  p2.assault.phy <- pnorm(combined_mod_phy_cp_means[1] - mu.assault.phy.mean,0, sqrt(2)) 
                    - pnorm(-mu.assault.phy.mean,0, sqrt(2)) # probability of level 2
  p3.assault.phy <- pnorm(combined_mod_phy_cp_means[2] - mu.assault.phy.mean, 0, sqrt(2)) - pnorm(combined_mod_phy_cp_means[1] - mu.assault.phy.mean, 0, sqrt(2)) # probability of level 3
  p4.assault.phy <- pnorm(combined_mod_phy_cp_means[3] - mu.assault.phy.mean, 0, sqrt(2)) - pnorm(combined_mod_phy_cp_means[2] - mu.assault.phy.mean, 0, sqrt(2)) # probability of level 4
  p5.assault.phy <- 1 - pnorm(combined_mod_phy_cp_means[3] - mu.assault.phy.mean, 0, sqrt(2)) # probability of level 5
  
  ## individual spaghetti plots
  tiff("Figure S3.tif", compression = "lzw", 
       height = 6.0, width = 9.0, units = "cm", res = 600, pointsize = 5)
  par(mfrow = c(2, 3))
  
  plot(v860 ~ v1666,
       data.multi,
       main = "Monogamy prescribed", 
       col = "white",
       ylim = c(0, 0.7), 
       xlab = "Assault frequency (v1666)", 
       ylab = "Probability of cultural rules (v860)")
  for (i in sample(1:nrow(mu.assault.phy), 100)) {
    lines(pnorm(-mu.assault.phy[i,], 0, sqrt(2)) ~ v1666, col = col.alpha("#f0f9e8", 0.2), lwd = 0.5)
  }
    lines(p1.assault.phy ~ v1666)
	
  plot(v860 ~ v1666, 
       data.multi,
       main = "Monogamy preferred",
       col = "white",
       ylim = c(0, 0.7),
       xlab = "Assault frequency (v1666)",
       ylab = "Probability of cultural rules (v860)")
  for (i in sample(1:nrow(mu.assault.phy), 100)) {
    lines(pnorm(combined_mod_phy$CP[i, 1] -mu.assault.phy[i, ], 0, sqrt(2)) - pnorm(-mu.assault.phy[i, ], 0, sqrt(2)) ~ v1666, col = col.alpha("#ccebc5", 0.2), lwd = 0.5)
  }
  lines(p2.assault.phy ~ v1666)
  
  plot(v860 ~ v1666, 
       data.multi,
       main = "Polygyny for leaders",
       col = "white",
       ylim = c(0, 0.7), 
       xlab = "Assault frequency (v1666)", 
       ylab = "Probability of cultural rules (v860)")
  for (i in sample(1:nrow(mu.assault.phy), 100)) {
    lines(pnorm(combined_mod_phy$CP[i, 2] - mu.assault.phy[i, ], 0, sqrt(2)) - pnorm(combined_mod_phy$CP[i, 1] - mu.assault.phy[i, ], 0, sqrt(2)) ~ v1666, col = col.alpha("#7bccc4", 0.2), lwd = 0.5)
  }
  lines(p3.assault.phy ~ v1666)
  
  plot(v860 ~ v1666, 
       data.multi,
       main = "Polygyny for upper class", 
       col = "white", 
       ylim = c(0, 0.7), 
       xlab = "Assault frequency (v1666)",
       ylab = "Probability of cultural rules (v860)")
  for (i in sample(1:nrow(mu.assault.phy), 100)) {
    lines(pnorm(combined_mod_phy$CP[i, 3] - mu.assault.phy[i, ], 0, sqrt(2)) - pnorm(combined_mod_phy$CP[i, 2] - mu.assault.phy[i, ], 0, sqrt(2)) ~ v1666, col = col.alpha("#43a2ca", 0.2), lwd = 0.5)
  }
  lines(p4.assault.phy ~ v1666)
  
  plot(v860 ~ v1666,
       data.multi,
       main = "Polygyny prevalent", 
       col = "white",
       ylim = c(0, 0.7), 
       xlab = "Assault frequency (v1666)",
       ylab = "Probability of cultural rules (v860)")
  for (i in sample(1:nrow(mu.assault.phy), 100)) {
    lines(1 - pnorm(combined_mod_phy$CP[i,3] - mu.assault.phy[i,],0, sqrt(2)) ~ v1666, col = col.alpha("#0868ac", 0.2), lwd = 0.5)
  }
  lines(p5.assault.phy ~ v1666)
  dev.off()
  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         