# ------------------------------------- #
# Example R-INLA code for a BYM CAR model
#  
# Author: Ian Buller (@idblr)
# Date created: February 15, 2019
# 
# Most recently modified on: 02/16/2019
# Most recently modified by: Ian Buller
# Modifications:
# A) Added overall label for spplot using the grid package (2/15/2019)
# B) Added code for merging two data.frames into one spatial polygon data frame (SPDF) (2/15/2019)
# C) Fixed section labels (2/15/2019)
# D) Replaced "Y" in INLA models with "Case_Count" to match variable in spdf (and dat_df) (2/15/2019)
# E) Added package "sp" and "grid" (2/15/2019)
# F) Major fixes to the interpretation section (UH and SH residuals added; accurate interpretation of residual relative risk and excess relative risk) (2/15/2019)
# G) Added area specific zero-probability calculation (2/16/2019)
# ------------------------------------- #

####################
##### PACKAGES #####
####################

library(sp) # for spatial data manipulation
library(INLA) # for R-INLA
library(spdep) # for adjacency matrix creation
library(tigris) # for TIGRIS county shapefile
library(colorspace) # for HCL color palette
library(grid) # for spplot customization

################
##### DATA #####
################

# County Shapefiles (Based on 2010 County Boundaries)
tigris_spdf <- tigris::counties()

# load your health data
health_df <- XXXXXXXXX # with a "FIPS" column and "Case_Counts" column

# Merge case data with county shapefile
spdf <- tigris::geo_join(tigris_spdf, health_df, "GEOID", "FIPS")
tigris_spdf <- NULL
spdf$Case_Count[is.na(spdf$Case_Count)] <- 0 # Set NA counties to 0 

# Alders Coordinate Reference System for United States
cref1 <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +7_0=0 +a=6570997 +b=6370997 +units=m +no_defs"
spdf <- sp::spTransform(spdf, sp::CRS(cref1)) # transform US counties to alders CRS

# Create adjacency matrix
spdf_adj <- spdep::poly2nb(spdf, queen = T) # set queen = F for rook
adj <- spdep::nb2mat(spdf_adj, style = "B") # create binary matrix
adj <- as(adj, "dgTMatrix") # create sparse matrix, memory conservation

# Create ID variable (if not already in SPDF)
spdf$idx <- 1:nrow(spdf)

# Extract data.frame from the SDPF
dat_df <- as.data.frame(spdf)

#################################################################################
##### R-INLA Conditional Autoregressive (CAR) Besag–York–Mollié (BYM) Model #####
#################################################################################

## See Chapter 6 (especially 6.1.2 for BYM) in the "Spatial and Spatio-Temporal Bayesian Models in R-INLA" by Blangiardo and Cameletti (2015)
# Section 6.3 is about ZIP models

####################
##### POISSON ######
####################

form_null <- Case_Count ~ 1 + f(idx, model = "bym", graph = adj) # intercept and bym random effects; graph is set to sparse adjacency matrix
# Can add any variable for multivariate analysis: Case_Count ~ 1 + VAR1 + f(idx, model = "bym", graph = adj)
m_null_pois <- INLA::inla(form_null,
                          data = dat_df, 
                          family = "poisson", # for Poisson model
                          control.predictor = list(compute = T), # to compute predicted values
                          control.compute = list(dic = T, waic = T), # to calculate model performance measurements; here, DIC and WAIC
                          offset = log(population), # offset here is optional, if expected counts use argument E
                          verbose = FALSE # set to TRUE if want a more detailed print out, helps debugging
)
summary(m_null_pois) # summary of model

###################################
##### Zero-inflated Poisson 0 ##### 
###################################

# ZIP0 assumes all zeros are structural

form_null <- Case_Count ~ 1 + f(idx, model = "bym", graph = adj)  # repeat of above
m_null_zip0 <- INLA::inla(form_null, 
                          data = dat_df, 
                          family = "zeroinflatedpoisson0", # for ZIP0 model
                          control.predictor = list(compute = T), 
                          control.compute = list(dic = T, waic = T),
                          offset = log(population),  # offset here is optional, if expected counts use argument E
                          verbose = FALSE
)
summary(m_null_zip0) # summary of model

###################################
##### Zero-inflated Poisson 1 ##### 
###################################

# ZIP1 assumes some zeros are structural and some are sampling zeroes

form_null <- Case_Count ~ 1 + f(idx, model = "bym", graph = adj) # repeat of above
m_null_zip1 <- INLA::inla(form_null,
                          data = dat_df, 
                          family = "zeroinflatedpoisson1",  # for ZIP1 model
                          control.predictor = list(compute = T), 
                          control.compute = list(dic = T, waic = T),
                          offset = log(population),  # offset here is optional, if expected counts use argument E
                          verbose = FALSE
)
summary(m_null_zip1) # summary of model

###################################
##### Hyperparameter Controls #####
###################################

# BYM has spatially unstructured and spatially structured random effects
# Can control the priors on both

form_null_control <- Case_Count ~ 1 + f(idx, model = "bym", graph = adj,
                                        hyper = list(prec.unstruct= list(prior = "loggamma", # spatially unstructured random effect (v)
                                                                         param = c(0.001,0.001) # modify as desired, default = logGamma(1,0.0005)
                                                                         ),
                                        prec.spatial = list(prior = "loggamma", # spatially structured random effect (u)
                                                            param = c(0.01,0.01) # modify as desired, default = logGamma(1,0.0005)
                                                            )
                                                    )

                                      )

m_null_control <- INLA::inla(form_null_control, 
                             family = "poisson", 
                             data = dat_df, 
                             control.predictor = list(compute = T), 
                             control.compute = list(dic = T, waic = T),
                             control.inla = list(tolerance = 1e-20, # can control precision
                                                 h = 1e-08 # can control h (search step parameter)
                                                 ),
                             offset = log(population),  # offset here is optional, if expected counts use argument E
                             verbose = FALSE
)
summary(m_null_control) # summary of model


##################
##### EXPORT #####
##################

# Export
save(m_null_pois, m_null_zip0, m_null_zip1, m_null_control,
     file = "inla_bym.Rdata"
     )

############################
###### INTERPRETATION ######
############################

# Check if any fixed effect 95% credibility interval span the null value (0)
# If it does not, statistically significant variable
# Intercept is the average outcome rate in the entire study region
# Exponentiate fixed effect coefficient means and 95% credibility interval
# Examine uncertainty of model by looking at the random effects and proportion of variance

## Exponentiate fixed effect coefficients and random effects
# Either exp() of the summary values (mean, 0.025quant, 0.975quant) or using the marginals (for some reason they can be slightly different)
exp(m_null_pois$summary.fixed[1,1]) # mean beta_0
exp(m_null_pois$summary.fixed[1,3]) # lower 95%CI beta_0
exp(m_null_pois$summary.fixed[1,5]) # upper 95%CI beta_0

exp.b0.mean <- INLA::inla.emarginal(exp,m_null_pois$marginals.fixed[[1]]) # mean beta_0
exp.b0.95CI <- INLA::inla.qmarginal(c(0.025,0.975), # 95%CI beta_0
                                    INLA::inla.tmarginal(exp,m_null_pois$marginals.fixed[[1]])
)
# REPEAT FOR ANY OTHER VARIABLES IN THE MODEL #

# If no other risk factors in the model (Disease Map), the county-specific relative risks are the exp(u + v)
# If other risk factors in the model (Ecological Regression), the county-specific residuals relative risks are the exp(u + v) compared to the entire extent after the risk factors are taken into account

nAreas <- nrow(dat_df) # number of areal units
UHres <-  m_null_pois$summary.random$idx$mean[1:nAreas] # area specific spatially unstructured residuals
SHres <-  m_null_pois$summary.random$idx$mean[(nAreas+1):(2*nAreas)] # area specific spatially structured residuals
csi <- m_null_pois$marginals.random$idx[1:nAreas] # reminder the first set of random effect marginals are u + v
zeta <- lapply(UH_coef, function(x) INLA::inla.emarginal(exp,x)) # posterior mean for random effects (relative risk, or if offset then counts)

a <- 0
prob_excess <- lapply(csi, function(x){1-INLA::inla.pmarginal(a,x)}) # posterior probability of (residual) relative risk exceeding 1 (excess risk)

## Examine uncertainty
# Proportion of variance
mat.marg <- matrix(NA, nrow=nAreas, ncol=100000) # create matrix to sample marginals
m <- m_null_pois$marginals.random$idx
for (i in 1:nAreas){
  #first nAreas values of the random effects are UH (u+v)
  #second nAreas values of the random effects are SH (u)
  u <- m[[nAreas+i]]
  mat.marg[i,] <- INLA::inla.rmarginal(100000, u)
  }
var.SH <- apply(mat.marg, 2, var) # emperical variance of spatially structured random effect

var.UH <- INLA::inla.rmarginal(100000, # expected variance of spatially unstructured random effect
                               INLA::inla.tmarginal(function(x) 1/x, m_null_pois$marginals.hyper$`Precision for ID (iid component)`)
                               )

perc.var.SH <- mean(var.SH/(var.SH+var.UH)) # calculate proportion of variance
perc.var.SH # View proportion of variance
# If closer to 1 then a large part of the variability is explained by the spatial structure
# Separate calculations if model was scaled (option scale.model = TRUE) see Blangiardo and Cameletti pp 186

## For ZIP models, pull out and interpret the zero-probability parameter
# inverse logit transformation of pi_0
# probability of pi_0 should be higher in ZIP0 model than ZIP1 models
# default = Normal(-2,1)
m_null_zip1$summary.hyperpar[1,] 
prob_zero1 <- m_null_zip1$summary.hyper[1,1]+(1-m_null_zip1$summary.hyper[1,1])*exp(-exp(m_null_zip1$summary.linear.predictor$mean)*log(spdf$population)) # area specific zero-probability 

# Any area specific predicted values have a 95% CI that crosses the null (0)
ci95 <- m_null_pois$summary.fitted.values[,c(3,5)]
ci95$sig <- with(ci95, ci95$'0.025quant' < 0 & ci95$'0.975quant' > 0) # can plot these, too (FALSE = Significant b.c. TRUE = cross null value)
table(is.na(ci95$sig)) # global count of CIs that cross null (0)

# size of 95% CI
ci95range <- abs(m_null_pois$summary.fitted.values[,5] - m_null_pois$summary.fitted.values[,3]) # upper 95% CI - lower 95% CI

##############################
##### DATA VISUALIZATION #####
##############################

# Append calculated values to SPDF
spdf$m_null_pois_fit <- m_null_pois$summary.fitted.values[,1] # area specific predicted mean values (RR if no offset, counts if offset)
spdf$m_null_pois_sd <- m_null_pois$summary.fitted.values[,2] # area specific predicted sd values (RR if no offset, counts if offset)
spdf$m_null_pois_ci95range <- ci95range # area specific size of 95% CI of predicted values (RR if no offset, counts if offset)
spdf$m_null_pois_UH <- UHres # area specific spatially unstructured residuals
spdf$m_null_pois_SH <- SHres # area specific spatially structured residuals
spdf$m_null_pois_RR <- unlist(zeta)/log(population) # area specific (residual) relative risk (if using offset to convert back to risk)
spdf$m_null_pois_RR <- unlist(zeta) # area specific (residual) relative risk (if not using offset because already risk)
spdf$m_null_pois_eRR <- unlist(probUH) # area specific probability of excess risk UH
spdf$m_null_pois_DIC <- m_null_pois$dic$local.dic # area specific DIC
spdf$m_null_pois_WAIC <- m_null_pois$dic$local.dic # area specific WAIC
spdf$m_null_zip1_probZ <- prob_zero1 # area specific probability of observing zero for ZIP1

# Plot
# col_lab <- XXXXXX # names of colorkey labels in vector
# at_lab <- XXXXX # location of colorkey lables in vector
# at_plot <- XXXXX # location of color breaks of plot in vector
# font_fam <- XXXXX # font family as a character
exfact <- 1 # PNG expansion factor (for higher res images)
grDevices::png("m_null_pois_fit.png", width = 2000*exfact, height = 2000*exfact) # create PNG
sp::spplot(spdf,  # SPDF
           "m_null_pois_fit", # column, variable
           col.regions=c("transparent", rev(colorspace::sequential_hcl(101))), # color scheme
           par.settings = list(axis.line = list(col =  'transparent')), # remove axess for maps
           colorkey=list(labels = list(labels = col_lab, # for custom colorkey labels
                                       at = at_lab, # for custom colorkey label positions, usually reduce the bottom and top number so they fit within window
                                       cex = 4*exfact, # for custom size of labels
                                       fontfamily= font_fam, # for custom font of labels
                                       fontface = 1), # for direction of labels
                         width = 1*exfact # for custom width of colorkey
           ),
           at = at_plot, # set custom color breaks
           lwd = 1*exfact  # custom line width
           )
grid::grid.text("INSERT LABEL FOR COLORKEY", # for overall label of colorkey
                x=unit(0.925, "npc"), # x position
                y=unit(0.50, "npc"), # y position
                rot=90, # for vertical colorkey
                gp = gpar(fontsize = 50*exfact, # for custom size of overall label
                          fontfamily= font_fam # for custom font of overall label
                          )
                )
dev.off() # turn off image device

# --------------- END OF CODE --------------- #