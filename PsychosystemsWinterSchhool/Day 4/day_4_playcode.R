# Time series data and stuffs

# This is a methodology that requires a lot of power unfortunately
# By summing the different predictors we will get more variance
# Could group them with Bayesian clustering off event boundaries
# 
# Contemporaneous networks are all effects that are faster than the
# time lag we employ
#
# Beta is lagged, Kappa is naur

# Psychonetrics uses Full Information Maximum Likelihood (FIML) for missing data handling

# We assume that properties such as mean, variance, autocorrelation etc. are all assumed to
# stay constant for the period we are measuring in

# kpss.test(data) will test if you have a linear trend in your data
# Can be detrended - but consider if that is what you want to do
# Detrending needs to be per variabel per person

# Also employs an equidistant mwasure assumption

# Continous time VAR is a way to handle missing data but is very hungry

# Individual network invariance test

#dplyr
if(!requireNamespace("dplyr", quietly= TRUE)){
  install.packages("dplyr")
}
library(dplyr)
#graphicalVAR
if(!requireNamespace("graphicalVAR", quietly= TRUE)){
  install.packages("graphicalVAR")
}
library(graphicalVAR)
#psychonetrics
if(!requireNamespace("psychonetrics", quietly= TRUE)){
  install.packages("psychonetrics")
}
library(psychonetrics)
#qgraph
if(!requireNamespace("qgraph", quietly= TRUE)){
  install.packages("qgraph")
}
library(qgraph)
#patchwork
if(!requireNamespace("patchwork",quietly= TRUE)){
  install.packages("patchwork")
}
library(patchwork)
# ggplot2
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)


# Renaming my columns
names(Data2)[names(Data2) %in% vars] <- varLabs

# Remove items:
Data2 <- Data2 %>% select(-Hungry,-Angry,-Music,-Procrastinate)
varLabs <- varLabs[!varLabs %in% c("Hungry","Angry","Music","Procrastinate")]

# Extracting single participant data
birth_date <- 07021996 # Enter date of birth here
set.seed(birth_date)
subject <- sample(Data2$id, 1)
my_data <- Data2[which(Data2$id == subject),]

# We then pick select variables because it can't estimate that much with pur power
vars <- c("Alone","Relax","Irritable", "Nervous")

varLabs <- c("Relax","Irritable","Worry","Nervous","Future","Anhedonia",
             "Tired","Hungry","Alone","Angry","Social_offline","Social_online",
             "Music","Procrastinate","Outdoors","C19_occupied","C19_worry",
             "Home")




# Trying to fit the mufuckery
library(graphicalVAR)
library(qgraph)
?graphicalVAR
var_fit <- graphicalVAR(my_data, gamma = 0.5, vars = vars)
qgraph(var_fit$PCC)

# Try to do the unregularized version
var_fit <- var_fit %>% runmodel() # Dinnae work

# We then move over to the single subject analysis
birth_date <- 666
# Enter date of birth here
set.seed(birth_date)
subject <-sample(Data2$id,1)
subject2 <- sample(Data2$id,1)
my_data<-Data2[which(Data2$id==subject),]
my_data2<-Data2[which(Data2$id==subject2),]

#
ind_var_fit <- graphicalVAR(my_data, gamma=0, vars = vars)
ind_var_fit2 <- graphicalVAR(my_data2, gamma=0, vars = vars)

# Comparing two participants
g1 <- qgraph(ind_var_fit$PCC)
g2 <- qgraph(ind_var_fit2$PCC)

avg_lay = averageLayout(g1, g2)
qgraph(ind_var_fit$PCC, layout=avg_lay)
qgraph(ind_var_fit$PCC, layout=avg_lay)


library(psychonetrics)
ind_var_fit <- gvar(my_data, vars = vars, estimator = "FIML")
ind_var_fit2 <- gvar(my_data2,vars = vars, estimator = "FIML")
ind_var_fit <- ind_var_fit %>% prune(alpha = 0.05)  %>% runmodel()
ind_var_fit2 <- ind_var_fit %>% prune(alpha = 0.05)  %>% runmodel()

ind_var_fit %>% fit()

ind_var_fit %>% getmatrix("omega_zeta")
ind_var_fit2 %>% getmatrix("PDC")

qgraph(ind_var_fit$PCC, layout=avg_lay)
qgraph(ind_var_fit2$PCC, layout=avg_lay)




# Multilevel Extension to Time-Series
# Combining power - and allows for separating effects

# Can look into individual differences by looking at standard deviations
# Time-series with grouped responses - and a poisson mixed level model
# They dont have this in the package but if you open it up I might be able
# to represent the model and response as count data


# VAR for panel data
# There does exist cross-lagged panel data approaches wchih can be used on our item
# responses captured during testing

# mVAR estimation


# Time varying var model
# What we are interested in is the lag varying var
# Actually we can use this potentially - no
# Bandwidth will impact in a tradeoff manner between long bandwidth(high stability to low
# sensitivity) to smnall bandwidth (less stable estimates but higher sensitivity)
# bwselect() could work as a datafriven approach

# Lunansky et al. (2020) for theoretical estimate of time scales

