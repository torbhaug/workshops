# Setting libraries
library(bootnet)
library(qgraph)
library(ggm)
library(graphicalVAR)

# Loading in the data
predictor_data <- read.csv("C:/Users/thaug/psychosys_winter_schhool/behav_ind_model_1.csv")
data <- read.csv("C:/Users/thaug/psychosys_winter_schhool/behav_model_1.csv")

# Calculating the correlations between the nodes
data.cor <- cor(predictor_data)
data.pcor <- cor2pcor(data.cor)

vars <- c("smoothed_response", "character_change", "object_change" "location_change", "time_of_day_change",
          "sound_presence_change")

data <- cbind(data, predictor_data)


# 
?graphicalVAR
model_data_fit <- graphicalVAR(data, gamma=0, vars = vars)

# Testing for significane like stuff in the thing
sin.ag(data.pcor, n = 1807)


# Visualizing
qgraph(data.pcor)

# We are then fitting a basic network model to the data and setting a threshold
# gamma sets the EBIC hyperparameter which is the default mode of the 
# model as with Lasso etc
# this has the ability to include participant ID in a vector is wanted with idvar = c()

# However, before we get this far we simply have to reduce the number of variables
# since the predictor stack is dummy thicc
#
# "Nooo, no thow complex eigenvalue error, your predictor stack is so fat"

# Either we manually select the variables we care the most about
include <- c("bymusa.absolute.change", "skogmusa.absolute.change",
             "Timepassage.absolute.change", "Relocation.absolute.change", "Causality.absolute.change",
             "Motorveien.absolute.change", "Bakeriet.absolute.change", "Furu.absolute.change",
             "Chapter.Break.absolute.change", "Foreground.Sound.absolute.change", "Background.Sound.absolute.change",
             "Speech.Sound.absolute.change", "Music.Sound.absolute.change"
             )
# Checking to find mismatches
setdiff(include, colnames(data))


data_subset <- predictor_data[, include]

# Or we can be dummy simple and just use the silly simple predictor set
# Here we have do some more checks though because we still break that thing - likely due to lacking variability
nrow(data)
ncol(data)


fit.pred <- graphicalVAR(
  data = predictor_data,
  lags = 0,
  nLambda = 10,
  lambda_min_beta = 0.01,
  lambda_min_kappa = 0.01,
  verbose = TRUE
  )


plot(network)

centralityPlot(network, scale = "raw0", include = c("Foreground.Sound.absolute.change", "Background.Sound.absolute.change", "Chapter.Break.absolute.change", "Relocation.absolute.change"))

centralityPlot()


if (!requireNamespace("vars", quietly = TRUE)) {
  install.packages("vars")
}

# Trying a basic VAR model
library(vars)

?VAR

lag_selection <- VARselect(data, lag.max = 10, type = "const")
lag_selection$criteria

p <- lag_selection$selection["AIC(n)"]

# Removing time so that it doesnt allahu akbar the entire regression
data <- data[, -1]

# Estimating VAR model
var_model <- VAR(
  data,
  p = 1,
  type = "trend"  # "none", "trend", or "both" also possible
)

summary(var_model)

# If one variable has a schock at one timepoint what is the 
irf(var_model,
    impulse = "sound_presence_change",
    response = "smoothed_response",
    n.ahead = 20,
    boot = TRUE)

