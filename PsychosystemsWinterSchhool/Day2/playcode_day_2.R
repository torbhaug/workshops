# Loading in the packages (not doing it right now because already have in environment)


# Loading the test data
dat <- read.csv("C:/Users/thaug/psychosys_winter_schhool/Day2/NA_2020_data.csv")

hist(dat$Q7)


# Selecting a subset fo variables to make the model easier to whatever
include <-  c("Q10", # I try to keep a regular sleep pattern
  "Q13", # I am worried about my current sleeping behavior
  "Q14", # My sleep interferes with my daily functioning
  "Q68", # I am happy with my physical health.
  "Q70", # I feel optimistic about the future.
  "Q75", # I am very happy
  "Q77", # I often feel alone
  "Q80" # I am happy with my love life
)

data_subset <- dat[, include]

# Renaming variables to make the network more clear
names(data_subset) <- c("regular_sleep",
                        "worried_sleep",
                        "sleep_interfere",
                        "happy_health",
                        "optimistic_future",
                        "very_happy",
                        "feel_alone",
                        "happy_love_life"
                        )

# Calculating the correlation for the variables selected
datsub.cor <- cor(data_subset, use = "pairwise.complete.obs")

# We then use that to calculate the partial correlations
datsub.pcor <- cor2pcor(datsub.cor)

# So we inverse this somehow - I dont understand how this works
datsub.Kappa <- solve(datsub.pcor)

#  Alternatively we can also do ti ourselves with the covariance
datsub.cov <- cov(data_subset, use = "pairwise.complete.obs")

)
Kappa <- solve(datsub.cov)
# Then we use the provided formula to calculate the partial correlation from the
# inverted matrix

-1 * Kappa[i, j] /
  (sqrt(Kappa[i,i]) *
     sqrt(Kappa[j, j])
   )

# I cannae be arsed so Imma just do them all
datsub.pcor <- -1 * cov2cor(Kappa)
diag(datsub.pcor) <- 1
round(datsub.pcor, 2)


# We can also get the partial correlations by creating regressions from our different nodes
fitFeelAlone <- lm(feel_alone ~ regular_sleep + worried_sleep + sleep_interfere + 
                     happy_health + optimistic_future + very_happy + happy_love_life,
                   data = data_subset)

fitOptimisticFuture <- lm(optimistic_future ~ regular_sleep + worried_sleep + sleep_interfere + 
                     happy_health + feel_alone + very_happy + happy_love_life,
                   data = data_subset)

# We then obtain the coefficients
coefFeelAlone <- coef(fitFeelAlone)
coefOptimisticFuture <- coef(fitOptimisticFuture)

# Calculating residual standard deviation
sdFeelAlone <- sigma(fitFeelAlone)
sdfOptimisticFuture <- sigma(fitOptimisticFuture)

# We can then use these to compute the partial correlation between the nodes
(coefFeelAlone["optimistic_future"] * sdfOptimisticFuture) / sdFeelAlone

# To skip all these individual steps we can simply estimate the network with bootnet
datsub.Network <- estimateNetwork(data_subset, default = "pcor", corMethod = "cor")
plot(datsub.Network, vsize = 12, layout = "spring")

centralityPlot(datsub.Network, include = c("Strength", "Closeness", "Betweenness"), scale = "z-scores")


# In order to calculate the reliability of the model we perform a non-parametric bootstrap
num_feet = 1000
#num_cores = paralell::detectCores()
feet_non_parametric <- bootnet(datsub.Network,
                               nBoots = num_feet,
                               default = "pcor")

plot(feet_non_parametric, order="sample")


# However, in order to calculate the metrics and the how stable they are we 
# drop that strap namsayin - AKA case drop strap

feet_case_drop <- bootnet(datsub.Network,
                          type = "case",
                          statistics = c("strength", "betweenness")
                          )
# We then plot what happens when the feds drop that case nmsyin
plot(
  feet_case_drop,
  statistics = c("strength", "betweenness")
)

# Finally we calculate the CS coefficient thingy
corStability(feet_case_drop,
              statistics = c("strength", "betweenness"))
