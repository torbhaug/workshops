# We are doing crazy Bayesian stuff this time
# Lets mail bombs to frequentists

install.packages("easybgm")
install.packages("bgms")
install.packages("BGGMgraph")
install.packages("gsl", type = "binary")


library(easybgm)
library(bgms)
library(BGGMgraph)
library(dplyr)

data_bae2 = read.csv("https://raw.githubusercontent.com/Bayesian-Graphical-Modelling-Lab/BGM_Workshops/6bc21f12553120bac42abebb232dafe55ba87426/data/exampledata.csv")[, -1]
head(data_bae2)

?easybgm

easybgm


bae_data <- bae_data %>% select(-c(1, 2))

set.seed(1)
fit = easybgm(data = data_bae2,
              type = "continuous",
              package = "bgms",
              iter = 1e5,
              burnin = 1e4,
              save = TRUE,
              centrality = TRUE,
              progress = TRUE)
summary(fit)       

plot_edgeevidence(fit, evidence_thresh = 10, legend = FALSE)

plot_edgeevidence(fit, evidence_thresh = 100, legend = FALSE, split = TRUE)

plot_network(fit, exec_prob = 0.5)

plot_parameterHDI(fit)

# Now for centrality stuff
plot_centrality(fit)


# Prior sensitivity check
# Write down / pre-register the prior you are intending to make
# THen do a prior sensitivity check

fit <- easybgm(
  data = data_bae2[, -c(1:3)],
  type = "ordinal",
  package = "bgms",
  iter = 2e3, # bgms performs 4 chains (runs) with "iter" iterations in parallel
  update_method = "adaptive-metropolis", # default "nuts" takes long!
)

summary(fit)
plot_edgeevidence(fit)

# We filter the data to get
# We redo the analyzis described in the paper to write down our results and whatnot
fit = easybgm(data = bae_data,
              type = "ordinal",
              package = "bgms",
              iter = 1e5,
              burnin = 1e4,
              update_method = "adaptive-metropolis",
              save = TRUE,
              centrality = TRUE,
              progress = TRUE)
summary(fit)


plot_edgeevidence(fit, evidence_thresh = 10, legend = FALSE, split = TRUE)

?bgms

# Comparing networks

library(bgms)
?ADHD

#dat <- ADHD
group_indicator <- ADHD[, 1]
data = ADHD[, -1]

# Check group sizes
table(group_indicator)

# Inspect symptom prevalences by group
prevalence_group0 = colMeans(data[group_indicator == 0, , drop = FALSE])
prevalence_group1 = colMeans(data[group_indicator == 1, , drop = FALSE])

prevalence_group0
prevalence_group1

??easybgm

fit_diag = bgm(
  x = data[group_indicator == 1, ],
  iter = 2000, # Number of MCMC iterations
  edge_selection = FALSE, #No edge selection
  seed = 1234)

fit_no_diag = bgm(
  x = data[group_indicator == 0, ],
  iter = 2000, # Number of MCMC iterations
  edge_selection = FALSE, #No edge selection
  seed = 1234)

layout = qgraph::qgraph(
  input = fit_diag$posterior_mean_pairwise,
  DoNotPlot = TRUE,
  layout = "spring"
)$layout

global_max <- max(
  (fit_diag$posterior_mean_pairwise),
  (fit_no_diag$posterior_mean_pairwise)
)

qgraph::qgraph(
  fit_diag$posterior_mean_pairwise,
  layout = layout,
  labels = colnames(ADHD[,-1]),
  theme = "TeamFortress",
  vsize = 10,
  esize = 20,
  legend = FALSE,
  maximum = global_max
)

qgraph::qgraph(
  fit_no_diag$posterior_mean_pairwise,
  layout = layout,
  labels = colnames(ADHD[,-1]),
  theme = "TeamFortress",
  vsize = 10,
  esize = 20,
  legend = FALSE,
  maximum = global_max
)


easy_fit = easybgm_compare(
  data = data,
  group_indicator = group_indicator,
  iter = 2000,
  package = "bgms",
  type = "ordinal",
  seed = 1234
)

summary(easy_fit)

beta_diag = bgm(
  x = data[group_indicator == 1, ],
  iter = 2000, # Number of MCMC iterations
  edge_selection = FALSE, #No edge selection
  edge_prior = "Beta-Bernoulli",
  seed = 1234)

beta_no_diag <- bgm(
  x = data[group_indicator == 1, ],
  iter = 2000, # Number of MCMC iterations
  edge_selection = FALSE, #No edge selection
  edge_prior = "Beta-Bernoulli",
  seed = 1234)



# Prior on the Network Structure & Clustering
qgraph::qgraph(
  beta_diag$posterior_mean_pairwise,
  layout = layout,
  labels = colnames(ADHD[,-1]),
  theme = "TeamFortress",
  vsize = 10,
  esize = 20,
  legend = FALSE,
  maximum = global_max
)

qgraph::qgraph(
  beta_no_diag$posterior_mean_pairwise,
  layout = layout,
  labels = colnames(ADHD[,-1]),
  theme = "TeamFortress",
  vsize = 10,
  esize = 20,
  legend = FALSE,
  maximum = global_max
)

easy_beta = easybgm_compare(
  data = data,
  group_indicator = group_indicator,
  iter = 2000,
  package = "bgms",
  type = "ordinal",
  seed = 1234
)

?easybgm_compare

easy_fit = easybgm_compare(
  data = data,
  group_indicator = group_indicator,
  iter = 2000,
  package = "bgms",
  type = "ordinal",
  edge_prior = "Beta-Bernoulli",
  seed = 1234
)
