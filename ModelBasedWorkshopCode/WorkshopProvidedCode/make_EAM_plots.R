rm(list = ls())
library(EMC2)

# This collection of functions was use to make the 4 EAM plots displayed in
# Figure 1 of the paper!


set.seed(1234)
matchfun <- function(d) d$S == d$lR


design_DDM <- design(formula = list(v ~ S, a ~ 1, t0 ~ 1),
                     Rlevels = 1:2,
                     factors = list(S = 1:2, subjects = 1:5, E = 1:3),
                     model = DDM)

# For each level of S
p_vector <- c(v=.5, v_S2 = 0,  a=log(1), t0=log(.2))

pdf(file = "DDMplot.pdf", width = 4, height = 3.5)
plot(design_DDM, p_vector, plot_legend = FALSE)
dev.off()

design_RDM <- design(formula = list(v ~ lM, B ~ 1, t0 ~ 1),
                     Rlevels = 1:2,
                     matchfun = matchfun,
                     factors = list(S = 1:2, subjects = 1:5),
                     model = RDM)
p_vector <- c(v=log(1), v_lMTRUE=log(2.2),
              B=log(.8), t0=log(.15))

pdf(file = "RDMplot.pdf", width = 4, height = 3.5)
plot(design_RDM, p_vector, factors = list(v = "lM"), plot_legend = FALSE)
dev.off()

design_LBA <- design(formula = list(v ~ lM, B ~ 1, t0 ~ 1),
                     Rlevels = 1:2,
                     matchfun = matchfun,
                     factors = list(S = 1:2, subjects = 1:5),
                     model = LBA)


p_vector <- c(v=.7, v_lMTRUE=1,
              B=log(1), t0=log(.15))

pdf(file = "LBAplot.pdf", width = 4, height = 3.5)
plot(design_LBA, p_vector, factors = list(v = "lM"), plot_legend = FALSE)
dev.off()

design_LNR <- design(formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                     Rlevels = 1:2,
                     matchfun = matchfun,
                     factors = list(S = 1:2, subjects = 1:5),
                     model = LNR)

p_vector <- c(m=-.8,m_lMTRUE=-.5, s = .5, t0=log(.15))


pdf(file = "LNRplot.pdf", width = 4, height = 3.5)
plot(design_LNR, p_vector, factors = list(m = "lM"), plot_legend = FALSE)
dev.off()




