################################################################################
# simulation study, v-1
################################################################################

setwd("C:/Users/wangs/OneDrive/codes/1- combine clinical and molecules. real data/codes/random_stacking_cox/codes/")

source("./functions/cox.fun.R")
source("./functions/helper.R")
source("./functions/drpca.R")
source("./functions/stk.cox_v6.R")
source("./functions/stk.glm_v5.R")
source("./functions/random.stk_v1.R")
source("./functions/ipflasso.R")
source("./functions/simulation_help.R")

# load libraries
library(ggplot2)
library(ggpubr)
library(glmnet)
library(survival)
library(dplyr)
library(tidyverse)
library(doParallel)
library(doRNG)
library(doSNOW)

# =============================================================================
# simulation study, v-1
# ==============================================================================

combs <- expand.grid(
    snr = c(1, 2.33),
    n = c(100, 500),
    r = c(0, 0.3, 0.7)
)

set.seed(123)

for (i in 1:nrow(combs)) {
    snr <- combs[i, 1]
    n <- combs[i, 2]
    r <- combs[i, 3]

    sd_sett <- getsd_sett(snr, r = r)
    combs[i, 4] <- sd_sett
}

colnames(combs) <- c("n", "snr", "r", "sd_sett")

# ======================================================================================
# Run Simulation
# ======================================================================================

nsim <- 100
nCores <- 10

whichRuns = c("stk", "random.stk", "competitors")

for(whichRun in whichRuns){

for (ind in 1:nrow(combs)) {
    
    snr <- combs[ind, 1]
    n <- combs[ind, 2]
    r <- combs[ind, 3]
    sd_sett <- combs[ind, 4]

    cat(
        as.character(Sys.time()), "\n",
        "sample size: ", n, "\n",
        "snr:", snr, "\n",
        "r:", r, "\n",
        "\n\n"
    )

    # generate test data

    rng <- RNGseq(nsim, 6534125)
    Cores <- pmin(nCores, detectCores() - 1)
    cl <- makeCluster(Cores)
    registerDoSNOW(cl)

    pb <- txtProgressBar(min = 1, max = nsim, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    foreach(
        i = 1:nsim,
        .packages = c("glmnet", "survival"),
        .options.snow = opts,
        .verbose = TRUE
    ) %dopar% {
        rngtools::setRNG(rng[[i]])
        # generate data
        data <- simData(
            n.obs = n,
            t.gene = 10000,
            r = r,
            sd_sett = sd_sett,
            weib_shape = 1.5,
            weib_scale = 0.1,
            censoring_target = 0.2
        )

        # plot(survival::Surv(data$time, data$status))

        y <- survival::Surv(data$time, data$status)
        times <- quantile(data$time[data$status == 1], probs = seq(0.05, 0.95, length = 10))
        cox.time <- quantile(data$time[data$status == 1], probs = 0.5)

        foldid <- sample(rep(1:10, length = n), n, replace = FALSE)

        argList <- list(
            X = data$X, y = y,
            times = times,
            cox.time = cox.time,
            blocks = data$blocks,
            foldid = foldid,
            intercept = TRUE,
            foldid.internal = NULL,
            cens.method = "KM"
        )

        # fit competitors
        if (whichRun == "competitors") {
            ncomp <- nrow(data$X) - 1
            argList$ncomp <- ncomp
            mod <- do.call(fit.others, argList)
        } else if (whichRun == "stk") {
            mod <- do.call(fit.stk, argList)
        } else if (whichRun == "random.stk") {
            mod <- do.call(fit.random.stk, argList)
        }

        nam <- paste0(n, "_", snr, "_", r, "_", i)
        save(
            list = "mod",
            file = paste0("../output/sim_", whichRun, "_", nam, ".RData")
        )
    }

    stopCluster(cl)
}
}

names(mod)
