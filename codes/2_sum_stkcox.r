setwd("C:/Users/wangs/OneDrive/codes/1- combine clinical and molecules. real data/codes/random_stacking_cox/codes/")

source("./functions/cox.fun.R")
source("./functions/helper.R")
source("./functions/drpca.R")
source("./functions/stk.cox_v6.R")
source("./functions/stk.glm_v5.R")
source("./functions/random.stk_v1.R")
source("./functions/ipflasso.R")

# load libraries
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(glmnet)
library(survival)
library(dplyr)
library(tidyverse)
library(doParallel)
library(doRNG)
library(doSNOW)

paraList <- list(
  nsim = 5,
  nfolds = 10,
  nCores = 10,
  datNames = c(
    "BLCA", "BRCA", "HNSC", "KIRC", "LAML",
    "LGG", "LIHC", "LUAD", "PAAD", "SKCM"
  )
)
nsim <- paraList$nsim
nfolds <- paraList$nfolds
nCores <- paraList$nCores
datNames <- paraList$datNames

# Summarize the results

#-------------------------------------------------------------------------------
#- prediction
#-------------------------------------------------------------------------------

# whichRun should be one of c("competitors", "stk", "random.stk")

whichRun <- "competitors"
ifCVMods <- FALSE
updateResults <- TRUE

###
ibss <- list()
aucs <- list()
logLiks <- list()
nCores <- 10

# set.seed(321)
# datNames <- c("BLCA", "HNSC")

for (nam in datNames) {
  print(nam)
  getData(nam, log2 = T, toBinary = F)

  n <- nrow(X)
  p <- ncol(X)

  # # try
  # X <- X[,1:500]
  # blocks <- blocks[1:500]

  #- get random splits, foldid, times
  getSplits(n, nsim, nfolds, seed = 1248372, foldid.internal = T)
  time <- y[, 1]
  status <- y[, 2]
  # times <- quantile(time[status == 1], probs = seq(0.1, 0.9, length = 10))
  timesRange <- quantile(time[status == 1], probs = c(0.1, 0.9))
  times <- seq(timesRange[1], timesRange[2], length = 20)

  Cores <- pmin(nCores, detectCores() - 1)
  cl <- makeCluster(Cores)
  registerDoSNOW(cl)

  pb <- txtProgressBar(min = 1, max = (nfolds * nsim), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  res <- foreach(
    i = 1:(nfolds * nsim),
    .packages = c("dplyr", "glmnet", "survival"),
    .options.snow = opts,
    .verbose = TRUE
  ) %dopar% {
    samID <- ceiling(i / nfolds)
    # j, which fold of the i.nsim split being the test
    j <- i - (samID - 1) * nfolds
    # foldid for i.nsim splits and j fold
    sam <- sampleSplits[[samID]]
    load(paste0("../output/", whichRun, "-", nam, "-", samID, "-", j, ".RData"))

    allStks <- grep("stk", names(mod), value = TRUE)

    # for(modName in allStks) {
    #
    #   obj = mod[[modName]]
    #   class(obj) <- "stk"
    #   mod[[modName]] <-
    #     s <- StkCoxSL(obj,
    #              optimLoss = c("ibsLoss", "logLik", "PaLogLik"),
    #              optimFun = c("ridge"),
    #              moreOptimFun = NULL,
    #              intercept = TRUE,
    #              foldid = NULL, foldid.cox = NULL, nfold = 10,
    #              limits = 0,
    #              ifWeights = c(TRUE, FALSE, FALSE),
    #              truncValue = 0.01,
    #              transPFun = function(x) log(x / (1 - x)))
    # }

    if (length(mod) == 1) {
      mod <- list(mod)
    }

    # training and validation data
    X.v <- X[sam == j, ]
    y.v <- y[sam == j]

    X.t <- X[sam != j, ]
    y.t <- y[sam != j]
    probs <- predProbs(X.v, mod, times, ifCVMods = ifCVMods)

    # # faltten list
    flattenList <- function(obj, target) {
      new <- obj[[target]]
      names(new) <- paste0(target, ":", names(new))
      obj <- append(obj, new)
      obj[[target]] <- NULL
      obj
    }

    for (k in names(probs)[grep("stk", names(probs))]) {
      probs <- flattenList(probs, target = k)
    }

    # get Time
    ibsTime <- lapply(probs, getTimeIBS,
      Xclin.t = X.t[, blocks == 1],
      Xclin.v = X.v[, blocks == 1],
      y.t = y.t,
      y.v = y.v,
      times = times,
      cens.method = "clinicalAdjusted",
      sumRes = FALSE
    )

    aucTime <- lapply(probs, getTimeAUC,
      Xclin.t = X.t[, blocks == 1],
      Xclin.v = X.v[, blocks == 1],
      y.t = y.t,
      y.v = y.v,
      times = times,
      cens.method = "equal",
      sumRes = FALSE
    )

    logLikTime <- lapply(probs, getTimelogLik,
      y.v = y.v,
      times = times,
      sumRes = FALSE
    )

    return(list(
      ibsTime = ibsTime,
      aucTime = aucTime,
      logLikTime = logLikTime
    ))
  }
  # from res to ibss, aucs, logLiks
  ibss[[nam]] <- lapply(res, `[[`, "ibsTime")
  aucs[[nam]] <- lapply(res, `[[`, "aucTime")
  logLiks[[nam]] <- lapply(res, `[[`, "logLikTime")

  stopCluster(cl)
}

if (updateResults) {
  save(list = c("aucs", "ibss", "logLiks"), file = paste0("../output/survPredRes", whichRun, "-", ifelse(ifCVMods, "T", "F"), ".RData"))
}

# load(paste0("../output/survPredRes", whichRun, "-V26", ifelse(ifCVMods, "T", "F"), ".RData"))

#-------------------------------------------------------------------------------
# Summarize the prediction
#-------------------------------------------------------------------------------

# output Results to Tables
# reNames <- c("Clin",
#              "Stk(SS,IBS)", "Stk(SS,Rdg)", "Stk(SS,CoxRdg)",
#              "Stk(Rlx,IBS)", "Stk(Rlx,Rdg)", "Stk(Rlx,CoxRdg)")


# function to summarize results

sumFun <- function(objs, fun = colMeans) {
  # objs, aucs, ibss, or logLiks
  areaSum <- function(obj) {
    obj <- na.omit(obj)
    tRange <- abs(range(obj[, "time"])[1] - range(obj[, "time"])[2])
    trapezoidal(obj[, "time"], obj[, 2]) / tRange
    # sum(obj[,2])
  }

  sumRes <- function(obj, fun) {
    res <- lapply(obj, function(x) unlist(lapply(x, areaSum)))
    res <- Reduce(rbind, res)
    do.call(fun, list(res))
  }

  tab <- lapply(objs, sumRes, fun = fun)
  tab <- Reduce(rbind, tab)
  rownames(tab) <- names(objs)
  return(tab)
}

round(sumFun(aucs, fun = colMeans), 4)
round(sumFun(ibss, fun = colMeans), 4)

sumResults <- function(target, ifCVMods, whichRun, ifSE = FALSE) {
  # AUC
  reMat <- c()
  reMatSe <- c()

  for (i in seq_along(whichRun)) {
    load(paste0("../output/survPredRes", whichRun[i], "-", ifelse(ifCVMods, "T", "F"), ".RData"))
    n <- length(get(target)[[1]])

    # remove duplicated clin
    tmp <- sumFun(get(target), fun = colMeans)
    if (i != 1) {
      tmp <- tmp[, -1]
    }

    reMat <- cbind(reMat, tmp)

    if (ifSE) {
      tmp <- sumFun(get(target), fun = function(x) apply(x, 2, sd) / sqrt(n))
      if (i != 1) {
        tmp <- tmp[, -1]
      }

      reMatSe <- cbind(reMatSe, tmp)
    }
  }

  # remove duplicated clin model

  re <- as.data.frame(reMat) %>%
    # select(clin, mol, naive, fs.las, pcr.las, ipf.las, `stk.pf:ibsLoss`, `stk.pf:logLik:ridge`, `stk.pf:PaLogLik:ridge`, `stk.comp:ibsLoss`,
    #        `stk.comp:logLik:ridge`, `stk.comp:PaLogLik:ridge`) %>%
    # rename_all( ~ c("Clin", "Mol", "Naive", "FS(las)", "PCR(las)", "IPF(las)", "Stk(PF,IBS)", "Stk(PF,Lik)", "Stk(PF,PaLik)",
    #                 "Stk(CP,IBS)", "Stk(CP,Lik)", "Stk(CP,PaLik)")) %>%
    mutate(across(everything(), ~ formatC(.x, format = "f", digits = 3)))

  # if (ifSE) {
  # re.se <- as.data.frame(reMatSe) %>%
  #     select(clin, mol, naive, fs.las, pcr.las, ipf.las, `stk.pf:ibsLoss`, `stk.pf:logLik:ridge`, `stk.pf:PaLogLik:ridge`, `stk.comp:ibsLoss`,
  # #            `stk.comp:logLik:ridge`, `stk.comp:PaLogLik:ridge`) %>%
  # #     rename_all( ~ c("Clin", "Mol", "Naive", "FS(las)", "PCR(las)", "IPF(las)", "Stk(PF,IBS)", "Stk(PF,Lik)", "Stk(PF,PaLik)",
  # #                     "Stk(CP,IBS)", "Stk(CP,Lik)", "Stk(CP,PaLik)")) %>%
  # mutate(across(everything(), ~ formatC(.x, format = "f", digits = 3)))
  #
  #   re <- mapply(function(mean, se) {
  #     paste0(mean, "(", se, ")")
  #   }, mean = re, se = re.se)
  #   rownames(re) <- rownames(re.se)
  # }

  # rowname data to first col
  re <- tibble::rownames_to_column(as.data.frame(re), var = "Data")
  # openxlsx::write.xlsx(re, file = paste0("../results/", target, "V1.xlsx"), colNames = T)
  write.csv(re, file = paste0("../results/real_", target, ".csv"), row.names = F)

  # print(re)
  return(re)
}

sumResults(
  target = "aucs", ifCVMods = FALSE,
  whichRun = c("competitors", "stk", "random.stk"), ifSE = FALSE
)

sumResults(
  target = "ibss", ifCVMods = FALSE,
  whichRun = c("competitors", "stk", "random.stk"), ifSE = FALSE
)

sumResults(
  target = "logLiks", ifCVMods = FALSE,
  whichRun = c("competitors", "stk", "random.stk"), ifSE = FALSE
)
