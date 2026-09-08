# =============================================================================
# test NN superLearner
# ==============================================================================

nnSL <- function(obj, epochs = 100, learning_rate = 0.001, weight_decay = 0.01) {
    require("torch")

    X <- obj$CVpred[, -c(1, 2)]
    y <- obj$CVpred[, 2]
    # normalization
    X_scaled <- scale(X) # ignore truncation for now, since we are using nn with regularization
    p <- ncol(X_scaled)

    x_tensor <- torch_tensor(X_scaled, dtype = torch_float())
    y_tensor <- torch_tensor(y, dtype = torch_float())$unsqueeze(2)
    # model
    # model <- nn_sequential(
    #     nn_linear(p, 10),
    #     nn_batch_norm1d(10),
    #     nn_relu(),
    #     nn_dropout(0.2),
    #     nn_linear(10, 5),
    #     nn_batch_norm1d(5),
    #     nn_relu(),
    #     nn_dropout(0.2),
    #     nn_linear(5, 1)
    # )

    model <- nn_sequential(

        nn_linear(p, 4),
        nn_batch_norm1d(4),
        nn_sigmoid(), 
        nn_dropout(0.2),

        # nn_linear(8, 4),
        # nn_batch_norm1d(4),
        # nn_relu(),
        # nn_dropout(0.2),

        nn_linear(4, 1) 
    )

    # loss + optimizer
    loss_fn <- nn_bce_with_logits_loss()
    optimizer <- optim_adam(model$parameters, lr = learning_rate, weight_decay = weight_decay)

    # training + loop
    model$train()

    for (epoch in 1:epochs) {
        optimizer$zero_grad()
        y_pred <- model(x_tensor)
        loss <- loss_fn(y_pred, y_tensor)
        loss$backward()
        optimizer$step()

        if (epoch %% 50 == 0) {
            cat("Epoch:", epoch, "Loss:", loss$item(), "\n")
        }
    }

    # prediction
    model$eval()
    prob <- torch_sigmoid(model(x_tensor))
    pred_class <- ifelse(as.numeric(prob) > 0.5, 1, 0)

    # output
    re <- list(
        model = model,
        scaler_mean = attr(X_scaled, "scaled:center"),
        scaler_sd = attr(X_scaled, "scaled:scale"),
        prob = prob,
        pred = pred_class
    )
    return(re)
}

predict_nnSL <- function(obj, newdata) {
    X_new <- as.matrix(newdata)
    X_new_scaled <- scale(X_new,
        center = obj$scaler_mean,
        scale  = obj$scaler_sd
    )

    x_new_tensor <- torch_tensor(X_new_scaled, dtype = torch_float())

    obj$model$eval()
    with_no_grad({
        logits <- obj$model(x_new_tensor)
        prob <- torch_sigmoid(logits)
    })

    return(prob = as.numeric(prob))
}

predict.stack.nnSL <- function(obj, newdata, times = NULL, ifCVMods = FALSE) {
    # NOTE: alpha is a list in the new version
    # loop for the method
    # for logLik and ibsLoss
    BLSurvP <- predBLSurvP(obj, newdata, times, ifCVMods = ifCVMods)
    # get times
    times <- attr(BLSurvP, "times")

    #- get predprob super model
    metaSurv <- matrix(NA,
        nrow = nrow(BLSurvP[[1]]), ncol = length(times),
        dimnames = list(
            rownames(BLSurvP),
            as.character(times)
        )
    )

    for (k in seq_along(times)) {
        Z <- mapply(subset,
            x = BLSurvP,
            MoreArgs = list(
                subset = T,
                select = k
            ), SIMPLIFY = T
        )

        metaSurv[, k] <- predict_nnSL(obj, Z)
    }

    return(metaSurv)
}

# ================================================================================
# Test on all datasets
# =================================================================================

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

whichRun <- "random.stk"
ifCVMods <- FALSE

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

        nn_model <- nnSL(mod$Rstk.comp, epochs = 500, learning_rate = 0.001, weight_decay = 0)

        # training and validation data
        X.v <- X[sam == j, ]
        y.v <- y[sam == j]

        X.t <- X[sam != j, ]
        y.t <- y[sam != j]

        BLSurvP <- predBLSurvP(mod$Rstk.comp, newdata = X.v, times = times, ifCVMods = T)
        # get times
        times <- attr(BLSurvP, "times")

        #- get predprob super model
        metaSurv <- matrix(NA,
            nrow = nrow(BLSurvP[[1]]), ncol = length(times),
            dimnames = list(
                rownames(BLSurvP),
                as.character(times)
            )
        )

        for (k in seq_along(times)) {
            Z <- mapply(subset,
                x = BLSurvP,
                MoreArgs = list(
                    subset = T,
                    select = k
                ), SIMPLIFY = T
            )

            metaSurv[, k] <- predict_nnSL(nn_model, Z)
        }

        probs <- metaSurv

        ibsTime <- getTimeIBS(probs,
            Xclin.t = X.t[, blocks == 1],
            Xclin.v = X.v[, blocks == 1],
            y.t = y.t,
            y.v = y.v,
            times = times,
            cens.method = "clinicalAdjusted",
            sumRes = FALSE
        )

        aucTime <- getTimeAUC(probs,
            Xclin.t = X.t[, blocks == 1],
            Xclin.v = X.v[, blocks == 1],
            y.t = y.t,
            y.v = y.v,
            times = times,
            cens.method = "clinicalAdjusted",
            sumRes = FALSE
        )

        logLikTime <- getTimelogLik(probs,
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


sumFun <- function(obj) {
    # objs, aucs, ibss, or logLiks
    areaSum <- function(obj) {
        obj <- na.omit(obj)
        tRange <- abs(range(obj[, "time"])[1] - range(obj[, "time"])[2])
        trapezoidal(obj[, "time"], obj[, 2]) / tRange
        # sum(obj[,2])
    }

    res <- unlist(lapply(obj, areaSum))
    res
}

unlist(lapply(lapply(aucs, sumFun), mean))
