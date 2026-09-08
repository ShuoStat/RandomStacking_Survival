# =============================================================================
# test NN superLearner
# ==============================================================================

nnSL <- function(
    X, y,
    epochs = 100,
    learning_rate = 0.001,
    weight_decay = 0.01,
    hidden_units = 4,
    dropout = 0.2,
    verbose = FALSE) {
    require(torch)

    X_scaled <- scale(X)
    scaler_mean <- attr(X_scaled, "scaled:center")
    scaler_sd <- attr(X_scaled, "scaled:scale")

    ## Avoid zero variance variables
    scaler_sd[scaler_sd == 0 | is.na(scaler_sd)] <- 1
    X_scaled <- sweep(X, 2, scaler_mean, "-")
    X_scaled <- sweep(X_scaled, 2, scaler_sd, "/")
    p <- ncol(X_scaled)

    ## 3. Convert to torch tensor
    x_tensor <- torch_tensor(X_scaled, dtype = torch_float())
    y_tensor <- torch_tensor(
        matrix(y, ncol = 1),
        dtype = torch_float()
    )

    ## 4. Neural network
    model <- nn_sequential(
        nn_linear(p, hidden_units),
        nn_batch_norm1d(hidden_units),
        nn_sigmoid(),
        nn_dropout(dropout),
        nn_linear(hidden_units, 1)
    )

    ## 5. Loss + optimizer
    loss_fn <- nn_bce_with_logits_loss()
    optimizer <- optim_adam(
        model$parameters,
        lr = learning_rate,
        weight_decay = weight_decay
    )

    ## 6. Training
    model$train()
    for (epoch in seq_len(epochs)) {
        optimizer$zero_grad()
        logits <- model(x_tensor)
        loss <- loss_fn(logits, y_tensor)
        loss$backward()
        optimizer$step()
        if (verbose && epoch %% 50 == 0) {
            cat(
                "Epoch:", epoch,
                "Loss:", loss$item(),
                "\n"
            )
        }
    }

    ## 7. Prediction
    model$eval()
    with_no_grad({
        logits <- model(x_tensor)
        prob <- torch_sigmoid(logits)
    })

    prob <- as.numeric(prob)
    pred_class <- ifelse(prob > 0.5, 1, 0)

    ## 8. Return
    re <- list(
        model = model,
        scaler_mean = scaler_mean,
        scaler_sd = scaler_sd,
        prob = prob,
        pred = pred_class,
        epochs = epochs,
        learning_rate = learning_rate,
        weight_decay = weight_decay,
        hidden_units = hidden_units,
        dropout = dropout
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


nnSL.cv <- function(
    X, y,
    epochs = c(100, 200),
    learning_rate = c(0.001, 0.01),
    weight_decay = c(0, 0.001, 0.01),
    hidden_units = c(2, 4, 8),
    dropout = c(0, 0.2),
    nfolds = 5,
    seed = 123,
    verbose = FALSE) {
    require(torch)
    set.seed(seed)

    X <- as.matrix(X)
    y <- as.numeric(y)
    n <- nrow(X)

    ## 1. Create CV folds
    foldid <- sample(rep(seq_len(nfolds), length.out = n))

    ## 2. Tuning grid
    tuning_grid <- expand.grid(
        epochs = epochs,
        learning_rate = learning_rate,
        weight_decay = weight_decay,
        hidden_units = hidden_units,
        dropout = dropout
    )

    n_grid <- nrow(tuning_grid)
    cv_nll <- numeric(n_grid)

    ## 3. Cross-validation
    for (g in seq_len(n_grid)) {
        if (verbose) {
            cat("\nCombination", g, "/", n_grid, "\n")
            print(tuning_grid[g, ])
        }

        fold_nll <- numeric(nfolds)
        for (k in seq_len(nfolds)) {
            ## Training / validation indices
            train_idx <- which(foldid != k)
            valid_idx <- which(foldid == k)
            ## Train neural network
            fit <- nnSL(
                X = X[train_idx, , drop = FALSE],
                y = y[train_idx],
                epochs = tuning_grid$epochs[g],
                learning_rate = tuning_grid$learning_rate[g],
                weight_decay = tuning_grid$weight_decay[g],
                hidden_units = tuning_grid$hidden_units[g],
                dropout = tuning_grid$dropout[g],
                verbose = FALSE
            )

            ## Validation prediction
            prob <- predict_nnSL(
                obj = fit,
                newdata = X[valid_idx, , drop = FALSE]
            )

            ## Validation NLL
            eps <- 1e-7
            prob <- pmin(pmax(prob, eps), 1 - eps)
            fold_nll[k] <- -mean(y[valid_idx] * log(prob) + (1 - y[valid_idx]) * log(1 - prob))
        }

        ## Mean CV NLL
        cv_nll[g] <- mean(fold_nll)

        if (verbose) {
            cat(
                "CV NLL:",
                round(cv_nll[g], 6),
                "\n"
            )
        }
    }

    ## 4. CV results
    cv_results <- cbind(tuning_grid, CV_NLL = cv_nll)

    ## Best tuning parameters
    best_id <- which.min(cv_nll)
    best_params <- tuning_grid[best_id, , drop = FALSE]

    ## 6. Return
    re <- list(
        epochs = best_params$epochs,
        learning_rate = best_params$learning_rate,
        weight_decay = best_params$weight_decay,
        hidden_units = best_params$hidden_units,
        dropout = best_params$dropout,
        cv_results = cv_results,
        foldid = foldid,
        verbose = verbose
    )

    return(re)
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

###
ibss <- list()
aucs <- list()
logLiks <- list()
nCores <- 20

# set.seed(321)
# datNames <- c("BLCA", "HNSC")

# to get the tuning parameters for the neural network super learner, we can use the following code:
# use top10 replicates to approximate the best tuning parameters

tuning_params_list <- list()

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

    pb <- txtProgressBar(min = 1, max = 10, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    res <- foreach(
        i = 1:10,
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
        obj <- mod$Rstk.comp

        # Cross-validation for neural network super learner
        nn_model <- nnSL.cv(
            X = as.matrix(obj$CVpred[, -c(1, 2)]), y = as.numeric(obj$CVpred[, 2]),
            epochs = c(500),
            learning_rate = c(0.001, 0.01),
            weight_decay = c(0, 0.01),
            hidden_units = c(2, 4, 8),
            dropout = c(0, 0.2),
            nfolds = 5,
            seed = 123,
            verbose = FALSE
        )

        # just return the tuning parameters, not the model itself
        tuning_params <- setNames(c(
            nn_model$epochs,
            nn_model$learning_rate,
            nn_model$weight_decay,
            nn_model$hidden_units,
            nn_model$dropout
        ), c("epochs", "learning_rate", "weight_decay", "hidden_units", "dropout"))

        return(tuning_params)
    }

    tuning_params_list[[nam]] <- res
    stopCluster(cl)
}

save(tuning_params_list, file = paste0("../output/", whichRun, "-tuning_params_list.RData"))

# get the best tuning parameters for each dataset
load(paste0("../output/", whichRun, "-tuning_params_list.RData"))

getOptim <- function(obj) {

    tuning_df <- as.data.frame(do.call(rbind, obj))

    param_freq <- lapply(names(tuning_df), function(p) {
        tab <- table(tuning_df[[p]])
        data.frame(
            parameter = p,
            value = names(tab),
            freq = as.integer(tab)
        ) %>% arrange(desc(freq))
    })
    
    optim_params <- lapply(param_freq, function(x) {
        x[1, 1:2]
    })

    Reduce(rbind, optim_params)
}

optimal_tuning_params <- lapply(tuning_params_list, getOptim)


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
        obj <- mod$Rstk.comp

        # Cross-validation for neural network super learner

        nn_model <- nnSL(
            X = as.matrix(obj$CVpred[, -c(1, 2)]),
            y = as.numeric(obj$CVpred[, 2]),
            epochs = as.numeric(optimal_tuning_params[[nam]][1, "value"]),
            learning_rate = as.numeric(optimal_tuning_params[[nam]][2, "value"]),
            hidden_units = as.numeric(optimal_tuning_params[[nam]][4, "value"]),
            dropout = as.numeric(optimal_tuning_params[[nam]][5, "value"]),
            weight_decay = as.numeric(optimal_tuning_params[[nam]][3, "value"])
        )

        # training and validation data
        X.v <- X[sam == j, ]
        y.v <- y[sam == j]

        X.t <- X[sam != j, ]
        y.t <- y[sam != j]

        BLSurvP <- predBLSurvP(obj, newdata = X.v, times = times, ifCVMods = T)
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
