
#-------------------------------------------------------------------------------
# source("../simulation/codes/ipflasso.R")
#-------------------------------------------------------------------------------

fit.others <- function(X, y,
                       ncomp = NULL,
                       blocks = rep(1, ncol(X)),
                       foldid,
                       times,
                       cox.time,
                       intercept = T,
                       foldid.internal = NULL,
                       ...) {
  require(glmnet)

  t1 <- Sys.time()
  s <- blocks == 1
  cvs <- cv.glmnet(X[, s], y, family = "cox", foldid = foldid)
  clin <- glmnet(X[, s], y, family = "cox", lambda = cvs$lambda.min)
  t2 <- Sys.time()

  beta <- coef(clin)[, 1]
  beta <- beta[beta != 0]
  lp <- getLp(X, beta)

  clin <- list(
    beta = beta,
    lp = lp,
    y = y,
    time = t2 - t1
  )

  #- molecular model
  t1 <- Sys.time()
  s <- blocks == 2
  cvs <- cv.glmnet(X[, s], y, family = "cox", foldid = foldid)
  mol <- glmnet(X[, s], y, family = "cox", lambda = cvs$lambda.min)
  t2 <- Sys.time()

  beta <- coef(mol)[, 1]
  beta <- beta[beta != 0]
  lp <- getLp(X, beta)

  mol <- list(
    beta = beta,
    lp = lp,
    y = y,
    time = t2 - t1
  )

  #- naive
  t1 <- Sys.time()
  cvs <- cv.glmnet(X, y, family = "cox", foldid = foldid)
  naive <- glmnet(X, y, family = "cox", lambda = cvs$lambda.min)
  t2 <- Sys.time()

  lp <- predict(naive, X)
  beta <- coef(naive)[, 1]
  beta <- beta[beta != 0]
  lp <- getLp(X, beta)

  naive <- list(
    beta = beta,
    lp = lp,
    y = y,
    time = t2 - t1
  )

  #- Dimension reduction ---
  #- PCA
  if (is.null(ncomp)) {
    ncomp <- nrow(X)
  }

  t1 <- Sys.time()
  cv.pcr <- cv.drpcr(X, y,
    alpha = 1,
    family = "cox",
    ncomp = ncomp,
    blocks = blocks,
    blocks.pc = 2,
    penalty.factors = list(c(1, 1)),
    foldid = foldid
  )

  pcr.las <- drpcr(X, y,
    alpha = 1,
    ncomp = ncomp,
    family = "cox",
    blocks = blocks,
    blocks.pc = 2,
    pf = cv.pcr$pf.min,
    lambda = cv.pcr$lambda.min
  )
  t2 <- Sys.time()

  beta <- pcr.las$back.mod$beta
  lp <- getLp(X, beta)

  pcr.las <- list(
    beta = beta,
    lp = lp,
    y = y,
    beta.pcr = c(pcr.las$a0, pcr.las$beta[, 1]),
    centers = pcr.las$centers,
    scales = pcr.las$scales,
    blocks = pcr.las$blocks,
    time = t2 - t1
  )

  # - Component-wise weighting ----
  t1 <- Sys.time()
  cv.ipf <- cv.ipflas(X, y,
    family = "cox",
    foldid = foldid,
    blocks = blocks,
    alpha = 1,
    penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8))
  )

  ipf.las <- ipflas(X, y,
    family = "cox",
    alpha = 1,
    blocks = blocks,
    pf = cv.ipf$pf.min,
    lambda = cv.ipf$lambda.min
  )
  t2 <- Sys.time()

  beta <- coef(ipf.las)[, 1]
  beta <- beta[beta != 0]
  lp <- getLp(X, beta)

  ipf.las <- list(
    beta = beta,
    lp = lp,
    y = y,
    pf = cv.ipf$pf.min,
    time = t2 - t1
  )

  list(
    clin = clin,
    mol = mol,
    naive = naive,
    pcr.las = pcr.las,
    ipf.las = ipf.las
  )
}


# stacking ---------------------------------------------------------------------

fit.stk <- function(X, y,
                    blocks = rep(1, ncol(X)),
                    foldid,
                    times,
                    cox.time,
                    intercept = TRUE,
                    foldid.internal = NULL,
                    cens.method = "clinicalAdjusted",
                    ...) {

  require(glmnet)

  #- clinical models, for validation
  t1 <- Sys.time()
  s <- blocks == 1
  cvs <- cv.glmnet(X[, s], y, family = "cox", foldid = foldid)
  clin <- glmnet(X[, s], y, family = "cox", lambda = cvs$lambda.min)
  t2 <- Sys.time()

  beta <- coef(clin)[, 1]
  beta <- beta[beta != 0]
  lp <- getLp(X, beta)

  clin <- list(
    beta = beta,
    lp = lp,
    y = y,
    time = t2 - t1
  )

  # # stacking
  # t1 <- Sys.time()
  # stk.pf <- stack.cox(X, y,
  #   times = times,
  #   foldid = foldid,
  #   foldid.internal = foldid.internal,
  #   nfold = 10,
  #   blocks = blocks,
  #   block.clin = 1,
  #   ifCox = TRUE,
  #   cox.time = cox.time,
  #   method = "bypf",
  #   cens.method = cens.method,
  #   penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
  #   track = TRUE,
  #   optimLoss = c("ibsLoss", "logLik", "PaLogLik"),
  #   optimFun = c("ridge"),
  #   moreOptimFun = NULL,
  #   intercept = TRUE,
  #   foldidSL = NULL,
  #   foldidSL.cox = NULL,
  #   nfoldSL = 10,
  #   limits = 0,
  #   ifWeights = TRUE,
  #   truncValue = 0.01,
  #   transPFun = function(x) log(x / (1 - x))
  # )
  # t2 <- Sys.time()
  # stk.pf$time <- t2 - t1

  #- 2
  t1 <- Sys.time()
  stk.comp <- stack.cox(X, y,
    times = times,
    foldid = foldid,
    foldid.internal = foldid.internal,
    nfold = 10,
    blocks = blocks,
    block.clin = 1,
    ifCox = TRUE,
    cox.time = cox.time,
    method = "bycomp",
    cens.method = cens.method,
    penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
    track = TRUE,
    optimLoss = c("ibsLoss", "logLik", "PaLogLik"),
    optimFun = c("lasso", "ridge"),
    moreOptimFun = NULL,
    intercept = TRUE,
    foldidSL = NULL,
    foldidSL.cox = NULL,
    nfoldSL = 10,
    limits = 0,
    ifWeights = TRUE,
    truncValue = 0.01,
    transPFun = function(x) log(x / (1 - x))
  )

  t2 <- Sys.time()
  stk.comp$time <- t2 - t1

  list(
    clin = clin,
    # stk.pf = stk.pf,
    stk.comp = stk.comp
  )
}

#
# fit.stkPF.pro <- function(X, y,
#                           blocks = rep(1, ncol(X)),
#                           foldid,
#                           times,
#                           cox.time,
#                           intercept = TRUE,
#                           foldid.internal = NULL, ...) {
#   require(glmnet)

#   #- clinical models, for validation
#   t1 <- Sys.time()
#   s <- blocks == 1
#   cvs <- cv.glmnet(X[, s], y, family = "cox", foldid = foldid)
#   clin <- glmnet(X[, s], y, family = "cox", lambda = cvs$lambda.min)

#   beta <- coef(clin)[, 1]
#   beta <- beta[beta != 0]
#   lp <- getLp(X, beta)

#   t2 <- Sys.time()

#   clin <- list(
#     beta = beta,
#     lp = lp,
#     y = y,
#     time = t2 - t1
#   )


#   #- subset stacking
#   t1 <- Sys.time()

#   indClin <- blocks == 1
#   blocks[!indClin] <- sort(rep(2:11, length = sum(!indClin)))

#   for (k in 2:11) {
#     ind <- blocks %in% c(1, k)
#     b <- as.numeric(blocks[ind] != 1) + 1
#     stkTmpt <- StkCoxBL(X[, ind], y,
#       times = times,
#       foldid = foldid,
#       foldid.internal = foldid.internal,
#       blocks = b,
#       block.clin = 1,
#       nfold = 10,
#       ifCox = TRUE,
#       cox.time = cox.time,
#       method = "bypf",
#       cens.method = "clinicalAdjusted",
#       penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
#       track = TRUE
#     )

#     if (k == 2) {
#       stk.subset <- stkTmpt
#     } else {
#       stk.subset <- merge.stk_cox(stk.subset, stkTmpt)
#     }
#   }
#   #
#   subNames <- paste0("sub", rep(1:10, each = 4), "_", c("1-1", "1-2", "1-4", "1-8"))

#   names(stk.subset$subMod$beta) <- subNames
#   names(stk.subset$subMod$lp) <- subNames
#   colnames(stk.subset$CVpredCox) <- subNames
#   colnames(stk.subset$CVpred) <- c("weight", "Zind", subNames)


#   for (j in seq_along(stk.subset$CVMods)) {
#     names(stk.subset$CVMods[[j]]$mods) <- subNames
#   }

#   # super learner
#   stk.subset <- StkCoxSL(stk.subset,
#     optimLoss = c("ibsLoss", "logLik", "PaLogLik"),
#     optimFun = c("ridge"),
#     moreOptimFun = NULL,
#     intercept = TRUE,
#     foldid = NULL, foldid.cox = NULL, nfold = 10,
#     limits = 0,
#     ifWeights = c(TRUE, FALSE, FALSE),
#     truncValue = 0.01,
#     transPFun = function(x) log(x / (1 - x))
#   )

#   t2 <- Sys.time()
#   stk.subset$time <- t2 - t1

#   # - relaxed
#   fast.relax.mod <- function(X, y,
#                              family,
#                              pfs = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
#                              blocks = rep(1, ncol(X)),
#                              foldid = NULL,
#                              ...) {
#     # lasso model
#     modList <- list()

#     for (pf in pfs) {
#       # get name for list
#       name <- paste0(pf, collapse = "-")

#       pf <- pf[as.numeric(factor(blocks))]
#       cvs <- cv.glmnet(X, y,
#         family = family,
#         penalty.factor = pf,
#         foldid = foldid, ...
#       )

#       las <- glmnet(X, y,
#         family = family,
#         lambda = cvs$lambda.min,
#         penalty.factor = pf, ...
#       )

#       beta <- las$beta[, 1]
#       selected_vars <- which(beta != 0)

#       if (length(selected_vars) < 2) {
#         # random add one
#         selected_vars <- c(selected_vars, sample(setdiff(1:ncol(X), selected_vars), 2 - length(selected_vars)))
#       }

#       las0 <- glmnet(X[, selected_vars, drop = FALSE], y,
#         family = family,
#         lambda = 0
#       )

#       # add lp
#       las0$lp <- predict(las0, X[, selected_vars, drop = FALSE])
#       las$lp <- predict(las, X)

#       modList[[paste0("las-", name)]] <- las
#       modList[[paste0("relax-", name)]] <- las0
#     }

#     return(modList)
#   }

#   t1 <- Sys.time()
#   # fit sub-models
#   stk.pf.relax <- stack.cox(X, y,
#     times = times,
#     foldid = foldid,
#     foldid.internal = foldid.internal,
#     blocks = blocks,
#     block.clin = 1,
#     nfold = 10,
#     ifCox = TRUE,
#     cox.time = cox.time,
#     method = function(...) fast.relax.mod(..., family = "cox", blocks = blocks, pfs = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8))),
#     cens.method = "clinicalAdjusted",
#     penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
#     track = TRUE,
#     # super learner argument
#     optimLoss = c("ibsLoss", "logLik", "PaLogLik"),
#     optimFun = "ridge",
#     moreOptimFun = NULL,
#     limits = 0,
#     ifWeights = TRUE,
#     truncValue = 0.01,
#     transPFun = function(x) log(x / (1 - x))
#   )

#   t2 <- Sys.time()
#   stk.pf.relax$time <- t2 - t1

#   list(
#     clin = clin,
#     stk.subset = stk.subset,
#     stk.pf.relax = stk.pf.relax
#   )
# }

#-------------------------------------------------------------------------------


fit.random.stk <- function(X, y,
                           blocks = rep(1, ncol(X)),
                           foldid,
                           times,
                           cox.time,
                           intercept = TRUE,
                           cens.method = "clinicalAdjusted",
                           ...) {
  #
  require(glmnet)

  #- clinical models, for validation
  t1 <- Sys.time()
  s <- blocks == 1
  cvs <- cv.glmnet(X[, s], y, family = "cox", foldid = foldid)
  clin <- glmnet(X[, s], y, family = "cox", lambda = cvs$lambda.min)

  beta <- coef(clin)[, 1]
  beta <- beta[beta != 0]
  lp <- getLp(X, beta)

  t2 <- Sys.time()

  clin <- list(
    beta = beta,
    lp = lp,
    y = y,
    time = t2 - t1
  )

  # random stacking comp
  t1 <- Sys.time()
  Rstk.comp <- random.stack.cox(X, y,
    times = times,
    bootID = NULL,
    randomP = c(table(blocks)[1], nrow(X)),
    nboot = 100,
    fastTune = 20,
    foldid = foldid,
    nfold = 10,
    blocks = blocks,
    block.clin = 1,
    ifCox = TRUE,
    cox.time = cox.time,
    method = "bycomp",
    cens.method = cens.method,
    penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
    optimLoss = c("ibsLoss", "logLik", "PaLogLik"),
    optimFun = "ridge",
    moreOptimFun = NULL,
    intercept = intercept,
    foldidSL = NULL,
    foldidSL.cox = NULL,
    nfoldSL = 10,
    limits = 0,
    ifWeights = c(TRUE, FALSE, FALSE),
    truncValue = 0.01,
    transPFun = function(x) log(x / (1 - x))
  )
  t2 <- Sys.time()
  Rstk.comp$time <- t2 - t1


  # # random stacking pf
  # t1 <- Sys.time()
  # Rstk.pf <- random.stack.cox(X, y,
  #   times = times,
  #   bootID = NULL,
  #   randomP = c(table(blocks)[1], nrow(X)),
  #   nboot = 100,
  #   fastTune = 20,
  #   foldid = foldid,
  #   nfold = 10,
  #   blocks = blocks,
  #   block.clin = 1,
  #   ifCox = TRUE,
  #   cox.time = cox.time,
  #   method = "bypf",
  #   cens.method = cens.method,
  #   penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
  #   optimLoss = c("ibsLoss", "logLik", "PaLogLik"),
  #   optimFun = "ridge",
  #   moreOptimFun = NULL,
  #   intercept = TRUE,
  #   foldidSL = NULL,
  #   foldidSL.cox = NULL,
  #   nfoldSL = 10,
  #   limits = 0,
  #   ifWeights = c(TRUE, FALSE, FALSE),
  #   truncValue = 0.01,
  #   transPFun = function(x) log(x / (1 - x))
  # )
  # t2 <- Sys.time()
  # Rstk.pf$time <- t2 - t1

  # output
  list(
    clin = clin,
    Rstk.comp = Rstk.comp
    # Rstk.pf = Rstk.pf
  )
}
