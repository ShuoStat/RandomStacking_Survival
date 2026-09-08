
################################################################################
# cox-v3, From cox-v2
################################################################################

# V26
# updates, bootstrap used in stacking
# using fixed lambda

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

#============================================================================== 
# Parameters
#============================================================================== 

paraList <- list(nsim   = 5,
                 nfolds = 10,
                 nCores = 25,
                 datNames = c("BLCA", "BRCA", "HNSC", "KIRC", "LAML", 
                              "LGG",  "LIHC", "LUAD", "PAAD", "SKCM"))
nsim   <- paraList$nsim
nfolds <- paraList$nfolds
nCores <- paraList$nCores
datNames <- paraList$datNames

#-------------------------------------------------------------------------------
#- Get clinical information
#-------------------------------------------------------------------------------

sumDat <- c()

for(nam in datNames){

  cat(nam, ":", as.character(Sys.time()), "\n\n")
  getData(nam, log2 = T, toBinary = F)

  print(sum(y[,2]))
  # get N clin
  N_clin <- sum(grepl("^clin_", colnames(X)))
  N_mol  <- sum(grepl("^m", colnames(X)))
  
  # get sampel size
  N <- nrow(X)
  
  # get N events
  N_event <- sum(y[,2])
  
  # get median survival time
  surv_fit <- survfit(y ~ 1)
  Meidan_S <- round(quantile(surv_fit$time, 0.5))
  
  inF <- setNames(c(nam, N, N_event, N_clin, N_mol, Meidan_S),
                  c("Data", "Sample Size", "Event", "N clinical variables", 
                    "N molecular variables", "Median Survival"))

  sumDat <- rbind(sumDat, inF)
}

print(sumDat)
write.csv(sumDat, file ="../results/basicInfo.csv")

#-------------------------------------------------------------------------------

# whichRun should be one of c("competitors", "stk", "stkPF.pro", "random.stk")
# "stkPF.pro" is not considered here

whichRun = "stk" 

# datNames = c("HNSC", "KIRC", "LAML", "LGG",  "LIHC", "LUAD", "PAAD", "SKCM")
# datNames = c("KIRC", "LAML", "LGG",  "LIHC", "LUAD", "PAAD", "SKCM")

for(nam in datNames){
  
  cat(nam, ":", as.character(Sys.time()), "\n\n")
  getData(nam, log2 = T, toBinary = F)
  
  n <- nrow(X)
  p <- ncol(X)
  
  getSplits(n, nsim, nfolds, seed = 1248372, foldid.internal = T)
  time   <- y[,1]
  status <- y[,2]
  
  times <- quantile(time[status == 1], probs = seq(0.1, 0.9, length = 10))
  cox.time <- quantile(time[status == 1], probs = 0.5)
  
  # # # try
#   X <- X[,1:300]
#   blocks <- blocks[1:300]
  
  Cores <- pmin(nCores, detectCores() - 1)
  cl <- makeCluster(Cores)
  registerDoSNOW(cl)
  
  rng <- RNGseq(nfolds * nsim, 6534125)
  pb <- txtProgressBar(min = 1, max = (nfolds * nsim), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  foreach(i = 1 : (nfolds * nsim),
          .packages = c("glmnet", "survival"),
          .options.snow = opts,
          .verbose = TRUE) %dopar% {
            
            rngtools::setRNG(rng[[i]])
            
            samID <- ceiling(i / nfolds)
            j <- i - (samID - 1) * nfolds
            
            sam <- sampleSplits[[samID]]
            foldid <- foldids[[samID]][[j]]
            foldid.internal <- foldid.internals[[samID]][[j]]
            X.t <- X[sam != j,]
            y.t <- y[sam != j]
            X.v <- X[sam == j,]
            y.v <- y[sam == j]
            
            argList <- list(X = X.t, y = y.t, 
                            times = times,
                            cox.time = cox.time, 
                            blocks = blocks,
                            foldid = foldid,
                            intercept = TRUE,
                            foldid.internal = foldid.internal)
            
            if (whichRun == "competitors") {
              ncomp <- nrow(X.t) - 1
              argList$ncomp <- ncomp
              mod <- do.call(fit.others, argList)
            } else if(whichRun == "stk"){
              mod <- do.call(fit.stk, argList)
            } else if (whichRun == "stkPF.pro") {
              mod <- do.call(fit.stkPF.pro, argList)
            } else if (whichRun == "random.stk") {
              mod <- do.call(fit.random.stk, argList)
            }
            
            save(list = "mod", 
                 file = paste0("../output/", whichRun, "-", nam, "-", samID, "-", j, ".RData"))
          }
  
  stopCluster(cl)
}


# summary data
# # V25 and V24 to stkPF.pro
# # whichRun = "competitors"
# 
# for(nam in datNames){
#   
#   cat(nam, ":", as.character(Sys.time()), "\n\n")
#   getData(nam, log2 = T, toBinary = F)
#   n <- nrow(X)
#   p <- ncol(X)
#   
#   getSplits(n, nsim, nfolds, seed = 1248372, foldid.internal = T)
#   
#   for(i in 1 : (nfolds * nsim)){
#     
#     samID <- ceiling(i / nfolds)
#     #- j, which fold of the i.nsim split being the test
#     j <- i - (samID - 1) * nfolds
#     #- foldid for i.nsim splits and j fold
#     
#     sam <- sampleSplits[[samID]]
#     load(paste0("../output/", nam, "-", samID, "-", j, "-coxOptimV16.RData"))
#     mod <- mod[!grepl("stk", names(mod))]
#     save(list = "mod", file = paste0("../output/", whichRun, "-", nam, "-", samID, "-", j, "-coxOptimV26.RData"))
#   }
# }
# 
# #-
# 
# whichRun = "stk"
# 
# for(nam in datNames){
#   
#   cat(nam, ":", as.character(Sys.time()), "\n\n")
#   getData(nam, log2 = T, toBinary = F)
#   n <- nrow(X)
#   p <- ncol(X)
#   
#   getSplits(n, nsim, nfolds, seed = 1248372, foldid.internal = T)
#   
#   for(i in 1 : (nfolds * nsim)){
#     
#     samID <- ceiling(i / nfolds)
#     #- j, which fold of the i.nsim split being the test
#     j <- i - (samID - 1) * nfolds
#     #- foldid for i.nsim splits and j fold
#     
#     sam <- sampleSplits[[samID]]
#     load(paste0("../output/", nam, "-", samID, "-", j, "-coxOptimV24.RData"))
#     # mod <- mod[!grepl("stk", names(mod))]
#     save(list = "mod", file = paste0("../output/", whichRun, "-", nam, "-", samID, "-", j, "-coxOptimV26.RData"))
#   }
# }
# 
# #- stkPF.pro

#-------------------------------------------------------------------------------
# Compare Bagging and CV base learner
#-------------------------------------------------------------------------------

calAUC <- function(marker, events){
  
  ordMarker <- order(marker, decreasing = T)
  events    <- events[ordMarker]
  marker    <- marker[ordMarker]
  
  # TP
  TP   <- c(0, cumsum(events)) /  sum(events)
  TP   <- TP[!duplicated(marker[ordMarker])]
  
  # FP
  FP   <- c(0, cumsum((1 - events))) /  sum((1 - events))
  FP   <- FP[!duplicated(marker[ordMarker])]
  
  AUC <- trapezoidal(FP, TP)
  AUC
}

#---

whichRun = "random.stk"
# ifCVMods = TRUE
nCores = 25

# set.seed(321)
# datNames <- c("BLCA", "HNSC")
res <- list()
for(nam in datNames) {
  
  print(nam)
  getData(nam, log2 = T, toBinary = F)
  getSplits(n, nsim, nfolds, seed = 1248372, foldid.internal = T)
  
  Cores <- pmin(nCores, detectCores() - 1)
  cl <- makeCluster(Cores)
  registerDoSNOW(cl)
  
  res[[nam]] <- foreach(i = 1 : (nfolds * nsim), 
                        .combine = "rbind",
                        .options.snow = opts,
                        .verbose = TRUE) %dopar% {
                          
                          samID <- ceiling(i / nfolds)
                          # j, which fold of the i.nsim split being the test
                          j <- i - (samID - 1) * nfolds
                          # foldid for i.nsim splits and j fold
                          sam <- sampleSplits[[samID]]
                          load(paste0("../output/", "random.stk", "-", nam, "-", samID, "-", j, "-coxRandom2.RData"))
                          cv1 <- mod$Rstk.comp$CVpred
                          load(paste0("../output/", "stk", "-", nam, "-", samID, "-", j, "-coxOptimV26.RData"))
                          cv2 <- mod$stk.comp$CVpred
                          
                          a1 <- calAUC(cv1[,3], cv1[,2])
                          a2 <- calAUC(cv2[,3], cv2[,2])
                          
                          b1 <- calAUC(cv1[,4], cv1[,2])
                          b2 <- calAUC(cv2[,4], cv2[,2])
                          
                          c(a1, a2, b1, b2)
                        }
  
}


lapply(res, colMeans)
  
#-------------------------------------------------------------------------------

sumResults <- function(target, ifCVMods, whichRun, ifSE = FALSE,){
  # AUC
  reMat <- c()
  reMatSe <- c()
  
  for(i in seq_along(whichRun)) {
    
    load(paste0("../output/survPredRes", whichRun[i], "-V1", ifelse(ifCVMods, "T", "F"), ".RData"))
    n <- length(get(target)[[1]])
    
    # remove duplicated clin
    tmp <- sumFun(get(target), fun = colMeans)
    if (i != 1)
      tmp <- tmp[,-1]
    
    reMat <- cbind(reMat, tmp)
    
    if (ifSE) {
      tmp <- sumFun(get(target), fun = function(x) apply(x, 2, sd) / sqrt(n))
      if (i != 1)
        tmp <- tmp[,-1]
      
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
  #   re.se <- as.data.frame(reMatSe) %>% 
  #     select(clin, mol, naive, fs.las, pcr.las, ipf.las, `stk.pf:ibsLoss`, `stk.pf:logLik:ridge`, `stk.pf:PaLogLik:ridge`, `stk.comp:ibsLoss`,
  #            `stk.comp:logLik:ridge`, `stk.comp:PaLogLik:ridge`) %>%
  #     rename_all( ~ c("Clin", "Mol", "Naive", "FS(las)", "PCR(las)", "IPF(las)", "Stk(PF,IBS)", "Stk(PF,Lik)", "Stk(PF,PaLik)",
  #                     "Stk(CP,IBS)", "Stk(CP,Lik)", "Stk(CP,PaLik)")) %>%
  #     mutate(across(everything(), ~ formatC(.x, format = "f", digits = 3)))
  #   
  #   re <- mapply(function(mean, se) {
  #     paste0(mean, "(", se, ")")
  #   }, mean = re, se = re.se)
  #   rownames(re) <- rownames(re.se)
  # }
  
  # rowname data to first col
  re <- tibble::rownames_to_column(as.data.frame(re), var = "Data")
  # openxlsx::write.xlsx(re, file = paste0("../results/", target, "V1.xlsx"), colNames = T)
  
  # print(re)
  return(re)
}

sumResults(target = "aucs", ifCVMods = TRUE, ifSE = TRUE)
sumResults(target = "ibss", ifCVMods = TRUE, ifSE = TRUE)
sumResults(target = "logLiks", ifCVMods = TRUE, ifSE = TRUE)

# output Results to Boxplots
# load("../output/survPredResV16.RData")

makePointRangePlot <- function(obj, title, xlab, ylab) {
  
  areaSum <- function(obj) {
    obj <- na.omit(obj)
    tRange <- abs(range(obj[ ,"time"])[1] - range(obj[ ,"time"])[2])
    trapezoidal(obj[,"time"], obj[,2]) / tRange
  }
  
  # data frame
  
  sumRes <- lapply(obj, function(x) unlist(lapply(x, areaSum)))
  sumRes <- Reduce(rbind, sumRes)
  
  # Tab
  sumRes <- sumRes %>% 
    as.data.frame() %>%
    rename_at(colnames(sumRes), ~ reNames) %>% 
    select(-starts_with("Stk(SS,")) %>%
    pivot_longer(cols = everything()) %>%
    group_by(name) %>%
    summarise(mean = mean(value),
              sd = sd(value) / sqrt(n())) %>%
    mutate(color = case_when(name %in% c("Clin", "Mol", "Naive") ~ "Basic",
                             name %in% c("FS(las)", "PCR(las)") ~ "DR",
                             name %in% c("IPF(las)") ~ "CW",
                             grepl("^Stk\\(PF", name) ~ "Stk(PF)",
                             grepl("^Stk\\(CP", name) ~ "Stk(CP)"),
           color = factor(color, levels = c("Basic", "DR", "CW", "Stk(PF)", "Stk(CP)")),
           name = factor(name, levels = reNames))
  
  # 
  ggplot(sumRes, aes(x = name, y = mean, color = color)) + 
    geom_line() +
    geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd)) +
    theme_bw() + 
    labs(title = title,
         x = xlab,
         y = ylab) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.title = element_blank())
}


# Auc
aucPointRangePlot <- list()
for(i in names(aucs)){
  aucPointRangePlot[[i]] <- makePointRangePlot(aucs[[i]], title = i, xlab = "Methods", ylab = "AUC")
}

g <- ggpubr::ggarrange(plotlist = aucPointRangePlot, common.legend = TRUE, nrow = 2,
                  ncol = 5)
ggsave("../results/aucBar.pdf", g, width = 12, height = 8)


# logLik
logLikPointRangePlot <- list()
for(i in names(logLiks)){
  logLikPointRangePlot[[i]] <- makePointRangePlot(logLiks[[i]], title = i, 
                                               xlab = "Methods", ylab = "logLik")
}

g <- ggpubr::ggarrange(plotlist = logLikPointRangePlot, common.legend = TRUE, nrow = 2,
                  ncol = 5)

ggsave("../results/logLikBar.pdf", g, width = 12, height = 8)

#-------------------------------------------------------------------------------
# Subset Stacking 
#-------------------------------------------------------------------------------

vsSS <- function(target, ifCVMods, ifSE = FALSE, to = c("table", "figure")) {
  
  # target: one of aucs, ibss, and logLiks
  # get all data
  
  reMat <- c()
  reMatSe <- c()
  
  whichRun <- c("competitors", "stk", "stkPF.pro")
  for(i in seq_along(whichRun)) {
    
    load(paste0("../output/survPredRes", whichRun[i], "-V26", ifelse(ifCVMods, "T", "F"), ".RData"))
    n <- length(get(target)[[1]])
    
    # remove duplicated clin
    tmp <- sumFun(get(target), fun = colMeans)
    if (i != 1)
      tmp <- tmp[,-1]
    
    reMat <- cbind(reMat, tmp)
    
    # for ifSE
    tmp <- sumFun(get(target), fun = function(x) apply(x, 2, sd) / sqrt(n))
    if (i != 1)
      tmp <- tmp[,-1]
    
    reMatSe <- cbind(reMatSe, tmp)
  }
  
  # selected moehods for comparison
  # to table
  sel <- c("stk.pf:ibsLoss", "stk.subset:ibsLoss", "stk.pf:logLik:ridge", "stk.subset:logLik:ridge", "stk.pf:PaLogLik:ridge", "stk.subset:PaLogLik:ridge")
  reNames <- c("Stk(PF,IBS)", "Stk(SS,IBS)", "Stk(PF,Lik)", "Stk(SS,Lik)", "Stk(PF,PaLik)", "Stk(SS,PaLik)")
  # I cannot use select because of `:`
  re <- as.data.frame(reMat) %>%
    select(all_of(sel)) %>%
    rename_all( ~ reNames)
  
  re.se <- as.data.frame(reMatSe) %>%
    select(all_of(sel)) %>%
    rename_all( ~ reNames)
  
  
  if (to == "table") {
    
    if (ifSE) {
      reTab <- mapply(function(mean, se) {
        mean <- formatC(mean, format = "f", digits = 3)
        se   <- formatC(se, format = "f", digits = 3)
        paste0(mean, "(", se, ")")
      }, mean = re, se = re.se)
      rownames(reTab) <- rownames(re.se)
    } else {
      reTab <- re
    }
    
    reTab <- tibble::rownames_to_column(as.data.frame(reTab), var = "Data")
    openxlsx::write.xlsx(reTab, file = paste0("../results/ss_", target, "V26.xlsx"), colNames = T)
    print(reTab)
  }
  
  if (to == "figure") {
    
    # to figures
    plotDat <- as.data.frame(re) %>% 
      rownames_to_column(var = "Data")  %>%
      pivot_longer(-1, values_to = "Mean") %>%
      left_join(as.data.frame(re.se) %>% 
                  rownames_to_column(var = "Data")  %>%
                  pivot_longer(-1, values_to = "SE"), by = c("Data", "name")) %>%
      mutate(method = ifelse(grepl("Stk\\(SS", name), "Stk(SS)", "Stk(PF)"),
             name = factor(name, levels = colnames(re), labels = colnames(re))) %>%
      # add filter
      filter(grepl(",Lik", name))
    
    ggplot(plotDat, aes(x = name, y = Mean, color = method))+ 
      geom_line(group = 1, color = 1) + 
      geom_pointrange(aes(ymin = Mean - SE, ymax = Mean + SE)) +
      # geom_text(aes(label = round(Mean, 3)), vjust = -0.5) +
      theme_bw() + 
      labs(title = "",
           x = "",
           y = case_when(target == "aucs" ~ "iAUC",
                         target == "ibss" ~ "IBS",
                         target == "logLiks"  ~ "iNNL")) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.title = element_blank(),
            legend.position = "none") +
      facet_grid( ~ Data, scales = "free")
  }
}

vsSS(target = "aucs", ifCVMods = TRUE, ifSE = TRUE, to = "figure") 
vsSS(target = "ibss", ifCVMods = TRUE, ifSE = TRUE, to = "figure") 
vsSS(target = "logLiks", ifCVMods = TRUE, ifSE = TRUE, to = "figure") 

#-------------------------------------------------------------------------------
# 
# # compare to the Best
# getStkSSTab <- function(objs, fun, name, reNames, ifMax = TRUE,
#                         compareTo = "best") {
#   
#   # areaSum <- function(obj) {
#   #   obj <- na.omit(obj)
#   #   tRange <- abs(range(obj[ ,"time"])[1] - range(obj[ ,"time"])[2])
#   #   trapezoidal(obj[,"time"], obj[,2]) / tRange
#   # }
#   # 
#   # # data frame
#   # sumRes <- function(obj, fun){
#   #   res <- lapply(obj, function(x) unlist(lapply(x, areaSum)))
#   #   res <- Reduce(rbind, res)
#   #   do.call(fun, list(res))
#   # }
#   
#   tab <- lapply(objs, sumRes, fun = fun)
#   tab <- Reduce(rbind, tab)
#   rownames(tab) <- names(objs)
#   
#   tab <- tab %>% 
#     as.data.frame() %>%
#     rename_at(colnames(tab), ~ reNames) %>% 
#     tibble::rownames_to_column(var = "Data")
#   
#   if (compareTo == "best") {
#     
#       # select the best AUC
#     tab <- tab %>% 
#       rowwise() %>%
#       mutate(best = ifelse(ifMax, max(across(-c(1, grep("Stk\\(SS", colnames(.))))),
#                            min(across(-c(1, grep("Stk\\(SS", colnames(.)))))),
#              bestMod = colnames(.)[-c(1, grep("Stk\\(SS", colnames(.)))][
#                ifelse(ifMax, which.max(across(-c(1, grep("Stk\\(SS", colnames(.))))),
#                       which.min(across(-c(1, grep("Stk\\(SS", colnames(.))))))]) %>%
#       ungroup() %>%
#       mutate(across(-grep("bestMod", colnames(.)), ~ format(.x, digits = 3, nsmall = 3L))) %>%
#       mutate(across(everything(), as.character)) %>%
#       select(c("Data", "bestMod", "best", "Stk(SS,IBS)", "Stk(SS,Rdg)", "Stk(SS,CoxRdg)")) 
#     
#   } else {
#     
#     tab <- tab %>%
#       select(c("Data", compareTo, "Stk(SS,IBS)", "Stk(SS,Rdg)", "Stk(SS,CoxRdg)")) %>%
#       mutate(across(-1, format, digits = 3, nsmall = 3L))
#   }
#   
#   openxlsx::write.xlsx(tab, paste0("../results/subset", name, "V16.xlsx"), colNames = T)
#   
#   return(tab)
# }
# 
# getStkSSTab(aucs, fun = colMeans, name = "auc", reNames = reNames, ifMax = TRUE,
#             compareTo = "Stk(PF,Rdg)")
# getStkSSTab(ibss, fun = colMeans, name = "ibs", reNames = reNames, ifMax = FALSE,
#             compareTo = "Stk(PF,Rdg)")
# 
# getStkSSTab(logLiks, fun = colMeans, name = "logLik", reNames = reNames, ifMax = FALSE,
#             compareTo = "Stk(PF,Rdg)")
# 
# #- plot
# 
# makePointRangePlot_ss <- function(obj, title, xlab, ylab) {
#   
#   areaSum <- function(obj) {
#     obj <- na.omit(obj)
#     tRange <- abs(range(obj[ ,"time"])[1] - range(obj[ ,"time"])[2])
#     trapezoidal(obj[,"time"], obj[,2]) / tRange
#   }
#   
#   # data frame
#   
#   sumRes <- lapply(obj, function(x) unlist(lapply(x, areaSum)))
#   sumRes <- Reduce(rbind, sumRes)
#   
#   # Tab
#   sumRes <- sumRes %>% 
#     as.data.frame() %>%
#     rename_at(colnames(sumRes), ~ reNames) %>% 
#     select(-starts_with("Stk(SS,")) %>%
#     pivot_longer(cols = everything()) %>%
#     group_by(name) %>%
#     summarise(mean = mean(value),
#               sd = sd(value) / sqrt(n())) %>%
#     mutate(color = case_when(name %in% c("Clin", "Mol", "Naive") ~ "Basic",
#                              name %in% c("FS(las)", "PCR(las)") ~ "DR",
#                              name %in% c("IPF(las)") ~ "CW",
#                              grepl("^Stk\\(PF", name) ~ "Stk(PF)",
#                              grepl("^Stk\\(CP", name) ~ "Stk(CP)"),
#            color = factor(color, levels = c("Basic", "DR", "CW", "Stk(PF)", "Stk(CP)")),
#            name = factor(name, levels = reNames))
#   
#   # 
#   ggplot(sumRes, aes(x = name, y = mean, color = color)) + 
#     geom_line() +
#     geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd)) +
#     theme_bw() + 
#     labs(title = title,
#          x = xlab,
#          y = ylab) + 
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#           legend.title = element_blank())
# }
# 
#-------------------------------------------------------------------------------
# compare to Stk(PF, Lik)
#-------------------------------------------------------------------------------
# 
# getStkSSTab <- function(objs, fun, name, reNames, ifMax = TRUE) {
#   
#   areaSum <- function(obj) {
#     obj <- na.omit(obj)
#     tRange <- abs(range(obj[ ,"time"])[1] - range(obj[ ,"time"])[2])
#     trapezoidal(obj[,"time"], obj[,2]) / tRange
#   }
#   
#   # data frame
#   sumRes <- function(obj, fun){
#     res <- lapply(obj, function(x) unlist(lapply(x, areaSum)))
#     res <- Reduce(rbind, res)
#     do.call(fun, list(res))
#   }
#   
#   tab <- lapply(objs, sumRes, fun = fun)
#   tab <- Reduce(rbind, tab)
#   rownames(tab) <- names(objs)
#   
#   tab <- tab %>% 
#     as.data.frame() %>%
#     rename_at(colnames(tab), ~ reNames) %>% 
#     tibble::rownames_to_column(var = "Data") %>%
#     # select the best AUC
#     rowwise() %>%
#     mutate(best = ifelse(ifMax, max(across(-c(1, grep("Stk\\(SS", colnames(.))))),
#                          min(across(-c(1, grep("Stk\\(SS", colnames(.)))))),
#            bestMod = colnames(.)[-c(1, grep("Stk\\(SS", colnames(.)))][
#              ifelse(ifMax, which.max(across(-c(1, grep("Stk\\(SS", colnames(.))))),
#                     which.min(across(-c(1, grep("Stk\\(SS", colnames(.))))))]) %>%
#     ungroup() %>%
#     mutate(across(-grep("bestMod", colnames(.)), ~ format(.x, digits = 3, nsmall = 3L))) %>%
#     mutate(across(everything(), as.character)) %>%
#     select(c("Data", "bestMod", "best", "Stk(SS,IBS)", "Stk(SS,Rdg)", "Stk(SS,CoxRdg)")) 
#   
#   openxlsx::write.xlsx(tab, paste0("../results/subset", name, "V16.xlsx"), colNames = T)
#   
#   return(tab)
# }
#-------------------------------------------------------------------------------
#- Compare CV models and single models.
#-------------------------------------------------------------------------------

whichRun = "stk"

load(paste0("../output/survPredRes", whichRun, "-V26T.RData"))
objT <- sumFun(logLiks, fun = colMeans)[,-1]
load(paste0("../output/survPredRes", whichRun, "-V26F.RData"))
objF <- sumFun(logLiks, fun = colMeans)[,-1]

plotDat <- as.data.frame(objT) %>% 
  rownames_to_column(var = "Data") %>%
  pivot_longer(col = c(-1)) %>%
  mutate(cvMods = T) %>%
  bind_rows({
    as.data.frame(objF) %>% 
      rownames_to_column(var = "Data") %>%
      pivot_longer(col = c(-1)) %>%
      mutate(cvMods = F)
  })

# rename and factorize
plotDat <- plotDat %>%
  mutate(name = factor(name, 
                       labels = c("Stk(CP, IBS)", "Stk(CP, Lik)", "Stk(CP, PaLik)", "Stk(PF, IBS)", "Stk(PF, Lik)", "Stk(PF, PaLik)"))) 

ggplot(plotDat, aes(x = name, y = value, color = cvMods)) + 
  geom_point() + 
  facet_grid(~ Data) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  xlab(NULL) + ylab("iAUC") + 
  theme(legend.position = "")

#-------------------------------------------------------------------------------
# Weights: alpha
#-------------------------------------------------------------------------------
# IPFLasso, PF distribution

pfs <- list()
for(nam in datNames) {
  
  pf <- c()
  for(i in 1:(nfolds * nsim)){
    samID <- ceiling(i / nfolds)
    j <- i - (samID - 1) * nfolds
    whichRun = "competitors"
    load(paste0("../output/", whichRun, "-", nam, "-", samID, "-", j, "-coxOptimV26.RData"))
    pf <- c(pf, paste0(mod$ipf.las$pf, collapse = "-"))
  }
  
  tmp <- setNames(rep(0, 4), c("1-1", "1-2", "1-4", "1-8"))
  pf  <- table(pf) / length(pf)
  tmp[names(pf)] <- pf
  pfs[[nam]] <- tmp
}

# data processing
pfs <- as.data.frame(pfs) %>%
  rownames_to_column(var = "PF") %>%
  pivot_longer(cols = -1,
               names_to = "Data", 
               values_to = "PF(%)")

gBar <- ggplot(pfs, aes(x = factor(PF), y = `PF(%)`)) + 
  geom_bar(stat="identity") +
  theme_bw() + 
  labs(x = NULL, title = "PF in IPF(las)") +
  facet_grid(~ Data) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# Alpha for Stk(PF, Rdg)

getAlpha <- function(datName, whichMethod, whichLoss, paraList) {
  
  nfolds = paraList$nfolds
  nsim   = paraList$nsim
  
  alpha <- vector("list", length(whichLoss))
  for(i in 1:(nfolds * nsim)){
    
    samID <- ceiling(i / nfolds)
    j <- i - (samID - 1) * nfolds
    
    # merge data
    # get stk
    whichRun = "stk"
    load(paste0("../output/", whichRun, "-", datName, "-", samID, "-", j, "-coxOptimV26.RData"))
    mod0 <- mod
    # get stk.subset 
    whichRun = "stkPF.pro"
    load(paste0("../output/", whichRun, "-", datName, "-", samID, "-", j, "-coxOptimV26.RData"))
    mod <- c(mod0, mod)
    rm(mod0)
    
    a <- mod[[whichMethod]]$alpha[whichLoss]
    
    for(k in seq_along(a)) {
      if (names(a)[k] %in% c("logLik:lasso", "logLik:ridge"))
        sub_a <- a[[k]][-1]
      else
        sub_a <- a[[k]]
      
      alpha[[k]] <- rbind(alpha[[k]], sub_a)
    }
  }
  names(alpha) = whichLoss
  return(alpha)
}

alphas <- list()
for(datName in datNames) {
  alphas[[datName]] <- getAlpha(datName, 
                                whichMethod = "stk.pf", 
                                whichLoss = c("ibsLoss", "logLik:ridge"),
                                paraList = paraList)
}


# Boxplot 

toggPlotData <- function(alpha, alphaName) {
  
  re <- c()
  for(i in seq_along(alpha)){
   
    a <- alpha[[i]]
    a <- a %>% 
      as.data.frame() %>%
      setNames(alphaName) %>%
      pivot_longer(cols = everything(.), 
                   names_to = "subModels",
                   values_to  = "w") %>%
      mutate(lossType = names(alpha)[i])
    
    re <- rbind(re, a)
  }
  return(re)
}


# get data
ggData <- c()
for(nam in names(alphas)) {
  alpha <- alphas[[nam]]
  tmp <- toggPlotData(alpha, c("PF(1-1)", "PF(1-2)", "PF(1-4)", "PF(1-8)"))
  tmp$data <- nam
  ggData <- rbind(ggData, tmp)
}

# boxplot and barPlot

gBox1 <- ggplot(ggData[ggData$lossType == "ibsLoss",], aes(x = factor(subModels), y = w)) + 
  geom_boxplot(outliers = FALSE) + 
  theme_bw() + 
  labs(y = "Weights", x = NULL, title = "PF weights in Stk(IBS)") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid( ~ data) 
  # theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "mm"))

gBox2 <- ggplot(ggData[ggData$lossType == "logLik:ridge",], aes(x = factor(subModels), y = w)) + 
  geom_boxplot(outliers = FALSE) + 
  theme_bw() + 
  labs(y = "Weights", x = NULL, title = "PF weights in Stk(PF,logLik)") + 
  facet_grid( ~ data) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

gg <- ggpubr::ggarrange(gBar, gBox1, gBox2, nrow = 3)
ggsave("../results/pf_alpha.pdf", plot = gg, width = 10, height = 8)
  
ggsave("../results/pf_alpha.jpeg", plot = gg, width = 10, height = 8, 
       units = "in", dpi = 600)


################################################################################
# Alpha in Subset Stacking
################################################################################

alphas <- list()

for(datName in datNames) {
  alphas[[datName]] <- getAlpha(datName, 
                                whichMethod = "stk.subset", 
                                whichLoss = c("ibsLoss", "logLik:ridge"),
                                paraList = paraList)
}

# Boxplot 
toggPlotData <- function(alpha, alphaName) {
  
  re <- c()
  for(i in seq_along(alpha)){
    
    a <- alpha[[i]]
    a <- a %>% 
      as.data.frame() %>%
      setNames(alphaName) %>%
      pivot_longer(cols = everything(.), 
                   names_to = "subModels",
                   values_to  = "w") %>%
      mutate(lossType = names(alpha)[i])
    
    re <- rbind(re, a)
  }
  return(re)
}

# get data
ggData <- c()
for(nam in names(alphas)) {
  alpha <- alphas[[nam]]
  subModName <- paste0("S", rep(1:10, each = 4), ",",
                       rep(c("1-1", "1-2", "1-4", "1-8"), times = 10)) 
  tmp <- toggPlotData(alpha, subModName) %>%
    mutate(data = nam,
           PF = sub(".*\\,(.*)", "\\1", subModels),
           SS = sub("\\,.*$", "", subModels),
           SS = factor(SS, levels = paste0("S", 1:10)),
           subModels = factor(subModels, levels = subModName)) 
  ggData <- rbind(ggData, tmp)
}

# boxplot and barPlot

gBox1 <- ggplot(ggData[ggData$lossType == "ibsLoss", ], 
                aes(x = factor(subModels), y = w, fill = SS)) + 
  geom_boxplot(outliers = TRUE,  outlier.size = 0.5) + 
  theme_bw() + 
  labs(y = "Weights", x = NULL) + 
  theme(axis.text.x = element_text(
    angle = 90, vjust = 0.5, hjust = 1, size = 6),
    legend.title = element_blank()) + 
  facet_wrap(vars(data), nrow = 5)

ggsave("../results/ss_alpha_ibs.jpeg", gBox1, width = 8, height = 8, units = "in",
       dpi = 300)

# end

gBox2 <- ggplot(ggData[ggData$lossType == "logLik:ridge",], 
                aes(x = factor(subModels), y = w, fill = SS)) + 
  geom_boxplot(outliers = TRUE, outlier.size = 0.5) + 
  theme_bw() + 
  labs(y = "Weights", x = NULL) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   size = 6),
        legend.title = element_blank()) + 
  facet_wrap(vars(data), nrow = 5) 

ggsave("../results/ss_alpha_logLik.jpeg", gBox2, width = 8, height = 8, units = "in",
       dpi = 300)




#-------------------------------------------------------------------------------
# Computational complexity
#-------------------------------------------------------------------------------

# get time Data

getTime <- function(datName, paraList){

  nfolds = paraList$nfolds
  nsim = paraList$nsim
  
  methodNames <- c("Clin", "Mol", "Naive", "FS(las)", "PCR(las)","IPF(las)", 
                   "Stk(PF)", "Stk(CP)", "Stk(SS)")
  
  time <- c()
  for(i in 1:(nfolds * nsim)){
    samID <- ceiling(i / nfolds)
    j <- i - (samID - 1) * nfolds
    load(paste0("../output/", datName, "-", samID, "-", j, "-coxOptimV16.RData"))
    mod0 <- mod
    load(paste0("../output/", datName, "-", samID, "-", j, "-coxOptimV23.RData"))
    
    mod0[["stk.subset"]] <- mod[["stk.subset"]]
    mod <- mod0 
    rm(mod0)
    
    tmp <- lapply(mod, `[[`, "time")
    time <- rbind(time, unlist(lapply(tmp, as.numeric, units = "mins")))
  }
  
  timeggData <- time %>% 
    as.data.frame() %>%
    setNames(methodNames) %>%
    pivot_longer(everything(),
                 names_to = "Method",
                 values_to = "Time") %>%
    group_by(Method) %>%
    summarise(Time = mean(Time)) %>% 
    mutate(Method = factor(Method, levels = methodNames))
  
  return(timeggData)
}

  
#- using BLCA as an example, ggplot

timeggData <- getTime("BLCA", paraList)
p1 <- ggplot(timeggData, aes(x = Method, y = Time)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 180)) +
  labs(x = NULL, title = "A, BLCA") + 
  coord_flip()

timeggData <- getTime("BRCA", paraList)
p2 <- ggplot(timeggData, aes(x = Method, y = Time)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 600)) +
  labs(x = NULL, title = "B, BRCA") + 
  coord_flip()

gg <- ggpubr::ggarrange(p1, p2, nrow = 1)
ggsave("../results/timeV16.pdf", plot = gg, width = 7, height = 3)

ggsave("../results/timeV16.jpeg", plot = gg, width = 7, height = 3,
       units = "in", dpi = 600)

#-------------------------------------------------------------------------------







pList <- list()

for(nam in names(aucsDatList)) {
  
  aucsDat <- aucsDatList[[nam]]
  names(aucsDat) <- reNames
  
  aucsPlot <- c()
  for(i in names(aucsDat)) {
    tmp<- cbind(as.data.frame(aucsDat[[i]]), method = i)
    aucsPlot <- rbind(aucsPlot, tmp)
  } 
  
  aucsPlot$method <- factor(aucsPlot$method, reNames)
  
  aucsPlot <- aucsPlot %>% 
    na.omit(aucsPlot) %>%
    mutate(time = round(time)) %>%
    group_by(method, time) %>%
    summarise(auc = mean(auc)) %>%
    filter(method %in% c("IPF(las)", "Stk(SS,IBS)", "Stk(SS,Rdg)", "Stk(SS,CoxRdg)"))
    # filter(!grepl("Las", method)) %>%
    # filter(!grepl("add1", method)) %>%
    # filter(!grepl("Stk\\(SS", method))
  
  nColor <- length(unique(aucsPlot$method))
  set.seed(1)
  colors <- randomcoloR::distinctColorPalette(k = nColor)
  
  pList[[nam]] <- ggplot(aucsPlot, aes(x = time, y = auc, colour = method)) + 
    geom_line() +
    scale_color_manual(values = colors) + 
    theme_bw() +
    scale_x_continuous(breaks = unique(aucsPlot$time)) +
    ggtitle(nam)
  
}

pdf("../results/PlotTimeAUCV16.pdf", width = 10, height = 6)

for (i in seq_along(pList)) {
  print(pList[[i]])
}

dev.off()

# save(list = c("aucs", "ibss", "logLiks"), file = "../output/survPredResV16.RData")

# load("../output/survPredResV16.RData")

reNames <- c("Clin", "Mol", "Naive", "FS(las)", "PCR(las)", "IPF(las)",
             "Stk(PF,IBS)", "Stk(PF,Las)", "Stk(PF,Rdg)",
             "Stk(PF,CoxLas)", "Stk(PF,CoxRdg)", "Stk(PF,add1)",
             "Stk(CP,IBS)", "Stk(CP,Las)", "Stk(CP,Rdg)",
             "Stk(CP,CoxLas)", "Stk(CP,CoxRdg)", "Stk(CP,add1)",
             "Stk(SS,IBS)", "Stk(SS,Las)", "Stk(SS,Rdg)",
             "Stk(SS,CoxLas)", "Stk(SS,CoxRdg)", "Stk(SS,add1)")

## output IBS
## remove logLik method based on Lasso

sumRes <- function(obj, fun = colMeans){
  res <- Reduce(rbind, obj)
  do.call(fun, list(res))
}

tab <- lapply(ibss, sumRes)
tab <- Reduce(rbind, tab)
rownames(tab) <- names(ibss)

tab <- tab %>% 
  as.data.frame() %>%
  rename_at(colnames(tab), ~ reNames) %>% 
  select(-c("Stk(PF,Las)", "Stk(PF,CoxLas)", 
            "Stk(CP,Las)", "Stk(CP,CoxLas)", 
            "Stk(SS,Las)", "Stk(SS,CoxLas)")) %>%
  filter(!rownames(.) %in% c("LUSC", "OV", "STAD")) %>%
  format(digits = 3, nsmall = 3L) %>%
  mutate(across(everything(), as.character)) %>%
  tibble::rownames_to_column() %>%
  openxlsx::write.xlsx(., "../results/survIBSV16.xlsx", colNames = T)

## output logLiks
tab <- lapply(logLiks, sumRes)
tab <- Reduce(rbind, tab)
rownames(tab) <- names(logLiks)

tab <- tab %>% 
  as.data.frame() %>%
  rename_at(colnames(tab), ~ reNames) %>% 
  select(-c("Stk(PF,Las)", "Stk(PF,CoxLas)", 
            "Stk(CP,Las)", "Stk(CP,CoxLas)", 
            "Stk(SS,Las)", "Stk(SS,CoxLas)")) %>%
  filter(!rownames(.) %in% c("LUSC", "OV", "STAD")) %>%
  format(digits = 1, nsmall = 1L) %>%
  mutate(across(everything(), as.character)) %>%
  tibble::rownames_to_column() %>%
  openxlsx::write.xlsx(., "../results/survlogLikV16.xlsx", colNames = T)

## output AUC
tab <- lapply(aucs, sumRes)
tab <- Reduce(rbind, tab)
rownames(tab) <- names(aucs)

tab <- tab %>% 
  as.data.frame() %>%
  rename_at(colnames(tab), ~ reNames) %>% 
  select(-c("Stk(PF,Las)", "Stk(PF,CoxLas)", 
            "Stk(CP,Las)", "Stk(CP,CoxLas)", 
            "Stk(SS,Las)", "Stk(SS,CoxLas)")) %>%
  filter(!rownames(.) %in% c("LUSC", "OV", "STAD")) %>%
  format(digits = 3, nsmall = 3L) %>%
  mutate(across(everything(), as.character)) %>%
  tibble::rownames_to_column() %>%
  openxlsx::write.xlsx(., "../results/survAUCV16.xlsx", colNames = T)

## Compare to the new proposal Stk(SS)

# tab <- t(cbind.data.frame(lapply(ibss, colMeans)))
# 
# tab <- tab %>% 
#   as.data.frame() %>%
#   rename_at(colnames(tab), ~ reNames) %>% 
#   select(-c("Stk(PF,Las)", "Stk(PF,CoxLas)", 
#             "Stk(CP,Las)", "Stk(CP,CoxLas)", 
#             "Stk(SS,Las)", "Stk(SS,CoxLas)")) %>%
#   filter(!rownames(.) %in% c("LUSC", "OV", "STAD")) %>% 
#   format(digits = 3, nsmall = 3L) %>%
#   mutate(across(everything(), as.character)) %>%
#   tibble::rownames_to_column() %>% 
  
# select the best AUC 
tab <- lapply(aucs, sumRes)
tab <- Reduce(rbind, tab)
rownames(tab) <- names(aucs)

tab <- tab %>% 
  as.data.frame() %>%
  rename_at(colnames(tab), ~ reNames) %>% 
  select(-c("Stk(PF,Las)", "Stk(PF,CoxLas)", "Stk(PF,add1)",
            "Stk(CP,Las)", "Stk(CP,CoxLas)", "Stk(CP,add1)",
            "Stk(SS,Las)", "Stk(SS,CoxLas)", "Stk(SS,add1)")) %>%
  filter(!rownames(.) %in% c("LUSC", "OV", "STAD")) %>% 
  tibble::rownames_to_column(var = "Data") %>%
  # select the best AUC
  rowwise() %>%
  mutate(bestAUC = max(across(-c(1, grep("Stk\\(SS", colnames(.))))),
         bestMod = colnames(.)[-c(1, grep("Stk\\(SS", colnames(.)))][which.max(across(-c(1, grep("Stk\\(SS", colnames(.)))))]) %>%
  ungroup() %>%
  mutate(across(-grep("bestMod", colnames(.)), ~ format(.x, digits = 3, nsmall = 3L))) %>%
  mutate(across(everything(), as.character)) %>%
  select(c("Data", "bestMod", "bestAUC", "Stk(SS,IBS)", "Stk(SS,Rdg)", "Stk(SS,CoxRdg)")) %>%
  openxlsx::write.xlsx(., "../results/subsetAUCV16.xlsx", colNames = T)

# select the best IBS

tab <- lapply(ibss, sumRes)
tab <- Reduce(rbind, tab)
rownames(tab) <- names(aucs)

tab <- tab %>% 
  as.data.frame() %>%
  rename_at(colnames(tab), ~ reNames) %>% 
  select(-c("Stk(PF,Las)", "Stk(PF,CoxLas)", "Stk(PF,add1)",
            "Stk(CP,Las)", "Stk(CP,CoxLas)", "Stk(CP,add1)",
            "Stk(SS,Las)", "Stk(SS,CoxLas)", "Stk(SS,add1)")) %>%
  filter(!rownames(.) %in% c("LUSC", "OV", "STAD")) %>% 
  tibble::rownames_to_column(var = "Data") %>%
  # select the best AUC
  rowwise() %>%
  mutate(bestIBS = min(across(-c(1, grep("Stk\\(SS", colnames(.))))),
         bestMod = colnames(.)[-c(1, grep("Stk\\(SS", colnames(.)))][which.min(across(-c(1, grep("Stk\\(SS", colnames(.)))))]) %>%
  ungroup() %>%
  mutate(across(-grep("bestMod", colnames(.)), ~ format(.x, digits = 3, nsmall = 3L))) %>%
  mutate(across(everything(), as.character)) %>%
  select(c("Data", "bestMod", "bestIBS", "Stk(SS,IBS)", "Stk(SS,Rdg)", "Stk(SS,CoxRdg)")) %>%
  openxlsx::write.xlsx(., "../results/subsetIBS.xlsx", colNames = T)


#- Boxplot ---------------------------------------------------------------------

barPlotList <- list()

for(i in datNames) {
  
  barData <- aucs[[i]] %>%
    Reduce(rbind, .) %>%
    as.data.frame() %>%
    rename_all( ~ reNames) %>%
    select(-c("Stk(PF,Las)", "Stk(PF,CoxLas)", "Stk(PF,add1)",
              "Stk(CP,Las)", "Stk(CP,CoxLas)", "Stk(CP,add1)",
              "Stk(SS,Las)", "Stk(SS,CoxLas)", "Stk(SS,add1)")) %>% 
    select(-starts_with("Stk(SS,")) %>%
    pivot_longer(cols = colnames(.)) %>%
    mutate(color = case_when(name %in% c("Clin", "Mol", "Naive") ~ "1",
                             name %in% c("FS(las)", "PCR(las)") ~ "2",
                             name %in% c("IPF(las)") ~ "3",
                             grepl("^Stk\\(PF", name) ~ "4",
                             grepl("^Stk\\(CP", name) ~ "5",
                             grepl("^Stk\\(SS", name) ~ "6"),
           name = factor(name, levels = rev(reNames)))
  
  barPlotList[[i]] <- ggplot(barData, aes(x = name, y = value, fill = color)) +
    geom_boxplot(color = "black") +
    labs(title = i,
         x = "Methods",
         y = "AUC") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    # coord_flip() 
  
}

pdf("../results/barPlotAUCV16.pdf", width = 8, height = 6)

for (i in seq_along(barPlotList)) {
  print(barPlotList[[i]])
}

dev.off()

#-------------------------------------------------------------------------------
# Tmpe explore worse results in IBS measure
#-------------------------------------------------------------------------------

ibss <- list()
aucs <- list()

nam = "BLCA"

getData(nam, log2 = T, toBinary = F)

n = nrow(X)
p = ncol(X)

#- get random splits, foldid, times
getSplits(n, nsim, nfolds, seed = 1248372, foldid.internal = T)
time   <- y[ ,1]
status <- y[ ,2]
times <- quantile(time[status == 1], probs = seq(0.1, 0.9, length = 10))

#- 
ibs <- c()
auc <- c()

pdf("../results/survP.pdf", width = 5, height = 5)

for(i in 1:(nfolds * nsim)){
  
  samID <- ceiling(i / nfolds)
  #- j, which fold of the i.nsim split being the test
  j <- i - (samID - 1) * nfolds
  #- foldid for i.nsim splits and j fold
  sam <- sampleSplits[[samID]]
  load(paste0("../output/", nam, "-", samID, "-", j, "-coxOptimV13.RData"))
  
  if (length(mod) == 1)
    print(i)
  
  X.v <- X[sam == j,]
  y.v <- y[sam == j]
  
  X.t <- X[sam != j,]
  y.t <- y[sam != j]
  
  # prob1 
  probs <- predProbs(X.v, mod[-length(mod)], times)
  
  # ipflass
  survProb <- probs$ipf.las
  ibsDat <- matrix(NA, length(times), 2, 
                   dimnames = list(seq_along(times), 
                                   c("ibs", "time")))
  
  Xclin.t = X.t[ ,blocks == 1]
  Xclin.v = X.v[ ,blocks == 1]
  
  # get BS data
  
  getBS <- function(survProb) {
    
    bsDat <- data.frame()
    for(i in seq_along(times)){
      
      bs <- brierScore(Xclin.t, Xclin.v, y.t, y.v, 
                       survProb[,i], times[i], cens.method = "equal")
      
      # ibsDat[i,"time"] <- times[i]
      # ibsDat[i,"ibs"]  <- sum(bs[,"weight"] * (bs[,"Zind"] - bs[,"survProb"]) ^ 2) / sum(bs[,"weight"])
      
      bs[, "time"] <- i
      bsDat <- rbind(bsDat, bs)
      
    }
    
    return(bsDat)
  }
  
  ipfBS <- getBS(probs$ipf.las)
  stkBS <- getBS(probs$stk.pf$`clinicalAdjusted:rdg`)
  
  par(mar = c(2, 2, 2, 2))
  
  plot(ipfBS[,"survProb"], stkBS[ ,"survProb"], 
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "IPF", ylab = "Stk",
       cex = 0.6)
  
  for(i in seq_len(nrow(ipfBS))){
    xval = ipfBS[i, "survProb"]
    yval = stkBS[i, "survProb"]
    lines(c(xval, xval), c(yval, stkBS[i,"Zind"]))
  }
  
  for(i in seq_len(nrow(ipfBS))){
    xval = ipfBS[i, "survProb"]
    yval = stkBS[i, "survProb"]
    lines(c(xval, ipfBS[i,"Zind"]), c(yval, yval), col = "red")
  }
  
  # ibsDat <- na.omit(ibsDat)
  # tRange <- abs(range(ibsDat[,"time"])[1] - range(ibsDat[,"time"])[2])
  # trapezoidal(ibsDat[,"time"], ibsDat[,"ibs"]) / tRange
}
dev.off()


#-------------------------------------------------------------------------------
# Alpha in subset methods
#-------------------------------------------------------------------------------

alphas <- list()

for(nam in datNames) {
  
  alpha <- vector("list", 5)
  for(i in 1:(nfolds * nsim)){
    
    samID <- ceiling(i / nfolds)
    j <- i - (samID - 1) * nfolds
    load(paste0("../output/", nam, "-", samID, "-", j, "-coxOptimV16.RData"))
    a <- mod$stk.subset$alpha
    
    for(k in seq_along(a)) {
      if (names(a)[k] %in% c("logLik:lasso", "logLik:ridge"))
        sub_a <- a[[k]][-1]
      else
        sub_a <- a[[k]]
      
      alpha[[k]] <- rbind(alpha[[k]], sub_a)
    }
    alphas[[nam]] <- alpha
  }
}
 

pdf("../results/subsetAlphaV16.pdf", width = 8, height = 2)

pnam <- c("ibsLoss", "logLik:lasso", "logLik:ridge", "PaLogLik:lasso", "PaLogLik:ridge")

for(i in names(alphas)) {
  
  par(mfrow = c(1, 5), mar = c(2, 2, 2, 2))
  # remove 2, logistic lasso
  for(j in seq_along(pnam)){
    # remove intercept except for solnp
    # note, here add abs for cox
    datPlot <- abs(alphas[[i]][[j]])
    boxplot(datPlot, xlab = "", xaxt = 'n')
    axis(1, at = 1:11, labels = c("Clin", paste0("M", 1:10)))
    mtext(paste0(i, ", (", pnam[j], ")"), side = 3, line = 0.2, adj = 0)
  }
}

dev.off()

#-------------------------------------------------------------------------------
# *PF* in Stk PF
#-------------------------------------------------------------------------------

alphas <- list()

for(nam in datNames) {
  
  alpha <- vector("list", 5)
  for(i in 1:(nfolds * nsim)){
    
    samID <- ceiling(i / nfolds)
    j <- i - (samID - 1) * nfolds
    load(paste0("../output/", nam, "-", samID, "-", j, "-coxOptimV13.RData"))
    a <- mod$stk.pf$alpha
    
    for(k in seq_along(a)) {
      if (names(a)[k] %in% c("clinicalAdjusted:las", "clinicalAdjusted:rdg"))
        sub_a <- a[[k]][-1]
      else
        sub_a <- a[[k]]
      
      alpha[[k]] <- rbind(alpha[[k]], sub_a)
    }
    alphas[[nam]] <- alpha
  }
}


pdf("../results/pfAlphaV13.pdf", width = 8, height = 2)

pnam <- c("IBS", "Las", "Rdg", "coxLas", "coxRdg")

for(i in names(alphas)) {
  
  par(mfrow = c(1, 5), mar = c(2, 2, 2, 2))
  # remove 2, logistic lasso
  for(j in seq_along(pnam)){
    # remove intercept except for solnp
    # note, here add abs for cox
    datPlot <- abs(alphas[[i]][[j]])
    boxplot(datPlot, xlab = "", xaxt = 'n')
    axis(1, at = 1:4, labels = c("1-1", "1-2", "1-4", "1-8"))
    mtext(paste0(i, ", (", pnam[j], ")"), side = 3, line = 0.2, adj = 0)
  }
}

dev.off()

#-------------------------------------------------------------------------------
# *Components* in Stk Comp
#-------------------------------------------------------------------------------

alphas <- list()

for(nam in datNames) {
  
  alpha <- vector("list", 5)
  for(i in 1:(nfolds * nsim)){
    
    samID <- ceiling(i / nfolds)
    j <- i - (samID - 1) * nfolds
    load(paste0("../output/", nam, "-", samID, "-", j, "-coxOptimV13.RData"))
    a <- mod$stk.comp$alpha
    
    for(k in seq_along(a)) {
      if (names(a)[k] %in% c("clinicalAdjusted:las", "clinicalAdjusted:rdg"))
        sub_a <- a[[k]][-1]
      else
        sub_a <- a[[k]]
      
      alpha[[k]] <- rbind(alpha[[k]], sub_a)
    }
    alphas[[nam]] <- alpha
  }
}


pdf("../results/compAlphaV13.pdf", width = 8, height = 2)

pnam <- c("IBS", "Las", "Rdg", "coxLas", "coxRdg")

for(i in names(alphas)) {
  
  par(mfrow = c(1, 5), mar = c(2, 2, 2, 2))
  # remove 2, logistic lasso
  for(j in seq_along(pnam)){
    # remove intercept except for solnp
    # note, here add abs for cox
    datPlot <- abs(alphas[[i]][[j]])
    boxplot(datPlot, xlab = "", xaxt = 'n')
    axis(1, at = 1:2, labels = c("clin", "mol"))
    mtext(paste0(i, ", (", pnam[j], ")"), side = 3, line = 0.2, adj = 0)
  }
}

dev.off()

#-------------------------------------------------------------------------------
#- Clinical Variables being selected
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#- Variables or Sparsity
#-------------------------------------------------------------------------------

getnvars <- function(method, m){
  
  #= m is the mod
  if (method %in% c("stk.pf", "stk.comp")){
    length(unique(unlist(m$beta[m$alpha[-1] != 0]))) - sum(m$alpha[-1] != 0)
  } else if (method %in% c("pcr.las", "pcr.ipf")){
    
    length(m$beta.pcr[m$beta.pcr != 0]) - 1
    
  } else if  (method %in% c("stk.super")){
    NA
  }else{
    length(m$beta[m$beta != 0]) - 1
  }
}

nfolds <- 10
nsim <- 5

d <- dir("../output/")
d <- lapply(strsplit(d[grepl("super", d)], "-"), `[[`, 1)
pass <- table(unlist(d))
pass <- names(pass[pass  == 50 ])
setdiff(nams, pass)

#- tmp

nams = pass
clinvars <- list()
nvars    <- list()

for(nam in nams){
  
  re.nvar <- c()
  re.clin <- c()
  
  for(i in 1:(nfolds * nsim)){
    
    i.nsim <- ceiling(i / nfolds)
    j <- i - (i.nsim - 1) * nfolds
    
    load(paste0("../output/", nam, "-", i.nsim, "-", j, "-bi.addsuper.RData"))
    re.nvar <- rbind(re.nvar, mapply(getnvars, m = mod, method = names(mod)))
    re.clin <- c(re.clin, names(mod$clin$beta)[-1])
  }
  
  nvars[[nam]] <- re.nvar
  clinvars[[nam]] <- table(re.clin)
  # output matrix: dev, ibs, cindex
}

predtable(nvars, digit = 1)

#-------------------------------------------------------------------------------
#- time complexity
#-------------------------------------------------------------------------------

gettime <- function(m){
  
  times <- lapply(m, `[[`, "time")
  times <- lapply(times, as.numeric, units = "secs")
  unlist(times)
}

times <- list()

for(nam in nams){
  
  tim <- c()
  
  for(i in 1:(nfolds * nsim)){
    
    i.nsim <- ceiling(i / nfolds)
    j <- i - (i.nsim - 1) * nfolds
    
    load(paste0("../output/", nam, "-", i.nsim, "-", j, "-bi.addsuper.RData"))
    tim <- rbind(tim, gettime(mod))
  }
  
  times[[nam]] <- tim
}

#- using BLCA as an example

d <- colMeans(times$BLCA)
d <- data.frame(Time = d, 
                Method = c("Clin", "Mol", "Las", 
                           "FS(las)", "FS(offset)", "Pcr(Las)",
                           "IPF(Las)", "FS(IPF)", "Pcr(IPF)", 
                           "Stk(PF)", "Stk(comp)", "Stk(multi)"))
d$Method <- factor(d$Method, d$Method)

library(ggplot2)
p <- ggplot(d, aes(x = Method, y = Time)) +
  geom_segment(aes(x = Method, xend = Method, y = 0, yend = Time),
               color = "black") +
  geom_point(aes(x = Method, y = Time), size = 2) + 
  theme_bw() + 
  # coord_flip() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1100)) + xlab("") +
  theme(axis.text.x = element_text(size = 7)) +
  ylab("Time, Seconds")

ggsave("../output/time.pdf", plot = p, width = 7, height = 3)

#- weights used in alpha -------------------------------------------------------

getweights <- function(mod){
  stkInd <- grepl("stk", names(mod))
  out <- lapply(mod[stkInd], `[[`, "alpha")
  lapply(out, as.matrix)
}

weights <- list()
for(nam in nams){
  
  weight <- c()
  for(i in 1:(nfolds * nsim)){
    
    i.nsim <- ceiling(i / nfolds)
    j <- i - (i.nsim - 1) * nfolds
    load(paste0("../output/", nam, "-", i.nsim, "-", j, "-bi.addsuper.RData"))
    if (i == 1)
      weight = getweights(mod)
    else
      weight = mapply(cbind.data.frame, weight, getweights(mod), SIMPLIFY = F)
  }
  
  weights[[nam]] <- weight
}

#- stkPf
stkPf <- lapply(weights, `[[`, "stk.pf")

#- rename
stkPf <- lapply(stkPf, function(x){
  
  #- remove intercept
  x <- x[-1, ]
  colnames(x) <- paste0("s", seq_len(ncol(x)))
  rownames(x) <- c("1-1", "1-2", "1-4", "1-8", "1-16")
  x
})

#- stkPf
stkComp <- lapply(weights, `[[`, "stk.comp")
stkComp <- lapply(stkComp, function(x){
  
  #- remove intercept
  x <- x[-1, ]
  colnames(x) <- paste0("s", seq_len(ncol(x)))
  rownames(x) <- c("Clin", "Mol")
  x
})

#- stkPf
stkMulti <- lapply(weights, `[[`, "stk.super")
stkMulti <- lapply(stkMulti, function(x){
  
  #- remove intercept
  x <- x[-1, ]
  colnames(x) <- paste0("s", seq_len(ncol(x)))
  rownames(x) <- c("Clin", "Mol", "Las", "IPF", "PCR(Las)")
  x
})


plot.box <- function(obj){
  
  obj <- lapply(obj, function(x){
    
    x %>% as.data.frame %>% rownames_to_column(var = "mods") %>%
      mutate(mods = factor(mods, mods)) %>%
      #- remove intercept
      pivot_longer(starts_with("s"), 
                   names_to = "N",
                   values_to = "alpha")
  })
  
  obj <- bind_rows(obj, .id = "Data")
  
  #- ggplot
  ggplot(obj, aes(x = mods, y = alpha)) + 
    geom_boxplot() + 
    stat_summary(fun = mean, geom = "point", 
                 shape = 18, size = 3, color = "brown1") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~ Data, ncol = 6) + 
    xlab("")
  
}

p <- plot.box(stkPf)
ggsave("../results/alpha.stkPf.pdf", plot = p, width = 6, height = 4)

p <- plot.box(stkComp)
ggsave("../results/alpha.stkComp.pdf", plot = p, width = 6, height = 4)

p <- plot.box(stkMulti)
ggsave("../results/alpha.stkMulti.pdf", plot = p, width = 6, height = 4)


#- plot for prediction performance ---------------------------------------------
#- cindex

#- "hnsc", "lusc", "ov", "blca", "luad"

library(ggplot2)
library(dplyr)

nam = "luad"
d <- reshape2::melt(cindexs[[nam]])
d$Var2 <- gsub(".index", "", d$Var2)

d <- d %>% mutate(Method = ifelse(grepl("las", d$Var2), "Las", "Rdg"),
                  Strategy = case_match(Var2,
                                        c("fs.las", "fs.rdg") ~ "FS",
                                        c("ipf.las", "ipf.rdg") ~ "IPF",
                                        c("pcr.las", "pcr.rdg") ~ "PCA",
                                        c("pls.las", "pls.rdg") ~ "PLS", 
                                        c("vae.las", "vae.rdg") ~ "AE",
                                        c("las", "rdg") ~ "Naive"),
                  Cindex = value) %>% 
  select(Method, Strategy, Cindex)

d$Strategy <- factor(d$Strategy, 
                     levels = c("Naive", "FS", "PCA", "PLS", "AE", "IPF"))

abl <- d %>% filter(Strategy == "Naive") %>% group_by(Method) %>% 
  summarise(median = median(Cindex))

p <- ggplot(d, aes(x = Strategy, y = Cindex, fill = Method)) +
  geom_boxplot() +
  theme_bw() + 
  geom_hline(yintercept = abl$median, linetype = "dashed") + 
  geom_hline(yintercept = 0.5, linetype = "longdash") +
  ggtitle(toupper(nam))

ggsave(paste0("../output/cindex.", nam, ".jpeg"), plot = p, width = 7, 
       height = 4, units = "in", dpi = 400)


#- prediction table

tab <- lapply(cindexs, colMeans)
tab <- round(Reduce(rbind, tab), 2) 
rownames(tab) <- toupper(names(cindexs))
colnames(tab) <- gsub(".index", "", colnames(tab))

write.csv(tab, file = "../output/cindex.csv")

#-------------------------------------------------------------------------------
#- computational time
#-------------------------------------------------------------------------------

s <- c("hnsc", "lusc", "ov", "blca", "luad")

plist <- list()

for(i in s){
  
  d <- times[[i]]
  d <- colSums(d)
  
  d <- data.frame(Time = d, method = names(d))
  d <- d %>% mutate(Method = ifelse(grepl("las", d$method), "Las", "Rdg"),
                    Strategy = case_match(method,
                                          c("fs.las", "fs.rdg") ~ "FS",
                                          c("ipf.las", "ipf.rdg") ~ "IPF",
                                          c("pcr.las", "pcr.rdg") ~ "PCA",
                                          c("pls.las", "pls.rdg") ~ "PLS", 
                                          c("vae.las", "vae.rdg") ~ "AE",
                                          c("las", "rdg") ~ "Naive")) %>% 
    select(Method, Strategy, Time)
  
  d$Strategy <- factor(d$Strategy, 
                       levels = rev(c("FS", "PCA", "PLS", "AE", "IPF")))
  
  plist[[i]] <- ggplot(d, aes(x = Strategy, y = Time)) +
    geom_segment(aes(x = Strategy, xend = Strategy, y = 0, yend = Time),
                 color = "black") +
    geom_point(aes(x = Strategy, y = Time, color = Method), size = 2) + 
    theme_bw() + 
    coord_flip() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5)) + xlab("") +
    ggtitle(i)
}

library(ggpubr)
p <- ggarrange(plotlist = plist, nrow = 1, common.legend = TRUE, legend="bottom")

ggsave("../output/time.jpeg", p, width = 9, height = 3, units = "in",
       dpi = 400)

#-------------------------------------------------------------------------------
#- Cuts for FS
#-------------------------------------------------------------------------------

tab <- function(x){
  
  s <- setNames(rep(0, 7), 
                c(0.001, 0.003, 0.005, 0.01, 0.03, 0.1, 0.3))
  
  f <- prop.table(table(x))
  s[names(f)] <- f
  s
}


jpeg("../output/cuts.jpeg", width = 9, height = 5, units= "in", res = 300)
par(mar = c(1, 3, 1, 1), oma = c(3, 1 ,1, 1))
mat <- matrix(1:10, nrow = 2, byrow = F)
layout(mat, heights = c(1, 1))

for(i in names(cuts)){
  
  freq <- apply(cuts[[i]], 2, tab, simplify = F)
  barplot(freq$Las, ylim = c(0, 0.6), col = "#00BFC4", xaxt='n')
  box()
  mtext("Frequency", side = 2, line = 2)
  mtext(toupper(i), side = 3, line = 0.2, adj = 0)
  
  
  barplot(freq$Rdg, ylim = c(0, 0.6), col = "#F8766D", las =2)
  box()
  mtext("Frequency", side = 2, line = 2)
}
dev.off()



#-------------------------------------------------------------------------------
#- repeated cross-validation
#-------------------------------------------------------------------------------

library(doParallel)
library(doRNG)
library(doSNOW)

# - 

for(nam in c("hnsc", "lusc", "ov", "blca", "luad")){
  
  # nsim <- 10
  nsim <- 2
  nfolds <- 10
  mol <- "rna" #-rna
  set.seed(1243872)
  get.data(nam, nfolds = nfolds, type = "cv", nsim = nsim, mol = mol)
  
  cores <- detectCores() - 2
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(min = 1, max = (nfolds * nsim), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cv.mods <- foreach(i = 1 : (nfolds * nsim), .packages = c("glmnet", "fav"), 
                     .options.snow = opts) %dopar% {
                       
                       i.nsim <- ceiling(i / nfolds)
                       
                       #- j, which fold of the i.nsim split being the test
                       j <- i - (i.nsim - 1) * nfolds
                       sam <- split.sample[[i.nsim]]
                       #- foldid for i.nsim splits and j fold
                       f <- foldids[[i.nsim]][[j]]

                       X.t <- X[sam != j,]
                       y.t <- y[sam != j]
                       
                       sel <- 1:ncol(X.t)
                       re1 <- tryCatch(fitfun(X.t[,sel], y.t, foldid = f, 
                                              blocks = blocks[sel], alpha = 1),
                                       error = function(e) "error")
                       
                       re0 <- tryCatch(fitfun(X.t[,sel], y.t, foldid = f, 
                                              blocks = blocks[sel], alpha = 0),
                                       error = function(e) "error")
                       
                       return(list(las.type = re1, rdg.type = re0))
                     }
  
  close(pb)
  cv.mods <- list(cv.mods = cv.mods,
                  split.sample = split.sample,
                  nsim = nsim,
                  foldids = foldids,
                  mol = mol)
  
  save(list = "cv.mods", 
       file = paste0("./output/lasrdg-", mol, "-", nam, ".RData"))
  stopCluster(cl)
}

#-------------------------------------------------------------------------------
#- PF in Fav
#-------------------------------------------------------------------------------

mol <- "rna"
pf.dat <- list()
clin.dat <- list()
#- get lam and pf of favoring
for(nam in c("hnsc", "lusc", "ov", "blca", "luad")){
  load(paste0("./output/cv-", mol, "-", nam, "1.RData"))
    mods <- cv.mods$cv.mods
    
    pf <- c()
    nclin <- c()
    for(i in seq_along(mods)){
      pf  <- c(pf, paste0(mods[[i]]$fav$pf, collapse = "-"))
      m <- mods[[i]][["fav"]][["ipf"]][[pf[i]]]$beta
      nclin <- c(nclin, sum(grepl("clinical", names(m))))
    }
    
    pf.dat[[nam]] <- pf
    clin.dat[[nam]] <- nclin
}


jpeg("./output/pf.jpeg", width = 8, height = 3, units = "in", res = 300)
par(mar = c(2, 2, 2, 1), mfrow = c(1, 5), oma = c(1, 2, 1, 1))

for(nam in c("hnsc", "lusc", "ov", "blca", "luad")){
  
  tab <- prop.table(table(pf.dat[[nam]]))
  tab.t <- setNames(c(0, 0, 0, 0), c("1-1", "1-2", "1-4", "1-8"))
  tab.t[names(tab)] <- tab
  
  m <- pf.dat[[nam]]
  p <- clin.dat[[nam]]

  tab <- tapply(p, m, "mean")
  tab.c <- setNames(c(0, 0, 0, 0), c("1-1", "1-2", "1-4", "1-8"))
  tab.c[names(tab)] <- tab
  
  barplot(tab.t, ylim = c(0, 1.0), las = 2)
  box()
  mtext(nam, side = 3, line = 0.2)
  s <- tab.t != 0
  text(c(0.5, 1.75, 3, 4.25)[s], rep(1, 4)[s], 
       labels = format(round(tab.c[s], 1), nsmall = 1), pos = 1)
  
  if (nam == "hnsc")
    mtext("Proportion of PF, %", side = 2, line = 2.2)
}
dev.off()

#-------------------------------------------------------------------------------
#- lambda distribution in real data
#-------------------------------------------------------------------------------

mol <- "rna"
lambda.dat <- list()

for(nam in c("hnsc", "lusc", "ov", "blca", "luad")){
  
  load(paste0("./output/cv-", mol, "-", nam, ".RData"))
  mods <- cv.mods$cv.mods
  
  lam <- c()
  for(i in seq_along(mods)){
    lam <- rbind(lam, unlist(lapply(mods[[i]], `[[`, "lambda")))
  }
  lambda.dat[[nam]] <- lam
}

jpeg("./output/lambda.jpeg", width = 7, height = 9, res = 300, units = "in")
par(mfrow = c(5, 3), mar = c(1.5, 1.5, 1.5, 1.5), oma = c(2, 3, 2, 2))

for(i in 1:5){
  
  lam <- lambda.dat[[i]]
  l <- c(0, quantile(lam, 0.975))
  breaks = c(quantile(l, seq(0, 1, 0.05)), max(lam))

  p <- hist(lam[,1], plot = F, breaks = breaks)
  p$density <- p$counts / nrow(lam)
  plot(p, xlim = l, ylim = c(0, 1),
       col = NULL, main = "", border = "green")
  
  p <- hist(lam[,2], plot = F, breaks = breaks)
  p$density <- p$counts / nrow(lam)
  plot(p, xlim = l, ylim = c(0, 1),
       col = NULL, main = "", border = "red", add = T)
  box()
  
  mtext(c("HNSC", "LUSC", "OV", "BLCA", "LUAD")[i], line = 2, side = 2)
  
  if (i == 1){
    mtext("Clin & mol", side = 3, line = 0.2)
  }
  
  p <- hist(lam[, 3], plot = F, breaks = breaks)
  p$density <- p$counts / 100
  plot(p, xlim = l, ylim = c(0, 1),
       col = NULL, main = "", border = "black")
  box()
  if (i == 1){
    mtext("Naive", side = 3, line = 0.2)
  }
  
  p <- hist(lam[, 4], plot = F, breaks = breaks)
  p$density <- p$counts / 100
  plot(p, xlim = l, ylim = c(0, 1),
       col = NULL, main = "",  border = "green")
  
  p <- hist(lam[, 5], plot = F, breaks = breaks)
  p$density <- p$counts / 100
  plot(p, xlim = l, ylim = c(0, 1),
       col = NULL, main = "", border = "red", add = T)
  box()
  if (i == 1){
    mtext("Favoring", side = 3, line = 0.2)
  }
}
dev.off()

#-------------------------------------------------------------------------------
#- selection of clin and ipf
#-------------------------------------------------------------------------------
#- S.clin, S.mol, IBS

get.n.var <- function(beta){
  
  beta <- names(beta)
  ind <- grepl("clinical",beta)
  c(sum(ind), sum(!ind))
}

mol <- "rna"
dat.vars <- list()
for(nam in c("hnsc", "lusc", "ov", "blca", "luad")){
  
  load(paste0("./output/cv-", mol, "-", nam, ".RData"))
  mods <- cv.mods$cv.mods
  tmp <- c()
  for(i in seq_along(mods)){
    betas <- lapply(mods[[i]], `[[`, "beta")
    n.vars <- unlist(lapply(betas, get.n.var))
    tmp <- rbind(tmp, n.vars)
  }  
  dat.vars[[nam]] <- tmp
}

dat <- c()
for(i in 1:5){
  
  re <- dat.vars[[i]]
  m  <- colMeans(re)
  sd <- apply(re, 2, sd)
  out <- paste0(round(m, 1), "(", round(sd, 1), ")")
  dat <- rbind(dat, out)
}

rownames(dat) <- c("hnsc", "lusc", "ov", "blca", "luad")
write.csv(dat, file = "./output/nvar.csv")

#-------------------------------------------------------------------------------
#- deviance changes
#-------------------------------------------------------------------------------

ibs.list <- list()
dev.list <- list()
for(nam in c("hnsc", "lusc", "ov", "blca", "luad")){
  
  print(nam)
  mol  <- "rna"
  nsim <- 10
  nfolds = 10
  #- para
  load(paste0("./output/cv-", mol, "-", nam, ".RData"))
  set.seed(1243872)
  get.data(nam, nfolds = 10, type = "cv", nsim = nsim, mol = mol)

  re.ibs <- c()
  re.dev <- c()
  for(i in 1:(nfolds * nsim)){ 
    
    obj <- cv.mods[["cv.mods"]][[i]]
    i.nsim <- ceiling(i / nfolds)
    #- used fold as validation
    j <- i - (i.nsim - 1) * nfolds
    sam <- cv.mods$split.sample[[i.nsim]]

    betas <- lapply(obj, `[[`, "beta")
    # betas <- c(null = list(NULL), betas)
    X.t <- X[sam != j, ]  
    y.t <- y[sam != j]
    
    X.v <- X[sam == j,]
    y.v <- y[sam == j]
    re.ibs <- rbind(re.ibs, performance(betas, X.t, y.t, X.v, y.v, 
                                        method = c("ibs")))
    re.dev <- rbind(re.dev, performance.dev(betas, X.v, y.v))
    
  }
  ibs.list[[nam]] <- re.ibs
  dev.list[[nam]] <- re.dev
}

#- plot
jpeg("./output/ibs.dev.jpeg", width = 6, height = 5, res = 300, units = "in")

nam <- c("HNSC", "LUSC", "OV", "BLCA", "LUAD")
layout(matrix(1:10, nrow = 2, byrow = F))
par(mar = c(2, 2, 2, 1), oma = c(1, 2, 1, 1))

for(i in 1:5){
  
  dev <- dev.list[[i]]
  
  boxplot(dev, las = 2)
  abline(h = 0)
  abline(h = median(dev[,1]), lty = "dashed")
  mtext(nam[i], side = 3, line = 0.2)

  
  if (i == 1)
    mtext("Dev", side = 2, line = 2.6)
  
  ibs <- ibs.list[[i]]
  
  colnames(ibs) <- c("clin", "mol", "nav", "fav")
  boxplot(ibs, las = 2)
  abline(h = median(ibs[,1]), lty = "dashed")
  if (i == 1)
    mtext("IBS", side = 2, line = 2.6)
}

dev.off()

#- Prove, CV needs to run several times

g <- rep(1:10, each =  10)
m <- list()
for(i in c("hnsc", "lusc", "ov", "blca", "luad")){
  m[[i]] <- round(apply(dev.list[[i]], 2, function(x) tapply(x, g, median)), 2)
}

#-------------------------------------------------------------------------------
















