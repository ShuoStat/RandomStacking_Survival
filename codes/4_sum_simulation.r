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
library(ggthemes)
library(glmnet)
library(survival)
library(dplyr)
library(tidyverse)
library(doParallel)
library(doRNG)
library(doSNOW)

# ===============================================================================


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

whichRuns <- c("competitors", "stk", "random.stk")
# whichRuns = c("competitors")
#

nCores <- 10
nsim <- 100
updateResults <- TRUE

for (whichRun in whichRuns) {
    ###
    ibss <- list()
    aucs <- list()
    logLiks <- list()

    # set.seed(321)
    # datNames <- c("BLCA", "HNSC")

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

        set.seed(1248372)
        testData <- simData(
            n.obs = 5000,
            t.gene = 10000,
            r = r,
            sd_sett = sd_sett,
            weib_shape = 1.5,
            weib_scale = 0.1,
            censoring_target = 0.2
        )

        timesRange <- quantile(testData$time[testData$status == 1], probs = c(0.05, 0.95))
        times <- seq(timesRange[1], timesRange[2], length = 20)
        y.v <- survival::Surv(testData$time, testData$status)

        # generate test data3
        rng <- RNGseq(nsim, 6534125)
        Cores <- pmin(nCores, detectCores() - 1)
        cl <- makeCluster(Cores)
        registerDoSNOW(cl)

        pb <- txtProgressBar(min = 1, max = nsim, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)

        res <- foreach(
            i = 1:nsim,
            .packages = c("glmnet", "survival", "dplyr"),
            .options.snow = opts,
            .verbose = TRUE
        ) %dopar% {
            rngtools::setRNG(rng[[i]])
            data <- simData(
                n.obs = n,
                t.gene = 10000,
                r = r,
                sd_sett = sd_sett,
                weib_shape = 1.5,
                weib_scale = 0.1,
                censoring_target = 0.2
            )

            blocks <- data$blocks
            y.t <- survival::Surv(data$time, data$status)

            # import data
            nam <- paste0(n, "_", snr, "_", r, "_", i)
            load(paste0("../output/sim_", whichRun, "_", nam, ".RData"))
            allStks <- grep("stk", names(mod), value = TRUE)
            probs <- predProbs(testData$X, mod, times, ifCVMods = FALSE)

            # faltten list
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
                Xclin.t = data$X[, blocks == 1],
                Xclin.v = testData$X[, blocks == 1],
                y.t = y.t,
                y.v = y.v,
                times = times,
                cens.method = "equal",
                sumRes = FALSE
            )

            aucTime <- lapply(probs, getTimeAUC,
                Xclin.t = data$X[, blocks == 1],
                Xclin.v = testData$X[, blocks == 1],
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

        nam <- paste0(n, "_", snr, "_", r)

        # from res to ibss, aucs, logLiks
        ibss[[nam]] <- lapply(res, `[[`, "ibsTime")
        aucs[[nam]] <- lapply(res, `[[`, "aucTime")
        logLiks[[nam]] <- lapply(res, `[[`, "logLikTime")

        stopCluster(cl)
    }

    if (updateResults) {
        save(list = c("aucs", "ibss", "logLiks"), file = paste0("../output/res_sim_", whichRun, ".RData"))
    }
}


sumFun <- function(objs, fun = colMeans) {
    # objs, aucs, ibss, or logLiks
    areaSum <- function(obj) {
        obj <- na.omit(obj)
        tRange <- abs(range(obj[, "time"])[1] - range(obj[, "time"])[2])
        trapezoidal(obj[, "time"], obj[, 2]) / tRange
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
round(sumFun(logLiks, fun = colMeans), 4)

sumResults <- function(target, whichRun) {
    reMat <- c()
    reMatSe <- c()

    for (i in seq_along(whichRun)) {
        load(paste0("../output/res_sim_", whichRun[i], ".RData"))
        n <- length(get(target)[[1]])

        # remove duplicated clin
        tmp <- sumFun(get(target), fun = colMeans)
        if (i != 1) {
            tmp <- tmp[, -1]
        }

        # calcualte SE
        reMat <- cbind(reMat, tmp)

        tmp <- sumFun(get(target), fun = function(x) apply(x, 2, sd) / sqrt(n))
        if (i != 1) {
            tmp <- tmp[, -1]
        }

        reMatSe <- cbind(reMatSe, tmp)
    }

    list(re = reMat, re.se = reMatSe)
}


re_auc <- sumResults(target = "aucs", whichRun = c("competitors", "stk", "random.stk"))
re_auc

re_ibs <- sumResults(target = "ibss", whichRun = c("competitors", "stk", "random.stk"))
re_ibs

re_logLik <- sumResults(target = "logLiks", whichRun = c("competitors", "stk", "random.stk"))
re_logLik


# Visualization
# Please make a plot to compare the results
# Point Plot to compare the mean of each model

plot_point <- function(df,
                       x_var = "mean",
                       y_var = "model",
                       ci_low = "lower", # column for CI lower bound
                       ci_high = "upper", # column for CI upper bound
                       group_var = NULL, # column for grouping (color mapping)
                       xlab = "Mean C-index",
                       title = "Model Performance Comparison",
                       point_size = 1,
                       line_width = 0.6,
                       brewer_palette = "Set1",
                       base_size = 14,
                       show_legend = FALSE) {
    # Base aesthetic mapping
    aes_mapping <- aes(y = .data[[x_var]], x = .data[[y_var]])

    # If group variable is provided, map color to it
    if (!is.null(group_var)) {
        aes_mapping <- modifyList(
            aes_mapping,
            aes(colour = .data[[group_var]])
        )
    }

    p <- ggplot(df, aes_mapping) +

        # Point range: point + confidence interval
        geom_pointrange(
            aes(ymin = .data[[ci_low]], ymax = .data[[ci_high]]),
            size = point_size, linewidth = line_width
        ) +

        # Scales & theme
        labs(y = xlab, x = "", title = title) +
        theme_bw(base_size = base_size) +
        theme(
            plot.title = element_text(hjust = 0, face = "bold"),
            axis.text.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, margin = margin(t = -2)),
            legend.title = element_blank(),
            legend.position = if (show_legend) "bottom" else "none",
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank()
        )

    # Add brewer color scale if group variable is present
    if (!is.null(group_var)) {
        p <- p + scale_color_brewer(palette = brewer_palette)
    }

    return(p)
}


plotData <- function(re, n, snr, r, sel, sel_title) {
    # re, object of class list, with elements re and re.se
    # n, number of samples
    # snr, signal-to-noise ratio
    # sel, vector of models to select
    # r, random seed
    nam <- paste0(n, "_", snr, "_", r)
    re_mean <- re$re[nam, sel]
    re_se <- re$re.se[nam, sel]

    df <- data.frame(mean = re_mean, lower = re_mean - 1.96 * re_se, upper = re_mean + 1.96 * re_se, model = sel_title)
    df$model <- factor(df$model, levels = sel_title)

    return(df)
}

sel <- c("clin", "mol", "naive", "pcr.las", "ipf.las", "stk.comp:ibsLoss", "stk.comp:logLik:ridge", "Rstk.comp:ibsLoss", "Rstk.comp:logLik:ridge")
sel_title <- c("Clin", "Mol", "Naive", "PCR(Las)", "IPFLasso", "Stk(IBS-Constrained)", "Stk(NLL-Ridge)", "Rstk(IBS-Constrained)", "Rstk(NLL-Rridge)")


# The effects of Sample Size on Model Performance

summary_df <- plotData(re_auc, n = 100, snr = 1, r = 0.3, sel = sel, sel_title = sel_title)
summary_df$group <- case_when(
    grepl("^Stk", summary_df$model) ~ "Stk",
    grepl("^Rstk", summary_df$model) ~ "Rstk",
    TRUE ~ summary_df$model
)

p1 <- plot_point(
    df = summary_df,
    x_var = "mean",
    y_var = "model",
    ci_low = "lower",
    ci_high = "upper",
    group_var = "group",
    show_legend = FALSE,
    brewer_palette = "Paired",
    xlab = "AUC",
    title = "A, N = 100"
)

#
summary_df <- plotData(re_auc, n = 500, snr = 1, r = 0.3, sel = sel, sel_title = sel_title)
summary_df$group <- case_when(
    grepl("^Stk", summary_df$model) ~ "Stk",
    grepl("^Rstk", summary_df$model) ~ "Rstk",
    TRUE ~ summary_df$model
)

p2 <- plot_point(
    df = summary_df,
    x_var = "mean",
    y_var = "model",
    ci_low = "lower",
    ci_high = "upper",
    group_var = "group",
    show_legend = FALSE,
    brewer_palette = "Paired",
    xlab = "AUC",
    title = "B, N = 500"
)

p <- ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1)
p

ggsave(p, file = "../results/sample_size.pdf", width = 8, height = 5)
ggsave(p, file = "../results/sample_size.tiff", width = 8, height = 5, dpi = 600, compress = "lzw")

# The effects of Signal-to-Noise Ratio on Model Performance


# The effects of Sample Size on Model Performance

summary_df <- plotData(re_auc, n = 100, snr = 1, r = 0.3, sel = sel, sel_title = sel_title)
summary_df$group <- case_when(
    grepl("^Stk", summary_df$model) ~ "Stk",
    grepl("^Rstk", summary_df$model) ~ "Rstk",
    TRUE ~ summary_df$model
)

p1 <- plot_point(
    df = summary_df,
    x_var = "mean",
    y_var = "model",
    ci_low = "lower",
    ci_high = "upper",
    group_var = "group",
    show_legend = FALSE,
    brewer_palette = "Paired",
    xlab = "AUC",
    title = "A, SNR = 1"
)

#
summary_df <- plotData(re_auc, n = 100, snr = 2.33, r = 0.3, sel = sel, sel_title = sel_title)
summary_df$group <- case_when(
    grepl("^Stk", summary_df$model) ~ "Stk",
    grepl("^Rstk", summary_df$model) ~ "Rstk",
    TRUE ~ summary_df$model
)

p2 <- plot_point(
    df = summary_df,
    x_var = "mean",
    y_var = "model",
    ci_low = "lower",
    ci_high = "upper",
    group_var = "group",
    show_legend = FALSE,
    brewer_palette = "Paired",
    xlab = "AUC",
    title = "B, SNR = 2.33"
)

p <- ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1)
p

ggsave(p, file = "../results/snr.pdf", width = 8, height = 5)
ggsave(p, file = "../results/snr.tiff", width = 8, height = 5, dpi = 600, compress = "lzw")


# The effects of correlation on Model Performance

summary_df <- plotData(re_auc, n = 100, snr = 1, r = 0, sel = sel, sel_title = sel_title)
summary_df$group <- case_when(
    grepl("^Stk", summary_df$model) ~ "Stk",
    grepl("^Rstk", summary_df$model) ~ "Rstk",
    TRUE ~ summary_df$model
)

p1 <- plot_point(
    df = summary_df,
    x_var = "mean",
    y_var = "model",
    ci_low = "lower",
    ci_high = "upper",
    group_var = "group",
    show_legend = FALSE,
    brewer_palette = "Paired",
    xlab = "AUC",
    title = "A, r = 0"
)

#
summary_df <- plotData(re_auc, n = 100, snr = 1, r = 0.3, sel = sel, sel_title = sel_title)
summary_df$group <- case_when(
    grepl("^Stk", summary_df$model) ~ "Stk",
    grepl("^Rstk", summary_df$model) ~ "Rstk",
    TRUE ~ summary_df$model
)

p2 <- plot_point(
    df = summary_df,
    x_var = "mean",
    y_var = "model",
    ci_low = "lower",
    ci_high = "upper",
    group_var = "group",
    show_legend = FALSE,
    brewer_palette = "Paired",
    xlab = "AUC",
    title = "B, r = 0.3"
)


summary_df <- plotData(re_auc, n = 100, snr = 1, r = 0.7, sel = sel, sel_title = sel_title)
summary_df$group <- case_when(
    grepl("^Stk", summary_df$model) ~ "Stk",
    grepl("^Rstk", summary_df$model) ~ "Rstk",
    TRUE ~ summary_df$model
)

p3 <- plot_point(
    df = summary_df,
    x_var = "mean",
    y_var = "model",
    ci_low = "lower",
    ci_high = "upper",
    group_var = "group",
    show_legend = FALSE,
    brewer_palette = "Paired",
    xlab = "AUC",
    title = "C, r = 0.7"
)

p <- ggpubr::ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
p

ggsave(p, file = "../results/cor.pdf", width = 12, height = 5)
ggsave(p, file = "../results/cor.tiff", width = 12, height = 5, dpi = 600, compress = "lzw")



