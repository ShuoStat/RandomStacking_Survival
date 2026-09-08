# Random Survival Stacking

**Random stacking for high-dimensional survival data**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains the code, data, and results associated with the development and evaluation of **Random Survival Stacking**, a randomization-based stacking framework for high-dimensional time-to-event data.

The proposed framework combines **bootstrap resampling, random feature selection, prediction aggregation, and survival-specific meta-learning** to improve the predictive accuracy and stability of stacking ensembles in high-dimensional survival settings.

---

## Overview

Stacking is a widely used ensemble learning framework that combines predictions from multiple base learners through a meta-learner. While traditional stacking has been successfully applied in many prediction settings, it can become vulnerable to **overfitting and high prediction variance in high-dimensional data**.

These challenges are particularly important for survival prediction, where the outcome is a time-to-event variable subject to censoring and the prediction target evolves over time.

To address these challenges, we developed **Random Survival Stacking**, a novel ensemble framework specifically designed for high-dimensional survival data.

The framework introduces two major innovations:

1. **Randomization and aggregation at the base-learning layer**

   * Bootstrap resampling is used to generate more diverse training samples.
   * Random feature selection is introduced to reduce dependence among base learners.
   * Predictions from multiple randomized models are aggregated to obtain more stable base-layer predictions.

2. **A survival-specific meta-learning framework**

   * The stacking framework is extended to time-to-event outcomes.
   * New meta-learners based on the **negative log-likelihood (NLL)** are developed.
   * **Ridge regularization** is incorporated into the meta-layer to improve stability in high-dimensional prediction settings.

The framework is designed particularly for applications involving **complex biomedical data**, including the integration of clinical and omics data for prognostic prediction.

---

## Method

### Random Survival Stacking

The overall framework can be summarized as:

```text
                    Training Data
                         │
              ┌──────────┴──────────┐
              │                     │
       Bootstrap Sampling     Random Feature Selection
              │                     │
              └──────────┬──────────┘
                         │
                Multiple Base Learners
                         │
              ┌──────────┴──────────┐
              │                     │
        Base Prediction 1    Base Prediction B
              │                     │
              └──────────┬──────────┘
                         │
                 Prediction Aggregation
                         │
                  Meta-Learning Layer
                         │
                  NLL-Ridge Meta-Learner
                         │
                         ▼
                Survival Prediction
```

Compared with conventional stacking, Random Survival Stacking introduces additional sources of randomness before the meta-learning stage. This is intended to generate more diverse and less correlated base predictions and thereby reduce ensemble variance.

### NLL-Ridge Meta-Learner

For survival outcomes, the meta-layer is based on the negative log-likelihood loss with ridge regularization.

The ridge penalty provides additional control of the meta-layer coefficients and helps prevent overfitting when the base predictions are highly correlated or when the number of candidate predictors is relatively large.

---

## Why Randomization?

Traditional stacking typically obtains base predictions using cross-validation. However, the training datasets used in different folds overlap substantially, which can result in highly correlated fitted models and therefore highly correlated base predictions.

Random Survival Stacking introduces additional randomization through:

* **Bootstrap resampling**
* **Random feature selection**
* **Prediction aggregation**

This creates greater diversity among the base learners.

The underlying motivation is:

> **More diverse base learners → less correlated predictions → improved ensemble stability → reduced prediction variance.**

However, randomization is not universally beneficial. When the underlying ensemble is already highly stable—for example, with large sample sizes—additional randomization may introduce unnecessary bias and can result in performance inferior to conventional stacking.

---

## Data Applications

The proposed method was evaluated using both **real-world cancer datasets** and **simulation studies**.

### Clinical and Omics Data Integration

A major application is the integration of heterogeneous biomedical information for survival prediction, including:

* Clinical variables
* Molecular/omics features
* High-dimensional genomic information
* Other potentially complementary predictors

The framework allows different data sources to be modeled using different base learners and subsequently integrated through the stacking framework.

---

## Simulation Studies

Extensive simulation studies were conducted to investigate the behavior of Random Survival Stacking under different data-generating conditions.

The simulations examine factors including:

* Sample size
* Signal-to-noise ratio
* Correlation among predictors/base learners

The simulation framework was designed to identify both the **advantages and limitations** of random survival stacking.

In general, Random Survival Stacking showed the greatest advantage in settings where conventional stacking was relatively unstable, such as:

* Small sample sizes
* Weak signals
* Greater prediction variance

In contrast, when the ensemble was already highly stable, particularly under large-sample settings, traditional stacking could perform better.

---

## Results

### Real-world datasets

Across the evaluated real-world cancer datasets, Random Survival Stacking generally demonstrated:

* Improved predictive accuracy
* Greater prediction stability
* More robust performance across datasets

Compared with traditional stacking, the proposed method showed consistent improvements in survival prediction performance across the evaluated datasets.

### Simulation studies

The simulation results further demonstrated that Random Survival Stacking outperformed traditional stacking and other competing approaches in the majority of evaluated scenarios.

However, the advantage of randomization depended on the stability of the underlying ensemble.

In particular:

> Random Survival Stacking may be less advantageous when the sample size is large and the conventional ensemble is already highly stable.

This finding highlights an important practical consideration: **randomization should be viewed as a strategy for controlling instability rather than as a universally superior replacement for conventional stacking.**

---

## Main Contributions

The main contributions of this work are:

### 1. Randomized base-layer stacking

We introduce bootstrap resampling and random feature selection into the base-learning layer of stacking to generate more diverse base predictions.

### 2. Prediction aggregation

Predictions from randomized base learners are aggregated to reduce prediction variance and improve stability.

### 3. Survival-specific meta-learning

We extend the stacking framework to time-to-event outcomes using survival-specific meta-learning based on negative log-likelihood loss.

### 4. NLL-Ridge meta-learner

We develop a ridge-regularized NLL-based meta-learner that provides additional control against overfitting and correlated base predictions.

### 5. Extensive empirical evaluation

The proposed framework is evaluated using both real-world cancer datasets and extensive simulation studies covering a range of sample sizes, signal strengths, correlations, and ensemble stability conditions.

---

## Repository Structure

The repository is organized as follows:

```text
RandomStacking_Survival/
│
├── codes/
│   ├── ...
│   └── ...
│
├── data/
│   ├── ...
│   └── ...
│
├── output/
│   ├── ...
│   └── ...
│
├── results/
│   ├── ...
│   └── ...
│
├── LICENSE
└── README.md
```

### `codes/`

Contains the analysis and simulation scripts used to implement Random Survival Stacking and competing methods.

### `data/`

Contains datasets or processed data required for the analyses, where redistribution is permitted.

> **Note:** Raw clinical or genomic datasets subject to third-party access restrictions are not redistributed in this repository. Users should obtain such datasets directly from the corresponding data providers.

### `output/`

Contains intermediate outputs generated during the computational analyses.

### `results/`

Contains summarized results, tables, figures, and other outputs used to evaluate the proposed framework.


---

## When Should Random Survival Stacking Be Used?

Random Survival Stacking is particularly motivated for:

* High-dimensional survival data
* Omics-based survival prediction
* Clinical + molecular data integration
* Small-to-moderate sample sizes
* Unstable base learners
* Settings with substantial prediction variance
* Biomedical prediction problems involving heterogeneous data sources

It may provide less benefit when:

* The sample size is very large;
* The base learners are already highly stable;
* Base predictions are strongly correlated;
* Additional randomization introduces substantial bias.

Therefore, the proposed method should be considered a **variance-control strategy for stacking**, rather than a universally superior stacking algorithm.

---


## License

This project is released under the **MIT License**. See [LICENSE](LICENSE) for details. The repository currently identifies itself as MIT-licensed.

---

## Contact

For questions, suggestions, or issues related to Random Survival Stacking, please open an issue in this repository or contact the corresponding author (wangsures@foxmail.com).

---

## Acknowledgements

We thank the researchers and data providers who contributed to the publicly available datasets used in this study.

---

**By the way, we already have an R implementation of Random Stacking, but unfortunately, we don’t currently have enough time to transfer the methods to Python. Please feel free to contact me if you’re interested in this work. **

---

**Random Survival Stacking**
*A randomized ensemble framework for stable and accurate high-dimensional survival prediction.*
