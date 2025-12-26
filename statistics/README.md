\# Disease Prevalence Estimation with Imperfect Diagnostic Tests



This project studies the estimation of disease prevalence in a population using

diagnostic tests with imperfect sensitivity and specificity. The analysis

compares estimation and inference methods under perfect and imperfect testing

scenarios, with a strong emphasis on uncertainty quantification.



\## Methodology

\- Binomial modeling of diagnostic test outcomes

\- Maximum Likelihood Estimation under perfect testing

\- Method of Moments estimation under imperfect sensitivity and specificity

\- Monte Carlo simulation to assess bias, variance and mean squared error

\- Bootstrap methods for confidence interval construction

\- Hypothesis testing for prevalence differences in pre–post intervention studies



\## Analysis

The performance of different estimators is evaluated through theoretical

derivations and simulation studies. Special attention is given to small sample

behavior, estimator truncation, and coverage properties of confidence intervals.

Scenarios with and without measurement error are compared.



\## Results

Results show that ignoring measurement error leads to biased prevalence

estimates. The method of moments estimator corrects for this bias but may produce

inadmissible values in small samples, motivating the use of a truncated

estimator. Bootstrap confidence intervals provide improved uncertainty

quantification for small sample sizes.



\## Repository Structure

\- `data/`: directory reserved for datasets (not used in this project, as all

&nbsp;  data are generated through simulation)

\- `figures/`: simulation results and diagnostic plots

\- `analysis.R`: R script implementing the estimators and simulations

\- `report/`: original assignment PDF (Spanish) for reference; 

&nbsp; all analysis and figures are self-contained in this repository





