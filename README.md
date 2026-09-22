# Parameter Sensitivity and Shock Transmission in a Stochastic RBC Model

This project develops a stochastic real business cycle (RBC) model in Dynare and
studies how household preference parameters and shock persistence shape
macroeconomic transmission. The baseline analysis is calibrated and is extended
in two directions: a fiscal block with stochastic government purchases and a
separate Bayesian estimation exercise using quarterly euro-area data.

The project is designed around five research questions:

1. How does the intertemporal elasticity of substitution affect consumption
   smoothing and capital accumulation?
2. How does the Frisch elasticity affect the labour-market response to shocks?
3. How does shock persistence change the magnitude and duration of business-cycle
   fluctuations?
4. How do productivity and labour-preference shocks generate different
   macroeconomic dynamics?
5. Which shocks account for fluctuations in output, consumption, employment and
   investment inside the calibrated model?

## Provenance and scope

The baseline RBC specification originated in Advanced Macroeconomics coursework.
The consolidated implementation, systematic sensitivity analysis, response
metrics, two-dimensional parameter grid, reproducible stochastic simulation and
government-purchases extension were developed as an independent follow-up.

The core and fiscal exercises are calibrated rather than estimated. Their
variance decompositions therefore describe the model under the assumed
calibration and must not be interpreted as empirical estimates of observed
business cycles. A separate estimation extension confronts the RBC model with
euro-area data and estimates selected shock-process parameters using Bayesian
methods.

## Economic environment

The representative household chooses consumption, employment and capital. The
firm has Cobb-Douglas production and capital depreciates geometrically. The core
model contains two AR(1) disturbances:

- a productivity shock;
- a labour-disutility shock.

Dynare uses the end-of-period convention for capital: production at date `t`
uses `k(-1)`, while current investment determines `k`.

### Preference parameters and fixed steady-state employment

`sigma` is the coefficient of relative risk aversion and the inverse of the
intertemporal elasticity of substitution. `frisch` (reported as varphi in the
figures and tables) is the inverse Frisch elasticity of labour supply.

A labour-disutility scale parameter, `chi`, is included in the intratemporal
condition:

```text
w*c^(-sigma) = chi*phi_t*n^frisch.
```

The key comparative-static design choice is that `chi` is recalibrated whenever
`sigma` or `frisch` changes so that deterministic steady-state employment remains
fixed at

```text
n_target = 1/3.
```

This prevents the sensitivity exercises from mixing elasticity effects with
changes in the steady-state employment level. The generated table
`steady_state_employment_control.csv` verifies this numerically.

## Parameter experiments

The core model compares:

```text
sigma   = {1, 2, 5, 10}
varphi  = {0.5, 1, 2, 5}
rho_z   = {0.5, 0.7, 0.9, 0.98}
rho_phi = {0.5, 0.7, 0.9, 0.98}
```

All impulse responses use a 40-quarter horizon. Quantitative sensitivity tables
report:

- impact response;
- signed peak response;
- time to peak;
- half-life after the peak;
- cumulative 20-quarter response;
- cumulative 40-quarter response.

## Fiscal extension

The fiscal version adds exogenous government purchases:

```text
c_t + i_t + g_t = y_t,
```

with government purchases following a persistent AR(1) process. A balanced
budget is represented explicitly by

```text
tax_t = g_t.
```

Taxes are lump sum and therefore non-distortionary. The aggregate resource
constraint summarizes the private and public budget constraints for the purpose
of equilibrium dynamics. Government purchases do not enter household utility.

The baseline government-purchases share is 20 percent of steady-state output.
The fiscal experiment varies

```text
rho_g = {0, 0.5, 0.9, 0.98}
```

and reports both impact and 20-quarter cumulative output multipliers.


## Bayesian estimation extension

A separate empirical extension estimates selected parameters of the RBC model
using quarterly data for the euro area (EA20) from 2000Q1 to 2025Q4. The two
observables are real GDP and total hours worked from Eurostat quarterly national
accounts. Both series are transformed into percentage log deviations from trend
using the quarterly Hodrick-Prescott filter with `lambda = 1600`.

The estimation keeps the preference and technology parameters calibrated and
estimates:

```text
rho_z
rho_phi
stderr eps_z
stderr eps_phi
```

Dynare combines the specified priors with the likelihood implied by the model
and the two observed series. The extension reports posterior parameter estimates,
smoothed structural innovations and posterior-distribution impulse responses. It
is implemented separately from the calibrated sensitivity analysis so that the
original exercises remain reproducible and directly comparable.

## Main outputs

Running the project produces the following figures:

1. `01_baseline_tfp_irfs.png`
2. `02_baseline_labour_irfs.png`
3. `03_sigma_sensitivity.png`
4. `04_labour_elasticity_sensitivity.png`
5. `05_tfp_persistence_sensitivity.png`
6. `06_labour_persistence_sensitivity.png`
7. `07_tfp_vs_labour_shock.png`
8. `08_sigma_varphi_heatmaps.png`
9. `09_baseline_variance_decomposition.png`
10. `10_stochastic_simulation.png`
11. `11_government_spending_irfs.png`
12. `12_fiscal_persistence_sensitivity.png`
13. `13_fiscal_variance_decomposition.png`
14. `14_bayesian_parameter_estimates.png`
15. `15_smoothed_structural_shocks.png`
16. `16_posterior_irfs_eps_z.png`
17. `17_posterior_irfs_eps_phi.png`

The main generated tables are:

- `baseline_irf_metrics.csv`
- `sigma_sensitivity_metrics.csv`
- `frisch_sensitivity_metrics.csv`
- `rho_z_sensitivity_metrics.csv`
- `rho_phi_sensitivity_metrics.csv`
- `steady_state_employment_control.csv`
- `baseline_hp_filtered_moments.csv`
- `baseline_autocorrelations.csv`
- `baseline_correlations.csv`
- `baseline_variance_decomposition.csv`
- `fiscal_rho_g_sensitivity_metrics.csv`
- `fiscal_multipliers.csv`
- `fiscal_variance_decomposition.csv`
- `research_question_map.csv`
- `bayesian_posterior_summary.csv`

The research-question map links each research question to its principal figure,
quantitative statistic and output table.

## Business-cycle moments

The baseline theoretical moments are computed with the quarterly Hodrick-Prescott
filter using `hp_filter=1600`. The corresponding file is therefore deliberately
named `baseline_hp_filtered_moments.csv`.

The output includes standard deviations and autocorrelations through lag 5 (or
all lags available from Dynare if fewer are returned). Contemporaneous
correlations are written separately.

## Variance decomposition

The variance-decomposition tables and stacked charts report unconditional model
variance shares associated with each structural innovation. These shares depend
on both the model transmission mechanism and the assumed innovation variances.
In the baseline calibration all structural innovations have standard deviation
0.01.

Accordingly, statements in a report should be phrased as, for example:

> Under the baseline calibration and assumed shock variances, productivity
> shocks account for X percent of model-generated output variance.

They should not be presented as estimates of the share of observed historical
output fluctuations unless the model is separately estimated or disciplined by
data.

## Reproducible simulation

The stochastic simulation uses 2,200 generated observations and discards the
first 200 as burn-in before plotting the sample. The seed is fixed at
`20260922`. Deterministic impulse responses and theoretical moments do not depend
on this seed.

## Files

- `rbc_core.mod`: core calibrated RBC model and sensitivity experiments.
- `rbc_fiscal.mod`: government-purchases extension.
- `run_rbc_project.m`: MATLAB entry point for the calibrated and fiscal analysis.
- `rbc_estimation.mod`: standalone Bayesian estimation extension.
- `run_rbc_estimation.m`: MATLAB entry point for the estimation and estimation figures.
- `euro_area_rbc_dynare.csv`: transformed EA20 observables used by Dynare.
- `figures/`: exported figures from the calibrated and estimated exercises.
- `tables/`: exported quantitative results.
- `results/`: generated intermediate output and not tracked in the repository.

To reproduce the calibrated project, run `run_rbc_project`. To reproduce the
Bayesian estimation extension, run `run_rbc_estimation`.

## Interpretation conventions

- IRFs are responses to one-standard-deviation innovations.
- Every baseline innovation has standard deviation 0.01.
- IRFs are multiplied by 100 and plotted as percentage deviations from steady
  state.
- `sigma` is the inverse intertemporal elasticity of substitution.
- `varphi` is the inverse Frisch elasticity; a larger value means a less elastic
  labour supply.
- The impact fiscal multiplier is the impact change in output divided by the
  impact change in government purchases, using the steady-state output-to-
  spending ratio to translate log responses into first-order level changes.
- The cumulative fiscal multiplier uses the first 20 quarters.
