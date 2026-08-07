This repository contains the Supplementary materials associated with the article:

**"Spatiotemporal modeling of chikungunya in Brazil through the interplay of climate, mobility, and socioeconomic factors"**  
**Authors:** Raquel Gardini Sanches Palasio, Patricia Marques Moralejo Bermudi, Gerson Laurindo Barbosa, Pedro Silva Peixoto, Sergio Muniz Oliva, Camila Lorenz, Marta Blangiardo, Monica Pirani, and Francisco Chiaravalloti-Neto.

Includes SPDE-based INLA modeling, visualization scripts, and mapping outputs.
---

This repository accompanies the manuscript by providing the R scripts and reduced model objects needed to reproduce key results.

## Repository Structure

### 1. Folders

- **`data/`**  
  Contains the datasets and auxiliary files used in the study:  
  - `*.RData`: model outputs.  
  - `borda_102033.*`: Shapefile defining Brazilian administrative borders for spatial modeling.

---
### 2. Scripts

- **`Chik_SPDE_nbi_model.R`**  
  Implements the full spatio-temporal model using R. 
  The model employs the **Stochastic Partial Differential Equations (SPDE)** approach with a **Negative Binomial likelihood** to capture the spatial and temporal dynamics of chikungunya incidence in Brazil,  within the Diffusion-based Extension of the Matérn Field (DEMF) (1,2,1) family, as proposed by Lindgren et al. 2024. 
 

  **Key features of the script:**
  1. Load and preprocess chikungunya incidence data and covariates (climate, mobility, socioeconomic, and environmental factors).  
  2. Load and transform shapefiles for municipalities and regions of Brazil.  
  3. Construct a spatial mesh using the SPDE approach.  
  4. Specify a spatio-temporal model with:
  5. Fit the model using the **INLA** framework and `inlabru` interface.  
  6. Extract posterior summaries (fixed effects, hyperparameters, random effects, fitted values).  
  7. Generate figures:
     - Plots of incidence rate ratios (IRR) with 95% credible intervals  
     - Non-linear climate effect plots  
     - Maps of predicted chikungunya incidence, including posterior mean and 95% CI
  8. Perform posterior variance partitioning analysis to quantify the relative contribution of each model component to the total variance of the linear predictor.
     - Posterior variance partitioning plot with 95% credible intervals

Lindgren F, Bakka H, Bolin D, Krainski E, Rue H. A diffusion-based spatiotemporal extension of Gaussian Matéern fields. SORT-Statistics and Operations Research Transactions. 2024;48(1):3–66. https://doi.org/10.57645/20.8080.02.13.

---

## Data Availability

The datasets and variable dictionary are publicly available at:  
[**Figshare (https://doi.org/10.6084/m9.figshare.30410647)**]  

The Figshare includes:
- Chikungunya incidence per municipality and time point (seasonal aggregation)  
- Covariates: climate (temperature, precipitation), mobility, socioeconomic indicators (Katz index, piped water, sewerage), environmental (forest cover, altitude, biome), ENSO classification  
- Spatial units: municipalities (`COD_MUN6`) and Brazilian regions  
- Temporal coverage: 2015–2024 (seasonal)

The complete model objects, posterior samples, and intermediate results required for full computational reproducibility are archived in Zenodo:
[**Zenodo archive: https://doi.org/10.5281/zenodo.17604821**]  

The Zenodo archive includes:

Posterior sample batches:
- amostras_blocos1_chik_original.RData
- amostras_blocos2_chik_original.RData
- amostras_blocos3_chik_original.RData

These correspond to three posterior sampling batches (200 samples each; 600 samples total) used for posterior variance partitioning analysis.

---

## Requirements

This project was developed in **R (version 4.3.3 — "Angel Food Cake")** with the following key packages:

```R
library(INLA)          # Bayesian inference via Integrated Nested Laplace Approximation
library(INLAspacetime) # Spatio-temporal extensions for INLA
library(inlabru)       # Bayesian spatial modeling using INLA
library(INLAOutputs)   # Utilities for INLA posterior summaries
library(geobr)         # Brazilian shapefiles (municipality and region)
library(ggspatial)     # Scale bars and north symbols in ggplot2
library(ggsn)          # North symbols and scale bars
library(fmesher)       # Mesh construction for SPDE

Package versions used:

INLA: 24.3.29
INLAspacetime: 0.1.8
inlabru: 2.10.1
fmesher: 0.1.7

Additional R packages
library(INLAOutputs)
library(geobr)
library(ggplot2)
library(ggspatial)
library(ggsn)
library(dplyr)
library(sf)
```
## Operating Systems 
Windows 10 / 11 (RStudio)

## Hardware Requirements
The full model (not included in this repository due to size constraints) was originally fitted on a workstation with:
32–64 GB RAM
Multicore CPU
However, the reduced objects provided in this repository (fit121_reduced.RData and fit121_reduced_map.RData) require less than 2 GB RAM, and load normally on standard systems.

### Installation Guide
```R
remotes::install_version("INLA", version="24.03.29",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
remotes::install_version("fmesher", "0.1.7")
install.packages("inlabru")
install.packages(c("geobr", "sf", "ggplot2", "dplyr", "ggspatial", "ggsn"))
```
## Typical installation time

On a standard laptop (8–16 GB RAM):
R packages: 5–10 minutes
INLA installation: 2–5 minutes

## Demo

A demonstration script is provided in:
Chik_SPDE_nbi_model.R

How to run the demo
### Load model output (reduced size)
```R
load("fit121_reduced_income.RData")       # contains fixed effects, hyperparameters, WAIC, DIC, CPO
load("fit121_reduced_income_map.RData")   # contains subset of fitted values for mapping
```
### Load input data
```R
data <- read_excel("data/Supplementary data.xlsx", sheet = "data")
```
### Expected Output

Posterior summaries of:
 - fixed effects
 - hyperparameters
 - selected random effects
 - Model fit statistics:
 - WAIC, DIC, CPO

Figures:
 - Relative risk plots
 - Fitted trend curves
 - Maps of predicted incidence and 95% CI

### Run Time
- On a normal computer (8–16GB RAM):
- Loading data: < 10 seconds
- Running demo plots: 10–30 seconds
  
(Note: The full model fitting is not included, as it takes ~24 hours and requires high-memory systems.)

### Instructions for Use
Users may apply the model code to their own datasets following these steps:

Step 1 — Prepare Data
 - Your data must include:
 - coordinates (xUTM, yUTM)
 - time index
 - population offset
 - covariates (Prec, Tmin.grp, logKatz, Forest, etc.)

Step 2 — Build Mesh (SPDE)
```R
mesh <- inla.mesh.2d(loc = cbind(coords.cv[, 1], coords.cv[, 2]),
         boundary = st_as_sfc(border))
```
Step 3 — Define Model Components
```R

M <- ~ -1 + Intercept(1) +
  myoffset(log(pop10), model = "offset") +
  sew(sewerage.fc, model ="factor_contrast") +  # Socioeconomic: sewerage
  pip(piped.fc, model ="factor_contrast") +     # Socioeconomic: piped water
  Tmin.grp(Tmin.grp, model= 'rw2') +            # Climate: min temperature (nonlinear)
  Prec +                                        # Climate: precipitation
  logKatz +                                     # Mobility: Katz index (Log)
  Forest +                                      # Environmental: forest cover
  Alt(ALT.fc, model = "factor_contrast") +      # Elevation factor
  Biome(BIO.fc, model ="factor_contrast") +     # Biome categorical factor
  field(list(
    space = cbind(xUTM, yUTM),               # Spatial coordinates
    time = ID.time.1),                       # Temporal index
    model = stmodel)                         # Spatiotemporal model defined earlier
```
Step 4 — Fit Model (optional, advanced)
Only recommended for users with high-performance hardware.

```R
lhood <- like(
  formula = y ~ .,          # Response variable is count of cases (y)
  family = "nbinomial",     # Negative binomial model for overdispersed counts
  data = dataf
)

fit <- bru(M,
           lhood,
            options = list(control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)))
```
## Reproduction Instructions

The repository includes:

   - Reduced model object (`fit121_reduced_income.RData`)
   - Subset of fitted values for mapping (`fit121_reduced_income_map.RData`)
   - Alternative reduced model object (fit121_reduced_sew_piped.RData)
   - Subset of fitted values for mapping from the alternative model (`fit121_reduced_map_sew_piped.RData`)
   - Full plotting and mapping scripts
   - Posterior sample batches (amostras_blocos) used for posterior variance partitioning.
     
   - These allow reproduction of all:
   - fixed-effect tables
   - hyperparameter summaries
   - relative risk curves
   - spatial maps (posterior mean, 95% CI)
   - Posterior variance partitioning summaries and figures.

without refitting the full model.

## Citation
If you use or adapt this code, please cite the original article:
Palasio, R.G.S., Bermudi, P.M.M., Barbosa, G.L., Peixoto, P.S., Oliva, S.M., Lorenz, C., Blangiardo, M., Pirani, M., & Chiaravalloti-Neto, F.
Spatiotemporal modeling of chikungunya in Brazil through the interplay of climate, mobility, and socioeconomic factors.
[Journal name, 2026, DOI – to be added upon publication]

Palasio, R., Patricia, B., Pirani, M., & Chiaravalloti-Neto, F. (2025). Code for article ”Spatio-temporal Modeling of chikungunya in Brazil though interplay of climate
with mobility and Socioeconomic Factors”. Zenodo. https://doi.org/10.5281/zenodo.17604821

Contact: Raquel Gardini Sanches Palasio
Thank you for your interest in this project!
