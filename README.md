This repository contains the complementary materials associated with the article:

**"Spatio-temporal modeling of chikungunya in Brazil: The Interplay of Climate, Mobility, and Socioeconomic Factors"**  
**Authors:** Raquel Gardini Sanches Palasio, Patricia Marques Moralejo Bermudi, Gerson Laurindo Barbosa, Pedro Silva Peixoto, Sergio Muniz Oliva, Camila Lorenz, Marta Blangiardo, Monica Pirani, and Francisco Chiaravalloti-Neto.

Includes SPDE-based INLA modeling, visualization scripts, and mapping outputs.
---

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

---

## Data Availability

The datasets and variable dictionary are publicly available at:  
[**Figshare (link to be inserted)**]  

The dataset includes:
- Chikungunya incidence per municipality and time point (seasonal aggregation)  
- Covariates: climate (temperature, precipitation), mobility, socioeconomic indicators (Katz index, piped water, sewerage), environmental (forest cover, altitude, biome), ENSO classification  
- Spatial units: municipalities (`COD_MUN6`) and Brazilian regions  
- Temporal coverage: 2015–2024 (seasonal)

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
```

## Citation
If you use or adapt this code, please cite the original article:
Palasio, R.G.S., Bermudi, P.M.M., Barbosa, G.L., Peixoto, P.S., Oliva, S.M., Lorenz, C., Blangiardo, M., Pirani, M., & Chiaravalloti-Neto, F.
Spatio-temporal modeling of chikungunya in Brazil: The Interplay of Climate, Mobility, and Socioeconomic Factors.
[Journal name, Year, DOI – to be added upon publication]

Contact: Raquel Gardini Sanches Palasio
Thank you for your interest in this project!
