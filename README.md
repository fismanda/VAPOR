# VAPOR
Ventilation-Aware Pandemic and Outbreak Risk Model

[![DOI](https://zenodo.org/badge/1026304551.svg)](https://doi.org/10.5281/zenodo.16422619)

**Authors:** David N. Fisman, Natalie J. Wilson, Callandra Moore, Clara Eunyoung Lee, Ashleigh R. Tuite  
**Affiliations:** Dalla Lana School of Public Health, University of Toronto  
**Version:** v1.0 | **Release date:** July 25, 2025

---

## Overview

**VAPOR** is a hybrid infectious disease transmission model that integrates:

- **Reed-Frost** (close contact transmission)
- **Wells-Riley** (airborne transmission)

It allows deterministic and stochastic simulations, including:
- R₀ and k estimation
- Route attribution
- Sensitivity to aerosol fraction, ventilation (ACH), and room volume
- Multi-patch modeling of outbreak spread across heterogeneous indoor spaces

---

## Repository Contents

- `Annotated deterministic single patch.R`: Deterministic single-compartment simulations with route attribution
- `Annotated stochastic including comparison with deterministic.R`: Stochastic simulations and comparisons to deterministic R₀
- `Multi-patch stochastic VAPOR.R`: Multi-room (meta-population) stochastic outbreak simulations
- `VAPOR4.docx`: Full manuscript

---

## Running the Model

The model is written in **R 4.5.0** and requires:

```r
library(MASS)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
```

## How to Cite

Fisman DN et al. Ventilation-Aware Pandemic and Outbreak Risk (VAPOR) Model. University of Toronto, 2025. DOI: 10.5281/zenodo.16422619
