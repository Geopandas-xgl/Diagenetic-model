# Diagenetic-model

## Numerical Model of Fluid–Rock Interaction in Carbonate Rocks  

**(Xiong et al., 2025, _Science Advances_)**

### Version: GLX  

**Last updated: 2025-09-03**

---

## 🔍 Introduction

This MATLAB-based model (developed in R2024a) simulates early phreatic meteoric diagenesis in carbonate sediments. It is adapted from Ahm et al. (2018, _GCA_) and Murphy et al. (2022, _GCA_), with significant enhancements for modern environments.

### 🔧 Key Features

1. **Simplified and Accessible**: Focuses on the most common isotopic signals—δ¹³C and δ¹⁸O—making the model easier to understand and apply.
2. **Modern Meteoric Conditions**: Incorporates depth-dependent linear decreases in both recrystallization rate (R) and vertical fluid flow velocity (u), more accurately reflecting modern phreatic meteoric systems.
3. **Isotopic Discrimination**: Distinguishes the influence of meteoric diagenesis from primary carbon isotope excursions (CIEs).

---

## ⚙️ Model Components

| Script                    | Description                                                  |
| ------------------------- | ------------------------------------------------------------ |
| `DiagModel_GLX_initial.m` | Main execution script: runs the model under Pre-CIE and Peak-CIE conditions, saves output, and calls the plot function. |
| `DiagModel_GLX_base.m`    | Defines initial conditions, model parameters, and constants. |
| `DiagModel_GLX_dJdt.m`    | ODE solver that tracks changes in mass and isotopic composition over time. |
| `DiagModel_GLX_fluxes.m`  | Computes the input/output fluxes of each model box at each time step. |
| `DiagModel_GLX_plot.m`    | Generates publication-quality visualizations of model results. |
| `Baseline condition.cocx` | Summary of baseline condition for primary LMC to second LMC recrystallization in phreatic meteoric zone. |
| `Data`                    | The C and O isotope from Tingri section. |

---

## 📧 Acknowledgments

Feedback and contributions are welcome!  
**Contact**: [xiongguolin@smail.nju.edu.cn](mailto:xiongguolin@smail.nju.edu.cn)

### 🔗 Project Repository

For code updates and related scripts, visit:  
[GitHub – Geopandas-xgl Repository](https://github.com/Geopandas-xgl/Diagenetic-model)
