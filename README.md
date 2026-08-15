# Automated XRD Crystallinity Index (CI) and Crystallite Size Analyzer

![MATLAB](https://img.shields.io/badge/MATLAB-R2020a%2B-blue.svg)
![Method](https://img.shields.io/badge/Algorithm-ALS%20%2B%20Scherrer-orange.svg)

## 📌 Overview
This repository provides a MATLAB pipeline for automated X-Ray Diffraction (XRD) pattern processing. The algorithm separates the **amorphous halo** from **crystalline diffraction peaks** using **Asymmetric Weighted Least Squares (ALS)** baseline estimation, calculates the **Crystallinity Index ($CI\%$)** via numerical integration, and computes the **apparent crystallite size ($D$)** of individual diffraction peaks using the **Scherrer Equation**.

---

## ✨ Key Features
* **Asymmetric Weighted Least Squares (ALS) Baseline Extraction:** Flexible and objective extraction of the broad amorphous background without manual baseline drawing.
* **Light Gaussian Pre-Smoothing:** Reduces high-frequency instrument noise while preserving true peak shapes.
* **Crystallinity Index ($CI\%$) Quantification:** Calculates phase area fractions using trapezoidal numerical integration (`trapz`).
* **Automated Peak & FWHM Detection:** Identifies local intensity maxima and calculates Full Width at Half Maximum (FWHM) values directly from the crystalline profile.
* **Scherrer Crystallite Size Estimation:** Computes crystallite size ($D$, in nm) for each identified peak based on Cu-$K\alpha$ radiation ($\lambda = 0.15406\text{ nm}$).
* **Two-Panel Publication-Ready Plots:** Visualizes raw vs. baseline data and isolates crystalline peaks with custom text annotations.
* **Automated Excel Export:** Saves comprehensive peak parameters and summary metrics across multiple sheets in an Excel workbook (`.xlsx`).

---

## 🔬 Computational Methodology

### 1. Amorphous Baseline Modeling (ALS)
The background $z$ is estimated by minimizing the penalized asymmetric least-squares objective function:

$$S = \sum_{i} w_i (y_i - z_i)^2 + \lambda \sum_{i} (\Delta^2 z_i)^2$$

Where:
* $\lambda = 10^5$: Smoothness parameter adjusting baseline flexibility.
* $p = 0.002$: Asymmetry weighting factor ensuring the baseline hugs the lower envelope ($w_i = p$ for $y_i > z_i$, and $w_i = 1 - p$ otherwise).

### 2. Crystallinity Index ($CI$)
The relative degree of crystallinity is determined from integrated peak areas:

$$CI = \left( \frac{A_{\text{crystalline}}}{A_{\text{total}}} \right) \times 100\%$$

Where:
* $A_{\text{crystalline}} = \int I_{\text{crystalline}}(\theta) \, d\theta$
* $A_{\text{total}} = \int I_{\text{smooth}}(\theta) \, d\theta$

### 3. Crystallite Size ($D$) — Scherrer Equation
The mean crystallite dimension for each crystalline reflection is calculated as:

$$D = \frac{K \lambda}{\beta \cos \theta}$$

Where:
* $K = 0.9$ (crystallite shape factor).
* $\lambda = 0.15406\text{ nm}$ (X-ray wavelength for Cu-$K\alpha$).
* $\beta$: Line broadening at half the maximum intensity (FWHM in radians).
* $\theta$: Bragg angle (half of the $2\theta$ peak center, in radians).

---

## 🚀 Getting Started

### Prerequisites
* **MATLAB** (R2018a or later recommended).
* **Signal Processing Toolbox** (optional, recommended for enhanced filtering).

### Input Data Format
Prepare your raw XRD data in an Excel file (`.xlsx`) with two numerical columns without header mismatch:
* **Column 1:** Diffraction angle $2\theta$ (degrees).
* **Column 2:** Measured intensity $I$ (arbitrary units, a.u.).

### Usage
1. Clone this repository:
   ```bash
   git clone [https://github.com/your-username/xrd-als-crystallinity-analyzer.git](https://github.com/your-username/xrd-als-crystallinity-analyzer.git)
