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

## 📊 Output Visualization Examples

The algorithm generates a 2-panel visualization:
1.  **Top Panel:** Displays the raw data, ALS background, the cumulative Multi-Gaussian fit, and color-coded filled areas separating the amorphous ($X_a$) and crystalline ($X_c$) phases.
2.  **Bottom Panel:** Isolates the pure crystalline peaks and annotates each peak with its exact $2\theta$ position and estimated crystallite size ($D$ in nm) derived from the Scherrer equation.

*(You can add an image here later by uploading a screenshot of your MATLAB plot to your repository and linking it like this: `![Output Example](path_to_image.png)`)*

---

## 🔬 Scientific Background
The methodology implemented in this code addresses the subjectivity and reproducibility issues often found in manual XRD profile fitting. 

The algorithm classifies peaks based on their Full Width at Half Maximum (FWHM):
*   **Crystalline Phase:** $\text{FWHM} \le 2.5^\circ$
*   **Amorphous Phase:** $\text{FWHM} > 2.5^\circ$

Crystallinity ($X_c$) is calculated as:
$$X_c = \left( \frac{A_c}{A_c + A_a} \right) \times 100\%$$

Crystallite Size ($D$) is estimated using the Scherrer equation:
$$D = \frac{K \lambda}{\beta \cos \theta}$$
Where $K = 0.9$ (shape factor) and $\lambda = 0.15406$ nm (Cu-K$\alpha$ radiation).

---

## 📄 Citation
If you use this algorithm or code in your research, please cite our paper:

---

## 👨‍💻 Authors
*   **Junianto Sesa** (Universitas Papua)
*   **Bayu Harnadi Nasrul** (Universitas Hasanuddin)
*   **Andi Tessiwoja Tenri Ola** (Universitas Lambung Mangkurat)
*   **Pryandi M. Tabaika** (Universitas Sebelas September)
*   **Nurul Fajri Ramadhani Tang** (Universitas Papua)

## 📜 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
