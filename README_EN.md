# IROHA_Basic
A quick hyperspectral data analyzer based on rule-based approaches

The software "IROHA_Basic" is an open-source/free desktop application designed for researchers and engineers to intuitively load, visualize, and segment (classify) their hyperspectral image data.
It was developed as a fundamental "Iroha" (ABC/basics) tool to unlock the invisible wavelength (spectral) information of light.

---

## Key Features

* High-Speed Loading of ENVI Formats
  * Supports the ENVI format and loading via header files (`.hdr`).
* ROI (Region of Interest) Selection and Spectral Visualization
  * Extract and display the average spectrum of selected areas through intuitive mouse operations while viewing pseudo-RGB images.
* Flexible and Easy Analysis Switching Between 5 Algorithms
  * Easily perform segmentation on hyperspectral data based on the spectral features of the selected ROIs.

---

## System Requirements

* Supported OS: Windows 10 / 11 (64-bit)
* Memory: 8 GB or more recommended (depends on the size of the hyperspectral data being handled).
* Display Resolution: 1280x720 or higher (Because text size changes dynamically according to display resolution, high UI scaling settings may degrade visibility).

---

## Installation Instructions

1. Download the latest `IROHA_Basic_Setup.msi` from the Releases page ( https://github.com/litandus-yako/IROHA_Basic/releases ).
2. Double-click the downloaded `.msi` file to launch the installer. If a Microsoft Defender warning appears, please bypass it by clicking "More info" -> "Run anyway" (at your own risk).
3. Follow the on-screen instructions to complete the installation.
4. After installation, you can launch the application from the "IROHA_Basic" desktop shortcut or the `.exe` file in the installation folder.

---

## Basic Usage

1. Open Data
   * Click the [Load HS Data (.hdr)] button at the top left of the window and select your target `.hdr` file. Once loaded, a pseudo-RGB image will appear on the right side of the window.
2. Image Adjustment
   * If you need to adjust the visibility of the pseudo-RGB image, use the "Brightness Multiplier" slider.
3. Specify ROI (Analysis Region)
   * Drag over the image to select a rectangular area and check the spectral characteristics of each region. The average spectrum of the selected region is displayed in the graph window at the bottom left.
4. Execute Classification
   * Select the algorithm for classification from the dropdown menu. A brief description of the algorithm can be viewed by clicking the [Details] button. Click the [Run] button to execute classification using the set ROIs as training data.
5. Save Data
   * The classification result will be displayed in a new window. Check the result and click the [Save Screen] button to save it as an image.

---

## For Those Needing Advanced Features (About IROHA_Pro)

`IROHA_Basic` is specialized for the interactive exploration and visualization of single images.
For research or business purposes requiring batch processing or quantitative analysis, please consider `IROHA_Pro` (paid version), which is scheduled for future release.

[Key features planned for IROHA_Pro (development roadmap is subject to change without notice)]
* Save/Load classification models and export feature quantities as CSV
* Bulk batch processing for multiple hyperspectral datasets
* Quantitative evaluation using PLSR (Partial Least Squares Regression)
* Vector format image export

---

## Third-Party Licenses

This application is built using open-source software such as Python, NumPy, scikit-learn, Spectral, and Tkinter. For full copyright and license details of each library, please refer to the `THIRD_PARTY_LICENSES.md` file in your installation directory.

---

## Inquiries / Bug Reports

For bug reports and feature requests, please submit them via GitHub [Issues] or the LITANDUS Website contact form ( https://sites.google.com/litandus.com/litandus/contact ).