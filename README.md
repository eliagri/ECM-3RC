# 3RC Equivalent Circuit Modeling (ECM)
## Dataset
Source: https://calce.umd.edu/battery-data

* Capacity Rating: 1100 mAh
* Cell Chemistry: LiFePO4
* Diameter: 25.4 mm
* Length: 65 mm
* Special Notes	Tab length not included in dimensions

# Battery ECM and Heat-Generation Analysis
This repository contains MATLAB scripts developed for equivalent circuit modeling (ECM), validation, and heat-generation estimation of lithium iron phosphate (LFP) battery cells.

The project focuses on:
* OCV characterization,
* 3RC ECM parameter identification,
* ECM validation against dynamic current profiles,
* and scaling of experimentally derived ECM parameters to a larger-format LFP cell for heat-generation studies.

# Project Overview

The workflow implemented in this repository consists of four main steps:

* Extraction of OCV(SOC) relationship from low-current discharge data.
* 3RC ECM Parameter Identification
* Validation of the ECM using independent dynamic load profiles.
* Scaling of ECM parameters from an experimentally characterized A123 cell to a larger Hithium LFP cell.
* Estimation of irreversible heat generation and simplified thermal response.

# Main Results

Key findings from the project include:

* Successful extraction of smooth OCV(SOC) curves for LFP chemistry.
* ECM validation RMSE values around 50–55 mV.
* Predicted cumulative heat generation of approximately 73 Wh during a 4-hour discharge.
* Estimated temperature rise of approximately 5°C in the scaled large-format cell.

# Future Work
Potential future improvements include:
* Temperature-dependent ECM parameters
* SOC-dependent resistance scaling
* Coupled electro-thermal modeling
* Validation using large-format cell measurements
* Integration with Simscape Battery or BattMo
