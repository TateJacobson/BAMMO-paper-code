# Code to reproduce results from "Predicting 5G Throughput with BAMMO, a Boosted Additive Model for data with Missing Observations"

## Contents:
- **mboost**: Contains a modified version of the mboost R package (Hothorn et al., 2010), implementing BAMMO, our proposed method of handling missing values.
- **Lumos5G-v1.0**: Contains the Lumos 5G data from Narayanan et al (2020) analyzed in the paper.
- **cleaning_lumos_cv_output.R**: Generates plots from the output of lumos_cv.R, mcar_etc_sims.R, and missingness_scenarios.R
- **lumos_cleaning.R**: Cleans the Lumos loop data for analysis.
- **lumos_cv.R**: Compares the prediction performance of BAMMO with that of alternative methods on the Lumos loop data (original and with groups of variables set to be missing completely at random).
- **lumos_gam.R**: Fits and tunes BAMMO on the Lumos loop data. Checks variable selections and computes importance in terms of MSE reduction. 
- **lumos_soil.R**: Computes the proposed SOIL variable importance measure for the BAMMO model on the Lumos loop data.
- **mcar_etc_sims.R**: Compares BAMMO with other methods of handling missing data in simulated data with values missing completely at random (MCAR), missing at random (MAR), and missing not at random (MNAR).
- **missingness_scenarios.R**: Compares BAMMO with other methods of handling missing data on augmented versions of the Lumos loop data with values MCAR, MAR, and MNAR.
- **plotting_effects.R**: Plots estimated partial effects for variables in the BAMMO model fit to the Lumos loop data.
- **timing_comparison.R**: Compares the computation time of BAMMO with that of alternative methods for handling missing data.

## References:
- Narayanan, A., Ramadan, E., Mehta, R., Hu, X., Liu, Q., Fezeu, R. A., Dayalan, U. K., Verma, S., Ji,
P., Li, T., et al. (2020). Lumos5G: Mapping and Predicting Commercial mmWave 5G Throughput.
In Proceedings of the ACM Internet Measurement Conference, pages 176-193.
- Torsten Hothorn, Peter Bühlmann, Thomas Kneib, Matthias Schmid, and Benjamin Hofner. *mboost:
Model-Based Boosting*, 2010. URL http://CRAN.R-project.org/package=mboost. R package
version 2.9-8.

## License: 
This work is licensed under [GPL v2](https://cran.r-project.org/web/licenses/GPL-2).

The [*mboost*](https://cran.r-project.org/web/packages/mboost/index.html) R package is licensed under [GPL v2](https://cran.r-project.org/web/licenses/GPL-2) and has been modified for this work. 
The [*Lumos5G 1.0*](https://ieee-dataport.org/open-access/lumos5g-dataset) dataset is licensed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/) and is being shared unmodified.

