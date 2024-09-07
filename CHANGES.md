## Change log

### Version 1.0.2

This release is that used to run the OSPAR 2024 CEMP assessment.

#### Weighted sums for TEQs etc.

Weighted sums of concentrations are now calculated using `determinand.link.sum`. The weights are supplied using the `weights` argument. This function supersedes `determinand.link.TEQDFP` which is now deprecated and will be removed at the next release.

Updated and corrected World Health Organisation TEFs for dioxins, furans and planar PCBs are now available in the data object `info_TEF`: there are four versions available:  

* DFP_2005; the 2005 values  
* DFP_2022; the 2022 values  
* DFP_HOLAS3; the values used in the HOLAS3 assessment; these are the 2005 values excluding three PCBs and are included for backward compatibility  
* DFP_CEMP; the values used in CEMP assessments <= 2024; these are the 2005 values excluding three PCBs and with the TEF for CDFO ten times too small; they are included for backward compatibility  

DFP_2005 and DFP_2022 use determinand OCDF (the correct code) rather than CDFO (which is a grouped determinand code). 

#### Auxiliary variables

Key auxiliary variables can now be plotted in `plot_assessment`. The default variables are:  

* biota: concentration, LNMEA, DRYWT%, LIPIDWT%  
* sediment: non-normalised concentration, normalised concentration, AL, CORG  
* water: no plots are currently produced  

The choice of auxiliary variables can be altered using the `auxiliary` argument, although the options here are still limited.

The merging of auxiliary variables with the data is now controlled using the `control` argument of `read_data`. `control$auxiliary$by_matrix` is a list which determines which auxiliary variables are merged by sample and matrix and which are just merged by sample. The default values are `c("DRYWT%", "LIPIDWT%)` for biota and `"all"` for sediment and water. Thus, by default, dry weight and lipid weight measurements are matched with chemical concentrations in the same tissue (matrix). but all other auxiliary variables in biota are matched at the sample level. For sediment (and water) all auxiliary variables (e.g. aluminium and organic carbon measurements) are matched with chemical concentrations in the same grain size fraction.

#### Reporting

`report_assessment` now has the full functionality required for the OSPAR CEMP assessment. This includes:  

* scatterplot matrices of concentrations of related compounds using the non-exported function `plot_multidata`  
* plots of assessments of related compounds using the non-exported function `plot_multiassessment`  
* plots of contaminant ratios using the non-exported function `plot_ratio`  

There is still some work required to make `report_assessment` suitable for all purposes.

#### Minor bug fixes  

* ensures that if `uncertainty` column is present in external data then so is `unit_uncertainty` and vice versa
* `plot_assessment` now correctly plots the 90% two-sided confidence intervals on VDSI estimates from imposex assessments 
* correct treatment of censoring data in `determinand.link.sum`
* `ctsm.check.sex.biota` now works for any auxiliary variable
* `get_timeseries` now always shows the series identifier for each timeseries
* `estimate_uncertainties` now traps for the case then `DRYWT%` and `LIPIDWT%` are not specified as auxiliary variables


### Version 1.0.1

Updates (mostly) required to run the OSPAR 2024 CEMP assessment. 

#### Data import

For OSPAR and HELCOM style assessments, data from Germany are now matched to stations by name for 2023 onwards. This applies to biota, sediment and water. Note that for HELCOM, biota data from Germany are already matched by name for all years.  

#### Uncertainty processing

harsat 1.0.0 replaced implausibly large relative uncertainties ($>=$ 100%) and replaced them with imputed values. However, implausibly small relative uncertainties were not changed. The code now replaces relative uncertainties $<=$ 1% with imputed values. 

The defaults can be changed using `control$relative_uncertainty` in `read_data`. To replicate the defaults in harsat 1.0.0, set `control$relative_uncertainty = c(0, 100)`. To keep all uncertainties, regardless of how ridiculous they are, set `control$relative_uncertainty = c(0, Inf)`.

Two minor bug fixes:  

* relative uncertainties were being filtered for all distributional types, but this is only a reliable procedure for determinands with `distribution == "lognormal"`; the checks are now only applied to lognormal data   
* some biological effect data with distributions other than normal or lognormal were being incorrectly deleted; this has now been corrected  

The oddity files have been updated to show:  

* implausible_uncertainties_reported.csv - all reported uncertainties that are replaced by imputed values  
* missing_uncertainties.csv - all uncertainties (normal or lognormal data) that are not reported and can't be imputed  
* implausible_uncertaintes_calculated.csv - all uncertainties that are calculated during the data processing (e.g. during normalisation) that are implausible and are set to missing  

#### Uncertainty coefficients

The function `ctsm_uncrt_workup` and related supporting functions are used in OSPAR assessments to update the fixed and proportional standard deviations which are subsequently used to impute missing uncertainties. These functions were ignored during the initial development of harsat and are now harsat compatible.  

#### Biological effect assessments

Imposex assessments: these are now fully reproducible with seeds for random number generation provided in the calls to `ctsm.VDS.cl` and `assess_imposex` 

Assessment functions for negative binomial data have been added. Negative binomial data includes MNC - the number of micronucleated cells.

#### Reporting

`report_assessment` generates default file names. These are based on the series identifier with additional station information. It is now possible to override this behaviour for a single report by providing a different file name using the `output_file` argument.

#### Reference tables  

* new values added to method_extraction table  

#### Minor bug fixes  

* correct behaviour of argument `return_early` in `create_timeseries`  
* pass `info` component of the harsat object to `determinand.link.sum`, `determinand.link.replace`, and `determinand.link.imposex`  
* ensure early return from `ctsm_convert_basis` when there is nothing to convert (avoids issues e.g. when all the data are biological effects)  
* ensure SURVT (in pargroup B-BIO) is recognised as a biological effect in `ctsm_get_datatype` (SURVT is the only determinand in this pargroup that isn't an auxiliary variable)  
* pass `good_status` to assessment functions for data with distributions other than normal and lognormal
* trap pathological case in estimation of `prtrend`; see #436
* ensure `ctsm_OHAT_legends` uses the symbology as specified in `write_summary_table` 


### Version 1.0.0

- Initial public release

### Version 0.1.3

- Various fixes

### Version 0.1.2

- Fixed issues when packaged; see #326, #328
- Updated AMAP data and packaging; see #329

### Version 0.1.1

- Fixed issue with auxiliary variables: see #289
- Small documentation improvements
- Added build processes for package bundles

### Version 0.1.0

- Initial release
