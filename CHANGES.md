## Change log

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
