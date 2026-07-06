# NA

## Change log

### Version 1.0.6

#### Major bug fix - p_overall_trend

For smooth models, the degrees of freedom used in the likelihood ratio
test of an overall temporal trend are now correct. Previously, they were
too small (by one) resulting in values of `p_overall_trend` that were
too significant. Fortunately, this had no effect on model selection as
the chosen model is based on AICc or AIC (depending on the distribution
of the response) and the p-values are calculated after the model is
chosen.

The statistical interpretation of smooth models in `report_assessment`
is now more nuanced. Smooth models chosen by AIC or AICc are not
necessarily significant at the conventional 5% level. The significance
of the final model is given by `p_overall_trend` which is based on a
likelihood ratio test that compares the smooth model and the mean model
(and essentially tests for any evidence of a change in concentrations
over time). This degree of significance is now characterised as weak (p
\>= 0.05), moderate (0.05 \< p \<= 0.01) and strong (p \< 0.01). In
theory, a smooth model chosen by AICc could have a `p_overall_trend` as
high as 0.135, but this can only happen when there are many years of
data and the improvement in AICc between the smooth model and the mean
model is marginal.

Note that, for smooth models, the significance of the overall trend can
be split into the significance of the nonlinear component
(`p_nonlinear_trend` based on a likelihood ratio test that compares the
smooth model with the linear model) and the linear component
(`p_linear_trend` based on a likelihood ratio test that compares the
linear model with the mean model). Both these p-values were calculated
correctly in previous releases.

### Version 1.0.5

This release is used to run the OSPAR 2026 assessment.

#### Graphics - year range

The year range on the x-axis in `plot_assessment` can now be specified
by the user. This could be used to make the year range consistent across
a series of related plots. The year-range is adjusted using the `xlim`
element of the new `control` argument. For example by specifying
`control = list(xlim = c(1999, 2026))`. It will affect all the images
produced by `plot_assessment` in that call.

#### Summary tables - timestamp

A timestamp can now be added to the summary table produced by
`write_summary_table` if argument `timestamp` is set to TRUE. The
default is not to produce a timestamp.

#### report_assessment

- the degrees of freedom used to add smoothers to the ratio plots have
  been reduced to avoid overfitting and potential crashing  
- a warning is now printed when a time series assessment has not
  converged and no statistical interpretation is provided  
- up to 49 related assessments can now be plotted in a 7x7 array
  (increased from 25); more than 49 and a message says there are too
  many related compounds to produce a meaningful plot  
- the scatter plot data matrix is now restricted to at most 20 related
  compounds; more than 20 and a message says there are too many related
  compounds to produce a meaningful plot

#### Minor bug fixes

- `create_timeseries` now recognises all ICES parameter groups which
  relate to chemicals in the checking function
  `ctsm.check.species_group.biota`  
- `create_timeseries` now recognises `"BL"` (blood) as a valid matrix
  for mammals in the checking function `ctsm.check.matrix.biota`  
- `read_data`: if geometry errors in the ICES station dictionary can not
  be fixed, the relevant stations are now deleted (currently only four
  stations affected)

### Version 1.0.4

This is an intermediate release to coincide with the 2025 harsat user
group meeting.

#### Summary tables

There have been several enhancements to `write_summary_table`. The most
important of these involves the construction of symbologies to
characterise summary features of each time series, such as its trend and
status using a shape and a colour.

The default symbology, used by OSPAR and HELCOM, has been made more
flexible through the argument `symbology_control`. As well as specifying
the colours for each threshold, the user can specify:

- the shape used for each type of trend  
- the p-value that a trend must attain to be treated as significant  
- whether the trends (shapes) are based on the recent change or the
  overall change

The formulation for specifying colours has changed slightly. It is now a
list of thresholds each with a ‘below’ colour and an ‘above’ colour
(rather than a list with ‘below’ colours for each threshold and a list
of ‘above’ colours for each threshold). The colour used for timeseries
when there are no thresholds is now black by default, but this can be
changed using `symbology_control`

In addition:

- custom symbologies can be applied with a user-supplied function  
- multiple symbologies can be applied, for example to assess against
  both environmental and human health standards  
- mutliple symbologies can be a mixture of the default symbologies and
  custom symbologies  
- names of the symbology columns can be changed / specified

There are lots of examples in the function documentation.

The other enhancements involve normalisation (see a later topic) and the
grouping of thresholds (e.g. where several thresholds such as the EAC
and EQS represent the same boundary between good and poor status). The
latter now has more checks and is implemented before any symbology is
applied. This means that any symbology should be defined using the
threshold groups, rather than their constituent thresholds. The argument
`collapse_AC` has been deprecated and replaced with the argument
`threshold_groups` which is more intuitive and consistent with the use
of threshold in most places through the code.

#### Units

Previous versions of `plot_assessment` and `report_assessment` only
worked with a pre-defined set of units, some of which were formatted to
look nice (e.g. with greek characters or superscripts). This restriction
has now been relaxed so that any unit supplied in the determinand
reference table can be used. A pre-defined set of units will still
undergo special formatting with all other units reported without special
formatting.

It is now easier (coding-wise) to add units to the pre-defined set.
However, doing so will require an issue to be raised and a subsequent
update to the harsat code.

#### Normalisation

Contaminant time series can be normalised in `create_timeseries`, for
example to 5% aluminium or 5% lipid. The way in which each time series
is normalised is now passed forwards in the timeseries object so that it
can be accessed in `plot_assessment`, `report_assessment`, and
`write_summary_table`. Specifically, if there is any normalisation, then
info\$normalise is set to `TRUE` and three extra variables are created
in the timeseries object:

- `normaliser` - the name of the normaliser for each timeseries;
  e.g. “AL”  
- `normaliser_value` - the value that concentrations are normalised to;
  e.g. 5  
- `normalise_unit` - the unit that `normaliser_value` is expressed as;
  e.g. “%”

These variables are also output in `write_summary_table`.

#### HAT

Experimental code for sending output to the XHAT has been uploaded. This
still needs a lot of work so please hold your breath until a later
release!

#### Minor bug fixes

- French pivot values used in sediment normalisation for OSPAR are now
  picked up correctly
- `read_data` now repairs invalid station geometries in ICES station
  dictionary extractions

### Version 1.0.3

This release is that used to run the OSPAR 2025 CEMP assessment.

#### ICES data extractions

`read_data` and `read_contaminants` now expect to find the variables
`amap_arctic_lme` and `casenumber` in the contaminants data file from an
ICES extraction (when argument `data_format` is set to `"ICES"`). These
variables were introduced in ICES extractions in mid-June 2024. A
warning is printed if the variables are missing. `read_stations` also
expects to fins `amap_arctic_lme` in the ICES station dictionary.

- `amap_arctic_lme` contains AMAP large marine ecosystem regional
  information
- `casenumber` is the accession identification id in the AIMS (Accession
  Information Management System) system introduced in 2024; it replaces
  accessionid from the now discontinued DAD system

#### AMAP assessments

There are two changes that affect data processing for AMAP assessments
(when argument `purpose` in `read_data` is set to `"AMAP"`).

First, the default behaviour by which `read_data` matches stations to
data from an ICES extraction is now to restrict eligible stations in
those in the AMAP area. Previously it didn’t matter where the stations
were. The default behaviour can be changed using `control$add_stations`
in `read_data`. See the help file of `add_stations` for more information
about station matching.

Second, the default values `control$region` are now

- `id = "amap_arctic_lme"`
- `names = "AMAP_arctic_lme"`
- `all = FALSE` (because of discrepancies between shape files for the
  AMAP area and the AMAP large marine ecosystems)

Users of external data might need to rename the regional variable in the
station dictionary to `amap_arctic_lme` or use `control$region` to pick
up the existing variable name.

#### Summations involving censored values

An argument `sum_censored` has been added to `determinand.link.sum`
which gives control over how censored values are treated when computing
a summed concentration.

`sum_censored = TRUE` (default) maintains the previous behaviour where
the sum is computed by adding all the non-censored measurements and the
censored values of the censored measurements. If all the measurements
are censored and all the censored values are equal to the limit of
detection, then the sum will be the sum of the limits of detection.

`sum_censored = FALSE` typically treats each censored value as zero when
calculating the sum. Usually, the output is the sum of the non-censored
measurements. There are two exceptions:

- when all measurements are censored, the output is the largest censored
  value (and is flagged as a censored measurement)
- when the sum of the non-censored measurements is less than the largest
  censored value, the output is the largest censored value (and is
  flagged as a censored measurement); this will be unusual and is most
  likely to occur when using weights (and small weights are applied to
  non-censored measurements and large weights are applied to the
  censored measurements)

The same argument has also been added to customised link functions such
as `determinand.link.BBKF` which call `determinand.link.sum`.

It is likely that most users will use `sum_censored = FALSE` since this
is compatible with Marine Strategy advice on the treatment of censored
values. If so, this might become the default option in a future release.

#### Summary tables - variable renaming

Several variables in the output from `write_summary_table` have been
renamed:

- `p_nonlinear` –\> `p_nonlinear_trend`
- `p_linear` –\> `p_linear_trend`
- `p_overall` –\> `p_overall_trend`
- `p_linear_trend` –\> `p_overall_change`
- `p_recent_trend` –\> `p_recent change`
- `linear_trend` –\> `overall_change`
- `recent_trend` –\> `recent_change`

The first three variables report evidence of systematic trends in the
data. The remaining variables all relate to changes in concentration
between the start and end of either the whole time series or the
‘recent’ part of the time series (typically the last 20 monitoring
years). The renaming was prompted by confusion about the original names
`p_linear_trend` and `linear_trend` which suggest linearity that is not
always the case. The harsat team struggled to come up with suitable
alternatives and hopes that the new names are at least less confusing,
if not crystal clear.

#### Calculation of recent_change

`overall_change` is the change in concentration over the whole
monitoring period. It is only calculated if there are at least five
years of data each with at least one non-censored measurement.

`recent_change` is the change in concentration in recent years,
typically taken to be the last twenty monitoring years. In previous
releases, `recent_change` was calculated if `overall_change` had been
calculated and there were at least five years of data in the recent
period. However there was no requirement on the number of years with
non-censored measurements in the recent period. This meant that the
evidence base for `recent_change` could be weak (e.g. if only one or two
years in the recent period had non-censored measurements). Occasionally
`recent_change` could be undefined for long time series with infrequent
monitoring and many censored measurements.

The calculation of `recent_change` can now be controlled using the
`control` argument of `run_assessment`. Specifically,
`control$recent_change` has two components:

- `n_year_fit` - default 5L
- `n_year_positive` - default 5L

where `n_year_fit` is the required number of years of data in the recent
period and `n_year_positive` is the required number of years of data
with at least one non-censored measurement in the recent period.

By default, `recent_change` will now only be calculated if there are at
least 5 years of data with non-censored measurements in the recent
period.

The behaviour of previous releases can be replicated (almost) by seting
`n_year_fit` to 5L and `n_year_positive` to 2L (the smallest value
allowed to avoid pathological behaviour).

#### Imposex assessments - annual indices equal to zero

Imposex assessments (based on individual measurements) involve the
estimation of cut-points that measure the transition from one imposex
stage to the next on the latent odds scale. This estimation pools the
data from several (often many) timeseries to improve the precision of
the cut-point estimates and is done before the call to `run_assessment`.

Data from timeseries / year combinations where all the individual
measurements are zero (equivalently, the annual index is zero) contain
no information on cut-points and are now removed from the estimation
procedure. This improves convergence. The estimation routine has been
renamed as `ordinal_theta_est` (previously `ctsm.VDS.index.opt`)

Estimation of the cut-points is followed by estimation of confidence
intervals for the annual indices. This is possible because the model
output also includes estimates of each index. The confidence intervals
are estimated from the posterior distribution of the parameter
estimates. However, this does not work for zero indices whose fitted
values would be infinite on the latent scale. A good solution would be
to calculate likelihood intervals, but this is a challenging numerical
problem and has been left for a future release. Instead, the previous
very ad-hoc approach has been replaced by a moderately ad-hoc approach.
Specifically, an infinite value for a zero index is replaced by a value
larger than all the fitted values for positive indices. This is done by
taking all the positive indices, fitting a linear model of fitted value
against square root index, and predicting what the fitted value should
be when the index is zero. The square root scale is based on the typical
relationship between fitted values and indices observed in the OSPAR
2025 assessment. The associated standard error is taken to be the upper
90th quantile of the standard errors associated with positive indices
(with a suitable adjustment for the number of measurements used in
each), which is reasonable since standard errors increase as the number
of zero measurements increase but not in a very predictable manner (at
least not in the OSPAR 2025 assessment). The upper confidence limit on a
zero index is then estimated by simulating from the posterior
distribution of the parameter estimates as before. (The lower confidence
limit is zero by definition.) The estimation routine has been renamed
`ordinal_theta_cl` (previously `ctsm.VDS.index.cl`).

#### Error handling

Error handling for the `sample` variable in the contaminants data file
(input to `read_data`) has been tightened up. This is only likely to
affect external data.

An error is now thrown if there are any missing values.

There are now checks for non-unique sample identifiers, for example when
the same sample identifier has been used in different years, or for
different stations or species. If these are found, a warning is printed
exhorting the user to sort out their data. However, the code also
attempts to create unique sample identifiers by pasting together `year`,
`station_code`, `species` (biota only) and `sample`.

#### Deleted functions

- `determinand.link.TEQDFP` was deprecated in release 1.0.2 and is now
  deleted; use `determinand.link.sum` instead.

#### Minor bug fixes

- `write_summary_table` now works for biological effects assessments
  with non-standard summary variables (e.g. imposex_class)
- `ctsm_symbology_OSPAR` now works when there are no non-parametric
  tests of status (e.g. when only imposex is assessed)
- `determinand.link.sum`: the uncertainty of the summed concentration is
  now independent of the order of the data; testing suggests that data
  order previously affected the uncertainty of a tiny proportion of
  samples, typically where all the measurements were censored
- `plot_assessment`: no extra page is produced in pdf plots of
  assessments with thresholds
- `write_summary_tabe`: a header is always produced when a new summary
  table is written to file; previously failed if the `append` argument
  was set to `TRUE` and there was no existing summary file to append to

### Version 1.0.2

This release is that used to run the OSPAR 2024 CEMP assessment.

#### Weighted sums for TEQs etc.

Weighted sums of concentrations are now calculated using
`determinand.link.sum`. The weights are supplied using the `weights`
argument. This function supersedes `determinand.link.TEQDFP` which is
now deprecated and will be removed at the next release.

Updated and corrected World Health Organisation TEFs for dioxins, furans
and planar PCBs are now available in the data object `info_TEF`: there
are four versions available:

- DFP_2005; the 2005 values  
- DFP_2022; the 2022 values  
- DFP_HOLAS3; the values used in the HOLAS3 assessment; these are the
  2005 values excluding three PCBs and are included for backward
  compatibility  
- DFP_CEMP; the values used in CEMP assessments \<= 2024; these are the
  2005 values excluding three PCBs and with the TEF for CDFO ten times
  too small; they are included for backward compatibility

DFP_2005 and DFP_2022 use determinand OCDF (the correct code) rather
than CDFO (which is a grouped determinand code).

#### Auxiliary variables

Key auxiliary variables can now be plotted in `plot_assessment`. The
default variables are:

- biota: concentration, LNMEA, DRYWT%, LIPIDWT%  
- sediment: non-normalised concentration, normalised concentration, AL,
  CORG  
- water: no plots are currently produced

The choice of auxiliary variables can be altered using the `auxiliary`
argument, although the options here are still limited.

The merging of auxiliary variables with the data is now controlled using
the `control` argument of `read_data`. `control$auxiliary$by_matrix` is
a list which determines which auxiliary variables are merged by sample
and matrix and which are just merged by sample. The default values are
`c("DRYWT%", "LIPIDWT%)` for biota and `"all"` for sediment and water.
Thus, by default, dry weight and lipid weight measurements are matched
with chemical concentrations in the same tissue (matrix). but all other
auxiliary variables in biota are matched at the sample level. For
sediment (and water) all auxiliary variables (e.g. aluminium and organic
carbon measurements) are matched with chemical concentrations in the
same grain size fraction.

#### Reporting

`report_assessment` now has the full functionality required for the
OSPAR CEMP assessment. This includes:

- scatterplot matrices of concentrations of related compounds using the
  non-exported function `plot_multidata`  
- plots of assessments of related compounds using the non-exported
  function `plot_multiassessment`  
- plots of contaminant ratios using the non-exported function
  `plot_ratio`

There is still some work required to make `report_assessment` suitable
for all purposes.

#### Minor bug fixes

- ensures that if `uncertainty` column is present in external data then
  so is `unit_uncertainty` and vice versa
- `plot_assessment` now correctly plots the 90% two-sided confidence
  intervals on VDSI estimates from imposex assessments
- correct treatment of censoring data in `determinand.link.sum`
- `ctsm.check.sex.biota` now works for any auxiliary variable
- `get_timeseries` now always shows the series identifier for each
  timeseries
- `estimate_uncertainties` now traps for the case then `DRYWT%` and
  `LIPIDWT%` are not specified as auxiliary variables

### Version 1.0.1

Updates (mostly) required to run the OSPAR 2024 CEMP assessment.

#### Data import

For OSPAR and HELCOM style assessments, data from Germany are now
matched to stations by name for 2023 onwards. This applies to biota,
sediment and water. Note that for HELCOM, biota data from Germany are
already matched by name for all years.

#### Uncertainty processing

harsat 1.0.0 replaced implausibly large relative uncertainties ($`>=`$
100%) and replaced them with imputed values. However, implausibly small
relative uncertainties were not changed. The code now replaces relative
uncertainties $`<=`$ 1% with imputed values.

The defaults can be changed using `control$relative_uncertainty` in
`read_data`. To replicate the defaults in harsat 1.0.0, set
`control$relative_uncertainty = c(0, 100)`. To keep all uncertainties,
regardless of how ridiculous they are, set
`control$relative_uncertainty = c(0, Inf)`.

Two minor bug fixes:

- relative uncertainties were being filtered for all distributional
  types, but this is only a reliable procedure for determinands with
  `distribution == "lognormal"`; the checks are now only applied to
  lognormal data  
- some biological effect data with distributions other than normal or
  lognormal were being incorrectly deleted; this has now been corrected

The oddity files have been updated to show:

- implausible_uncertainties_reported.csv - all reported uncertainties
  that are replaced by imputed values  
- missing_uncertainties.csv - all uncertainties (normal or lognormal
  data) that are not reported and can’t be imputed  
- implausible_uncertaintes_calculated.csv - all uncertainties that are
  calculated during the data processing (e.g. during normalisation) that
  are implausible and are set to missing

#### Uncertainty coefficients

The function `ctsm_uncrt_workup` and related supporting functions are
used in OSPAR assessments to update the fixed and proportional standard
deviations which are subsequently used to impute missing uncertainties.
These functions were ignored during the initial development of harsat
and are now harsat compatible.

#### Biological effect assessments

Imposex assessments: these are now fully reproducible with seeds for
random number generation provided in the calls to `ctsm.VDS.cl` and
`assess_imposex`

Assessment functions for negative binomial data have been added.
Negative binomial data includes MNC - the number of micronucleated
cells.

#### Reporting

`report_assessment` generates default file names. These are based on the
series identifier with additional station information. It is now
possible to override this behaviour for a single report by providing a
different file name using the `output_file` argument.

#### Reference tables

- new values added to method_extraction table

#### Minor bug fixes

- correct behaviour of argument `return_early` in `create_timeseries`  
- pass `info` component of the harsat object to `determinand.link.sum`,
  `determinand.link.replace`, and `determinand.link.imposex`  
- ensure early return from `ctsm_convert_basis` when there is nothing to
  convert (avoids issues e.g. when all the data are biological
  effects)  
- ensure SURVT (in pargroup B-BIO) is recognised as a biological effect
  in `ctsm_get_datatype` (SURVT is the only determinand in this pargroup
  that isn’t an auxiliary variable)  
- pass `good_status` to assessment functions for data with distributions
  other than normal and lognormal
- trap pathological case in estimation of `prtrend`; see \#436
- ensure `ctsm_OHAT_legends` uses the symbology as specified in
  `write_summary_table`

### Version 1.0.0

- Initial public release

### Version 0.1.3

- Various fixes

### Version 0.1.2

- Fixed issues when packaged; see \#326, \#328
- Updated AMAP data and packaging; see \#329

### Version 0.1.1

- Fixed issue with auxiliary variables: see \#289
- Small documentation improvements
- Added build processes for package bundles

### Version 0.1.0

- Initial release
