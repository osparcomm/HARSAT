# Add stations to contaminant data from an ICES extraction

Adds the station name and station code to the contaminant data from an
ICES extraction. This is done by either matching the station names
submitted with the data to the station dictionary, or by matching the
sample coordinates to the station dictionary, or a combination of both.

## Usage

``` r
add_stations(data, stations, info)
```

## Arguments

- data:

  A data frame with the contaminant data from an ICES extraction

- stations:

  A data frame with the ICES station dictionary

- info:

  A HARSAT information list which must contain the elements `purpose`,
  `compartment`, and `add_stations`. The latter is a list of control
  parameters supplied through `control_default` or `control_modify`
  which control how the station matching is achieved. See details.

## Value

A data frame containing the contaminant data augmented by variables
containing the station code and the station name

## Details

`info$add_stations` is a list of control parameters that modify the
station matching process. The default values depend on `info$purpose`
and are given at the end of this section. The elements of
`info$add_stations` are:

- `method`: a string specifying whether the stations are matched by
  `"name"`, `"coordinates"`, or `"both"`. If `info$purpose` is
  `"custom"`, `method` is restricted to either `"name"` (the default) or
  `"coordinates"`. If `info$purpose` is `"OSPAR"`, `"HELCOM"` or
  `"AMAP"`, then method is set to `"both"` by default and stations are
  matched by name or coordinates according to rules specified by OSPAR,
  HELCOM or AMAP data providers. Currently, stations are matched by name
  for Denmark, France (biota and water - all years; sediment 2009
  onwards), Germany (biota HELCOM - all years; biota OSPAR, biota AMAP,
  sediment, water 2023 onwards), Ireland, Norway, Portugal, Spain (2005
  onwards), Sweden, The Netherlands (2007 onwards), United Kingdom. All
  other stations are matched by coordinates.

- `area`: a vector of strings containing one or more of `"OSPAR"`,
  `"HELCOM"` and `"AMAP"`; this restricts the stations to those in the
  corresponding convention area(s); `NULL` matches to all stations in
  the station dictionary

- `datatype`: a logical specifying whether the stations should be
  restricted to those with an appropriate datatype. If `TRUE`, a
  contaminant measurement in biota (for example) will only be matched to
  stations with `station_datatype` containing the string `"CF"`.
  Similarly, a biological effect measurement in biota will only be
  matched to stations with `station_datatype` containing the string
  `"EF"`

- `temporal`: a logical with `TRUE` indicating that stations should be
  restricted to those with `station_purpm` containing the string `"T"`

- `governance_type`: a string: `"none"`, `"data"`, `"stations"` or
  `"both"` which controls how data and station governance are used in
  station matching. Data governance information is found in
  `data$is_amap_monitoring`, `data$is_helcom_monitoring` and
  `data$is_ospar_monitoring` which are based on the monitoring programme
  (MPROG) information provided in submissions to the ICES data base.
  Station governance information is found in
  `stations$station_programgovernance` which is provided by data
  submitters to the ICES station dictionary. `governance_type` is used
  in conjunction with `governance_id` (see below) as follows:

  - `"none"` means that data and station governance are both ignored

  - `"data"` means that matching will be restricted by data governance
    but not station governance. For example, if
    `governance_id == "HELCOM"`, then data will only be matched to a
    station if `is_helcom_monitoring == TRUE` (with all stations
    considered regardless of station governance). If `governance_id`
    takes multiple values, e.g. `c("OSPAR", "AMAP")`, then data will
    only be matched to a station if e.g. either
    `is_ospar_monitoring == TRUE` or `is_amap_monitoring == TRUE`.

  - `"stations"` means that matching will be restricted by station
    governance but not by data governance. For example, if
    `governance_id == "HELCOM"`, then the stations will be restricted to
    those where `station_programgovernance` contains `"HELCOM"` (with
    all data considered regardless of data governance). If
    `governance_id` takes mulitple values, e.g. `c("OSPAR", "AMAP")`,
    then the stations will be restricted to those where
    `station_programgovernance` contains e.g. either `"OSPAR"` or
    `"AMAP"`.

  - `"both"` uses both data and station governance. For example, if
    `governance_id == "HELCOM"` then data will only be matched if
    `is_helcom_monitoring == TRUE` and the candidate stations will be
    restricted to those where `station_programmegovernance` contains
    `"HELCOM"`. If `governance_id` contains multiple values, e.g.
    `c("OSPAR", "AMAP")`, then data with `is_ospar_monitoring == TRUE`
    and `is_amap_monitoring == FALSE` are matched to stations where
    `station_programgovernance` contains `"OSPAR"`; data with
    `is_ospar_monitoring == FALSE` and `is_amap_monitoring == TRUE` are
    matched with stations where `station_programgovernance` contains
    `"AMAP"`; but measurements where `is_ospar_monitoring == TRUE` and
    `is_amap_monitoring == TRUE` are matched to stations where
    `station_programgovernance` contains either `"OSPAR"` or `"AMAP"`.

- `governance_id`: used in conjunction with `governance_type`. If
  `governance_type == "none"`, then `governance_id` should be set to
  `NULL`. Otherwise, `governance_id` must be a vector of strings
  containing one or more of `"AMAP"`, `"HELCOM"` and `"OSPAR"`.

- `grouping`: a logical with `TRUE` indicating that stations will be
  grouped into meta-stations as specified by
  `stations$station_asmtmimeparent` which is provided by data submitters
  to the ICES station dictionary. Defaults to `FALSE` apart from when
  `info$purpose == "OSPAR"`.

- `check_coordinates`: a logical with `TRUE` indicating that, when
  stations are matched by name, the sample coordinates must also be
  within the station geometry. Not implemented yet, so defaults to
  `FALSE`.

The default values of `info$add_stations` depend on `info$purpose` as
follows:

- `"AMAP"`:
  `list(method = "both", area = "AMAP", datatype = FALSE, temporal = FALSE, governance_type = "none", governance_id = NULL, group = FALSE, check_coordinates = FALSE)`

- `"HELCOM"`:
  `list(method = "both", area = "HELCOM", datatype = FALSE, temporal = FALSE, governance_type = "none", governance_id = NULL, group = FALSE, check_coordinates = FALSE)`

- `"OSPAR"`:
  `list(method = "both", area = "OSPAR", datatype = TRUE, temporal = TRUE, governance_type = "both", governance_id = c("OSPAR", "AMAP"), group = TRUE, check_coordinates = FALSE)`

- `"custom"`:
  `list(method = "name", area = NULL, datatype = FALSE, temporal = FALSE, governance_type = "none", governance_id = NULL, group = FALSE, check_coordinates = FALSE)`
