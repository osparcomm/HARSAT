library(rprojroot)
working_directory <- is_r_package$find_file()

data_dir = file.path(working_directory, "tests", "datasets", "external-1")
info_dir = file.path(working_directory, "information", "AMAP")

call_read_data <- function(compartment = "biota", data_format = "external", purpose = "AMAP",
                           contaminants = "EXTERNAL_FO_PW_DATA.csv") {
  read_data(
    compartment = compartment,
    purpose = purpose,
    contaminants = contaminants,
    stations = "EXTERNAL_AMAP_STATIONS.csv",
    data_dir = data_dir,
    data_format = data_format,
    info_dir = info_dir
  )
}

test_that('ctsm_get_info reads assessment values', {

  biota_data <- call_read_data()
  biota_data <- tidy_data(biota_data)
  expect_type(biota_data$info, 'list')

  ok <- ctsm_get_info(biota_data$info$species, biota_data$data$species, "assess")
  expect_equal(all(ok), TRUE)

  # Reload the data with a missing species
  biota_data <- call_read_data(contaminants = "EXTERNAL_FO_PW_DATA_missing_species.csv")
  biota_data <- tidy_data(biota_data)

  # Verify the data set contains our spurious species
  expect_equal("Anlobicephala melas" %in% biota_data$data$species, TRUE)

  caller <- function() {
    ctsm_get_info(biota_data$info$species, biota_data$data$species, "assess")
  }

  # Verify we get an error when calling ctsm_get_info
  expect_error(caller(),
    regexp = "^\nThe following values are not in the reference table\\.\nPlease add them to the reference table or edit your data to continue\\.\nAnlobicephala melas"
  )

  # Now, when in error mode, we should get the same message
  expect_error(ctsm_validate_reference(biota_data$info$species, biota_data$data, 'species', warn=FALSE),
    regexp = "^\nThe following values are not in the reference table\\.\nPlease add them to the reference table or edit your data to continue\\.\nAnlobicephala melas"
  )

  # However, in warning mode, we should get a new data frame
  modified <- ctsm_validate_reference(biota_data$info$species, biota_data$data, 'species', warn=TRUE)
  expect_type(modified, 'list')

  # And in this new data frame, the extra species should no longer be present
  expect_equal("Anlobicephala melas" %in% modified$species, FALSE)
})
