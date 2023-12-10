library(rprojroot)
working_directory <- is_testthat$find_file()

data_dir = file.path(working_directory, "..", "datasets", "external-2", "data")
info_dir = file.path(working_directory, "..", "datasets", "external-2", "information")

output_dir <- file.path(tempdir(), 'config-tests')
dir.create(output_dir, showWarnings = FALSE)
unlink(file.path(output_dir, '*'))

test_that('plotting works with configurations', {

    # set directives that specify run-type

    My_data_type <- "biota"
    My_purpose <- "AMAP" ## options [AMAP/OSPAR/HELCOM/custom]

    # set directives that control input

    My_data_dir <- data_dir
    My_referencetableset <- info_dir
    My_datafile <- "GL_RS_SAMBA_sw.csv"
    My_stationfile <- "EXTERNAL_AMAP_STATIONS.csv"
    My_dataformat <- "external" ## options [ICES/external]

    # set directives that control output

    My_output_dir <- output_dir
    My_output_file <- "GL_RS_SAMBA_NOSUBSERIES.csv"
    My_outfile_write <- "replace" ## options [replace/append]
    My_graphics_format <- "png" ## options [pdf/png]

    # set directives that control run

    My_alpha_power <- 0.05
    My_target_power <- 0.8
    My_trend_start_year <- ""
    My_trend_end_year <- 2020
    My_minimum_number_years <- 6 ## minimum number of years is not currently applied in harsat 
    My_allowed_gap <- 4
    # assign("My_plot_subset_choice", NULL)
    My_plot_subset_choice <- NULL ## expression(species %in% "Phoca hispida")

    biota_data <- read_data(
        compartment = c(My_data_type),
        purpose = c(My_purpose),
        contaminants = My_datafile,
        stations = My_stationfile,
        data_dir = My_data_dir,
        data_format = c(My_dataformat),
        info_files = list(),
        info_dir = My_referencetableset,
        extraction = NULL,
        max_year = My_trend_end_year,
        oddity_dir = "oddities",
        control = list()
        # control = list(use_stage = TRUE)
    )

    biota_data <- tidy_data(biota_data)

    biota_timeseries <- create_timeseries(
        biota_data,
        determinands = ctsm_get_determinands(biota_data$info),
        determinands.control = NULL,
        oddity_path = "oddities",
        return_early = FALSE,
        print_code_warnings = FALSE,
        get_basis = get_basis_most_common,
        normalise = FALSE,
        normalise.control = list()
    )

    biota_assessment <- run_assessment(
        biota_timeseries,
        subset = NULL,
        AC = NULL,
        get_AC_fn = NULL,
        recent_trend = 20L,
        parallel = FALSE
    )

    ## print("test1")

    check_assessment(biota_assessment, save_result = FALSE)

    write_summary_table(
        biota_assessment,
        output_file = My_output_file,
        output_dir = My_output_dir,
        export = TRUE,
        determinandGroups = NULL,
        classColour = NULL,
        collapse_AC = NULL
    )

    ## Check data files
    data.files <- list.files(My_output_dir, pattern = '\\.csv$')
    expect_equal(data.files, c('GL_RS_SAMBA_NOSUBSERIES.csv'))

    ## print("test2")

    ## If we pass NULL, all plots should be generated.
    plot_assessment(
        biota_assessment,
        subset = NULL,
        output_dir = My_output_dir,
        file_type = c("data", "index"),
        file_format = c(My_graphics_format)
    )
    plot.files <- list.files(My_output_dir, pattern = '\\.png$')
    expect_length(plot.files, 22)
    unlink(file.path(output_dir, '*.png'))

    ## print("test3")

    ## If we pass My_plot_subset_choice set to NULL, same
    My_plot_subset_choice <- NULL
    plot_assessment(
        biota_assessment,
        subset = My_plot_subset_choice,
        output_dir = My_output_dir,
        file_type = c("data", "index"),
        file_format = c(My_graphics_format )
    )
    plot.files <- list.files(My_output_dir, pattern = '\\.png$')
    expect_length(plot.files, 22)
    unlink(file.path(output_dir, '*.png'))

    ## print("test4")

    # If we pass My_plot_subset_choice set to a logical FALSE
    My_plot_subset_choice <- FALSE

    plot_assessment(
        biota_assessment,
        subset = My_plot_subset_choice,
        output_dir = My_output_dir,
        file_type = c("data", "index"),
        file_format = c(My_graphics_format )
    )
    plot.files <- list.files(My_output_dir, pattern = '\\.png$')
    expect_length(plot.files, 0)
    unlink(file.path(output_dir, '*.png'))

    ## print("test5")

    # If we pass My_plot_subset_choice set to a filtered expression
    My_plot_subset_choice <- expression(station_code %in% 'A90')

    plot_assessment(
        biota_assessment,
        subset = My_plot_subset_choice,
        output_dir = My_output_dir,
        file_type = c("data", "index"),
        file_format = c(My_graphics_format )
    )
    plot.files <- list.files(My_output_dir, pattern = '\\.png$')
    expect_length(plot.files, 8)
    expect_setequal(plot.files, c(
        "A90 Greenland Ittoqqortoormiit AMAP CD Phoca hispida LI data.png", 
        "A90 Greenland Ittoqqortoormiit AMAP CD Phoca hispida LI index.png",
        "A90 Greenland Ittoqqortoormiit AMAP HG Phoca hispida LI data.png", 
        "A90 Greenland Ittoqqortoormiit AMAP HG Phoca hispida LI index.png",
        "A90 Greenland Ittoqqortoormiit AMAP HG Phoca hispida MU data.png", 
        "A90 Greenland Ittoqqortoormiit AMAP HG Phoca hispida MU index.png",
        "A90 Greenland Ittoqqortoormiit AMAP SE Phoca hispida LI data.png", 
        "A90 Greenland Ittoqqortoormiit AMAP SE Phoca hispida LI index.png"
    ))
    unlink(file.path(output_dir, '*.png'))

    ## print("test6")

    ## And finally, if we pass an expression directly
    plot_assessment(
        biota_assessment,
        subset = station_code %in% 'A100',
        output_dir = My_output_dir,
        file_type = c("data", "index"),
        file_format = c(My_graphics_format )
    )
    plot.files <- list.files(My_output_dir, pattern = '\\.png$')
    expect_length(plot.files, 6)
    unlink(file.path(output_dir, '*.png'))
})
