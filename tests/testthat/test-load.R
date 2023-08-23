test_that('basic functions all exist', {
  expect_type(read_data, 'closure')
  expect_type(ctsm_tidy_data, 'closure')
  expect_type(ctsm_create_timeSeries, 'closure')
  expect_type(run_assessment, 'closure')
  expect_type(check_assessment, 'closure')
  expect_type(write_summary_table, 'closure')
})
