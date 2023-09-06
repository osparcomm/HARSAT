test_that('basic functions all exist', {
  expect_type(read_data, 'closure')
  expect_type(tidy_data, 'closure')
  expect_type(create_timeseries, 'closure')
  expect_type(run_assessment, 'closure')
  expect_type(check_assessment, 'closure')
  expect_type(write_summary_table, 'closure')
})
