test_that("ess_footnote_text works", {
  expect_equal(
    ess_footnote_text(width = 100),
    paste0(
      "An ESS reduction up to ~60% is not unexpected based on the 2021 survey of NICE's technology\n",
      "appraisals (https://onlinelibrary.wiley.com/doi/full/10.1002/jrsm.1511), whereas a reduction of\n",
      ">75% is less common and it may be considered suboptimal."
    )
  )
})
