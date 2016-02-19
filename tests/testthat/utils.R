expect_rng_unchanged <- function(f) {
  oldseed <- .GlobalEnv$.Random.seed
  f()
  expect_identical(oldseed, .GlobalEnv$.Random.seed)
}

expect_different <- function(x, y) {
  eval(bquote(expect_false( identical(.(x), .(y)) )))
}
