test_that("compact output tiene columnas correctas", {
  df <- annotate_lipid("PC 16:0/18:1")
  expect_true("Class" %in% names(df))
  expect_true("total_cl" %in% names(df))
  expect_equal(df$Class, "PC")
  expect_equal(df$total_cl, 34)
  expect_equal(df$total_cs, 1)
})

test_that("ether lipids se detectan correctamente", {
  df <- annotate_lipid("PC O-18:1/20:4", detail = "standard")
  expect_true(df$is_ether)
  expect_false(df$is_plasmalogen)
  df2 <- annotate_lipid("PE P-18:0/20:4", detail = "standard")
  expect_true(df2$is_plasmalogen)
})

test_that("sphingolipids no son marcados como ISTD", {
  df <- annotate_lipid("Cer d18:1/16:0", detail = "standard")
  expect_false(df$is_istd)
  expect_equal(df$sphingoid_prefix, "d")
})

test_that("ISTDs deuterados se detectan correctamente", {
  df <- annotate_lipid("LPE 18:1(d7)", detail = "standard")
  expect_true(df$is_istd)
})

test_that("15-MHDA se parsea como FA 17:0", {
  df <- annotate_lipid("15-MHDA", detail = "standard")
  expect_equal(df$Class, "FA")
  expect_equal(df$total_cl, 17)
})

test_that("shorthand canonico es correcto", {
  df <- annotate_lipid("PC 16:0/18:1", shorthand = TRUE)
  expect_equal(df$shorthand_lm, "PC(16:0/18:1)")
  df2 <- annotate_lipid("HexCer 18:1/16:0", shorthand = TRUE)
  expect_equal(df2$shorthand_lm, "HexCer(d18:1/16:0)")
})

test_that("tres niveles de detalle funcionan", {
  lipids <- c("PC 16:0/18:1", "TG 16:0/18:1/18:1")
  n_compact  <- ncol(annotate_lipid(lipids, detail = "compact"))
  n_standard <- ncol(annotate_lipid(lipids, detail = "standard"))
  n_full     <- ncol(annotate_lipid(lipids, detail = "full"))
  expect_true(n_compact < n_standard)
  expect_true(n_standard < n_full)
})
