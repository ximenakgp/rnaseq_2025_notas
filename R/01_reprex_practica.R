# Practica del uso de reprex::reprex()

<details>
  <summary>Detalles</summary>

  ```r
  1 + 3
  #> [1] 4
  stop("algun error con la suma")
  #> Error: algun error con la suma
  options(width = 120)
  sessioninfo::session_info()
  #> ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
  #>  setting  value
  #>  version  R version 4.4.2 (2024-10-31)
  #>  os       Ubuntu 20.04.3 LTS
  #>  system   x86_64, linux-gnu
  #>  ui       X11
  #>  language (EN)
  #>  collate  en_US.UTF-8
  #>  ctype    en_US.UTF-8
  #>  tz       America/Mexico_City
  #>  date     2025-01-29
  #>  pandoc   3.2 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
  #>
  #> ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
  #>  package     * version date (UTC) lib source
  #>  cli           3.6.3   2024-06-21 [2] CRAN (R 4.4.2)
  #>  digest        0.6.37  2024-08-19 [2] CRAN (R 4.4.2)
  #>  evaluate      1.0.3   2025-01-10 [2] CRAN (R 4.4.2)
  #>  fastmap       1.2.0   2024-05-15 [2] CRAN (R 4.4.2)
  #>  fs            1.6.5   2024-10-30 [2] CRAN (R 4.4.2)
  #>  glue          1.8.0   2024-09-30 [2] CRAN (R 4.4.2)
  #>  htmltools     0.5.8.1 2024-04-04 [2] CRAN (R 4.4.2)
  #>  knitr         1.49    2024-11-08 [2] CRAN (R 4.4.2)
  #>  lifecycle     1.0.4   2023-11-07 [2] CRAN (R 4.4.2)
  #>  reprex        2.1.1   2024-07-06 [2] CRAN (R 4.4.2)
  #>  rlang         1.1.5   2025-01-17 [2] CRAN (R 4.4.2)
  #>  rmarkdown     2.29    2024-11-04 [2] CRAN (R 4.4.2)
  #>  rstudioapi    0.17.1  2024-10-22 [2] CRAN (R 4.4.2)
  #>  sessioninfo   1.2.2   2021-12-06 [2] CRAN (R 4.1.2)
  #>  withr         3.0.2   2024-10-28 [2] CRAN (R 4.4.2)
  #>  xfun          0.50    2025-01-07 [2] CRAN (R 4.4.2)
  #>  yaml          2.3.10  2024-07-26 [2] CRAN (R 4.4.2)
  #>
  #>  [1] /home/ximenagp/R/x86_64-pc-linux-gnu-library/4.4
  #>  [2] /usr/local/lib/R/site-library
  #>  [3] /usr/lib/R/site-library
  #>  [4] /usr/lib/R/library
  #>
  #> ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
  ```
