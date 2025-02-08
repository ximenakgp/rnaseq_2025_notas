library("sessioninfo")

## Información de reproducibilidad
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
