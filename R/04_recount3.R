# Usar recount3
## ----'start', message=FALSE-----------------------------------
## Load recount3 R package
library("recount3")


## ----'quick_example'------------------------------------------
## Revisemos todos los proyectos con datos de humano en recount3
human_projects <- available_projects()

# Explorar que tipo de objeto tenemos
class(human_projects)
# No es lo mismo que usar typeof
typeof(human_projects)

# Para ver el tamaño
dim(human_projects)
# Para explorarlo
head(human_projects)
tail(human_projects)
head(human_pojects[order(human_projects$n_samples),])
head(human_pojects[order(human_projects$n_samples),])

summary(human_projects)

## Encuentra tu proyecto de interés. Aquí usaremos
## SRP009615 de ejemplo
proj_info <- subset(
  human_projects,
  project == "SRP009615" & project_type == "data_sources"
)

dim(proj_info)

## Crea un objeto de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
## ?create_rse
## Nos dice que solo requiere un renglon en el data
rse_gene_SRP009615 <- create_rse(proj_info)
## Explora el objeto RSE
rse_gene_SRP009615
# 63856 genes 12 muestras
# 175 variables de nuestras muestras
# 10 variables de info de los genes
# 8 variables del metadata

metadata(rse_gene_SRP009615)

rowRanges(rse_gene_SRP009615)

# La longitud del gen en recont3 se usa para normalizar en algunos casos
rowData(rse_gene_SRP009615)

sort(tablerowRanges(rse_gene_SRP009615)

# mcols y rowData es lo mismo
## ----"interactive_display", eval = FALSE----------------------
# Explora los proyectos disponibles de forma interactiva
proj_info_interactive <- interactiveDisplayBase::display(human_projects)

# Selecciona un solo renglón en la tabla y da click en "send".
# BigWigURL
# ## Aquí verificamos que solo seleccionaste un solo renglón.
stopifnot(nrow(proj_info_interactive) == 1)
# ## Crea el objeto RSE
rse_gene_interactive <- create_rse(proj_info_interactive)


## ----"tranform_counts"----------------------------------------
## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
## Para otras transformaciones como RPKM y TPM, revisa transform_counts().
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)


## ----"expand_attributes"--------------------------------------
## Para este estudio en específico, hagamos más fácil de usar la
## información del experimento
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)
colData(rse_gene_SRP009615)[,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP009615)))
]

#iSEE::iSEE(rse_gene_SRP009615)
