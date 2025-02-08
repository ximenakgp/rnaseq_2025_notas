## ----first_rse------------------------------------------------
## Lets build our first SummarizedExperiment object
library("SummarizedExperiment")
## ?SummarizedExperiment

## De los ejemplos en la ayuda oficial

## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6
## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
## Información de nuestros genes
rowRanges <- GRanges(
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))
## Información de nuestras muestras
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)
## Juntamos ahora toda la información en un solo objeto de R
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

## Exploremos el objeto resultante
rse

## Número de genes y muestras
dim(rse)

## IDs de nuestros genes y muestras
dimnames(rse)

## Nombres de tablas de cuentas que tenemos (RPKM, CPM, counts, logcounts, etc)
assayNames(rse)

## El inicio de nuestra tabla de cuentas
head(assay(rse))

## Información de los genes en un objeto de Bioconductor
rowRanges(rse)

## Tabla con información de los genes
rowData(rse) # es idéntico a 'mcols(rowRanges(rse))'

## Tabla con información de las muestras
colData(rse)


## ----rse_exercise---------------------------------------------
## Comando 1
rse[1:2, ]
## Comando 2
rse[, c("A", "D", "F")] # Nos quedamos con todos los genes y solo la columna A, D y F

stopifnot(identical(rse[, c(1,4,6)], rse[, c("A", "D", "F")]))

## ----isee_basic, eval = FALSE---------------------------------
# ## Explora el objeto rse de forma interactiva

initial <- list()

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new(
  "RowDataTable",
  Selected = "gene_1",
  Search = "",
  SearchColumns = c("", "", "", "", ""),
  HiddenColumns = character(0),
  VersionInfo = list(
    iSEE = structure(
      list(c(2L, 18L, 0L)),
      class = c("package_version", "numeric_version")
    )
  ),
  PanelId = c(RowDataTable = 1L),
  PanelHeight = 500L,
  PanelWidth = 4L,
  SelectionBoxOpen = FALSE,
  RowSelectionSource = "---",
  ColumnSelectionSource = "---",
  DataBoxOpen = FALSE,
  RowSelectionDynamicSource = FALSE,
  ColumnSelectionDynamicSource = FALSE,
  RowSelectionRestrict = FALSE,
  ColumnSelectionRestrict = FALSE,
  SelectionHistory = list()
)

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new(
  "FeatureAssayPlot",
  Assay = "counts",
  XAxis = "Column data",
  XAxisColumnData = "Treatment",
  XAxisFeatureName = "gene_1",
  XAxisFeatureSource = "---",
  XAxisFeatureDynamicSource = FALSE,
  YAxisFeatureName = "gene_1",
  YAxisFeatureSource = "RowDataTable1",
  YAxisFeatureDynamicSource = TRUE,
  FacetRowByColData = "Treatment",
  FacetColumnByColData = "Treatment",
  ColorByColumnData = "Treatment",
  ColorByFeatureNameAssay = "counts",
  ColorBySampleNameColor = "#FF0000",
  ShapeByColumnData = "Treatment",
  SizeByColumnData = NA_character_,
  TooltipColumnData = character(0),
  FacetRowBy = "None",
  FacetColumnBy = "None",
  ColorBy = "Column data",
  ColorByDefaultColor = "#000000",
  ColorByFeatureName = "gene_1",
  ColorByFeatureSource = "---",
  ColorByFeatureDynamicSource = FALSE,
  ColorBySampleName = "A",
  ColorBySampleSource = "---",
  ColorBySampleDynamicSource = FALSE,
  ShapeBy = "None",
  SizeBy = "None",
  SelectionAlpha = 0.1,
  ZoomData = numeric(0),
  BrushData = list(),
  VisualBoxOpen = FALSE,
  VisualChoices = "Color",
  ContourAdd = FALSE,
  ContourColor = "#0000FF",
  FixAspectRatio = FALSE,
  ViolinAdd = TRUE,
  PointSize = 1,
  PointAlpha = 1,
  Downsample = FALSE,
  DownsampleResolution = 200,
  CustomLabels = FALSE,
  CustomLabelsText = "A",
  FontSize = 1,
  LegendPointSize = 1,
  LegendPosition = "Bottom",
  HoverInfo = TRUE,
  LabelCenters = FALSE,
  LabelCentersBy = "Treatment",
  LabelCentersColor = "#000000",
  VersionInfo = list(
    iSEE = structure(
      list(c(2L, 18L, 0L)),
      class = c("package_version", "numeric_version")
    )
  ),
  PanelId = c(FeatureAssayPlot = 1L),
  PanelHeight = 500L,
  PanelWidth = 4L,
  SelectionBoxOpen = FALSE,
  RowSelectionSource = "---",
  ColumnSelectionSource = "---",
  DataBoxOpen = FALSE,
  RowSelectionDynamicSource = FALSE,
  ColumnSelectionDynamicSource = FALSE,
  RowSelectionRestrict = FALSE,
  ColumnSelectionRestrict = FALSE,
  SelectionHistory = list()
)

################################################################################
# Settings for Column data plot 1
################################################################################

initial[["ColumnDataPlot1"]] <- new(
  "ColumnDataPlot",
  XAxis = "None",
  YAxis = "Treatment",
  XAxisColumnData = "Treatment",
  FacetRowByColData = "Treatment",
  FacetColumnByColData = "Treatment",
  ColorByColumnData = "Treatment",
  ColorByFeatureNameAssay = "counts",
  ColorBySampleNameColor = "#FF0000",
  ShapeByColumnData = "Treatment",
  SizeByColumnData = NA_character_,
  TooltipColumnData = character(0),
  FacetRowBy = "None",
  FacetColumnBy = "None",
  ColorBy = "None",
  ColorByDefaultColor = "#000000",
  ColorByFeatureName = "gene_1",
  ColorByFeatureSource = "---",
  ColorByFeatureDynamicSource = FALSE,
  ColorBySampleName = "A",
  ColorBySampleSource = "---",
  ColorBySampleDynamicSource = FALSE,
  ShapeBy = "None",
  SizeBy = "None",
  SelectionAlpha = 0.1,
  ZoomData = numeric(0),
  BrushData = list(),
  VisualBoxOpen = FALSE,
  VisualChoices = "Color",
  ContourAdd = FALSE,
  ContourColor = "#0000FF",
  FixAspectRatio = FALSE,
  ViolinAdd = TRUE,
  PointSize = 1,
  PointAlpha = 1,
  Downsample = FALSE,
  DownsampleResolution = 200,
  CustomLabels = FALSE,
  CustomLabelsText = "A",
  FontSize = 1,
  LegendPointSize = 1,
  LegendPosition = "Bottom",
  HoverInfo = TRUE,
  LabelCenters = FALSE,
  LabelCentersBy = "Treatment",
  LabelCentersColor = "#000000",
  VersionInfo = list(
    iSEE = structure(
      list(c(2L, 18L, 0L)),
      class = c("package_version", "numeric_version")
    )
  ),
  PanelId = c(ColumnDataPlot = 1L),
  PanelHeight = 500L,
  PanelWidth = 4L,
  SelectionBoxOpen = FALSE,
  RowSelectionSource = "---",
  ColumnSelectionSource = "---",
  DataBoxOpen = FALSE,
  RowSelectionDynamicSource = FALSE,
  ColumnSelectionDynamicSource = FALSE,
  RowSelectionRestrict = FALSE,
  ColumnSelectionRestrict = FALSE,
  SelectionHistory = list()
)

################################################################################
# Settings for Row data plot 1
################################################################################

initial[["RowDataPlot1"]] <- new(
  "RowDataPlot",
  XAxis = "None",
  YAxis = "feature_id",
  XAxisRowData = "feature_id",
  FacetRowByRowData = "rowRanges_seqnames",
  FacetColumnByRowData = "rowRanges_seqnames",
  ColorByRowData = "feature_id",
  ColorBySampleNameAssay = "counts",
  ColorByFeatureNameColor = "#FF0000",
  ShapeByRowData = "rowRanges_seqnames",
  SizeByRowData = "rowRanges_start",
  TooltipRowData = character(0),
  FacetRowBy = "None",
  FacetColumnBy = "None",
  ColorBy = "None",
  ColorByDefaultColor = "#000000",
  ColorByFeatureName = "gene_1",
  ColorByFeatureSource = "---",
  ColorByFeatureDynamicSource = FALSE,
  ColorBySampleName = "A",
  ColorBySampleSource = "---",
  ColorBySampleDynamicSource = FALSE,
  ShapeBy = "None",
  SizeBy = "None",
  SelectionAlpha = 0.1,
  ZoomData = numeric(0),
  BrushData = list(),
  VisualBoxOpen = FALSE,
  VisualChoices = "Color",
  ContourAdd = FALSE,
  ContourColor = "#0000FF",
  FixAspectRatio = FALSE,
  ViolinAdd = TRUE,
  PointSize = 1,
  PointAlpha = 1,
  Downsample = FALSE,
  DownsampleResolution = 200,
  CustomLabels = FALSE,
  CustomLabelsText = "gene_1",
  FontSize = 1,
  LegendPointSize = 1,
  LegendPosition = "Bottom",
  HoverInfo = TRUE,
  LabelCenters = FALSE,
  LabelCentersBy = "rowRanges_seqnames",
  LabelCentersColor = "#000000",
  VersionInfo = list(
    iSEE = structure(
      list(
        c(2L, 18L, 0L)
      ),
      class = c("package_version", "numeric_version")
    )
  ),
  PanelId = c(RowDataPlot = 1L),
  PanelHeight = 500L,
  PanelWidth = 4L,
  SelectionBoxOpen = FALSE,
  RowSelectionSource = "---",
  ColumnSelectionSource = "---",
  DataBoxOpen = FALSE,
  RowSelectionDynamicSource = FALSE,
  ColumnSelectionDynamicSource = FALSE,
  RowSelectionRestrict = FALSE,
  ColumnSelectionRestrict = FALSE,
  SelectionHistory = list()
)

################################################################################
# Settings for Sample assay plot 1
################################################################################

initial[["SampleAssayPlot1"]] <- new(
  "SampleAssayPlot",
  Assay = "counts",
  XAxis = "Row data",
  XAxisRowData = "rowRanges_seqnames",
  XAxisSampleName = "A",
  XAxisSampleSource = "---",
  XAxisSampleDynamicSource = FALSE,
  YAxisSampleName = "A",
  YAxisSampleSource = "ColumnDataTable1",
  YAxisSampleDynamicSource = TRUE,
  FacetRowByRowData = "rowRanges_seqnames",
  FacetColumnByRowData = "rowRanges_seqnames",
  ColorByRowData = "rowRanges_strand",
  ColorBySampleNameAssay = "counts",
  ColorByFeatureNameColor = "#FF0000",
  ShapeByRowData = "rowRanges_seqnames",
  SizeByRowData = "rowRanges_start",
  TooltipRowData = character(0),
  FacetRowBy = "None",
  FacetColumnBy = "None",
  ColorBy = "Row data",
  ColorByDefaultColor = "#000000",
  ColorByFeatureName = "gene_1",
  ColorByFeatureSource = "---",
  ColorByFeatureDynamicSource = FALSE,
  ColorBySampleName = "A",
  ColorBySampleSource = "---",
  ColorBySampleDynamicSource = FALSE,
  ShapeBy = "None",
  SizeBy = "None",
  SelectionAlpha = 0.1,
  ZoomData = numeric(0),
  BrushData = list(),
  VisualBoxOpen = FALSE,
  VisualChoices = "Color",
  ContourAdd = FALSE,
  ContourColor = "#0000FF",
  FixAspectRatio = FALSE,
  ViolinAdd = TRUE,
  PointSize = 1,
  PointAlpha = 1,
  Downsample = FALSE,
  DownsampleResolution = 200,
  CustomLabels = FALSE,
  CustomLabelsText = "gene_1",
  FontSize = 1,
  LegendPointSize = 1,
  LegendPosition = "Bottom",
  HoverInfo = TRUE,
  LabelCenters = FALSE,
  LabelCentersBy = "rowRanges_seqnames",
  LabelCentersColor = "#000000",
  VersionInfo = list(
    iSEE = structure(
      list(
        c(2L, 18L, 0L)
      ),
      class = c("package_version", "numeric_version")
    )
  ),
  PanelId = c(SampleAssayPlot = 1L),
  PanelHeight = 500L,
  PanelWidth = 4L,
  SelectionBoxOpen = FALSE,
  RowSelectionSource = "---",
  ColumnSelectionSource = "---",
  DataBoxOpen = FALSE,
  RowSelectionDynamicSource = FALSE,
  ColumnSelectionDynamicSource = FALSE,
  RowSelectionRestrict = FALSE,
  ColumnSelectionRestrict = FALSE,
  SelectionHistory = list()
)

################################################################################
# Settings for Column data table 1
################################################################################

initial[["ColumnDataTable1"]] <- new(
  "ColumnDataTable",
  Selected = "A",
  Search = "",
  SearchColumns = "",
  HiddenColumns = character(0),
  VersionInfo = list(
    iSEE = structure(
      list(
        c(2L, 18L, 0L)
      ),
      class = c("package_version", "numeric_version")
    )
  ),
  PanelId = c(ColumnDataTable = 1L),
  PanelHeight = 500L,
  PanelWidth = 4L,
  SelectionBoxOpen = FALSE,
  RowSelectionSource = "---",
  ColumnSelectionSource = "---",
  DataBoxOpen = FALSE,
  RowSelectionDynamicSource = FALSE,
  ColumnSelectionDynamicSource = FALSE,
  RowSelectionRestrict = FALSE,
  ColumnSelectionRestrict = FALSE,
  SelectionHistory = list()
)

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new(
  "ComplexHeatmapPlot",
  Assay = "counts",
  CustomRows = TRUE,
  CustomRowsText = "gene_1\ngene_3\ngene_5",
  ClusterRows = TRUE,
  ClusterRowsDistance = "spearman",
  ClusterRowsMethod = "ward.D2",
  DataBoxOpen = FALSE,
  VisualChoices = "Annotations",
  ColumnData = character(0),
  RowData = character(0),
  CustomBounds = FALSE,
  LowerBound = NA_real_,
  UpperBound = NA_real_,
  AssayCenterRows = TRUE,
  AssayScaleRows = TRUE,
  DivergentColormap = "blue < white < orange",
  ShowDimNames = "Rows",
  LegendPosition = "Bottom",
  LegendDirection = "Horizontal",
  VisualBoxOpen = FALSE,
  NamesRowFontSize = 10,
  NamesColumnFontSize = 10,
  ShowColumnSelection = TRUE,
  OrderColumnSelection = TRUE,
  VersionInfo = list(
    iSEE = structure(
      list(c(2L, 18L, 0L)),
      class = c("package_version", "numeric_version")
    )
  ),
  PanelId = c(ComplexHeatmapPlot = 1L),
  PanelHeight = 500L,
  PanelWidth = 4L,
  SelectionBoxOpen = FALSE,
  RowSelectionSource = "---",
  ColumnSelectionSource = "---",
  RowSelectionDynamicSource = FALSE,
  ColumnSelectionDynamicSource = FALSE,
  RowSelectionRestrict = FALSE,
  ColumnSelectionRestrict = FALSE,
  SelectionHistory = list()
)
library("iSEE")
iSEE::iSEE(rse)


## ----download_sce_layer---------------------------------------
## Descarguemos unos datos de spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")
sce_layer

## Revisemos el tamaño de este objeto
lobstr::obj_size(sce_layer)


## ----explore_sce_layer, eval = FALSE--------------------------
iSEE::iSEE(sce_layer)
