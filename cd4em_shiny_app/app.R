# This script is used to explore the expression of Tbx21 across all cell types from CD4EM.
#  CD4 EM cell type, split samples based on top 20% and rest of samples for 8 different genes. 
#  For each gene, identify whether gene expressed highly in CD4 EM also expressed highly in other cell types, 
#  such as B, monocytic cell etc.
# TL;DR: It answers the question whether Tbx21 CD4EM outliers have truely higher expression across other cell types.

library(dplyr)
library(tibble)
library(ggplot2)
library(shiny)
library(tidySingleCellExperiment)
library(shinyWidgets)

friendly_cols <- dittoSeq::dittoColors()

# Set theme
custom_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )
  )

genes <- c("Tbx21", "Zeb2", "Gata3", "Rorc", "Bcl6", "Prmd1", "Tcf7", "Hif1a", "Lmo4")
genes_tbl <- tibble(
  symbol = c("Tbx21", "Zeb2", "Gata3", "Rorc", "Bcl6", "Prmd1", "Tcf7", "Hif1a", "Lmo4"),
  ensemble_id = c("ENSG00000073861", "ENSG00000169554", "ENSG00000107485", "ENSG00000143365", 
                  "ENSG00000113916", "ENSG00000057657", "ENSG00000081059", "ENSG00000100644",
                  "ENSG00000143013")
)

# Generate pseudobulk
# metadata = get_metadata() |>
#   join_census_table() |>
#   keep_quality_cells() |>
#   filter(feature_count > 500)
# 
# cd4_samples = metadata |>
#   filter(cell_type_unified_ensemble |> str_detect("cd4")) |>
#   select(sample_id) |> pull() |> unique()
# 
# se = metadata|>
# 
#   # Extract the same samples that all 8 cells present
#   filter(sample_id %in% cd4_samples) |>
#   get_pseudobulk(features = genes_tbl$ensemble_id)
# 
# # Filter colSums > 0
# se = se[, colSums(assay(se, "counts")) > 0]
# 
# # change id to symbol from gene tbl
# id_to_symbol <- genes_tbl$symbol
# names(id_to_symbol) <- genes_tbl$ensemble_id
# rownames(se) <- id_to_symbol[rownames(se)]
# 
# # Re-aggregate counts by (sample_id, cell_type_aggregated).
# # The HDF5 pseudobulk was computed at sample+cell_type_unified_ensemble level, so multiple
# # fine-grained subtypes that map to the same coarse label must be summed before CPM is calculated.
# cell_type_coarse <- tribble(
#   ~ cell_type_unified_ensemble, ~ cell_type_aggregated,
#   "b","B cell",
#   "b memory","B cell",
#   "b naive","B cell",
#   "blood","Blood",
#   "bone","Bone",
#   "cartilage","Cartilage",
#   "cd14 mono","Monocyte",
#   "cd16 mono","Monocyte",
#   "cd4 fh em","CD4 Fh EM",
#   "cd4 naive","CD4 Other",
#   "cd4 tcm","CD4 Other",
#   "cd4 tem","CD4 EM",
#   "cd4 th1 em","CD4 EM",
#   "cd4 th1/th17 em","CD4 EM",
#   "cd4 th17 em","CD4 EM",
#   "cd4 th2 em","CD4 EM",
#   "cd8 naive","CD8",
#   "cd8 tcm","CD8",
#   "cd8 tem","CD8",
#   "cdc","Dendritic",
#   "cytotoxic","Cytotoxic",
#   "dc","Dendritic",
#   "endocrine","Endocrine",
#   "endothelial","Endothelial",
#   "epidermal","Epithelial",
#   "epithelial","Epithelial",
#   "erythrocyte","Erythroid",
#   "fat","Fat",
#   "glial","Glial",
#   "granulocyte","Granulocyte",
#   "ilc","Innate Lymphoid",
#   "immune","Immune",
#   "lens","Lens",
#   "liver","Liver",
#   "macrophage","Macrophage",
#   "mait","MAIT",
#   "mast","Mast",
#   "mesothelial","Mesothelial",
#   "monocytic","Monocyte",
#   "muscle","Muscle",
#   "myoepithelial","Epithelial",
#   "neuron","Neuron",
#   "nk","NK",
#   "nkt","NKT",
#   "other","Other",
#   "pdc","Dendritic",
#   "pericyte","Pericyte",
#   "plasma","Plasma",
#   "pneumocyte","Pneumocyte",
#   "progenitor","Progenitor",
#   "renal","Renal",
#   "reproductive","Reproductive",
#   "secretory","Epithelial",
#   "sensory","Neuron",
#   "stromal","Stromal",
#   "t","T",
#   "t cd4","CD4 Other",
#   "t cd8","CD8",
#   "tgd","Treg",
#   "treg","CD4 Other"
# )
# 
# {
#   count_mat <- as.matrix(assay(se, "counts"))
#   cd        <- as.data.frame(colData(se))
#   agg_key   <- paste0(cd$sample_id, "___", cd$cell_type_aggregated)
#   unique_keys <- unique(agg_key)
# 
#   agg_counts <- vapply(unique_keys, function(k) {
#     idx <- which(agg_key == k)
#     if (length(idx) == 1L) count_mat[, idx] else rowSums(count_mat[, idx, drop = FALSE])
#   }, numeric(nrow(se)))
#   colnames(agg_counts) <- unique_keys
# 
#   agg_cd <- cd[match(unique_keys, agg_key), , drop = FALSE]
#   rownames(agg_cd) <- unique_keys
# 
#   se <- SummarizedExperiment::SummarizedExperiment(
#     assays  = list(counts = agg_counts),
#     rowData = rowData(se),
#     colData = S4Vectors::DataFrame(agg_cd)
#   )
#   se <- se[, colSums(assay(se, "counts")) > 0]
#   assay(se, "cpm") <- edgeR::cpm(se, log = TRUE, prior.count = 1)
# }

# Load data
data_dir <- "/Users/shen.m/Documents/GitHub/Immune_study/cd4em_shiny_app/data"
se <- HDF5Array::loadHDF5SummarizedExperiment("./data/cd4_se_h5")

# Prepare choices for dropdowns
cell_types <- unique(as.character(se$cell_type_aggregated))
disease_groups <- unique(as.character(se$disease_groups))
tissue_groups <- unique(as.character(se$tissue_groups))
ethnicity_groups <- unique(as.character(se$ethnicity_groups))
age_groups <- unique(as.character(se$age_groups))
sex <- unique(as.character(se$sex))
assay <- unique(as.character(se$assay))


# To show overview visualisation for tbx21 first, then users can select cell type and expression level of interest
ui <- fluidPage(
  titlePanel("CD4EM Gene Expression Explorer"),

  p("This app explores the expression of key genes across all cell types in 5286 samples, with a focus on CD4 EM cells.",
    "For each gene, samples within CD4 EM are divided into the top 20% (high expression) and the remaining 80% (low expression).",
    "The app visualises whether high expression in CD4 EM is also observed in other cell types.",
    "Raw counts are re-summed per (sample, aggregated cell type) before CPM normalisation."),

  p(tags$strong("Two views are available via the 'Stratify by CD4 EM high/low expression' checkbox:"),
    tags$br(),
    tags$b("Unchecked (default):"),
    "Samples are split into High (top 20%) and Low (remaining 80%) based on CD4 EM expression of the selected gene.",
    "Each cell type is shown as a separate facet, answering the question:",
    tags$em("Do samples with high CD4 EM expression also show elevated expression in other cell types?"),
    tags$br(),
    tags$b("Checked:"),
    "No High/Low splitting is applied.",
    "All CD4-containing samples are shown together, with cell type on the x-axis, answering the question:",
    tags$em("What is the baseline distribution of this gene across cell types within CD4-containing samples?")),

  tags$details(
    tags$summary(tags$strong("Cell type aggregation mapping (click to expand)")),
    tags$p(
      "Fine-grained cellNexus labels (", tags$code("cell_type_unified_ensemble"), ") are collapsed into coarse groups ",
      "(", tags$code("cell_type_aggregated"), ") as follows:"
    ),
    tags$ul(
      tags$li(tags$strong("CD4 EM"), ": cd4 tem, cd4 th1 em, cd4 th1/th17 em, cd4 th17 em, cd4 th2 em"),
      tags$li(tags$strong("CD4 Fh EM"), ": cd4 fh em"),
      tags$li(tags$strong("CD4 Other"), ": cd4 naive, cd4 tcm, t cd4, treg"),
      tags$li(tags$strong("CD8"), ": cd8 naive, cd8 tcm, cd8 tem, t cd8"),
      tags$li(tags$strong("B cell"), ": b, b memory, b naive"),
      tags$li(tags$strong("Monocyte"), ": cd14 mono, cd16 mono, monocytic"),
      tags$li(tags$strong("Dendritic"), ": cdc, dc, pdc"),
      tags$li(tags$strong("Epithelial"), ": epidermal, epithelial, myoepithelial, secretory"),
      tags$li(tags$strong("Neuron"), ": neuron, sensory"),
      tags$li(tags$strong("Treg"), ": tgd"),
      tags$li(tags$strong("T"), ": t"),
      tags$li("All remaining labels map 1-to-1: Blood, Bone, Cartilage, Cytotoxic, Endocrine,",
              "Endothelial, Erythroid, Fat, Glial, Granulocyte, Immune, Innate Lymphoid, Lens,",
              "Liver, Macrophage, MAIT, Mast, Mesothelial, Muscle, NK, NKT, Other, Pericyte,",
              "Plasma, Pneumocyte, Progenitor, Renal, Reproductive, Stromal.")
    )
  ),

  # Add gene symbol dropdown to the top row
  fluidRow(
    column(4,
      selectInput("gene_symbol", "Gene Symbol:", choices = genes_tbl$symbol, selected = genes_tbl$symbol[1])
    ),
    column(4,
      pickerInput("cell_type", "Cell Type:", choices = cell_types, selected = cell_types, multiple = TRUE, options = list('actions-box' = TRUE, 'live-search' = TRUE, size = 5))
    ),
    column(4,
      pickerInput("disease_group", "Disease Groups:", choices = disease_groups, selected = disease_groups, multiple = TRUE, options = list('actions-box' = TRUE, 'live-search' = TRUE, size = 5))
    ),
    column(4,
      pickerInput("tissue_group", "Tissue Groups:", choices = tissue_groups, selected = tissue_groups, multiple = TRUE, options = list('actions-box' = TRUE, 'live-search' = TRUE, size = 5))
    ),
    column(4,
      pickerInput("ethnicity_group", "Ethnicity Groups:", choices = ethnicity_groups, selected = ethnicity_groups, multiple = TRUE, options = list('actions-box' = TRUE, 'live-search' = TRUE, size = 5))
    ),
    column(4,
      pickerInput("age_group", "Age Groups:", choices = age_groups, selected = age_groups, multiple = TRUE, options = list('actions-box' = TRUE, 'live-search' = TRUE, size = 5))
    ),
    column(4,
      pickerInput("sex", "Sex:", choices = sex, selected = sex, multiple = TRUE, options = list('actions-box' = TRUE, 'live-search' = TRUE, size = 5))
    ),
    column(4,
      pickerInput("assay", "Assay:", choices = assay, selected = assay, multiple = TRUE, options = list('actions-box' = TRUE, 'live-search' = TRUE, size = 5))
    )
  ),
  fluidRow(
    column(4,
      checkboxInput("stratify_cd4em_expression", "Stratify by CD4 EM high/low expression", value = FALSE)
    ),
    column(4,
      checkboxInput("stratify_sex", "Stratify by sex", value = FALSE)
    ),
    column(4,
      checkboxInput("use_CPM_counts", "Use CPM counts", value = FALSE)
    )
  ),
  # Plot appears below the dropdowns
  plotOutput("exprPlot", height = "500px")
)

server <- function(input, output, session) {
  output$exprPlot <- renderPlot({
    gene_symbol <- input$gene_symbol  # Use selected gene symbol

    # Map gene symbol to ENSG ID with "_X"
    gene_info <- genes_tbl %>% filter(symbol == gene_symbol)
    gene_id <- gene_info$symbol

    selected_assay <- if (isTRUE(input$use_CPM_counts)) "cpm" else "counts"

    base_data <- se |> join_features(gene_id, shape = "wide", assay = selected_assay)

    if (isTRUE(input$stratify_cd4em_expression)) {
      label_tbl <- base_data |>
        filter(cell_type_aggregated == "CD4 EM") |>
        mutate(high_cd4 = .data[[gene_id]] >= quantile(.data[[gene_id]], 0.8, na.rm = TRUE)) |>
        select(sample_id, high_cd4)

      plot_data <- base_data |>
        inner_join(label_tbl) |>
        mutate(high_cd4 = ifelse(high_cd4 == TRUE, "High", "Low"))

      x_var   <- "high_cd4"
      x_label <- "CD4 EM expression group"
    } else {
      plot_data <- base_data
      x_var   <- "cell_type_aggregated"
      x_label <- "Cell type"
    }

    # Check if any filter is empty
    validate(
      need(!is.null(input$cell_type) && length(input$cell_type) > 0, 'Please select at least one Cell Type.'),
      need(!is.null(input$disease_group) && length(input$disease_group) > 0, 'Please select at least one Disease Group.'),
      need(!is.null(input$tissue_group) && length(input$tissue_group) > 0, 'Please select at least one Tissue Group.'),
      need(!is.null(input$ethnicity_group) && length(input$ethnicity_group) > 0, 'Please select at least one Ethnicity Group.'),
      need(!is.null(input$age_group) && length(input$age_group) > 0, 'Please select at least one Age Group.'),
      need(!is.null(input$sex) && length(input$sex) > 0, 'Please select at least one Sex.'),
      need(!is.null(input$assay) && length(input$assay) > 0, 'Please select at least one Assay.')
    )

    # Multi-selection filtering
    plot_data <- plot_data |> filter(cell_type_aggregated %in% input$cell_type)
    plot_data <- plot_data |> filter(disease_groups %in% input$disease_group)
    plot_data <- plot_data |> filter(tissue_groups %in% input$tissue_group)
    plot_data <- plot_data |> filter(ethnicity_groups %in% input$ethnicity_group)
    plot_data <- plot_data |> filter(age_groups %in% input$age_group)
    plot_data <- plot_data |> filter(sex %in% input$sex)
    plot_data <- plot_data |> filter(assay %in% input$assay)

    # Plot
    title_str <- if (length(input$cell_type) == length(cell_types)) {
      paste("Expression of", gene_symbol, "across all cell types and CD4EM groups")
    } else {
      paste("Expression of", gene_symbol,
            if (!is.null(input$cell_type) && length(input$cell_type) > 0 && length(input$cell_type) < length(cell_types)) paste("in", paste(input$cell_type, collapse = ", ")) else "")
    }

    expr_label <- paste0(gene_symbol, " Expression (", toupper(selected_assay), ")")

    sex_palette <- c("#0072B2", "#D55E00", "#999999")   # Wong colour-blind-safe blue / vermillion / grey (unknown)

    if (isTRUE(input$stratify_sex)) {
      p <- ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[gene_id]], fill = sex, color = sex)) +
        geom_boxplot(alpha = 0.5, outlier.size = 0.5, position = position_dodge(width = 0.8)) +
        scale_y_log10() +
        scale_fill_manual(values = sex_palette) +
        scale_color_manual(values = sex_palette) +
        labs(x = x_label, y = expr_label, title = title_str, fill = "Sex", color = "Sex") +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text = element_text(size = 7)
        )
    } else {
      p <- ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[gene_id]])) +
        geom_boxplot() +
        scale_y_log10() +
        labs(x = x_label, y = expr_label, title = title_str) +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text = element_text(size = 7)
        )
    }

    # When stratifying by expression, facet by cell type if multiple selected.
    # When not stratifying, cell type is already on the x-axis so no facet is needed.
    if (isTRUE(input$stratify_cd4em_expression) && !is.null(input$cell_type) && length(input$cell_type) > 1) {
      p <- p + facet_grid(~ cell_type_aggregated)
    }
    p
  })
} 

# Run app
shinyApp(ui, server)

# Run locally
# shiny::runApp("cd4em_shiny_app/")



# Before deployment, need to make all dependencies reproducible. For self-developed package, need to install from GitHub
#   renv::install("MangiolaLaboratory/cellNexus")
#   renv::snapshot()
#   rsconnect::deployApp('cd4em_shiny_app/')

