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

genes <- c("Tbx21", "Zeb2", "Gata3", "Rorc", "Bcl6", "Prmd1", "Tcf7", "Hif1a")
genes_tbl <- tibble(
  symbol = c("Tbx21", "Zeb2", "Gata3", "Rorc", "Bcl6", "Prmd1", "Tcf7", "Hif1a"),
  ensemble_id = c("ENSG00000073861", "ENSG00000169554", "ENSG00000107485", "ENSG00000143365", 
                  "ENSG00000113916", "ENSG00000057657", "ENSG00000081059", "ENSG00000100644")
) |> 
   mutate(ensemble_id_modified = paste0(ensemble_id, "_X"))

cell_type_coarse <- tribble(
  ~ cell_type_unified_ensemble, ~ cell_type_aggregated,
  "b","B cell",
  "b memory","B cell",
  "b naive","B cell",
  "blood","Blood",
  "bone","Bone",
  "cartilage","Cartilage",
  "cd14 mono","Monocyte",
  "cd16 mono","Monocyte",
  "cd4 fh em","CD4 EM",
  "cd4 naive","CD4 Other",
  "cd4 tcm","CD4 Other",
  "cd4 tem","CD4 EM",
  "cd4 th1 em","CD4 EM",
  "cd4 th1/th17 em","CD4 EM",
  "cd4 th17 em","CD4 EM",
  "cd4 th2 em","CD4 EM",
  "cd8 naive","CD8",
  "cd8 tcm","CD8",
  "cd8 tem","CD8",
  "cdc","Dendritic",
  "cytotoxic","Cytotoxic",
  "dc","Dendritic",
  "endocrine","Endocrine",
  "endothelial","Endothelial",
  "epidermal","Epithelial",
  "epithelial","Epithelial",
  "erythrocyte","Erythroid",
  "fat","Fat",
  "glial","Glial",
  "granulocyte","Granulocyte",
  "ilc","Innate Lymphoid",
  "immune","Immune",
  "lens","Lens",
  "liver","Liver",
  "macrophage","Macrophage",
  "mait","MAIT",
  "mast","Mast",
  "mesothelial","Mesothelial",
  "monocytic","Monocyte",
  "muscle","Muscle",
  "myoepithelial","Epithelial",
  "neuron","Neuron",
  "nk","NK",
  "nkt","NKT",
  "other","Other",
  "pdc","Dendritic",
  "pericyte","Pericyte",
  "plasma","Plasma",
  "pneumocyte","Pneumocyte",
  "progenitor","Progenitor",
  "renal","Renal",
  "reproductive","Reproductive",
  "secretory","Epithelial",
  "sensory","Neuron",
  "stromal","Stromal",
  "t","T",
  "t cd4","CD4 Other",
  "t cd8","CD8",
  "tgd","Treg",
  "treg","CD4 Other"
)

# Load data
data_dir <- "/Users/shen.m/Documents/GitHub/Immune_study/cd4em_shiny_app/data"

se <- readRDS("data/cd4em_se.rds")

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

  # Interface discription
  p("This app explores the expression of key genes across all cell types in 7513 samples, with a focus on CD4 EM cells.",
    "For each gene, samples within CD4 EM are divided into the top 20% (high expression) and the remaining 80% (low expression).",
    "The app visualises whether high expression in CD4 EM is also observed in other cell types"),

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
  # Plot appears below the dropdowns
  plotOutput("exprPlot", height = "500px")
)

server <- function(input, output, session) {
  output$exprPlot <- renderPlot({
    gene_symbol <- input$gene_symbol  # Use selected gene symbol

    # Map gene symbol to ENSG ID with "_X"
    gene_info <- genes_tbl %>% filter(symbol == gene_symbol)
    gene_id <- gene_info$ensemble_id_modified

    # Label CD4 EM high expression
    label_tbl <- se |>
      join_features(gene_id, shape = "wide") |>
      filter(cell_type_aggregated == "CD4 EM") |>
      mutate(high_cd4 = .data[[gene_id]] >= quantile(.data[[gene_id]], 0.8, na.rm = TRUE)) |>
      select(sample_id, high_cd4)

    # Merge and plot
    plot_data <- se |>
      join_features(gene_id, shape = "wide") |>
      inner_join(label_tbl) |>
      mutate(high_cd4 = ifelse(high_cd4 == TRUE, "High", "Low"))

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
    p <- ggplot(plot_data, aes(x = high_cd4, y = .data[[gene_id]])) +
      geom_boxplot() +
      scale_y_log10() +
      labs(
        x = NULL,
        y = paste0(gene_symbol, " Expression"),
        title = if (length(input$cell_type) == length(cell_types)) {
          paste("Expression of", gene_symbol, "across all cell types and CD4EM groups")
        } else {
          paste("Expression of", gene_symbol,
                if (!is.null(input$cell_type) && length(input$cell_type) > 0 && length(input$cell_type) < length(cell_types)) paste("in", paste(input$cell_type, collapse = ", ")) else "")
        }
      ) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text = element_text(size = 7)
      )

    # Facet if showing multiple cell types
    if (!is.null(input$cell_type) && length(input$cell_type) > 1) {
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

