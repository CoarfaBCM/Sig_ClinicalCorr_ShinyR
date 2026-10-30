# Signature Activity – Clinical Association Explorer
#
# Interactive front-end around run_clinical_association() (R/clinicalAssociation.R).
# Upload gene signatures + an expression matrix + clinical metadata, run the
# analysis, browse heatmaps / scatter / box plots, and download the complete
# result bundle (CSVs + PDFs) as a single zip file.
#
# Required packages: shiny, pheatmap, ggplot2, DT, commonmark
# Install if needed: install.packages(c("shiny", "pheatmap", "ggplot2", "DT", "commonmark"))
# Also requires a command-line "zip" utility on PATH (standard on macOS/Linux).

suppressPackageStartupMessages({
  library(shiny)
  library(pheatmap)
  library(ggplot2)
  library(DT)
})

source("R/clinicalAssociation.R")

options(shiny.maxRequestSize = 300 * 1024^2)  # allow up to 300 MB uploads

BWR_PAL    <- colorRampPalette(c("blue", "white", "red"))(101)
BWR_BREAKS <- seq(-1, 1, length.out = 102)
RAINBOW_PAL  <- c("grey80", colorRampPalette(c("white", "red"))(100))
RAINBOW_BRKS <- seq(0, 4, length.out = 102)

# ── Bundled example dataset (TB Discovery cohort, Lymphoid/Myeloid dysregulation signatures) ──
EXAMPLE_SIG_DIR     <- "data"
EXAMPLE_SIG_PATTERN <- "sig.*.txt"
EXAMPLE_EXPR        <- "data/TBDiscovery.signatureGenes.csv"  # signature genes only (full 171 MB matrix exhausts shinyapps.io memory)
EXAMPLE_CLINICAL    <- "data/SampleClusters_res1.csv"
EXAMPLE_DATA_DIR    <- "data"

# ── "How it works" flowchart, reused verbatim from the tutorial's Mermaid diagram ──
OVERVIEW_MERMAID <- r"-----(flowchart TD
    A["Signature Files\nfolder + glob pattern"] --> D
    B["Expression CSV\ngenes × samples"] --> D
    C["Clinical CSV\nsamples × variables"] --> D

    D["Load inputs & auto-detect variable types\nclassify columns → Numerical or Categorical"] --> E
    E["Filter expression to signature genes\napply optional sample whitelist"] --> F

    F(["↻ For each clinical variable"])

    F --> G["Subset samples\nkeep only samples with a valid value"]
    G --> H["Per-gene z-score\nz = (expr − μ) / σ  across this variable's samples only"]
    H --> I["Compute signature activity scores\nactivity = Σz(up-genes) − Σz(down-genes)"]

    I --> J{"Variable\ntype?"}
    J -->|"Numerical"| K["Pearson correlation\nr + p-value per signature"]
    J -->|"Categorical"| L["Kruskal-Wallis test\nBH FDR → q-values"]

    K --> M["CSV + Clustered Heatmap PDF\nblue–white–red palette"]
    L --> N["CSV + Clustered Heatmap PDF\ngrey–white–red palette"]

    style A fill:#EDF7F0,stroke:#2E8B57
    style B fill:#EDF7F0,stroke:#2E8B57
    style C fill:#EDF7F0,stroke:#2E8B57
    style D fill:#E3F0F7,stroke:#2D7DA8
    style E fill:#E3F0F7,stroke:#2D7DA8
    style F fill:#F5F9FD,stroke:#8AAABB
    style G fill:#E3F0F7,stroke:#2D7DA8
    style H fill:#E3F0F7,stroke:#2D7DA8
    style I fill:#E3F0F7,stroke:#2D7DA8
    style J fill:#FFF5E0,stroke:#C8900A
    style K fill:#E3F0F7,stroke:#2D7DA8
    style L fill:#F0EBF7,stroke:#7B5EA7
    style M fill:#2D7DA8,color:#FFFFFF,stroke:#2D7DA8
    style N fill:#7B5EA7,color:#FFFFFF,stroke:#7B5EA7
)-----"

# ── Tutorial content, pre-rendered once from the original R package's tutorial ──
# commonmark turns the ```mermaid fence into a dead <pre><code> block (plain text,
# not a diagram) since mermaid.js only renders elements with class="mermaid"; swap
# it for a live mermaid div so the flowchart actually renders in the Tutorial tab.
TUTORIAL_HTML <- commonmark::markdown_html(
  paste(readLines("docs/clinicalAssociation_tutorial.md", warn = FALSE), collapse = "\n"),
  extensions = TRUE
)
TUTORIAL_HTML <- sub(
  '(?s)<pre><code class="language-mermaid">.*?</code></pre>',
  paste0('<div class="mermaid">', OVERVIEW_MERMAID, '</div>'),
  TUTORIAL_HTML,
  perl = TRUE
)

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════

ui <- fluidPage(
  titlePanel("Signature Activity – Clinical Association Explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      wellPanel(
        h4("Quick Start"),
        actionButton("run_example_btn", "Run Example Data",
                     class = "btn-info", width = "100%"),
        br(), br(),
        downloadButton("download_example_btn", "Download Example Data", width = "100%"),
        helpText("Run the bundled example data to see how the analysis works.")
      ),
      hr(),
      h4("Run your own data"),
      h4("1. Signatures"),
      fileInput("sig_files", "Signature files (one per file, multiple allowed)",
                multiple = TRUE, accept = c(".txt", ".tsv", ".xls", ".csv")),
      helpText("Each file: line 1 = signature name; remaining lines = ",
               tags$code("GENE\\tlog2FC"), " (tab-separated)."),

      h4("2. Expression matrix"),
      fileInput("expr_file", "Expression CSV (genes × samples)", accept = ".csv"),
      helpText("First column = gene symbols; remaining columns = samples."),

      h4("3. Clinical metadata"),
      fileInput("clinical_file", "Clinical CSV (samples × variables)", accept = ".csv"),
      helpText("First column = sample IDs. Numeric vs. categorical variables ",
               "are auto-detected."),

      numericInput("min_samples", "Minimum valid samples per variable", value = 6, min = 2),
      numericInput("figure_size", "Heatmap size (inches)", value = 20, min = 4, max = 40),

      actionButton("run_btn", "Run Analysis", class = "btn-primary", width = "100%"),
      hr(),
      uiOutput("download_ui")
    ),

    mainPanel(
      width = 9,
      tags$head(tags$style(HTML("
        .tutorial-content table { border-collapse: collapse; margin: 12px 0; width: 100%; }
        .tutorial-content th, .tutorial-content td { border: 1px solid #ddd; padding: 6px 10px; text-align: left; }
        .tutorial-content th { background-color: #f5f5f5; }
        .tutorial-content code { background-color: #f5f5f5; padding: 1px 4px; border-radius: 3px; }
        .tutorial-content pre { background-color: #f5f5f5; padding: 10px; border-radius: 4px; overflow-x: auto; }
      "))),
      tabsetPanel(
        tabPanel("Overview", br(),
                 h3("Signature Activity – Clinical Association Explorer"),
                 p("A browser front-end around ", tags$code("run_clinical_association()"),
                   ": score gene signatures per sample and test whether that activity ",
                   "tracks with clinical variables – no coding required."),
                 p("Upload one or more gene signatures, an expression matrix, and a clinical ",
                   "metadata table. The app computes a per-sample ", tags$strong("activity score"),
                   " for each signature, auto-detects whether each clinical variable is ",
                   "numerical or categorical, and tests the association (Pearson correlation ",
                   "or Kruskal-Wallis, respectively) – then lets you browse heatmaps and ",
                   "per-signature plots, and download the full result bundle."),
                 tags$ul(
                   tags$li("Multiple signatures at once, auto up/down-gene parsing"),
                   tags$li("Numerical vs. categorical clinical variables auto-detected"),
                   tags$li("Clustered correlation / significance heatmaps"),
                   tags$li("Per-signature scatter (numerical) and violin (categorical) plots"),
                   tags$li("One-click download of every CSV and PDF produced")
                 ),
                 p(tags$em("New here? Click “Run Example Data” in the sidebar to see the ",
                           "app in action with a bundled test dataset, or check the Tutorial tab ",
                           "for full documentation.")),
                 h4("How it works"),
                 tags$script(src = "https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.min.js"),
                 tags$div(class = "mermaid", OVERVIEW_MERMAID),
                 tags$script(HTML("if (window.mermaid) { mermaid.initialize({startOnLoad: true}); }"))
        ),
        tabPanel("Run Log", br(), verbatimTextOutput("run_log")),
        tabPanel("Heatmaps", br(),
                 radioButtons("heatmap_type", "Clinical variable type",
                              choices = c("Numerical – Pearson correlation" = "numerical",
                                          "Categorical – Kruskal-Wallis significance" = "categorical"),
                              inline = TRUE),
                 plotOutput("heatmap_plot", height = "750px")),
        tabPanel("Association Plots", br(),
                 fluidRow(
                   column(6, selectInput("var_select", "Clinical variable", choices = NULL)),
                   column(6, selectInput("sig_select", "Signature", choices = NULL))
                 ),
                 plotOutput("assoc_plot", height = "500px")),
        tabPanel("Result Tables", br(),
                 h5("Significant Pearson correlations (p < 0.05; others = 0)"),
                 DTOutput("pcc_table"),
                 hr(),
                 h5("Kruskal-Wallis −log10(p-value)"),
                 DTOutput("pval_table")),
        tabPanel("Tutorial", br(),
                 div(class = "tutorial-content", HTML(TUTORIAL_HTML)))
      )
    )
  )
)

# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════

server <- function(input, output, session) {

  rv <- reactiveValues(
    has_run     = FALSE,
    out_dir     = NULL,
    base_name   = NULL,
    dir_support = NULL,
    var_files   = NULL,   # named: clinical variable -> support CSV path
    var_types   = NULL    # named: clinical variable -> "numerical" / "categorical"
  )

  # Shared pipeline: runs run_clinical_association() against whichever
  # signature folder / expression / clinical paths are passed in, then
  # populates rv with everything the plot/table/download outputs need.
  # Used by both "Run Analysis" (uploaded files) and "Run Example Data".
  run_pipeline <- function(sig_dir, sig_pattern, expr_path, clinical_path) {
    run_id      <- paste0("run_", as.integer(Sys.time()), "_", sample.int(1e6, 1))
    output_root <- file.path(tempdir(), run_id, "results")

    log_lines <- character(0)
    run_error <- NULL

    withProgress(message = "Running clinical association analysis...", value = 0.2, {
      withCallingHandlers(
        {
          tryCatch({
            run_clinical_association(
              sig_folder    = sig_dir,
              sig_pattern   = sig_pattern,
              expression    = expr_path,
              clinical      = clinical_path,
              output_root   = output_root,
              min_samples   = input$min_samples,
              support_plots = TRUE,
              figure_size   = input$figure_size,
              save_temp     = TRUE
            )
          }, error = function(e) {
            run_error <<- conditionMessage(e)
          })
        },
        message = function(m) {
          log_lines <<- c(log_lines, conditionMessage(m))
          invokeRestart("muffleMessage")
        },
        warning = function(w) {
          log_lines <<- c(log_lines, paste0("Warning: ", conditionMessage(w), "\n"))
          invokeRestart("muffleWarning")
        }
      )
      incProgress(0.8)
    })

    output$run_log <- renderText({
      msg <- paste(log_lines, collapse = "")
      if (!is.null(run_error)) msg <- paste0(msg, "\nERROR: ", run_error, "\n")
      msg
    })

    if (!is.null(run_error)) {
      showNotification(paste("Run failed:", run_error), type = "error", duration = NULL)
      rv$has_run <- FALSE
      return(invisible(NULL))
    }

    rv$out_dir     <- normalizePath(dirname(output_root))
    rv$base_name   <- basename(output_root)
    rv$dir_support <- file.path(rv$out_dir, paste0("support-clinical.", rv$base_name))

    clinical_split <- load_and_split_clinical(clinical_path)

    num_vars <- NULL
    if (!is.null(clinical_split$numerical)) {
      num_names <- setdiff(colnames(clinical_split$numerical), clinical_split$sample_col)
      num_files <- file.path(rv$dir_support,
                              paste0(rv$base_name, ".Numerical_variable.",
                                     make.names(num_names), ".activity_and_clinvar.csv"))
      names(num_files) <- num_names
      num_vars <- num_files[file.exists(num_files)]
    }

    cat_vars <- NULL
    if (!is.null(clinical_split$categorical)) {
      cat_names <- setdiff(colnames(clinical_split$categorical), clinical_split$sample_col)
      cat_files <- file.path(rv$dir_support,
                              paste0(rv$base_name, ".Categorical_variable.",
                                     make.names(cat_names), ".activity_and_clinvar.csv"))
      names(cat_files) <- cat_names
      cat_vars <- cat_files[file.exists(cat_files)]
    }

    rv$var_files <- c(num_vars, cat_vars)
    rv$var_types <- setNames(
      c(rep("numerical", length(num_vars)), rep("categorical", length(cat_vars))),
      c(names(num_vars), names(cat_vars))
    )

    updateSelectInput(session, "var_select", choices = names(rv$var_files))

    rv$has_run <- TRUE
    showNotification("Analysis complete.", type = "message")
  }

  observeEvent(input$run_btn, {
    if (is.null(input$sig_files)) {
      showNotification("Please upload at least one signature file.", type = "error"); return()
    }
    if (is.null(input$expr_file)) {
      showNotification("Please upload an expression matrix CSV.", type = "error"); return()
    }
    if (is.null(input$clinical_file)) {
      showNotification("Please upload a clinical metadata CSV.", type = "error"); return()
    }

    work_dir <- file.path(tempdir(), paste0("upload_", as.integer(Sys.time()), "_", sample.int(1e6, 1)))
    sig_dir  <- file.path(work_dir, "signatures")
    dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)

    for (i in seq_len(nrow(input$sig_files))) {
      file.copy(input$sig_files$datapath[i],
                file.path(sig_dir, input$sig_files$name[i]), overwrite = TRUE)
    }

    run_pipeline(sig_dir, "*", input$expr_file$datapath, input$clinical_file$datapath)
  })

  observeEvent(input$run_example_btn, {
    run_pipeline(EXAMPLE_SIG_DIR, EXAMPLE_SIG_PATTERN, EXAMPLE_EXPR, EXAMPLE_CLINICAL)
  })

  # ── Signature dropdown follows the selected clinical variable ──────────────

  observeEvent(input$var_select, {
    req(input$var_select, rv$var_files[[input$var_select]])
    df <- tryCatch(read.csv(rv$var_files[[input$var_select]], check.names = FALSE, row.names = 1),
                   error = function(e) NULL)
    req(df)
    updateSelectInput(session, "sig_select",
                       choices = setdiff(colnames(df), input$var_select))
  })

  # ── Heatmaps (one tab, type-switchable) ─────────────────────────────────

  output$heatmap_plot <- renderPlot({
    req(rv$has_run, input$heatmap_type)

    if (input$heatmap_type == "numerical") {
      f <- file.path(rv$out_dir, paste0(rv$base_name, ".numerical_var_correlations_significant_pcc.csv"))
      validate(need(file.exists(f), "No numerical variables passed the minimum-sample filter."))
      mat <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
      pheatmap(mat, color = BWR_PAL, breaks = BWR_BREAKS,
               cluster_rows = nrow(mat) >= 2, cluster_cols = ncol(mat) >= 2,
               clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               main = "Pearson correlation (r): signature activity vs. numerical clinical variables\n(significant only, p < 0.05; others set to 0)",
               fontsize = 9, border_color = "black")
    } else {
      f <- file.path(rv$out_dir, paste0(rv$base_name, ".categorical_variables.minusLog10_pvalues.csv"))
      validate(need(file.exists(f), "No categorical variables passed the minimum-sample/group filter."))
      mat <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
      pheatmap(mat, color = RAINBOW_PAL, breaks = RAINBOW_BRKS,
               cluster_rows = nrow(mat) >= 2, cluster_cols = ncol(mat) >= 2,
               clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               main = "Kruskal-Wallis significance (−log10 p-value): signature activity vs. categorical clinical variables",
               fontsize = 9, border_color = "black")
    }
  })

  # ── Association plots (one tab: scatter for numerical, boxplot for categorical) ──

  output$assoc_plot <- renderPlot({
    req(input$var_select, input$sig_select, rv$var_files[[input$var_select]])
    var_type <- rv$var_types[[input$var_select]]
    df <- read.csv(rv$var_files[[input$var_select]], check.names = FALSE, row.names = 1)
    validate(need(input$sig_select %in% colnames(df), "Select a signature."))

    if (identical(var_type, "numerical")) {
      x <- df[[input$var_select]]
      y <- df[[input$sig_select]]
      ct <- tryCatch(cor.test(x, y, method = "pearson"), error = function(e) NULL)
      subtitle <- if (!is.null(ct)) sprintf("Pearson r = %.3f    p = %.3g", ct$estimate, ct$p.value) else ""

      ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
        geom_point(color = "red", alpha = 0.7) +
        geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7, formula = y ~ x) +
        labs(x = input$var_select, y = paste(input$sig_select, "activity score"),
             title = sprintf("%s activity vs. %s (numerical – scatter plot)", input$sig_select, input$var_select),
             subtitle = subtitle) +
        theme_classic(base_size = 13)
    } else {
      plot_df <- data.frame(group = factor(df[[input$var_select]]), score = df[[input$sig_select]])
      kw_p    <- tryCatch(signif(kruskal.test(score ~ group, data = plot_df)$p.value, 3), error = function(e) NA)
      anova_p <- tryCatch(signif(summary(aov(score ~ group, data = plot_df))[[1]][["Pr(>F)"]][1], 3),
                           error = function(e) NA)

      ggplot(plot_df, aes(x = group, y = score)) +
        geom_violin(fill = "steelblue", alpha = 0.6, trim = FALSE, color = "black") +
        geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
        labs(x = input$var_select, y = paste(input$sig_select, "activity score"),
             title = sprintf("%s activity across %s groups (categorical – violin plot)", input$sig_select, input$var_select),
             subtitle = sprintf("Kruskal-Wallis p = %s    ANOVA p = %s", kw_p, anova_p)) +
        theme_classic(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  })

  # ── Result tables ────────────────────────────────────────────────────────

  # write.csv() leaves the row-name column header blank; read it back as real
  # row names and promote it to a properly named column before handing the
  # data to DT (a blank column header otherwise breaks DataTables' JS).
  read_result_table <- function(f) {
    mat <- read.csv(f, row.names = 1, check.names = FALSE)
    cbind(`Clinical variable` = rownames(mat), mat)
  }

  output$pcc_table <- renderDT({
    req(rv$has_run)
    f <- file.path(rv$out_dir, paste0(rv$base_name, ".numerical_var_correlations_significant_pcc.csv"))
    validate(need(file.exists(f), "No numerical variables were tested."))
    datatable(read_result_table(f), options = list(scrollX = TRUE), rownames = FALSE)
  })

  output$pval_table <- renderDT({
    req(rv$has_run)
    f <- file.path(rv$out_dir, paste0(rv$base_name, ".categorical_variables.minusLog10_pvalues.csv"))
    validate(need(file.exists(f), "No categorical variables were tested."))
    datatable(read_result_table(f), options = list(scrollX = TRUE), rownames = FALSE)
  })

  # ── Download everything ──────────────────────────────────────────────────

  output$download_ui <- renderUI({
    if (isTRUE(rv$has_run)) {
      downloadButton("download_all", "Download All Results (.zip)",
                      class = "btn-success", width = "100%")
    } else {
      helpText("Run the analysis to enable downloads.")
    }
  })

  output$download_all <- downloadHandler(
    filename = function() paste0("clinical_association_results_", format(Sys.Date(), "%Y%m%d"), ".zip"),
    content = function(file) {
      req(rv$has_run)
      items <- setdiff(list.files(rv$out_dir, full.names = FALSE), "signatures")
      old_wd <- setwd(rv$out_dir)
      on.exit(setwd(old_wd), add = TRUE)
      zip_path <- normalizePath(file, mustWork = FALSE)
      utils::zip(zipfile = zip_path, files = items, flags = "-r9Xq")
    },
    contentType = "application/zip"
  )

  output$download_example_btn <- downloadHandler(
    filename = function() "example_data_TBDiscovery.zip",
    content = function(file) {
      items <- list.files(EXAMPLE_DATA_DIR, full.names = FALSE)
      old_wd <- setwd(EXAMPLE_DATA_DIR)
      on.exit(setwd(old_wd), add = TRUE)
      zip_path <- normalizePath(file, mustWork = FALSE)
      utils::zip(zipfile = zip_path, files = items, flags = "-r9Xq")
    },
    contentType = "application/zip"
  )
}

shinyApp(ui, server)
