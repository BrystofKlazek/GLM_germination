##############################################################################
# GLM Germination Analysis — Shiny + bslib  (v4)
#
# Pipeline:
#   Upload → Transform (long→wide) → Data Prep (cum→increments)
#   → Models (quasibinomial F-tests) → Post-hoc (emmeans) → Plots → Diagnostics
##############################################################################

# ── Packages ────────────────────────────────────────────────────────────────
required_pkgs <- c("shiny", "bslib", "readxl", "readr", "emmeans", "ggplot2",
                   "dplyr", "tidyr", "broom", "DT", "multcomp", "multcompView", "car")
missing <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing)) {
  message("Instaluji balíčky do r: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(readxl); library(readr)
  library(emmeans); library(ggplot2); library(dplyr); library(tidyr)
  library(broom); library(DT); library(multcomp); library(car)
})

# ── Helper ──────────────────────────────────────────────────────────────────
read_upload <- function(path, name) {
  ext <- tolower(tools::file_ext(name))
  if (ext %in% c("xls", "xlsx")) read_excel(path)
  else if (ext == "csv") read_csv(path, show_col_types = FALSE)
  else if (ext == "tsv") read_tsv(path, show_col_types = FALSE)
  else read_csv(path, show_col_types = FALSE)
}

# ── UI ──────────────────────────────────────────────────────────────────────
ui <- page_navbar(
  title = "GLM - analýza vzcházivosti",
  fillable = FALSE,
  theme = bs_theme(
    bootswatch = "lumen",
    base_font  = font_google("Source Sans 3"),
    code_font  = font_google("JetBrains Mono")
  ),

  # ── Sidebar (shared across all tabs) ──
  sidebar = sidebar(
    width = 300,
    card(
      card_header("Nahrát data"),
      fileInput("file", NULL, accept = c(".csv", ".tsv", ".xls", ".xlsx")),
      p(class = "text-muted small",
        "V dlouhém formátu: varianta (postřik či podobné), datum sběru, kumulativní počet. Jeden řádek na květináč/jamu/podobné na datum.")
    ),
    card(
      card_header("Sloupce"),
      selectInput("col_variant", "Varianta", choices = "(nahrajte soubor)", selected = NULL),
      selectInput("col_date", "Datum", choices = "(nahrajte soubor)", selected = NULL),
      selectInput("col_count", "Počet", choices = "(nahrajte soubor)", selected = NULL)
    ),
    card(
      card_header("Nastavení"),
      numericInput("N_seeds", "Semínek ve sledované nádobě", value = 3, min = 1),
      selectInput("p_method", "Korekce P-hodnot pro porovnávání",
                  choices = c("sidak", "bonferroni", "tukey", "scheffe", "none")),
      numericInput("alpha", "Hladina testu", value = 0.05,
                   min = 0.001, max = 0.2, step = 0.01),
      selectInput("model_choice", "Volba modelu (interakce)",
                  choices = c("Automaticky (dle testu)" = "auto",
                              "Vždy s interakcí (varianta × čas)" = "M3",
                              "Vždy bez interakce (varianta + čas)" = "M2"),
                  selected = "auto"),
      selectInput("family_choice", "Rodina modelu",
                  choices = c("Kvazibinomický (doporučeno)" = "quasibinomial",
                              "Binomický" = "binomial"),
                  selected = "quasibinomial"),
      p(class = "text-muted small",
        "Kvazibinomický model zohledňuje overdisperzi (variabilitu navíc) ",
        "a je konzervativnější. Binomický předpokládá, že veškerá variabilita ",
        "je čistě binomická — vhodný pouze pokud je disperze blízká 1."),
      checkboxInput("drop_zero_times", "Vyřadit data, kde nic nevyklíčilo", value = TRUE),
      p(class = "text-muted small",
        "Vyřadí časové body, ve kterých mají všechny nádoby kumulativní počet = 0. ",
        "Tyto body nepřináší informaci o rozdílech mezi variantami a mohou ",
        "způsobit numerickou nestabilitu modelu (extrémní intervaly spolehlivosti).")
    )
  ),

  # ── Tab 1: Upload ──
  nav_panel("1) Nahrát tabulku",
    card(
      card_header("Nahraná tabulka"),
      DTOutput("raw_table"),
      uiOutput("upload_summary")
    )
  ),

  # ── Tab 2: Příprava dat ──
  nav_panel("2) Příprava dat",
    layout_columns(
      col_widths = c(12),

      card(
        card_header("Krok 1: Transformace (dlouhý na široký formát)"),
        card_body(class = "bg-light",
          p(tags$strong("Proč to tu je?"), " V datech se očekává, že je jen jeden sloupec pro datum a květináč se v něm tedy opakuje. ",
            "Tohle to přetvaruje, takže jeden řádek je právě jeden květináč a má víc dat. Tenhle krok je tady hlavně pro kontrolu, že data jsou zpracovávána správně. ",
            "pot_id je vnitřní mechanismus, aby zpracování dat rozpoznalo, který květináč je který."),
          p(tags$strong("Monotónnost -"), " pokud celkový počet semenáčků z jednoho data na druhé klesne (mortalita), ",
            "jsou zde označena. Binomická data očekávají, že bude přírůst pouze kladný. Jinak by se musel dělat mnohem složitější",
	    " multinomický model. Tady to tedy kompenzuji kvazibinomickým modelem, který alespoň trošku tuhle variabilitu navíc",
	    " vyrovnává. Zatím jsou zde zapsány i negativní hodnoty, v části na přípravu dat se však nastaví na nulu.")
        )
      ),

      card(
        card_header("Spustit transformaci dat"),
        layout_columns(
          col_widths = c(6, 6),
          card_body(
            actionButton("transform_btn", "Transformace z dlouhého na širokého",
                         class = "btn-success btn-lg"),
            uiOutput("transform_status")
          ),
          uiOutput("transform_summary")
        )
      ),

      uiOutput("transform_tables_ui"),

      card(
        card_header("Krok 2: Kondicionální přírůstky"),
        card_body(class = "bg-light",
          p(tags$strong("Proč tohle dělám?"), " Kumulativní přírůstky v jednom květináči jsou silně korelované. ",
            "Kondicionální přírůstky, které měří, kolik nových semenáčků vzešlo z doposud nevzešlých semen, ",
            "jsou mnohem méně závislé a tedy můžeme uvažovat kvazibinomický GLM."),
          p(tags$strong("Vzorce:"), tags$br(),
            tags$code("new_germ = max(0, cumulative \u2212 cum_prev)"),
            " S omezením na nulu pokud je v datech úmrtnost", tags$br(),
            tags$code("remaining = N \u2212 cum_prev"), tags$br(),
            "Odpověď/Response/Model: ", tags$code("cbind(new_germ, remaining \u2212 new_germ)"),
            " s kvazibinomickou rodinou.")
        )
      ),

      card(
        card_header("Převod dat"),
        layout_columns(
          col_widths = c(6, 6),
          card_body(
            actionButton("convert_btn", "Převod na kondicionální přírůstky",
                         class = "btn-success btn-lg"),
            uiOutput("convert_status")
          ),
          uiOutput("conversion_summary")
        )
      ),

      uiOutput("prep_tables_ui")
    )
  ),

  # ── Tab 4: Časové modely ──
  nav_panel("3) Časové modely",
    card(
      card_header("Vnořené GLM s časovou interakcí"),
      layout_columns(
        col_widths = c(8, 4),
        p(class = "text-muted",
          "Tři vnořené modely porovnány pomocí statistických testů. ",
          "Kvazibinomický model počítá s overdisperzí (dalo by se říct nepodchycenou složkou modelu) z mortality či podobných srand ",
	  "(přidává do výpočtu parametr \u03C6, který zvýší SE - směrodatnou odchylku - násobně o \u221A\u03C6)."),
        actionButton("run_btn", "Spustit analýzu",
                     class = "btn-success btn-lg", style = "width:100%;")
      )
    ),
    uiOutput("ftest_summary"),
    uiOutput("overdisp_summary"),
    uiOutput("model_choice_info"),
    card(card_header("Testy (porovnání vnořených modelů)"),
         verbatimTextOutput("lr_tests")),
    card(card_header("AIC / BIC porovnání"),
         uiOutput("aic_bic_note"),
         DTOutput("aic_bic_table")),
    card(card_header("Predikované průměry — podmíněný model (emmeans)"),
         uiOutput("emm_info"),
         DTOutput("emm_table")),
    card(card_header("Párová porovnání variant (post-hoc)"),
         uiOutput("posthoc_info"),
         DTOutput("posthoc_table")),
    card(card_header("Statistické skupiny (CLD)",
                     tags$span(class = "text-muted small ms-2",
                       "Mají-li dvě varianty společné písmenko, nejsou na zvolené hladině statisticky rozlišitelné.")),
         DTOutput("cld_table")),
    card(card_header("Shrnutí modelů (technické detaily)"),
         verbatimTextOutput("model_summaries")),

    # ── Kumulativní model (orientační) ──
    card(
      card_header("Orientační kumulativní model"),
      card_body(class = "bg-light",
        p(tags$strong("Pozor:"), " Tento model fituje GLM přímo na kumulativních počtech ",
          tags$code("cbind(cumulative, N - cumulative) ~ treatment * time"), ". ",
          "Pozorování ze stejného květináče v různých datech jsou silně korelovaná, ",
          "takže p-hodnoty a CI jsou pouze orientační a nelze z nich dělat formální závěry. ",
          "Pro formální testování použijte kondicionální model výše nebo jednoduchý test.")
      )
    ),
    uiOutput("cum_model_summary"),
    card(card_header("Kumulativní model — predikované průměry (emmeans)"),
         DTOutput("cum_emm_table")),
    card(card_header("Kumulativní model — párová porovnání"),
         DTOutput("cum_posthoc_table")),
    card(card_header("Kumulativní model — statistické skupiny (CLD)"),
         DTOutput("cum_cld_table"))
  ),

  # ── Tab 5: Jednoduchý test ──
  nav_panel("4) Jednoduchý model",
    card(
      card_header("Test bez času"),
      card_body(class = "bg-light",
        p(tags$strong("Jednoduchý přístup:"), " Nebere v potaz časovou strukturu experimentu a pouze se kouká na výsledné počty v posledním datu ",
          "a následně testuje, jestli se jednotlivé varianty liší v konečné klíčivosti pomocí ",
          "GLM (dle zvolené rodiny modelu). Porovnává ", tags$code("glm(~ varianta)"), " s ",
          tags$code("glm(~ 1)"), " (konstantní model) pomocí statistického testu."),
        p("Odpovídá na otázku: liší se varianty v celkovém počtu vzešlých semenáčků bez ohledu na časové rozložení?")
      )
    ),
    uiOutput("final_summary"),
    uiOutput("final_overdisp_summary"),
    card(card_header("Test: varianta vs konstantní model"),
         verbatimTextOutput("final_ftest")),
    card(card_header("Predikované průměry (emmeans) a statistické skupiny (CLD)"),
         DTOutput("final_cld_table")),
    card(card_header("Párová porovnání variant"),
         DTOutput("final_posthoc_table"))
  ),

  # ── Tab 6: Výsledky a grafy ──
  nav_panel("5) Výsledky a grafy",
    uiOutput("ftest_summary_results"),

    # ── Graf 1: Sloupcový graf celkové vzcházivosti + tabulky ──
    card(full_screen = TRUE,
      card_header("Celková vzcházivost dle varianty (jednoduchý test)"),
      p(class = "text-muted",
        "Modelové odhady (emmeans) z GLM na posledním datu. ",
        "Písmenka = statistické skupiny z jednoduchého testu (záložka 4)."),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel("Upravit popisky grafu",
          textInput("lbl_bar_title", "Název", ""),
          textInput("lbl_bar_subtitle", "Podnadpis", ""),
          textInput("lbl_bar_y", "Osa Y", "Vzcházivost (%)"),
          textInput("lbl_bar_x", "Osa X", "")
        )
      ),
      plotOutput("final_bar_plot", height = "550px")
    ),
    card(card_header("Predikované průměry a statistické skupiny (emmeans + CLD) — jednoduchý test"),
         DTOutput("final_emm_results")),
    card(card_header("Párová porovnání — jednoduchý test"),
         DTOutput("final_posthoc_results")),

    # ── Graf 2: Postup klíčivosti v čase + tabulky ──
    card(full_screen = TRUE,
      card_header("Postup klíčivosti v čase"),
      p(class = "text-muted",
        "Kumulativní vzcházivost v jednotlivých datech měření (orientační model)."),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel("Upravit popisky grafu",
          textInput("lbl_time_title", "Název", "Postup klíčivosti dle varianty"),
          textInput("lbl_time_subtitle", "Podnadpis", ""),
          textInput("lbl_time_y", "Osa Y", "Kumulativní klíčivost (%)"),
          textInput("lbl_time_x", "Osa X", ""),
          textInput("lbl_time_legend", "Legenda", "Datum")
        )
      ),
      plotOutput("timeline_plot", height = "550px")
    ),
    card(card_header("Souhrnná kumulativní vzcházivost varianty (průměr přes časy, orientační)"),
         p(class = "text-muted small", "Jedna hodnota na variantu. Jde o marginální modelový odhad zprůměrovaný přes data měření."),
         DTOutput("cum_emm_results")),
    card(card_header("Statistické skupiny po jednotlivých datech (CLD, orientační)"),
         DTOutput("cum_cld_results")),
    card(card_header("Párová porovnání — kumulativní model (orientační)"),
         DTOutput("cum_posthoc_results")),

    # ── Graf 3: Podmíněná pravděpodobnost + tabulky ──
    card(full_screen = TRUE,
      card_header("Podmíněná pravděpodobnost klíčení (odhady z modelu)"),
      p(class = "text-muted",
        "Modelové odhady podmíněné pravděpodobnosti klíčení (emmeans) s intervalem spolehlivosti. ",
        "Podmíněná = pravděpodobnost, že dosud nevzešlé semínko vzejde v daném intervalu."),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel("Upravit popisky grafu",
          textInput("lbl_emm_title", "Název", ""),
          textInput("lbl_emm_subtitle", "Podnadpis", ""),
          textInput("lbl_emm_y", "Osa Y", ""),
          textInput("lbl_emm_x", "Osa X", ""),
          textInput("lbl_emm_legend", "Legenda", "Datum")
        )
      ),
      plotOutput("main_plot", height = "600px")
    ),
    card(card_header("Souhrnný odhad varianty — podmíněný model (průměr přes časy)"),
         p(class = "text-muted small", "Jedna hodnota na variantu. Jde o marginální modelový odhad zprůměrovaný přes data měření."),
         DTOutput("emm_results")),
    card(card_header("Statistické skupiny po jednotlivých datech (CLD) — podmíněný model"),
         DTOutput("cld_results")),
    card(card_header("Párová porovnání — podmíněný model"),
         DTOutput("posthoc_results")),

    # ── Graf 4: Teplotní mapa ──
    card(full_screen = TRUE,
      card_header("Teplotní mapa vzcházivosti"),
      p(class = "text-muted",
        "Sytost barvy znázorňuje míru vzcházivosti: tmavší = větší vzcházivost. ",
        "Čísla v buňkách = průměrné procento vzcházivosti."),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel("Upravit popisky grafu",
          textInput("lbl_heat_title", "Název", "Teplotní mapa klíčivosti"),
          textInput("lbl_heat_subtitle", "Podnadpis", ""),
          textInput("lbl_heat_y", "Osa Y", ""),
          textInput("lbl_heat_x", "Osa X", "Datum měření"),
          textInput("lbl_heat_legend", "Legenda", "Klíčivost (%)")
        )
      ),
      plotOutput("heatmap_plot", height = "600px")
    )
  ),

  # ── Tab 8: Diagnostiky ──
  nav_panel("6) Diagnostiky",
    card(full_screen = TRUE,
      card_header("Reziduální diagnostika"),
      plotOutput("diag_plot", height = "700px")),
    card(card_header("Kontrola overdisperze"),
         verbatimTextOutput("overdisp_text"))
  ),

  # ── Tab 9: Log ──
  nav_panel("7) Log",
    card(card_header("Záznam analýzy"),
         verbatimTextOutput("log_text"))
  )
)


# ── Server ──────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Helper: use custom label if non-empty, otherwise default
  lbl <- function(input_id, default) {
    val <- input[[input_id]]
    if (is.null(val) || trimws(val) == "") default else val
  }

  rv <- reactiveValues(
    wide_data = NULL, transformed = FALSE, mono_issues = NULL,
    treatment_order = NULL,
    inc_data = NULL, full_data = NULL, converted = FALSE,
    M1 = NULL, M2 = NULL, M3 = NULL,
    lr12 = NULL, lr23 = NULL,
    best_model = NULL, posthoc = NULL, cld_df = NULL, emm_plot = NULL,
    final_M0 = NULL, final_M1 = NULL, final_ftest = NULL, final_cld = NULL, final_posthoc = NULL,
    cum_model = NULL, cum_emm = NULL, cum_posthoc = NULL, cum_cld = NULL,
    # Marginální emmeans (průměr přes časy) pro tab Výsledky
    emm_marginal = NULL, cum_emm_marginal = NULL,
    aic_bic = NULL,
    log = ""
  )

  add_log <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
    rv$log <- paste0(rv$log, msg, "\n")
  }

  reset_all <- function() {
    rv$wide_data <- NULL; rv$transformed <- FALSE; rv$mono_issues <- NULL
    rv$treatment_order <- NULL
    rv$inc_data <- NULL; rv$full_data <- NULL; rv$converted <- FALSE
    rv$M1 <- NULL; rv$M2 <- NULL; rv$M3 <- NULL
    rv$lr12 <- NULL; rv$lr23 <- NULL
    rv$best_model <- NULL; rv$posthoc <- NULL; rv$cld_df <- NULL; rv$emm_plot <- NULL
    rv$final_M0 <- NULL; rv$final_M1 <- NULL; rv$final_ftest <- NULL; rv$final_cld <- NULL; rv$final_posthoc <- NULL
    rv$cum_model <- NULL; rv$cum_emm <- NULL; rv$cum_posthoc <- NULL; rv$cum_cld <- NULL
    rv$emm_marginal <- NULL; rv$cum_emm_marginal <- NULL
    rv$aic_bic <- NULL
  }

  # ── File upload ──
  uploaded <- reactiveVal(NULL)

  observeEvent(input$file, {
    reset_all()
    df <- read_upload(input$file$datapath, input$file$name)
    add_log("Loaded: ", input$file$name, " (", nrow(df), "\u00D7", ncol(df), ")")
    uploaded(df)

    # Update column dropdowns
    cols <- names(df)
    updateSelectInput(session, "col_variant", choices = cols, selected = cols[1])
    updateSelectInput(session, "col_date", choices = cols,
                      selected = if (length(cols) >= 2) cols[2] else cols[1])
    updateSelectInput(session, "col_count", choices = cols,
                      selected = if (length(cols) >= 3) cols[3] else cols[1])
  })

  # (column update now happens inside observeEvent above)

  # ════════════════════════════════════════════════════════════
  # TAB 1: Upload
  # ════════════════════════════════════════════════════════════

  output$raw_table <- renderDT({
    req(uploaded())
    datatable(uploaded(), options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
  })

  output$upload_summary <- renderUI({
    req(uploaded(), input$col_variant, input$col_date)
    df <- uploaded()
    nv <- n_distinct(df[[input$col_variant]])
    nd <- n_distinct(df[[input$col_date]])
    div(class = "alert alert-success mt-3",
      tags$strong("\u2713 "), sprintf("%d varianty, %d data, %d řádků. ", nv, nd, nrow(df)),
      "Následuje část ", tags$strong("2) Transformace"), ". "
    )
  })

  # ════════════════════════════════════════════════════════════
  # TAB 2: Transform
  # ════════════════════════════════════════════════════════════

  observeEvent(input$transform_btn, {
    req(uploaded(), input$col_variant, input$col_date, input$col_count)
    tryCatch({
      df <- uploaded()
      long <- df[, c(input$col_variant, input$col_date, input$col_count)]
      names(long) <- c("treatment", "time", "cumulative")
      long <- long %>%
        mutate(treatment = as.character(treatment),
               time = as.character(time),
               cumulative = as.numeric(cumulative))

      date_order <- unique(long$time)
      long$time <- factor(long$time, levels = date_order)
      add_log("Data: ", paste(date_order, collapse = " \u2192 "))

      # Zachovat pořadí variant podle prvního výskytu ve vstupních datech
      treat_order <- unique(long$treatment)
      long$treatment <- factor(long$treatment, levels = treat_order)
      add_log("Varianty: ", paste(treat_order, collapse = ", "))

      long <- long %>%
        group_by(treatment, time) %>%
        mutate(pot_id = row_number()) %>%
        ungroup()

      wide <- long %>%
        pivot_wider(id_cols = c(treatment, pot_id),
                    names_from = time, values_from = cumulative) %>%
        arrange(treatment, pot_id)

      # Monotonicity check
      mono <- data.frame(treatment = character(), pot_id = integer(),
                         from_date = character(), to_date = character(),
                         from_val = numeric(), to_val = numeric(),
                         stringsAsFactors = FALSE)
      for (i in seq_len(nrow(wide))) {
        for (j in seq_along(date_order)[-1]) {
          v_prev <- wide[[date_order[j - 1]]][i]
          v_curr <- wide[[date_order[j]]][i]
          if (!is.na(v_prev) && !is.na(v_curr) && v_curr < v_prev) {
            mono <- rbind(mono, data.frame(
              treatment = wide$treatment[i], pot_id = wide$pot_id[i],
              from_date = date_order[j - 1], to_date = date_order[j],
              from_val = v_prev, to_val = v_curr, stringsAsFactors = FALSE))
          }
        }
      }

      rv$wide_data <- wide; rv$mono_issues <- mono; rv$transformed <- TRUE
      rv$treatment_order <- treat_order
      rv$inc_data <- NULL; rv$full_data <- NULL; rv$converted <- FALSE
      add_log("\u2713 Transformace: ", nrow(wide), " nádob, ", nrow(mono), " problémů s monotónností")
    }, error = function(e) {
      add_log("\u274C ", conditionMessage(e))
      showNotification(conditionMessage(e), type = "error")
    })
  })

  output$transform_status <- renderUI({
    req(rv$transformed)
    n <- nrow(rv$mono_issues)
    if (n > 0) {
      div(class = "alert alert-warning mt-3",
        sprintf("%d nádob má poznatelnou mortalitu. Záporné hodnoty jsou zatím zachovány, později budou omezeny na nulu.", n),
        tags$br(), "Následuje ", tags$strong("3) Příprava dat"), ".")
    } else {
      div(class = "alert alert-success mt-3",
        sprintf("%d nádob bez znatelných problémů s monotónností. ", nrow(rv$wide_data)),
        "Následuje ", tags$strong("3) Příprava dat"), ".")
    }
  })

  output$transform_summary <- renderUI({
    req(rv$transformed, rv$wide_data)
    wide <- rv$wide_data
    dcols <- setdiff(names(wide), c("treatment", "pot_id"))
    card(class = "bg-light",
      p(tags$strong("Varianty: "), n_distinct(wide$treatment)),
      p(tags$strong("Nádoby: "), nrow(wide)),
      p(tags$strong("Data: "), paste(dcols, collapse = ", ")),
      p(tags$strong("Problémy: "), if (nrow(rv$mono_issues) == 0) "\u2713 Žádné"
        else paste0("\u26A0 ", nrow(rv$mono_issues)))
    )
  })

  output$transform_tables_ui <- renderUI({
    req(rv$transformed)
    tables <- list(
      card(card_header("Široký formát (jeden řádek na květináč)"),
           DTOutput("wide_table"))
    )
    if (nrow(rv$mono_issues) > 0) {
      tables[[2]] <- card(card_header("Porušení monotónnosti"),
                          DTOutput("mono_table"))
    }
    tables
  })

  output$wide_table <- renderDT({
    req(rv$wide_data)
    datatable(rv$wide_data, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
  })

  output$mono_table <- renderDT({
    req(rv$mono_issues, nrow(rv$mono_issues) > 0)
    datatable(rv$mono_issues, options = list(scrollX = TRUE), rownames = FALSE)
  })

  # ════════════════════════════════════════════════════════════
  # TAB 3: Data Prep
  # ════════════════════════════════════════════════════════════

  observeEvent(input$convert_btn, {
    if (!rv$transformed) {
      showNotification("Je prvně potřeba transformace dat v předchozím kroku!", type = "warning"); return()
    }
    tryCatch({
      wide <- rv$wide_data
      N <- input$N_seeds
      dcols <- setdiff(names(wide), c("treatment", "pot_id"))

      long <- wide %>%
        pivot_longer(cols = all_of(dcols), names_to = "time", values_to = "cumulative") %>%
        mutate(time = factor(time, levels = dcols),
               treatment = factor(treatment, levels = rv$treatment_order)) %>%
        arrange(treatment, pot_id, time) %>%
        group_by(treatment, pot_id) %>%
        mutate(cum_prev = lag(cumulative, default = 0),
               new_germ = cumulative - cum_prev,
               remaining = N - cum_prev) %>%
        ungroup()

      n_neg <- sum(long$new_germ < 0, na.rm = TRUE)
      if (n_neg > 0) add_log("Oseknutí na nulu ", n_neg, " záporných přírůstků")

      long <- long %>%
        mutate(new_germ = pmax(0, new_germ))

      # Vyřadit časové body kde u VŠECH nádob kumulativní počet = 0
      if (isTRUE(input$drop_zero_times)) {
        zero_times <- long %>%
          group_by(time) %>%
          summarise(all_zero = all(cumulative == 0), .groups = "drop") %>%
          filter(all_zero) %>%
          pull(time)
        if (length(zero_times) > 0) {
          long <- long %>% filter(!time %in% zero_times)
          # Přepočítat cum_prev po vyřazení — první zbývající timepoint musí mít cum_prev = 0
          long <- long %>%
            mutate(time = droplevels(time)) %>%
            arrange(treatment, pot_id, time) %>%
            group_by(treatment, pot_id) %>%
            mutate(cum_prev = lag(cumulative, default = 0),
                   new_germ = pmax(0, cumulative - cum_prev),
                   remaining = N - cum_prev) %>%
            ungroup()
          add_log("\u2713 Vyřazeno ", length(zero_times), " časových bodů s nulovým klíčením: ",
                  paste(zero_times, collapse = ", "))
        }
      }

      # Uložit kompletní data PŘED filtrací (pro kumulativní grafy)
      rv$full_data <- long

      # Vyřadit pozorování kde remaining <= 0 (všechna semínka už vzešla) — jen pro model
      n_zero <- sum(long$remaining <= 0, na.rm = TRUE)
      if (n_zero > 0) add_log("Vyřazeno ", n_zero, " pozorování s remaining \u2264 0 pro model (všechna semínka vzešla)")
      long <- long %>% filter(remaining > 0)

      # Omezit new_germ na maximum remaining
      long <- long %>%
        mutate(new_germ = pmin(new_germ, remaining))

      rv$inc_data <- long; rv$converted <- TRUE
      rv$M1 <- NULL; rv$M2 <- NULL; rv$M3 <- NULL
      add_log("\u2713 Increments: ", nrow(long), " observations")
    }, error = function(e) {
      add_log("\u274C ", conditionMessage(e))
      showNotification(conditionMessage(e), type = "error")
    })
  })

  output$convert_status <- renderUI({
    req(rv$converted)
    div(class = "alert alert-success mt-3",
      sprintf("%d pozorování připraveno. ", nrow(rv$inc_data)),
      "Následuje část ", tags$strong("3) Časové modely"), ".")
  })

  output$conversion_summary <- renderUI({
    req(rv$converted, rv$inc_data, rv$full_data)
    dat <- rv$inc_data; full <- rv$full_data; N <- input$N_seeds
    mean_rate <- mean(dat$new_germ / dat$remaining, na.rm = TRUE)
    finals <- full %>% group_by(treatment, pot_id) %>%
      summarise(f = max(cumulative), .groups = "drop")
    card(class = "bg-light",
      p(tags$strong("Pozorování (pro model): "), nrow(dat)),
      p(tags$strong("Semínek na nádobu: "), N),
      p(tags$strong("Průměrná podmíněná šance vzrůstu: "), sprintf("%.1f%%", mean_rate * 100)),
      p(tags$strong("Průměrná celková vzcházivost: "), sprintf("%.1f%%", mean(finals$f / N) * 100))
    )
  })

  output$prep_tables_ui <- renderUI({
    req(rv$converted)
    tagList(
      card(card_header("Data s přírůstky"),
           DTOutput("inc_table")),
      card(card_header("Trajektorie jednotlivých nádob"),
           plotOutput("prep_plot", height = "380px"))
    )
  })

  output$inc_table <- renderDT({
    req(rv$inc_data)
    show_cols <- c("treatment", "pot_id", "time", "cumulative", "cum_prev", "new_germ", "remaining")
    datatable(rv$inc_data[, show_cols],
              options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = c("cumulative", "cum_prev", "new_germ", "remaining"), digits = 0)
  })

  output$prep_plot <- renderPlot({
    req(rv$full_data)
    dat <- rv$full_data; N <- input$N_seeds
    ggplot(dat, aes(x = time, y = cumulative / N, colour = treatment,
                    group = interaction(treatment, pot_id))) +
      geom_line(alpha = 0.25, linewidth = 0.5) +
      stat_summary(aes(group = treatment), fun = mean,
                   geom = "line", linewidth = 1.4) +
      stat_summary(aes(group = treatment), fun = mean,
                   geom = "point", size = 3) +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
      labs(x = "Datum", y = "Kumulativní vzcházivost (%)", colour = "Varianta") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  })

  # ════════════════════════════════════════════════════════════
  # TAB 4: Models
  # ════════════════════════════════════════════════════════════

  observeEvent(input$run_btn, {
    if (!rv$converted) {
      showNotification("Je potřeba převod dat!", type = "warning"); return()
    }
    dat <- rv$inc_data
    dat$treatment <- factor(dat$treatment, levels = rv$treatment_order)
    fam <- if (input$family_choice == "binomial") binomial else quasibinomial
    fam_name <- input$family_choice
    test_type <- if (fam_name == "quasibinomial") "F" else "Chisq"
    add_log("Spouští se ", fam_name, " GLM...")

    tryCatch({
      rv$M1 <- glm(cbind(new_germ, remaining - new_germ) ~ time,
                    family = fam, data = dat)
      rv$M2 <- glm(cbind(new_germ, remaining - new_germ) ~ treatment + time,
                    family = fam, data = dat)
      rv$M3 <- glm(cbind(new_germ, remaining - new_germ) ~ treatment * time,
                    family = fam, data = dat)

      add_log("M1 dev=", round(rv$M1$deviance, 1),
              "  M2 dev=", round(rv$M2$deviance, 1),
              "  M3 dev=", round(rv$M3$deviance, 1))

      # ── AIC/BIC — přímo z modelů pokud binomický, jinak z binomických protějšků ──
      if (fam_name == "binomial") {
        rv$aic_bic <- data.frame(
          Model = c("M1 (~time)", "M2 (~treatment + time)", "M3 (~treatment * time)"),
          AIC = c(AIC(rv$M1), AIC(rv$M2), AIC(rv$M3)),
          BIC = c(BIC(rv$M1), BIC(rv$M2), BIC(rv$M3)),
          stringsAsFactors = FALSE
        )
      } else {
        B1 <- glm(cbind(new_germ, remaining - new_germ) ~ time,
                  family = binomial, data = dat)
        B2 <- glm(cbind(new_germ, remaining - new_germ) ~ treatment + time,
                  family = binomial, data = dat)
        B3 <- glm(cbind(new_germ, remaining - new_germ) ~ treatment * time,
                  family = binomial, data = dat)
        rv$aic_bic <- data.frame(
          Model = c("M1 (~time)", "M2 (~treatment + time)", "M3 (~treatment * time)"),
          AIC = c(AIC(B1), AIC(B2), AIC(B3)),
          BIC = c(BIC(B1), BIC(B2), BIC(B3)),
          stringsAsFactors = FALSE
        )
      }
      add_log("AIC: ", paste(sprintf("%.1f", rv$aic_bic$AIC), collapse = ", "))

      rv$lr12 <- anova(rv$M1, rv$M2, test = test_type)
      rv$lr23 <- anova(rv$M2, rv$M3, test = test_type)

      # Správný sloupec pro p-hodnotu závisí na typu testu
      p_col <- if (test_type == "F") "Pr(>F)" else "Pr(>Chi)"
      p_treat <- rv$lr12[[p_col]][2]
      p_inter <- rv$lr23[[p_col]][2]
      add_log("Treatment: p=", format.pval(p_treat, 4),
              "  Interaction: p=", format.pval(p_inter, 4))

      # ── Model selection: manual override or automatic ──
      alpha <- input$alpha
      choice <- input$model_choice
      if (choice == "M3") {
        rv$best_model <- "M3"
        add_log("Model zvolen ručně: M3 (s interakcí)")
      } else if (choice == "M2") {
        rv$best_model <- "M2"
        add_log("Model zvolen ručně: M2 (bez interakce)")
      } else {
        # Automatic: F-test decides
        if (!is.na(p_inter) && p_inter < alpha) {
          rv$best_model <- "M3"
        } else {
          rv$best_model <- "M2"
        }
        add_log("Automatický výběr (dle testu): ", rv$best_model)
      }

      mod <- if (rv$best_model == "M3") rv$M3 else rv$M2

      conf_level <- 1 - input$alpha

      if (rv$best_model == "M3") {
        emm <- emmeans(mod, pairwise ~ treatment | time,
                       adjust = input$p_method, type = "response", level = conf_level)
        rv$posthoc <- as.data.frame(emm$contrasts)
        emm_main <- emmeans(mod, ~ treatment | time, type = "response", level = conf_level)
        rv$emm_plot <- as.data.frame(emm_main)
      } else {
        emm <- emmeans(mod, pairwise ~ treatment,
                       adjust = input$p_method, type = "response", level = conf_level)
        rv$posthoc <- as.data.frame(emm$contrasts)
        emm_main <- emmeans(mod, ~ treatment, type = "response", level = conf_level)
        emm_time <- emmeans(mod, ~ treatment + time, type = "response", level = conf_level)
        rv$emm_plot <- as.data.frame(emm_time)
      }

      # CLD: compute letter groupings
      rv$cld_df <- NULL
      tryCatch({
        if (rv$best_model == "M3") {
          emm_bt <- emmeans(mod, ~ treatment | time, type = "response", level = conf_level)
          cld_res <- multcomp::cld(emm_bt, adjust = input$p_method, Letters = letters, sort = FALSE, alpha = input$alpha)
          cld_out <- as.data.frame(cld_res)
          cld_out$.group <- trimws(cld_out$.group)
          rv$cld_df <- cld_out
          add_log("\u2713 CLD computed for ", length(unique(cld_out$time)), " timepoints")
        } else {
          cld_res <- multcomp::cld(emm_main, adjust = input$p_method, Letters = letters, sort = FALSE, alpha = input$alpha)
          cld_out <- as.data.frame(cld_res)
          cld_out$.group <- trimws(cld_out$.group)
          rv$cld_df <- cld_out
          add_log("\u2713 CLD computed")
        }
      }, error = function(e) {
        add_log("\u26A0 CLD error: ", conditionMessage(e))
        tryCatch({
          if (rv$best_model == "M3") {
            emm_bt <- emmeans(mod, ~ treatment | time, type = "response", level = conf_level)
            cld_res <- multcomp::cld(emm_bt, adjust = input$p_method, Letters = letters, alpha = input$alpha)
            cld_out <- as.data.frame(cld_res)
            if (".group" %in% names(cld_out)) cld_out$.group <- trimws(cld_out$.group)
            rv$cld_df <- cld_out
            add_log("\u2713 CLD computed (fallback)")
          } else {
            cld_res <- multcomp::cld(emm_main, adjust = input$p_method, Letters = letters, alpha = input$alpha)
            cld_out <- as.data.frame(cld_res)
            if (".group" %in% names(cld_out)) cld_out$.group <- trimws(cld_out$.group)
            rv$cld_df <- cld_out
            add_log("\u2713 CLD computed (fallback)")
          }
        }, error = function(e2) {
          add_log("\u26A0 CLD fallback also failed: ", conditionMessage(e2))
        })
      })

      add_log("\u2713 Temporal analysis complete")

      # ── Marginální emmeans (průměr přes časy) pro podmíněný model ──
      tryCatch({
        emm_marg <- emmeans(mod, ~ treatment, type = "response", level = conf_level)
        rv$emm_marginal <- as.data.frame(emm_marg)
        add_log("\u2713 Marginální emmeans (podmíněný model) computed")
      }, error = function(e) {
        add_log("\u26A0 Marginální emmeans selhaly: ", conditionMessage(e))
        rv$emm_marginal <- NULL
      })

      # ── Final Germination: simple GLM on last date only ──
      add_log("Running final germination test...")
      tryCatch({
        # Get last timepoint data — z KOMPLETNÍCH dat (ne filtrovaných)
        full <- rv$full_data
        last_time <- levels(full$time)[length(levels(full$time))]
        final_dat <- full[full$time == last_time, ]
        final_dat$treatment <- factor(final_dat$treatment, levels = rv$treatment_order)
        N <- input$N_seeds

        # M0: constant model (no treatment effect)
        rv$final_M0 <- glm(cbind(cumulative, N - cumulative) ~ 1,
                           family = fam, data = final_dat)
        # M1: treatment effect
        rv$final_M1 <- glm(cbind(cumulative, N - cumulative) ~ treatment,
                           family = fam, data = final_dat)

        rv$final_ftest <- anova(rv$final_M0, rv$final_M1, test = test_type)
        p_col_final <- if (test_type == "F") "Pr(>F)" else "Pr(>Chi)"
        p_final <- rv$final_ftest[[p_col_final]][2]
        add_log("Final germination ", test_type, "-test: p = ", format.pval(p_final, 4))

        # CLD and post-hoc for final germination
        tryCatch({
          emm_final <- emmeans(rv$final_M1, ~ treatment, type = "response", level = conf_level)
          
          # Párová porovnání
          emm_pairs <- emmeans(rv$final_M1, pairwise ~ treatment,
                               adjust = input$p_method, type = "response", level = conf_level)
          rv$final_posthoc <- as.data.frame(emm_pairs$contrasts)
          
          # CLD
          cld_final <- multcomp::cld(emm_final, adjust = input$p_method,
                                      Letters = letters, sort = FALSE, alpha = input$alpha)
          cld_final_df <- as.data.frame(cld_final)
          cld_final_df$.group <- trimws(cld_final_df$.group)
          rv$final_cld <- cld_final_df
          add_log("\u2713 Final CLD + post-hoc computed")
        }, error = function(e) {
          add_log("\u26A0 Final CLD/post-hoc failed: ", conditionMessage(e))
          rv$final_cld <- NULL
          rv$final_posthoc <- NULL
        })

        add_log("\u2713 Final germination analysis complete")
      }, error = function(e) {
        add_log("\u26A0 Final germination error: ", conditionMessage(e))
      })

      # ── Orientační kumulativní model ──
      add_log("Orientační kumulativní model...")
      tryCatch({
        full <- rv$full_data
        full$treatment <- factor(full$treatment, levels = rv$treatment_order)
        N <- input$N_seeds

        rv$cum_model <- glm(cbind(cumulative, N - cumulative) ~ treatment * time,
                            family = fam, data = full)

        cum_emm <- emmeans(rv$cum_model, ~ treatment | time, type = "response", level = conf_level)
        rv$cum_emm <- as.data.frame(cum_emm)

        cum_pairs <- emmeans(rv$cum_model, pairwise ~ treatment | time,
                             adjust = input$p_method, type = "response", level = conf_level)
        rv$cum_posthoc <- as.data.frame(cum_pairs$contrasts)

        tryCatch({
          cum_cld_res <- multcomp::cld(cum_emm, adjust = input$p_method,
                                        Letters = letters, sort = FALSE, alpha = input$alpha)
          cum_cld_df <- as.data.frame(cum_cld_res)
          cum_cld_df$.group <- trimws(cum_cld_df$.group)
          rv$cum_cld <- cum_cld_df
          add_log("\u2713 Kumulativní model: emmeans + CLD + post-hoc hotovo")
        }, error = function(e) {
          add_log("\u26A0 Kumulativní CLD selhalo: ", conditionMessage(e))
          rv$cum_cld <- NULL
        })

        # Marginální emmeans (průměr přes časy) pro kumulativní model
        tryCatch({
          cum_emm_marg <- emmeans(rv$cum_model, ~ treatment, type = "response", level = conf_level)
          rv$cum_emm_marginal <- as.data.frame(cum_emm_marg)
          add_log("\u2713 Marginální emmeans (kumulativní model) computed")
        }, error = function(e) {
          add_log("\u26A0 Kumulativní marginální emmeans selhaly: ", conditionMessage(e))
          rv$cum_emm_marginal <- NULL
        })
      }, error = function(e) {
        add_log("\u26A0 Kumulativní model selhal: ", conditionMessage(e))
        rv$cum_model <- NULL; rv$cum_emm <- NULL; rv$cum_posthoc <- NULL; rv$cum_cld <- NULL
      })

      add_log("\u2713 All analyses complete")
    }, error = function(e) {
      add_log("\u274C ", conditionMessage(e))
      showNotification(conditionMessage(e), type = "error")
    })
  })

  output$lr_tests <- renderPrint({
    req(rv$lr12)
    test_name <- if ("Pr(>F)" %in% names(rv$lr12)) "F-test" else "LR test (Chi-sq)"
    cat(test_name, ": M1 (~datum)  vs  M2 (~varianta + datum)\n")
    cat("H0: Kategorie varianta nemá efekt\n\n")
    print(rv$lr12)
    cat("\n\n", test_name, ": M2 (~varianta + datum)  vs  M3 (~varianta * datum)\n")
    cat("H0: Bez interakce varianta \u00D7 datum\n\n")
    print(rv$lr23)
  })

  output$aic_bic_note <- renderUI({
    if (input$family_choice == "quasibinomial") {
      p(class = "text-muted small",
        "AIC/BIC nelze přímo spočítat pro kvazibinomické modely, ",
        "proto se používají binomické protějšky jako doplňková informace. Nižší = lepší.")
    } else {
      p(class = "text-muted small", "AIC/BIC přímo z binomických modelů. Nižší = lepší.")
    }
  })

  output$aic_bic_table <- renderDT({
    req(rv$aic_bic)
    df <- rv$aic_bic
    # Mark best model
    df$AIC_best <- ifelse(df$AIC == min(df$AIC), "\u2713", "")
    df$BIC_best <- ifelse(df$BIC == min(df$BIC), "\u2713", "")
    datatable(df, options = list(dom = "t", ordering = FALSE), rownames = FALSE) %>%
      formatRound(columns = c("AIC", "BIC"), digits = 1)
  })

  output$model_choice_info <- renderUI({
    req(rv$best_model, rv$aic_bic)
    choice <- input$model_choice
    aic_best <- rv$aic_bic$Model[which.min(rv$aic_bic$AIC)]
    bic_best <- rv$aic_bic$Model[which.min(rv$aic_bic$BIC)]

    manual <- choice %in% c("M2", "M3")
    label <- if (rv$best_model == "M3") "M3 (s interakcí)" else "M2 (bez interakce)"

    div(class = if (manual) "alert alert-warning" else "alert alert-info",
      tags$strong("Zvolený model: "), label,
      if (manual) " (ruční volba)" else " (automaticky dle testu)",
      tags$br(),
      tags$small(
        sprintf("AIC preferuje: %s | BIC preferuje: %s", aic_best, bic_best),
        if (manual) tags$span(class = "text-muted",
          " — Ruční volba přepisuje automatický výběr.")
      )
    )
  })

  # Plain-language F-test summary (shared helper)
  make_ftest_summary <- function(show_simple = FALSE) {
    req(rv$lr12, rv$lr23, rv$best_model)
    alpha <- input$alpha

    # Extrahovat p-hodnotu — robustně najít sloupec s p-hodnotou
    get_p <- function(aov_table, row = 2) {
      p_cols <- grep("^Pr\\(", names(aov_table), value = TRUE)
      if (length(p_cols) == 0) return(NA_real_)
      val <- aov_table[[p_cols[1]]][row]
      if (is.null(val)) NA_real_ else val
    }
    p_treat <- get_p(rv$lr12)
    p_inter <- get_p(rv$lr23)

    treat_sig <- !is.na(p_treat) && p_treat < alpha
    inter_sig <- !is.na(p_inter) && p_inter < alpha
    manual <- input$model_choice %in% c("M2", "M3")
    using_m3 <- rv$best_model == "M3"

    # AIC/BIC summary line
    aic_line <- NULL
    if (!is.null(rv$aic_bic)) {
      aic_best <- rv$aic_bic$Model[which.min(rv$aic_bic$AIC)]
      bic_best <- rv$aic_bic$Model[which.min(rv$aic_bic$BIC)]
      aic_line <- tags$p(class = "text-muted small",
        sprintf("AIC preferuje %s, BIC preferuje %s.", aic_best, bic_best))
    }

    # Simple test result (only on results tab)
    simple_line <- NULL
    if (show_simple && !is.null(rv$final_ftest)) {
      p_final <- get_p(rv$final_ftest)
      final_sig <- !is.na(p_final) && p_final < alpha
      simple_line <- tags$p(
        tags$strong("4. Jednoduchý test (jen poslední datum): "),
        if (final_sig)
          tags$span(style = "color: #d32f2f;",
            sprintf("\u2714 Varianty se liší v celkové vzcházivosti (p = %s)", format.pval(p_final, digits = 3)))
        else
          tags$span(style = "color: #666;",
            sprintf("\u2718 Nebyl detekován rozdíl v celkové vzcházivosti (p = %s)", format.pval(p_final, digits = 3)))
      )
    }

    card(
      card_header("Shrnutí výsledků"),
      tags$div(style = "font-size: 1.05rem; line-height: 1.8;",
        tags$p(tags$strong("Časový model"), " — analyzuje průběh klíčení v čase (kondicionální přírůstky):",
               style = "margin-bottom: 2px; color: #555;"),
        tags$p(
          tags$strong("1. Existuje jednotný efekt varianty v průběhu klíčení? "),
          if (treat_sig)
            tags$span(style = "color: #d32f2f;",
              sprintf("\u2714 ANO (p = %s)", format.pval(p_treat, digits = 3)))
          else
            tags$span(style = "color: #666;",
              sprintf("\u2718 Statisticky nevýznamný rozdíl (p = %s)", format.pval(p_treat, digits = 3))),
          tags$br(),
          tags$small(class = "text-muted",
            "Testuje, zda mají jednotlivé varianty jednotný efekt v čase - při silném výsledku interakce může vyjít nevýznamné, i když celkově varianta 
	    data silně ovlivňuje.")
        ),
        tags$p(
          tags$strong("2. Mění se efekt variant v čase? "),
          if (inter_sig)
            tags$span(style = "color: #d32f2f;",
              sprintf("\u2714 ANO (p = %s)", format.pval(p_inter, digits = 3)))
          else
            tags$span(style = "color: #666;",
              sprintf("\u2718 Ne — efekt konzistentní napříč daty (p = %s)", format.pval(p_inter, digits = 3))),
          tags$br(),
          tags$small(class = "text-muted",
            "Testuje, zda se rozdíly mezi variantami mění v různých datech (např. jedna začne brzy, jiná pozdě). Celkově tedy interakce čas:varianta. Vyjde-li
	    interakce jako významná, hůře se interpretuje test předchozí.")
        ),
        tags$p(
          tags$strong("3. Použitý model: "),
          if (using_m3) "M3 (s interakcí) — varianty porovnány v každém čase zvlášť."
          else "M2 (bez interakce) — varianty porovnány celkově.",
          if (manual) tags$span(style = "color: #e65100;", " (ruční volba)")
        ),
        aic_line,
        if (show_simple) tags$hr(style = "margin: 8px 0;"),
        if (show_simple) tags$p(tags$strong("Jednoduchý test"), " — ignoruje čas, kouká jen na konečný výsledek:",
               style = "margin-bottom: 2px; color: #555;"),
        simple_line,
        if (show_simple && !is.null(simple_line)) tags$small(class = "text-muted",
          "Testuje, zda se varianty liší v celkovém počtu vzešlých semenáčků na konci experimentu. ",
          "Může se lišit od časového modelu — varianty mohou klíčit různě rychle, ale skončit na podobném výsledku.")
      )
    )
  }

  output$ftest_summary <- renderUI({ make_ftest_summary(show_simple = FALSE) })
  output$ftest_summary_results <- renderUI({ make_ftest_summary(show_simple = TRUE) })

  output$model_summaries <- renderPrint({
    req(rv$M1)
    cat("-- M1: ~datum --\n"); print(summary(rv$M1))
    cat("\n-- M2: ~varianta + datum --\n"); print(summary(rv$M2))
    cat("\n-- M3: ~varianta * datum --\n"); print(summary(rv$M3))
  })

  output$overdisp_summary <- renderUI({
    req(rv$best_model)
    mod <- if (rv$best_model == "M3") rv$M3 else rv$M2
    phi <- mod$deviance / mod$df.residual
    is_quasi <- input$family_choice == "quasibinomial"

    if (phi > 1.5) {
      div(class = if (is_quasi) "alert alert-warning" else "alert alert-danger",
        tags$strong("Overdisperze: "),
        sprintf("\u03C6 = %.2f (SE nafouklé \u00D7%.2f)", phi, sqrt(phi)),
        if (!is_quasi) tags$span(
          " — Binomický model ignoruje tuto overdisperzi! P-hodnoty mohou být příliš optimistické. Je k dispozici kvazibinomický model.")
        else tags$span(" — Kvazibinomický model toto zohledňuje v SE a p-hodnotách.")
      )
    } else if (phi < 0.5) {
      div(class = "alert alert-info",
        tags$strong("Underdisperze: "), sprintf("\u03C6 = %.2f", phi))
    } else {
      div(class = "alert alert-success",
        tags$strong("Disperze: "), sprintf("\u03C6 = %.2f — blízká 1, je přiměřený i binomický model.", phi))
    }
  })

  # ════════════════════════════════════════════════════════════
  # TAB 5: Post-hoc
  # ════════════════════════════════════════════════════════════

  output$emm_info <- renderUI({
    req(rv$best_model)
    txt <- if (rv$best_model == "M3")
      "Odhady podmíněné pravděpodobnosti klíčení pro každou variantu v každém datu (model s interakcí)."
    else
      "Odhady podmíněné pravděpodobnosti klíčení pro každou variantu (průměr přes data, model bez interakce)."
    div(class = "alert alert-light", tags$small(txt))
  })

  output$emm_table <- renderDT({
    req(rv$emm_plot)
    df <- rv$emm_plot; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$posthoc_info <- renderUI({
    req(rv$best_model)
    txt <- if (rv$best_model == "M3")
      "Interakce mezi datem a variantou byla uznána za důležitou"
    else
      "Interakce nezohledněna, varianta je porovnávána celkově jako průměr přes data"
    div(class = "alert alert-info",
        tags$strong("Přístup: "), txt,
        paste0(" (s uvažovanou kompenzí testů metodou ", input$p_method, ")"))
  })

  output$posthoc_table <- renderDT({
    req(rv$posthoc)
    df <- rv$posthoc; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$cld_table <- renderDT({
    req(rv$cld_df)
    df <- rv$cld_df; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  # ════════════════════════════════════════════════════════════
  # Orientační kumulativní model — renderery
  # ════════════════════════════════════════════════════════════

  output$cum_model_summary <- renderUI({
    req(rv$cum_model)
    phi <- rv$cum_model$deviance / rv$cum_model$df.residual
    is_quasi <- input$family_choice == "quasibinomial"
    div(class = "alert alert-secondary",
      tags$strong("Kumulativní model (orientační): "),
      sprintf("treatment * time, \u03C6 = %.2f", phi),
      if (is_quasi) sprintf(", SE inflace = \u00D7%.2f", sqrt(phi)) else "",
      tags$br(),
      tags$small(class = "text-muted",
        "Pozorování z téhož květináče jsou korelovaná — výsledky jsou pouze orientační.")
    )
  })

  output$cum_emm_table <- renderDT({
    req(rv$cum_emm)
    df <- rv$cum_emm; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$cum_posthoc_table <- renderDT({
    req(rv$cum_posthoc)
    df <- rv$cum_posthoc; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$cum_cld_table <- renderDT({
    req(rv$cum_cld)
    df <- rv$cum_cld; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  # ════════════════════════════════════════════════════════════
  # TAB 5b: Final Germination
  # ════════════════════════════════════════════════════════════

  output$final_summary <- renderUI({
    req(rv$final_ftest)
    alpha <- input$alpha
    p_cols <- grep("^Pr\\(", names(rv$final_ftest), value = TRUE)
    p_val <- if (length(p_cols) > 0) rv$final_ftest[[p_cols[1]]][2] else NA_real_
    sig <- !is.na(p_val) && p_val < alpha
    phi <- rv$final_M1$deviance / rv$final_M1$df.residual
    is_quasi <- input$family_choice == "quasibinomial"
    test_label <- if (is_quasi) "F-test" else "LR test"

    card(
      card_header("Celková vzcházivost bez časového rozložení"),
      tags$div(style = "font-size: 1.1rem; line-height: 2;",
        tags$p(
          tags$strong("Liší se varianty v celkové vzcházivosti? "),
          if (sig)
            tags$span(style = "color: #d32f2f;",
              sprintf("\u2714 ANO (%s p = %s)", test_label, format.pval(p_val, digits = 3)))
          else
            tags$span(style = "color: #666;",
              sprintf("\u2718 Nebyl detekován rozdíl (%s p = %s)", test_label, format.pval(p_val, digits = 3)))
        ),
        if (is_quasi)
          tags$p(style = "font-size: 0.9rem; color: #666;",
            sprintf("Disperze \u03C6 = %.2f, SE inflace = \u00D7%.2f", phi, sqrt(phi)))
        else
          tags$p(style = "font-size: 0.9rem; color: #666;",
            sprintf("Binomický model (bez korekce na overdisperzi). Disperze \u03C6 = %.2f", phi))
      )
    )
  })

  output$final_overdisp_summary <- renderUI({
    req(rv$final_M1)
    phi <- rv$final_M1$deviance / rv$final_M1$df.residual
    is_quasi <- input$family_choice == "quasibinomial"

    if (phi > 1.5) {
      div(class = if (is_quasi) "alert alert-warning" else "alert alert-danger",
        tags$strong("Overdisperze: "),
        sprintf("\u03C6 = %.2f (SE nafouklé \u00D7%.2f)", phi, sqrt(phi)),
        if (!is_quasi) tags$span(
          " — Binomický model ignoruje tuto overdisperzi! P-hodnoty mohou být příliš optimistické. Je k dispozici kvazibinomický model.")
        else tags$span(" — Kvazibinomický model toto zohledňuje v SE a p-hodnotách.")
      )
    } else if (phi < 0.5) {
      div(class = "alert alert-info",
        tags$strong("Underdisperze: "), sprintf("\u03C6 = %.2f", phi))
    } else {
      div(class = "alert alert-success",
        tags$strong("Disperze: "), sprintf("\u03C6 = %.2f — blízká 1, je přiměřený i binomický model.", phi))
    }
  })

  output$final_ftest <- renderPrint({
    req(rv$final_ftest)
    fam_label <- if (input$family_choice == "quasibinomial") "Kvazibinomický" else "Binomický"
    cat(fam_label, " GLM na celkové kumulativní vzcházivosti\n")
    cat("M0: ~ 1          (bez efektu varianty)\n")
    cat("M1: ~ varianta   (efekt varianty)\n\n")
    print(rv$final_ftest)
    cat("\n")
    cat("-- M1 shrnutí --\n")
    print(summary(rv$final_M1))
  })

  output$final_cld_table <- renderDT({
    req(rv$final_cld)
    df <- rv$final_cld; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$final_posthoc_table <- renderDT({
    req(rv$final_posthoc)
    df <- rv$final_posthoc; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$final_bar_plot <- renderPlot({
    req(rv$full_data, rv$final_cld)
    dat <- rv$full_data; N <- input$N_seeds
    cld_df <- rv$final_cld
    ci_pct <- round((1 - input$alpha) * 100)

    # Modelové odhady a CI přímo z emmeans (rv$final_cld)
    plot_data <- data.frame(
      treatment = cld_df$treatment,
      mean_pct = cld_df$prob * 100,
      ci_lo = cld_df$asymp.LCL * 100,
      ci_hi = cld_df$asymp.UCL * 100,
      letter = trimws(cld_df$.group),
      stringsAsFactors = FALSE
    )
    plot_data$treatment <- factor(plot_data$treatment, levels = rv$treatment_order)

    last_time <- levels(dat$time)[length(levels(dat$time))]

    ggplot(plot_data, aes(x = treatment, y = mean_pct, fill = treatment)) +
      geom_col(alpha = 0.85, width = 0.7) +
      geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                    width = 0.25, linewidth = 0.8) +
      geom_text(aes(label = letter, y = ci_hi + 2),
                size = 5, fontface = "bold", vjust = 0) +
      geom_text(aes(label = sprintf("%.1f%%", mean_pct)),
                vjust = -0.5, size = 3.5, colour = "grey30") +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
      labs(x = lbl("lbl_bar_x", NULL),
           y = lbl("lbl_bar_y", "Vzcházivost (%)"),
           title = lbl("lbl_bar_title", sprintf("Celková vzcházivost v posledním datu (%s)", last_time)),
           subtitle = lbl("lbl_bar_subtitle", sprintf("Modelové odhady (emmeans) \u00B1 %d%% CI. Písmenka = statistické skupiny.", ci_pct))) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", size = 18),
            plot.subtitle = element_text(colour = "grey50", size = 12),
            axis.text.x = element_text(size = 12, face = "bold"),
            panel.grid.major.x = element_blank())
  })

  # ════════════════════════════════════════════════════════════
  # TAB 6: Plots (model-based)
  # ════════════════════════════════════════════════════════════

  # emmeans conditional probability: treatments on x, faceted by time
  output$main_plot <- renderPlot({
    req(rv$emm_plot); df <- rv$emm_plot
    ci_pct <- round((1 - input$alpha) * 100)

    # Zachovat pořadí variant ze vstupního souboru
    df$treatment <- factor(df$treatment, levels = rv$treatment_order)

    if (!"time" %in% names(df)) {
      ggplot(df, aes(x = treatment, y = prob, fill = treatment)) +
        geom_col(width = 0.7, alpha = 0.85) +
        geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.25) +
        scale_y_continuous(labels = scales::percent_format()) +
        labs(x = lbl("lbl_emm_x", NULL),
             y = lbl("lbl_emm_y", "Pravděpodobnost klíčení"),
             title = lbl("lbl_emm_title", "Odhady pravděpodobností klíčení dle varianty"),
             subtitle = lbl("lbl_emm_subtitle", sprintf("Chybové úsečky = %d%% interval spolehlivosti", ci_pct))) +
        theme_minimal(base_size = 14) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
              plot.title = element_text(face = "bold", size = 16),
              plot.subtitle = element_text(colour = "grey50"),
              panel.grid.major.x = element_blank())
    } else {
      ggplot(df, aes(x = treatment, y = prob, fill = time)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
        geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                      position = position_dodge(width = 0.8), width = 0.25) +
        scale_y_continuous(labels = scales::percent_format()) +
        scale_fill_brewer(palette = "Set2") +
        labs(x = lbl("lbl_emm_x", NULL),
             y = lbl("lbl_emm_y", "Podmíněná pravděpodobnost klíčení"),
             fill = lbl("lbl_emm_legend", "Datum"),
             title = lbl("lbl_emm_title", "Podmíněná pravděpodobnost klíčení"),
             subtitle = lbl("lbl_emm_subtitle", sprintf("Odhady pravděpodobností klíčení dle varianty a data (%d%% CI)", ci_pct))) +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
              plot.title = element_text(face = "bold", size = 16),
              plot.subtitle = element_text(colour = "grey50"),
              legend.position = "top",
              panel.grid.major.x = element_blank())
    }
  })

  # Timeline: treatments on x, colored by date, grouped bars with SE
  output$timeline_plot <- renderPlot({
    req(rv$full_data)
    ci_pct <- round((1 - input$alpha) * 100)

    if (!is.null(rv$cum_emm)) {
      timeline <- data.frame(
        treatment = rv$cum_emm$treatment,
        time = rv$cum_emm$time,
        mean_pct = rv$cum_emm$prob * 100,
        ci_lo = rv$cum_emm$asymp.LCL * 100,
        ci_hi = rv$cum_emm$asymp.UCL * 100,
        stringsAsFactors = FALSE
      )
      sub_default <- sprintf("Modelové odhady (emmeans) \u00B1 %d%% CI — orientační (korelovaná data)", ci_pct)
    } else {
      dat <- rv$full_data; N <- input$N_seeds
      z_crit <- qnorm(1 - input$alpha / 2)
      timeline <- dat %>%
        group_by(treatment, time) %>%
        summarise(
          mean_pct = mean(cumulative / N * 100),
          se = sd(cumulative / N * 100) / sqrt(n()),
          .groups = "drop"
        ) %>%
        mutate(ci_lo = pmax(0, mean_pct - z_crit * se),
               ci_hi = pmin(100, mean_pct + z_crit * se))
      sub_default <- sprintf("Surová data, průměr \u00B1 %d%% CI", ci_pct)
    }

    timeline$treatment <- factor(timeline$treatment, levels = rv$treatment_order)

    ggplot(timeline, aes(x = treatment, y = mean_pct, fill = time)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
      geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                    position = position_dodge(width = 0.8), width = 0.25) +
      scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0))) +
      scale_fill_brewer(palette = "Set2") +
      labs(x = lbl("lbl_time_x", NULL),
           y = lbl("lbl_time_y", "Kumulativní klíčivost (%)"),
           fill = lbl("lbl_time_legend", "Datum"),
           title = lbl("lbl_time_title", "Postup klíčivosti dle varianty"),
           subtitle = lbl("lbl_time_subtitle", sub_default)) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold", size = 18),
            plot.subtitle = element_text(colour = "grey50", size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
            legend.position = "top",
            panel.grid.major.x = element_blank())
  })

  # Heatmap: treatments × dates
  output$heatmap_plot <- renderPlot({
    req(rv$full_data)
    dat <- rv$full_data; N <- input$N_seeds

    heat_data <- dat %>%
      group_by(treatment, time) %>%
      summarise(mean_pct = mean(cumulative / N * 100), .groups = "drop")

    heat_data$treatment <- factor(heat_data$treatment, levels = rv$treatment_order)

    ggplot(heat_data, aes(x = time, y = treatment, fill = mean_pct)) +
      geom_tile(colour = "white", linewidth = 1.5) +
      geom_text(aes(label = sprintf("%.1f%%", mean_pct)),
                size = 4.5, fontface = "bold",
                colour = ifelse(heat_data$mean_pct > 50, "white", "grey20")) +
      scale_fill_gradient(low = "#f7f7f7", high = "#1b7837",
                          limits = c(0, NA), name = lbl("lbl_heat_legend", "Klíčivost (%)")) +
      labs(x = lbl("lbl_heat_x", "Datum měření"),
           y = lbl("lbl_heat_y", NULL),
           title = lbl("lbl_heat_title", "Teplotní mapa klíčivosti"),
           subtitle = lbl("lbl_heat_subtitle", "Varianty seřazeny dle celkové vzcházivosti (nejvyšší nahoře)")) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold", size = 18),
            plot.subtitle = element_text(colour = "grey50", size = 12),
            axis.text.y = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(size = 12),
            panel.grid = element_blank(),
            legend.position = "right")
  })

  # ════════════════════════════════════════════════════════════
  # Duplikáty tabulek pro tab Výsledky a grafy
  # ════════════════════════════════════════════════════════════

  # Jednoduchý test — emmeans + CLD v jedné tabulce (jen jedno datum, není co marginalizovat)
  output$final_emm_results <- renderDT({
    req(rv$final_cld)
    df <- rv$final_cld; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$final_posthoc_results <- renderDT({
    req(rv$final_posthoc)
    df <- rv$final_posthoc; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  # Kumulativní model — emmeans = marginální průměr přes časy
  output$cum_emm_results <- renderDT({
    req(rv$cum_emm_marginal)
    df <- rv$cum_emm_marginal; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  # Kumulativní model — CLD = skupiny po jednotlivých datech
  output$cum_cld_results <- renderDT({
    req(rv$cum_cld)
    df <- rv$cum_cld; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$cum_posthoc_results <- renderDT({
    req(rv$cum_posthoc)
    df <- rv$cum_posthoc; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  # Podmíněný model — emmeans = marginální průměr přes časy
  output$emm_results <- renderDT({
    req(rv$emm_marginal)
    df <- rv$emm_marginal; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$cld_results <- renderDT({
    req(rv$cld_df)
    df <- rv$cld_df; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$posthoc_results <- renderDT({
    req(rv$posthoc)
    df <- rv$posthoc; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  # ════════════════════════════════════════════════════════════
  # TAB 8: Diagnostics (technical)
  # ════════════════════════════════════════════════════════════

  output$diag_plot <- renderPlot({
    req(rv$best_model)
    mod <- if (rv$best_model == "M3") rv$M3 else rv$M2
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)); plot(mod)
  })

  output$overdisp_text <- renderPrint({
    req(rv$best_model)
    mod <- if (rv$best_model == "M3") rv$M3 else rv$M2
    phi <- mod$deviance / mod$df.residual
    fam_label <- if (input$family_choice == "quasibinomial") "kvazibinomický" else "binomický"
    cat("Model:", rv$best_model, "(", fam_label, ")\n")
    cat("Reziduální deviance:", round(mod$deviance, 2), "\n")
    cat("Reziduální st. volnosti:", mod$df.residual, "\n")
    cat("Disperze \u03C6:             ", round(phi, 3), "\n")
    if (input$family_choice == "quasibinomial") {
      cat("Inflace SE:             \u00D7", round(sqrt(phi), 3), "\n\n")
      if (phi > 1.5) cat("\u26A0 Výrazná overdisperze. SE nafouklé \u00D7", round(sqrt(phi), 2), ".\n")
      else if (phi < 0.5) cat("\u2139 Underdisperze.\n")
      else cat("\u2713 Disperze blízká 1. Kvazibinomická korekce minimální.\n")
    } else {
      cat("\n")
      if (phi > 1.5) cat("\u26A0 Disperze = ", round(phi, 2), " — data vykazují overdisperzi. Zvažte kvazibinomický model.\n")
      else if (phi < 0.5) cat("\u2139 Underdisperze.\n")
      else cat("\u2713 Disperze blízká 1. Binomický model je přiměřený.\n")
    }
  })

  # ════════════════════════════════════════════════════════════
  # TAB 8: Log
  # ════════════════════════════════════════════════════════════

  output$log_text <- renderPrint({ cat(rv$log) })
}

shinyApp(ui, server)
