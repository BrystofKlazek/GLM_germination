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
  if (ext %in% c("xls", "xlsx")) {
    read_excel(path)
  } else if (ext == "csv") {
    read_csv(path, show_col_types = FALSE)
  } else if (ext == "tsv") {
    read_tsv(path, show_col_types = FALSE)
  } else {
    read_csv(path, show_col_types = FALSE)
  }
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
      radioButtons(
        "contrast_mode",
        "Rezim post-hoc porovnani",
        choices = c(
          "Proti kontrole (Dunnett)" = "dunnett",
          "Vsechny pary" = "allpairs"
        ),
        selected = "dunnett"
      ),
      conditionalPanel(
        condition = "input.contrast_mode == 'dunnett'",
        selectInput("control_variant", "Kontrolní varianta (pro Dunnett)",
                    choices = "(nahrajte soubor)", selected = NULL),
        p(class = "text-muted small",
          "Dunnett porovnava kazdou variantu jen proti zvolene kontrole. ",
          "Mene prisny nez all-pairs. Vhodny kdyz je hlavni otazka: lisi se osetreni od kontroly? ",
          "Písmenka v grafech: a = jako kontrola, b = lepší než kontrola, c = horší než kontrola.")
      ),
      conditionalPanel(
        condition = "input.contrast_mode == 'allpairs'",
        selectInput("p_method", "Korekce P-hodnot (pro všechny páry)",
                    choices = c("sidak", "bonferroni", "tukey", "scheffe", "none"))
      ),
      numericInput("alpha", "Hladina testu", value = 0.05,
                   min = 0.001, max = 0.2, step = 0.01),
      radioButtons(
        "infer_scale",
        "Škála pro contrasts a CLD",
        choices = c(
          "Model scale (doporučeno)" = "link",
          "Response scale (experimentální)" = "response"
        ),
        selected = "link"
      ),
      p(class = "text-muted small",
        "Model scale = statisticky standardní a nejlépe odpovídá inferenci modelu. ",
        "Response scale = interpretačně přirozenější, ale experimentální a méně statisticky přesné. Contrasts i CLD se mohou lišit od standardního modelového pohledu."
      ),
	conditionalPanel(
	  condition = "input.contrast_mode == 'allpairs'",
	  checkboxInput(
		  "show_group_bands",
		  "Zobrazit grouping bands (simultánní intervaly)",
		  value = TRUE
		),
		p(class = "text-muted small",
		  "Tenké přerušované úsečky = simultánní intervaly pro celou rodinu odhadů. ",
		  "Mají lépe odpovídat CLD než běžné CI, ale nejsou s CLD úplně ekvivalentní."
	)),
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
    card(card_header("Kontrasty variant (rozdíly + adjustované CI)"),
         uiOutput("posthoc_info"),
         conditionalPanel(
           condition = "input.contrast_mode == 'allpairs'",
           DTOutput("posthoc_table")
         ),
         conditionalPanel(
           condition = "input.contrast_mode == 'dunnett'",
           DTOutput("dunnett_table")
         )),
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(card_header("Statistické skupiny (CLD)",
                       tags$span(class = "text-muted small ms-2",
                         "Mají-li dvě varianty společné písmenko, nejsou na zvolené hladině statisticky rozlišitelné.")),
           DTOutput("cld_table"))
    ),
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
    card(card_header("Kumulativní model — kontrasty variant (rozdíly + adjustované CI)"),
         conditionalPanel(
           condition = "input.contrast_mode == 'allpairs'",
           DTOutput("cum_posthoc_table")
         ),
         conditionalPanel(
           condition = "input.contrast_mode == 'dunnett'",
           DTOutput("dunnett_cum_table")
         )),
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(card_header("Kumulativní model — statistické skupiny (CLD)"),
           DTOutput("cum_cld_table"))
    )
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
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(card_header("Predikované průměry (emmeans) a statistické skupiny (CLD)"),
           DTOutput("final_cld_table"))
    ),
    card(card_header("Kontrasty variant (rozdíly + adjustované CI)"),
         conditionalPanel(
           condition = "input.contrast_mode == 'allpairs'",
           DTOutput("final_posthoc_table")
         ),
         conditionalPanel(
           condition = "input.contrast_mode == 'dunnett'",
           DTOutput("dunnett_final_table")
         ))
  ),

  # ── Tab 6: Výsledky a grafy ──
  nav_panel("5) Výsledky a grafy",
    uiOutput("ftest_summary_results"),

    # ── Graf 1: Sloupcový graf celkové vzcházivosti + tabulky ──
    card(full_screen = TRUE,
      card_header("Celková vzcházivost dle varianty (jednoduchý test)"),
      p(class = "text-muted",
        "Modelové odhady (emmeans) z GLM na posledním datu. ",
        "Písmenka: v režimu CLD = statistické skupiny, v režimu Dunnett = a (jako kontrola), b (lepší), c (horší)."),
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
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(full_screen = TRUE,
          card_header("Jednoduchý test — matice signifikance kontrastů"),
          p(class = "text-muted", "Zelená = bez průkazného rozdílu po korekci, červená = průkazný rozdíl."),
          plotOutput("final_sig_heatmap", height = "650px")
      )
    ),
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(full_screen = TRUE,
          card_header("Jednoduchý test — forest plot kontrastů"),
          p(class = "text-muted", "Intervaly kontrastů po korekci; svislá čára v nule = žádný rozdíl."),
          plotOutput("final_contrast_plot", height = "650px")
      )
    ),
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(card_header("Predikované průměry a statistické skupiny (emmeans + CLD) — jednoduchý test"),
           DTOutput("final_emm_results"))
    ),
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(card_header("Kontrasty variant — jednoduchý test"),
           DTOutput("final_posthoc_results"))
    ),
    card(full_screen = TRUE,
        card_header("Jednoduchý test — Dunnett forest plot (proti kontrole)"),
        p(class = "text-muted", "Kontrasty kazde varianty jen proti kontrole s Dunnettovou korekci."),
        plotOutput("dunnett_final_contrast_plot", height = "550px")
    ),
    conditionalPanel(
      condition = "input.contrast_mode == 'dunnett'",
      card(card_header("Kontrasty proti kontrole (Dunnett) — jednoduchý test"),
           DTOutput("dunnett_final_results"))
    ),

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
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(full_screen = TRUE,
          card_header("Kumulativní model — matice signifikance kontrastů"),
          p(class = "text-muted", "Fasetová matice podle dat; zelená = bez rozdílu, červená = rozdíl."),
          plotOutput("cum_sig_heatmap", height = "650px")
      )
    ),
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(full_screen = TRUE,
          card_header("Kumulativní model — forest plot kontrastů"),
          p(class = "text-muted", "Kontrasty mezi variantami v jednotlivých datech s adjustovanými intervaly."),
          plotOutput("cum_contrast_plot", height = "650px")
      )
    ),
    card(card_header("Souhrnná kumulativní vzcházivost varianty (průměr přes časy, orientační)"),
         p(class = "text-muted small", "Jedna hodnota na variantu — průměr modelových odhadů přes data měření. ",
           "Sloupec n_times = přes kolik časů se průměrovalo."),
         DTOutput("cum_emm_results")),
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(card_header("Statistické skupiny po jednotlivých datech (CLD, orientační)"),
           DTOutput("cum_cld_results"))
    ),
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(card_header("Kontrasty variant — kumulativní model (orientační)"),
           DTOutput("cum_posthoc_results"))
    ),
    card(full_screen = TRUE,
        card_header("Kumulativní model — Dunnett forest plot (proti kontrole)"),
        p(class = "text-muted", "Kontrasty proti kontrole v jednotlivych datech s Dunnettovou korekci (orientacni)."),
        plotOutput("dunnett_cum_contrast_plot", height = "650px")
    ),
    conditionalPanel(
      condition = "input.contrast_mode == 'dunnett'",
      card(card_header("Kontrasty proti kontrole (Dunnett) — kumulativní model (orientační)"),
           DTOutput("dunnett_cum_results"))
    ),

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
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(full_screen = TRUE,
          card_header("Podmíněný model — matice signifikance kontrastů"),
          p(class = "text-muted", "Fasetová matice podle dat; vychází přímo z adjustovaných párových kontrastů."),
          plotOutput("cond_sig_heatmap", height = "650px")
      )
    ),
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(full_screen = TRUE,
          card_header("Podmíněný model — forest plot kontrastů"),
          p(class = "text-muted", "Kontrasty mezi variantami v jednotlivých datech s adjustovanými intervaly."),
          plotOutput("cond_contrast_plot", height = "650px")
      )
    ),
    card(card_header("Souhrnný odhad varianty — podmíněný model (průměr přes časy)"),
         p(class = "text-muted small", "Jedna hodnota na variantu — průměr modelových odhadů přes data měření. ",
           "Sloupec n_times = přes kolik časů se průměrovalo (saturované varianty mohou mít méně)."),
         DTOutput("emm_results")),
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(card_header("Statistické skupiny po jednotlivých datech (CLD) — podmíněný model"),
           DTOutput("cld_results"))
    ),
    conditionalPanel(
      condition = "input.contrast_mode == 'allpairs'",
      card(card_header("Kontrasty variant — podmíněný model"),
           DTOutput("posthoc_results"))
    ),
    card(full_screen = TRUE,
        card_header("Podmíněný model — Dunnett forest plot (proti kontrole)"),
        p(class = "text-muted", "Kontrasty kazde varianty proti kontrole s Dunnettovou korekci."),
        plotOutput("dunnett_cond_contrast_plot", height = "650px")
    ),
    conditionalPanel(
      condition = "input.contrast_mode == 'dunnett'",
      card(card_header("Kontrasty proti kontrole (Dunnett) — podmíněný model"),
           DTOutput("dunnett_cond_results"))
    ),

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

  standardize_interval_cols <- function(df) {
    low_name <- grep("(^lower\\.CL$|LCL$)", names(df), value = TRUE)[1]
    up_name  <- grep("(^upper\\.CL$|UCL$)", names(df), value = TRUE)[1]
    if (!is.na(low_name) && nzchar(low_name) && low_name != "lower.CL") {
      names(df)[names(df) == low_name] <- "lower.CL"
    }
    if (!is.na(up_name) && nzchar(up_name) && up_name != "upper.CL") {
      names(df)[names(df) == up_name] <- "upper.CL"
    }
    df
  }
  
  get_inference_emm <- function(emm_obj, scale = "link") {
    if (identical(scale, "response")) {
      regrid(emm_obj, transform = "response")
    } else {
      emm_obj
    }
  }

  make_contrast_table <- function(emm_obj, adjust_method, conf_level, scale = "link") {
    obj <- get_inference_emm(emm_obj, scale = scale)
    pw <- pairs(obj, adjust = adjust_method)
    out <- as.data.frame(summary(
      pw,
      infer = c(TRUE, TRUE),
      level = conf_level,
      adjust = adjust_method
    ))
    standardize_interval_cols(out)
  }

  make_dunnett_table <- function(emm_obj, control_name, conf_level, scale = "link") {
    obj <- get_inference_emm(emm_obj, scale = scale)
    levs <- levels(obj)$treatment
    if (is.null(levs)) levs <- levels(obj)[[1]]
    ref_idx <- which(levs == control_name)
    if (length(ref_idx) == 0) ref_idx <- 1
    pw <- contrast(obj, method = "trt.vs.ctrl", ref = ref_idx, adjust = "dunnettx")
    out <- as.data.frame(summary(
      pw,
      infer = c(TRUE, TRUE),
      level = conf_level,
      adjust = "dunnettx"
    ))
    standardize_interval_cols(out)
  }

  # Dunnett letter assignment: a = same as control, b = better, c = worse
  make_dunnett_letters <- function(dunnett_df, emm_df, control_name, alpha) {
    if (is.null(dunnett_df) || nrow(dunnett_df) == 0) return(NULL)

    # Get control estimate on response scale
    ctrl_row <- emm_df[emm_df$treatment == control_name, , drop = FALSE]
    if (nrow(ctrl_row) == 0) return(NULL)
    ctrl_prob <- ctrl_row$prob[1]

    has_time <- "time" %in% names(dunnett_df)
    
    # Parse contrast names to get treatment names
    results <- list()
    
    if (has_time) {
      times <- unique(dunnett_df$time)
      for (tt in times) {
        sub_d <- dunnett_df[dunnett_df$time == tt, , drop = FALSE]
        ctrl_prob_t <- emm_df$prob[emm_df$treatment == control_name & emm_df$time == tt][1]
        
        letter_rows <- data.frame(
          treatment = control_name,
          time = tt,
          .dunnett_letter = "a",
          stringsAsFactors = FALSE
        )
        
        for (i in seq_len(nrow(sub_d))) {
          cname <- as.character(sub_d$contrast[i])
          # extract treatment name from contrast like "TrtX - Control"
          parts <- strsplit(cname, " - ", fixed = TRUE)[[1]]
          trt_name <- trimws(parts[1])
          # check if this is "Control - TrtX" format instead
          if (trt_name == control_name && length(parts) >= 2) {
            trt_name <- trimws(parts[2])
          }
          
          p_val <- sub_d$p.value[i]
          trt_prob <- emm_df$prob[emm_df$treatment == trt_name & emm_df$time == tt]
          if (length(trt_prob) == 0) trt_prob <- NA
          trt_prob <- trt_prob[1]
          
          if (is.na(p_val) || p_val >= alpha) {
            ltr <- "a"
          } else if (!is.na(trt_prob) && trt_prob > ctrl_prob_t) {
            ltr <- "b"
          } else {
            ltr <- "c"
          }
          
          letter_rows <- rbind(letter_rows, data.frame(
            treatment = trt_name, time = tt,
            .dunnett_letter = ltr, stringsAsFactors = FALSE
          ))
        }
        results[[length(results) + 1]] <- letter_rows
      }
      do.call(rbind, results)
    } else {
      letter_rows <- data.frame(
        treatment = control_name,
        .dunnett_letter = "a",
        stringsAsFactors = FALSE
      )
      
      for (i in seq_len(nrow(dunnett_df))) {
        cname <- as.character(dunnett_df$contrast[i])
        parts <- strsplit(cname, " - ", fixed = TRUE)[[1]]
        trt_name <- trimws(parts[1])
        if (trt_name == control_name && length(parts) >= 2) {
          trt_name <- trimws(parts[2])
        }
        
        p_val <- dunnett_df$p.value[i]
        trt_prob <- emm_df$prob[emm_df$treatment == trt_name]
        if (length(trt_prob) == 0) trt_prob <- NA
        trt_prob <- trt_prob[1]
        
        if (is.na(p_val) || p_val >= alpha) {
          ltr <- "a"
        } else if (!is.na(trt_prob) && trt_prob > ctrl_prob) {
          ltr <- "b"
        } else {
          ltr <- "c"
        }
        
        letter_rows <- rbind(letter_rows, data.frame(
          treatment = trt_name,
          .dunnett_letter = ltr, stringsAsFactors = FALSE
        ))
      }
      letter_rows
    }
  }

  make_cld_table <- function(emm_obj, adjust_method, alpha, scale = "link") {
    obj <- get_inference_emm(emm_obj, scale = scale)
    cld_res <- multcomp::cld(
      obj,
      adjust = adjust_method,
      Letters = letters,
      sort = FALSE,
      alpha = alpha
    )
    out <- as.data.frame(cld_res)
    if (".group" %in% names(out)) out$.group <- trimws(out$.group)
    out
  }

	make_group_band_table <- function(emm_obj, conf_level, scale = "link") {
	  if (identical(scale, "response")) {
	    obj <- regrid(emm_obj, transform = "response")
	    out <- as.data.frame(confint(obj, adjust = "mvt", level = conf_level))
	  } else {
	    out <- as.data.frame(confint(
	      emm_obj,
	      adjust = "mvt",
	      level = conf_level,
	      type = "response"
	    ))
	  }

	  out <- standardize_interval_cols(out)
	  names(out)[names(out) == "lower.CL"] <- "band.LCL"
	  names(out)[names(out) == "upper.CL"] <- "band.UCL"
	  out
	}

  extract_contrast_pairs <- function(df) {
    if (!"contrast" %in% names(df) || nrow(df) == 0) return(df)
    pair_mat <- t(vapply(
      as.character(df$contrast),
      function(x) {
        parts <- strsplit(x, " - ", fixed = TRUE)[[1]]
        
        part1 <- if (length(parts) >= 1) trimws(parts[1]) else NA_character_
        part2 <- if (length(parts) >= 2) trimws(parts[2]) else NA_character_
        c(part1, part2)
      },
      character(2)
    ))
    df$treatment1 <- pair_mat[, 1]
    df$treatment2 <- pair_mat[, 2]
    df
  }

  contrast_heatmap_data <- function(df, treatment_order, alpha) {
    if (is.null(df) || nrow(df) == 0 || !"p.value" %in% names(df)) return(NULL)
    x <- extract_contrast_pairs(df)
    if (!all(c("treatment1", "treatment2") %in% names(x))) return(NULL)

    keep_cols <- intersect(c("time", "treatment1", "treatment2", "p.value"), names(x))
    base <- x[, keep_cols, drop = FALSE]
    base$status <- ifelse(base$p.value < alpha, "rozdíl", "bez rozdílu")

    swap <- base
    tmp <- swap$treatment1
    swap$treatment1 <- swap$treatment2
    swap$treatment2 <- tmp

    if ("time" %in% names(base)) {
      diag_df <- expand.grid(
        treatment1 = treatment_order,
        treatment2 = treatment_order,
        time = unique(as.character(base$time)),
        stringsAsFactors = FALSE
      )
    } else {
      diag_df <- expand.grid(
        treatment1 = treatment_order,
        treatment2 = treatment_order,
        stringsAsFactors = FALSE
      )
    }
    diag_df$status <- ifelse(diag_df$treatment1 == diag_df$treatment2, "stejná varianta", NA_character_)

    out <- bind_rows(base, swap, diag_df) %>%
      filter(!is.na(status)) %>%
      distinct()

    out$treatment1 <- factor(out$treatment1, levels = rev(treatment_order))
    out$treatment2 <- factor(out$treatment2, levels = treatment_order)
    out$status <- factor(out$status, levels = c("rozdíl", "bez rozdílu", "stejná varianta"))
    if ("time" %in% names(out)) {
      out$time <- factor(out$time, levels = unique(as.character(out$time)))
    }
    out
  }

  make_contrast_plot <- function(df, title_txt, subtitle_txt) {
    req(df)
    plot_df <- standardize_interval_cols(df)
    validate(
      need("contrast" %in% names(plot_df), "Chybí sloupec contrast."),
      need("estimate" %in% names(plot_df), "Chybí sloupec estimate."),
      need("lower.CL" %in% names(plot_df) && "upper.CL" %in% names(plot_df), "Chybí intervaly kontrastů.")
    )

    if ("p.value" %in% names(plot_df)) {
      plot_df$sig <- ifelse(plot_df$p.value < input$alpha, "rozdíl", "bez rozdílu")
    } else {
      plot_df$sig <- "kontrast"
    }
    
    plot_df$contrast <- factor(plot_df$contrast, levels = rev(unique(as.character(plot_df$contrast))))

    x_lab <- if (identical(input$infer_scale, "response")) "Rozdíl odhadů (response scale)" else "Rozdíl odhadů (model scale)"

    g <- ggplot(plot_df, aes(y = contrast, x = estimate, colour = sig)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
      geom_segment(aes(x = lower.CL, xend = upper.CL, yend = contrast), linewidth = 0.8) +
      geom_point(size = 2.4) +
      scale_colour_manual(
        values = c("rozdíl" = "#c0392b", "bez rozdílu" = "#1b7837", "kontrast" = "#2c3e50"),
        name = NULL
      ) +
      labs(
        x = x_lab,
        y = NULL,
        title = title_txt,
        subtitle = subtitle_txt
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(colour = "grey45"),
        panel.grid.major.y = element_blank()
      )

    if ("time" %in% names(plot_df)) {
      g <- g + facet_wrap(~ time, scales = "free_y")
    }
    g
  }

  make_sig_heatmap_plot <- function(df, title_txt, subtitle_txt) {
    heat_df <- contrast_heatmap_data(df, rv$treatment_order, input$alpha)
    validate(need(!is.null(heat_df) && nrow(heat_df) > 0, "Není k dispozici matice signifikance."))

    g <- ggplot(heat_df, aes(x = treatment2, y = treatment1, fill = status)) +
      geom_tile(colour = "white", linewidth = 1) +
      geom_text(aes(label = ifelse(status == "rozdíl", "×",
                                   ifelse(status == "bez rozdílu", "ns", "—"))),
                size = 4.8, fontface = "bold") +
      scale_fill_manual(
        values = c("rozdíl" = "#d73027", "bez rozdílu" = "#1a9850", "stejná varianta" = "#d9d9d9"),
        name = NULL,
        drop = FALSE
      ) +
      labs(
        x = "Varianta 2",
        y = "Varianta 1",
        title = title_txt,
        subtitle = subtitle_txt
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(colour = "grey45"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold")
      )

    if ("time" %in% names(heat_df)) {
      g <- g + facet_wrap(~ time)
    }
    g
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
    emm_marginal = NULL, cum_emm_marginal = NULL,
    aic_bic = NULL,
    log = "",
    emm_group_band = NULL,
    final_group_band = NULL,
    cum_group_band = NULL,
    dunnett_posthoc = NULL,
    dunnett_final_posthoc = NULL,
    dunnett_cum_posthoc = NULL
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
    rv$emm_group_band <- NULL
    rv$final_group_band <- NULL
    rv$cum_group_band <- NULL
    rv$dunnett_posthoc <- NULL
    rv$dunnett_final_posthoc <- NULL
    rv$dunnett_cum_posthoc <- NULL
  }

  uploaded <- reactiveVal(NULL)

  observeEvent(input$file, {
    reset_all()
    df <- read_upload(input$file$datapath, input$file$name)
    add_log("Loaded: ", input$file$name, " (", nrow(df), "\u00D7", ncol(df), ")")
    uploaded(df)

    cols <- names(df)
    updateSelectInput(session, "col_variant", choices = cols, selected = cols[1])
    
    sel_date <- if (length(cols) >= 2) cols[2] else cols[1]
    updateSelectInput(session, "col_date", choices = cols, selected = sel_date)
    
    sel_count <- if (length(cols) >= 3) cols[3] else cols[1]
    updateSelectInput(session, "col_count", choices = cols, selected = sel_count)
  })

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
      updateSelectInput(session, "control_variant", choices = treat_order, selected = treat_order[1])
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
    
    prob_txt <- if (nrow(rv$mono_issues) == 0) "\u2713 Žádné" else paste0("\u26A0 ", nrow(rv$mono_issues))

    card(class = "bg-light",
      p(tags$strong("Varianty: "), n_distinct(wide$treatment)),
      p(tags$strong("Nádoby: "), nrow(wide)),
      p(tags$strong("Data: "), paste(dcols, collapse = ", ")),
      p(tags$strong("Problémy: "), prob_txt)
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
    tagList(tables)
  })

  output$wide_table <- renderDT({
    req(rv$wide_data)
    datatable(rv$wide_data, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
  })

  output$mono_table <- renderDT({
    req(rv$mono_issues, nrow(rv$mono_issues) > 0)
    datatable(rv$mono_issues, options = list(scrollX = TRUE), rownames = FALSE)
  })

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

      long <- long %>% mutate(new_germ = pmax(0, new_germ))

      if (isTRUE(input$drop_zero_times)) {
        zero_times <- long %>%
          group_by(time) %>%
          summarise(all_zero = all(cumulative == 0), .groups = "drop") %>%
          filter(all_zero) %>%
          pull(time)
        if (length(zero_times) > 0) {
          long <- long %>% filter(!time %in% zero_times)
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

      rv$full_data <- long

      n_zero <- sum(long$remaining <= 0, na.rm = TRUE)
      if (n_zero > 0) add_log("Vyřazeno ", n_zero, " pozorování s remaining \u2264 0 pro model (všechna semínka vzešla)")
      long <- long %>% filter(remaining > 0)

      long <- long %>% mutate(new_germ = pmin(new_germ, remaining))

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
      card(card_header("Data s přírůstky"), DTOutput("inc_table")),
      card(card_header("Trajektorie jednotlivých nádob"), plotOutput("prep_plot", height = "380px"))
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
      stat_summary(aes(group = treatment), fun = mean, geom = "line", linewidth = 1.4) +
      stat_summary(aes(group = treatment), fun = mean, geom = "point", size = 3) +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
      labs(x = "Datum", y = "Kumulativní vzcházivost (%)", colour = "Varianta") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  })

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
      rv$M1 <- glm(cbind(new_germ, remaining - new_germ) ~ time, family = fam, data = dat)
      rv$M2 <- glm(cbind(new_germ, remaining - new_germ) ~ treatment + time, family = fam, data = dat)
      rv$M3 <- glm(cbind(new_germ, remaining - new_germ) ~ treatment * time, family = fam, data = dat)

      add_log("M1 dev=", round(rv$M1$deviance, 1),
              "  M2 dev=", round(rv$M2$deviance, 1),
              "  M3 dev=", round(rv$M3$deviance, 1))

      if (fam_name == "binomial") {
        rv$aic_bic <- data.frame(
          Model = c("M1 (~time)", "M2 (~treatment + time)", "M3 (~treatment * time)"),
          AIC = c(AIC(rv$M1), AIC(rv$M2), AIC(rv$M3)),
          BIC = c(BIC(rv$M1), BIC(rv$M2), BIC(rv$M3)),
          stringsAsFactors = FALSE
        )
      } else {
        B1 <- glm(cbind(new_germ, remaining - new_germ) ~ time, family = binomial, data = dat)
        B2 <- glm(cbind(new_germ, remaining - new_germ) ~ treatment + time, family = binomial, data = dat)
        B3 <- glm(cbind(new_germ, remaining - new_germ) ~ treatment * time, family = binomial, data = dat)
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

      p_col <- if (test_type == "F") "Pr(>F)" else "Pr(>Chi)"
      p_treat <- rv$lr12[[p_col]][2]
      p_inter <- rv$lr23[[p_col]][2]
      add_log("Treatment: p=", format.pval(p_treat, 4),
              "  Interaction: p=", format.pval(p_inter, 4))

      alpha <- input$alpha
      choice <- input$model_choice
      if (choice == "M3") {
        rv$best_model <- "M3"
        add_log("Model zvolen ručně: M3 (s interakcí)")
      } else if (choice == "M2") {
        rv$best_model <- "M2"
        add_log("Model zvolen ručně: M2 (bez interakce)")
      } else {
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
	  emm_main <- emmeans(mod, ~ treatment | time)
	  rv$posthoc <- make_contrast_table(emm_main, input$p_method, conf_level, input$infer_scale)
	  rv$emm_plot <- as.data.frame(summary(emm_main, type = "response", level = conf_level))
	  rv$emm_group_band <- make_group_band_table(emm_main, conf_level, input$infer_scale)
	} else {
	  emm_main <- emmeans(mod, ~ treatment)
	  rv$posthoc <- make_contrast_table(emm_main, input$p_method, conf_level, input$infer_scale)
	  emm_time <- emmeans(mod, ~ treatment + time)
	  rv$emm_plot <- as.data.frame(summary(emm_time, type = "response", level = conf_level))
	  rv$emm_group_band <- make_group_band_table(emm_main, conf_level, input$infer_scale)
	}

      rv$dunnett_posthoc <- NULL
      tryCatch({
        ctrl <- input$control_variant
        if (!is.null(ctrl) && ctrl %in% rv$treatment_order) {
          rv$dunnett_posthoc <- make_dunnett_table(emm_main, ctrl, conf_level, input$infer_scale)
          add_log("Dunnett contrasts (podminenych) vs ", ctrl)
        }
      }, error = function(e) {
        add_log("Dunnett (podmineny) selhal: ", conditionMessage(e))
      })

      rv$cld_df <- NULL
      tryCatch({
        if (rv$best_model == "M3") {
          emm_cld <- emmeans(mod, ~ treatment | time)
          rv$cld_df <- make_cld_table(emm_cld, input$p_method, input$alpha, input$infer_scale)
          add_log("\u2713 CLD computed for ", length(unique(rv$cld_df$time)),
                  " timepoints on ", input$infer_scale, " scale")
        } else {
          emm_cld <- emmeans(mod, ~ treatment)
          rv$cld_df <- make_cld_table(emm_cld, input$p_method, input$alpha, input$infer_scale)
          add_log("\u2713 CLD computed on ", input$infer_scale, " scale")
        }
      }, error = function(e) {
        add_log("\u26A0 CLD error: ", conditionMessage(e))
        tryCatch({
          if (rv$best_model == "M3") {
            emm_cld <- emmeans(mod, ~ treatment | time)
            obj <- get_inference_emm(emm_cld, input$infer_scale)
            cld_res <- multcomp::cld(obj, adjust = input$p_method, Letters = letters, alpha = input$alpha)
            cld_out <- as.data.frame(cld_res)
            if (".group" %in% names(cld_out)) cld_out$.group <- trimws(cld_out$.group)
            rv$cld_df <- cld_out
          } else {
            emm_cld <- emmeans(mod, ~ treatment)
            obj <- get_inference_emm(emm_cld, input$infer_scale)
            cld_res <- multcomp::cld(obj, adjust = input$p_method, Letters = letters, alpha = input$alpha)
            cld_out <- as.data.frame(cld_res)
            if (".group" %in% names(cld_out)) cld_out$.group <- trimws(cld_out$.group)
            rv$cld_df <- cld_out
          }
          add_log("\u2713 CLD computed (fallback) on ", input$infer_scale, " scale")
        }, error = function(e2) {
          add_log("\u26A0 CLD fallback also failed: ", conditionMessage(e2))
        })
      })
      add_log("\u2713 Temporal analysis complete")

      tryCatch({
        df_emm <- rv$emm_plot
        if ("time" %in% names(df_emm)) {
          rv$emm_marginal <- df_emm %>%
            filter(!is.na(prob)) %>%
            group_by(treatment) %>%
            summarise(
              prob = mean(prob),
              SE = sqrt(mean(SE^2)),
              asymp.LCL = mean(asymp.LCL),
              asymp.UCL = mean(asymp.UCL),
              n_times = n(),
              .groups = "drop"
            )
        } else {
          rv$emm_marginal <- df_emm
        }
        add_log("\u2713 Marginální odhad (podmíněný model) computed")
      }, error = function(e) {
        add_log("\u26A0 Marginální odhad selhal: ", conditionMessage(e))
        rv$emm_marginal <- NULL
      })

      add_log("Running final germination test...")
      tryCatch({
        full <- rv$full_data
        last_time <- levels(full$time)[length(levels(full$time))]
        final_dat <- full[full$time == last_time, ]
        final_dat$treatment <- factor(final_dat$treatment, levels = rv$treatment_order)
        N <- input$N_seeds

        rv$final_M0 <- glm(cbind(cumulative, N - cumulative) ~ 1, family = fam, data = final_dat)
        rv$final_M1 <- glm(cbind(cumulative, N - cumulative) ~ treatment, family = fam, data = final_dat)

        rv$final_ftest <- anova(rv$final_M0, rv$final_M1, test = test_type)
        p_col_final <- if (test_type == "F") "Pr(>F)" else "Pr(>Chi)"
        p_final <- rv$final_ftest[[p_col_final]][2]
        add_log("Final germination ", test_type, "-test: p = ", format.pval(p_final, 4))

        tryCatch({
          emm_final <- emmeans(rv$final_M1, ~ treatment)
          rv$final_group_band <- make_group_band_table(
  		emm_final, conf_level, input$infer_scale
	  )
	  rv$final_posthoc <- make_contrast_table(emm_final, input$p_method, conf_level, input$infer_scale)
          rv$final_cld <- make_cld_table(emm_final, input$p_method, input$alpha, input$infer_scale)
          add_log("\u2713 Final CLD + post-hoc computed on ", input$infer_scale, " scale")
          rv$dunnett_final_posthoc <- NULL
          tryCatch({
            ctrl <- input$control_variant
            if (!is.null(ctrl) && ctrl %in% rv$treatment_order) {
              rv$dunnett_final_posthoc <- make_dunnett_table(emm_final, ctrl, conf_level, input$infer_scale)
              add_log("Dunnett contrasts (jednoduchy) vs ", ctrl)
            }
          }, error = function(e2) {
            add_log("Dunnett (jednoduchy) selhal: ", conditionMessage(e2))
          })
        }, error = function(e) {
          add_log("\u26A0 Final CLD/post-hoc failed: ", conditionMessage(e))
          rv$final_cld <- NULL
          rv$final_posthoc <- NULL
        })

        add_log("\u2713 Final germination analysis complete")
      }, error = function(e) {
        add_log("\u26A0 Final germination error: ", conditionMessage(e))
      })

      add_log("Orientační kumulativní model...")
      tryCatch({
        full <- rv$full_data
        full$treatment <- factor(full$treatment, levels = rv$treatment_order)
        N <- input$N_seeds

        rv$cum_model <- glm(cbind(cumulative, N - cumulative) ~ treatment * time, family = fam, data = full)

        cum_emm <- emmeans(rv$cum_model, ~ treatment | time)
        rv$cum_group_band <- make_group_band_table(
 		 cum_emm, conf_level, input$infer_scale
	)
	rv$cum_emm <- as.data.frame(summary(cum_emm, type = "response", level = conf_level))
        rv$cum_posthoc <- make_contrast_table(cum_emm, input$p_method, conf_level, input$infer_scale)

        rv$dunnett_cum_posthoc <- NULL
        tryCatch({
          ctrl <- input$control_variant
          if (!is.null(ctrl) && ctrl %in% rv$treatment_order) {
            rv$dunnett_cum_posthoc <- make_dunnett_table(cum_emm, ctrl, conf_level, input$infer_scale)
            add_log("Dunnett contrasts (kumulativni) vs ", ctrl)
          }
        }, error = function(e) {
          add_log("Dunnett (kumulativni) selhal: ", conditionMessage(e))
        })

        tryCatch({
          rv$cum_cld <- make_cld_table(cum_emm, input$p_method, input$alpha, input$infer_scale)
          add_log("\u2713 Kumulativní model: emmeans + CLD + post-hoc hotovo na ", input$infer_scale, " scale")
        }, error = function(e) {
          add_log("\u26A0 Kumulativní CLD selhalo: ", conditionMessage(e))
          rv$cum_cld <- NULL
        })

        tryCatch({
          df_cum <- rv$cum_emm
          rv$cum_emm_marginal <- df_cum %>%
            filter(!is.na(prob)) %>%
            group_by(treatment) %>%
            summarise(
              prob = mean(prob),
              SE = sqrt(mean(SE^2)),
              asymp.LCL = mean(asymp.LCL),
              asymp.UCL = mean(asymp.UCL),
              n_times = n(),
              .groups = "drop"
            )
          add_log("\u2713 Marginální odhad (kumulativní model) computed")
        }, error = function(e) {
          add_log("\u26A0 Kumulativní marginální odhad selhal: ", conditionMessage(e))
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
    
    label_txt <- if (rv$best_model == "M3") "M3 (s interakcí)" else "M2 (bez interakce)"
    div_class <- if (manual) "alert alert-warning" else "alert alert-info"
    manual_txt <- if (manual) " (ruční volba)" else " (automaticky dle testu)"
    manual_span <- if (manual) tags$span(class = "text-muted", " — Ruční volba přepisuje automatický výběr.") else NULL

    div(class = div_class,
      tags$strong("Zvolený model: "), label_txt, manual_txt,
      tags$br(),
      tags$small(
        sprintf("AIC preferuje: %s | BIC preferuje: %s", aic_best, bic_best),
        manual_span
      )
    )
  })

  make_ftest_summary <- function(show_simple = FALSE) {
    req(rv$lr12, rv$lr23, rv$best_model)
    alpha <- input$alpha

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

    aic_line <- NULL
    if (!is.null(rv$aic_bic)) {
      aic_best <- rv$aic_bic$Model[which.min(rv$aic_bic$AIC)]
      bic_best <- rv$aic_bic$Model[which.min(rv$aic_bic$BIC)]
      aic_line <- tags$p(class = "text-muted small",
        sprintf("AIC preferuje %s, BIC preferuje %s.", aic_best, bic_best))
    }

    simple_line <- NULL
    if (show_simple && !is.null(rv$final_ftest)) {
      p_final <- get_p(rv$final_ftest)
      final_sig <- !is.na(p_final) && p_final < alpha
      
      final_span <- if (final_sig) {
        tags$span(style = "color: #d32f2f;", sprintf("\u2714 Varianty se liší v celkové vzcházivosti (p = %s)", format.pval(p_final, digits = 3)))
      } else {
        tags$span(style = "color: #666;", sprintf("\u2718 Nebyl detekován rozdíl v celkové vzcházivosti (p = %s)", format.pval(p_final, digits = 3)))
      }
      
      simple_line <- tags$p(
        tags$strong("4. Jednoduchý test (jen poslední datum): "),
        final_span
      )
    }

    treat_span <- if (treat_sig) {
      tags$span(style = "color: #d32f2f;", sprintf("\u2714 ANO (p = %s)", format.pval(p_treat, digits = 3)))
    } else {
      tags$span(style = "color: #666;", sprintf("\u2718 Statisticky nevýznamný rozdíl (p = %s)", format.pval(p_treat, digits = 3)))
    }
    
    inter_span <- if (inter_sig) {
      tags$span(style = "color: #d32f2f;", sprintf("\u2714 ANO (p = %s)", format.pval(p_inter, digits = 3)))
    } else {
      tags$span(style = "color: #666;", sprintf("\u2718 Ne — efekt konzistentní napříč daty (p = %s)", format.pval(p_inter, digits = 3)))
    }
    
    model_txt <- if (using_m3) "M3 (s interakcí) — varianty porovnány v každém čase zvlášť." else "M2 (bez interakce) — varianty porovnány celkově."
    manual_span <- if (manual) tags$span(style = "color: #e65100;", " (ruční volba)") else NULL

    card(
      card_header("Shrnutí výsledků"),
      tags$div(style = "font-size: 1.05rem; line-height: 1.8;",
        tags$p(tags$strong("Časový model"), " — analyzuje průběh klíčení v čase (kondicionální přírůstky):",
               style = "margin-bottom: 2px; color: #555;"),
        tags$p(
          tags$strong("1. Existuje jednotný efekt varianty v průběhu klíčení? "),
          treat_span,
          tags$br(),
          tags$small(class = "text-muted",
            "Testuje, zda mají jednotlivé varianty jednotný efekt v čase - při silném výsledku interakce může vyjít nevýznamné, i když celkově varianta data silně ovlivňuje.")
        ),
        tags$p(
          tags$strong("2. Mění se efekt variant v čase? "),
          inter_span,
          tags$br(),
          tags$small(class = "text-muted",
            "Testuje, zda se rozdíly mezi variantami mění v různých datech (např. jedna začne brzy, jiná pozdě). Celkově tedy interakce čas:varianta. Vyjde-li interakce jako významná, hůře se interpretuje test předchozí.")
        ),
        tags$p(
          tags$strong("3. Použitý model: "),
          model_txt,
          manual_span
        ),
        aic_line,
        if (show_simple) tags$hr(style = "margin: 8px 0;") else NULL,
        if (show_simple) tags$p(tags$strong("Jednoduchý test"), " — ignoruje čas, kouká jen na konečný výsledek:", style = "margin-bottom: 2px; color: #555;") else NULL,
        simple_line,
        if (show_simple && !is.null(simple_line)) tags$small(class = "text-muted", "Testuje, zda se varianty liší v celkovém počtu vzešlých semenáčků na konci experimentu. Může se lišit od časového modelu — varianty mohou klíčit různě rychle, ale skončit na podobném výsledku.") else NULL
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

    div_class <- if (is_quasi) "alert alert-warning" else "alert alert-danger"
    warn_span <- if (!is_quasi) tags$span(" — Binomický model ignoruje tuto overdisperzi! P-hodnoty mohou být příliš optimistické. Je k dispozici kvazibinomický model.") else tags$span(" — Kvazibinomický model toto zohledňuje v SE a p-hodnotách.")

    if (phi > 1.5) {
      div(class = div_class,
        tags$strong("Overdisperze: "),
        sprintf("\u03C6 = %.2f (SE nafouklé \u00D7%.2f)", phi, sqrt(phi)),
        warn_span
      )
    } else if (phi < 0.5) {
      div(class = "alert alert-info",
        tags$strong("Underdisperze: "), sprintf("\u03C6 = %.2f", phi))
    } else {
      div(class = "alert alert-success",
        tags$strong("Disperze: "), sprintf("\u03C6 = %.2f — blízká 1, je přiměřený i binomický model.", phi))
    }
  })

  output$emm_info <- renderUI({
    req(rv$best_model)
    txt <- if (rv$best_model == "M3") {
      "Odhady podmíněné pravděpodobnosti klíčení pro každou variantu v každém datu (model s interakcí)."
    } else {
      "Odhady podmíněné pravděpodobnosti klíčení pro každou variantu (průměr přes data, model bez interakce)."
    }
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
    is_dunnett <- identical(input$contrast_mode, "dunnett")
    
    txt <- if (rv$best_model == "M3") {
      "Interakce mezi datem a variantou byla uznána za důležitou; tabulka níže ukazuje kontrasty mezi variantami v každém datu."
    } else {
      "Interakce nezohledněna; tabulka níže ukazuje celkové kontrasty mezi variantami."
    }
    
    infer_txt <- if (input$infer_scale == "link") {
      "Inference běží na model scale (doporučeno)."
    } else {
      "Inference běží na response scale (experimentální)."
    }
    
    method_txt <- if (is_dunnett) {
      sprintf("Dunnett (kontrasty proti kontrole: %s).", input$control_variant)
    } else {
      paste0("Všechny páry (s korekcí metodou ", input$p_method, ").")
    }

    tagList(
      div(class = "alert alert-info",
          tags$strong("Přístup: "), txt, " ",
          tags$strong("Metoda: "), method_txt
      ),
      tags$br(),
      tags$small(infer_txt)
    )
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

  output$cum_model_summary <- renderUI({
    req(rv$cum_model)
    phi <- rv$cum_model$deviance / rv$cum_model$df.residual
    is_quasi <- input$family_choice == "quasibinomial"
    
    se_txt <- if (is_quasi) sprintf(", SE inflace = \u00D7%.2f", sqrt(phi)) else ""

    div(class = "alert alert-secondary",
      tags$strong("Kumulativní model (orientační): "),
      sprintf("treatment * time, \u03C6 = %.2f", phi),
      se_txt,
      tags$br(),
      tags$small(class = "text-muted", "Pozorování z téhož květináče jsou korelovaná — výsledky jsou pouze orientační.")
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

  output$final_summary <- renderUI({
    req(rv$final_ftest)
    alpha <- input$alpha
    p_cols <- grep("^Pr\\(", names(rv$final_ftest), value = TRUE)
    p_val <- if (length(p_cols) > 0) rv$final_ftest[[p_cols[1]]][2] else NA_real_
    sig <- !is.na(p_val) && p_val < alpha
    phi <- rv$final_M1$deviance / rv$final_M1$df.residual
    is_quasi <- input$family_choice == "quasibinomial"
    test_label <- if (is_quasi) "F-test" else "LR test"

    sig_span <- if (sig) {
      tags$span(style = "color: #d32f2f;", sprintf("\u2714 ANO (%s p = %s)", test_label, format.pval(p_val, digits = 3)))
    } else {
      tags$span(style = "color: #666;", sprintf("\u2718 Nebyl detekován rozdíl (%s p = %s)", test_label, format.pval(p_val, digits = 3)))
    }
    
    disp_p <- if (is_quasi) {
      tags$p(style = "font-size: 0.9rem; color: #666;", sprintf("Disperze \u03C6 = %.2f, SE inflace = \u00D7%.2f", phi, sqrt(phi)))
    } else {
      tags$p(style = "font-size: 0.9rem; color: #666;", sprintf("Binomický model (bez korekce na overdisperzi). Disperze \u03C6 = %.2f", phi))
    }

    card(
      card_header("Celková vzcházivost bez časového rozložení"),
      tags$div(style = "font-size: 1.1rem; line-height: 2;",
        tags$p(
          tags$strong("Liší se varianty v celkové vzcházivosti? "),
          sig_span
        ),
        disp_p
      )
    )
  })

  output$final_overdisp_summary <- renderUI({
    req(rv$final_M1)
    phi <- rv$final_M1$deviance / rv$final_M1$df.residual
    is_quasi <- input$family_choice == "quasibinomial"

    div_class <- if (is_quasi) "alert alert-warning" else "alert alert-danger"
    warn_span <- if (!is_quasi) tags$span(" — Binomický model ignoruje tuto overdisperzi! P-hodnoty mohou být příliš optimistické. Je k dispozici kvazibinomický model.") else tags$span(" — Kvazibinomický model toto zohledňuje v SE a p-hodnotách.")

    if (phi > 1.5) {
      div(class = div_class,
        tags$strong("Overdisperze: "),
        sprintf("\u03C6 = %.2f (SE nafouklé \u00D7%.2f)", phi, sqrt(phi)),
        warn_span
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
	  req(rv$full_data, rv$final_M1)
	  
	  is_dunnett <- identical(input$contrast_mode, "dunnett")
	  
	  # In allpairs mode, require CLD; in dunnett mode, require dunnett results
	  if (is_dunnett) {
	    req(rv$dunnett_final_posthoc)
	  } else {
	    req(rv$final_cld)
	  }

	  dat <- rv$full_data
	  band_df <- rv$final_group_band
	  ci_pct <- round((1 - input$alpha) * 100)

	  # always get displayed means on response scale
	  emm_disp <- as.data.frame(
	    summary(
	      emmeans(rv$final_M1, ~ treatment),
	      type = "response",
	      level = 1 - input$alpha
	    )
	  )

	  emm_disp <- standardize_interval_cols(emm_disp)

	  plot_data <- data.frame(
	    treatment = emm_disp$treatment,
	    mean_pct = emm_disp$prob * 100,
	    ci_lo = emm_disp$lower.CL * 100,
	    ci_hi = emm_disp$upper.CL * 100,
	    stringsAsFactors = FALSE
	  )

	  # add letters: CLD or Dunnett a/b/c
	  if (is_dunnett) {
	    dlett <- make_dunnett_letters(
	      rv$dunnett_final_posthoc, emm_disp,
	      input$control_variant, input$alpha
	    )
	    if (!is.null(dlett)) {
	      letter_df <- data.frame(
	        treatment = dlett$treatment,
	        letter = dlett$.dunnett_letter,
	        stringsAsFactors = FALSE
	      )
	    } else {
	      letter_df <- data.frame(treatment = plot_data$treatment, letter = "", stringsAsFactors = FALSE)
	    }
	  } else {
	    cld_df <- rv$final_cld
	    letter_df <- data.frame(
	      treatment = cld_df$treatment,
	      letter = trimws(cld_df$.group),
	      stringsAsFactors = FALSE
	    )
	  }

	  plot_data <- dplyr::left_join(plot_data, letter_df, by = "treatment")

	  # add grouping bands if available and not in dunnett mode
	  if (!is_dunnett && !is.null(band_df)) {
	    band_keep <- band_df[, intersect(c("treatment", "band.LCL", "band.UCL"), names(band_df)), drop = FALSE]

	    if (all(c("treatment", "band.LCL", "band.UCL") %in% names(band_keep))) {
	      plot_data <- dplyr::left_join(plot_data, band_keep, by = "treatment")
	      plot_data$band.LCL <- plot_data$band.LCL * 100
	      plot_data$band.UCL <- plot_data$band.UCL * 100
	    }
	  }

	  plot_data$treatment <- factor(plot_data$treatment, levels = rv$treatment_order)
	  last_time <- levels(dat$time)[length(levels(dat$time))]

	  plot_data$label_y <- plot_data$ci_hi + 2

	if (!is_dunnett && isTRUE(input$show_group_bands) &&
	    all(c("band.LCL", "band.UCL") %in% names(plot_data))) {
	  plot_data$label_y <- pmax(plot_data$ci_hi, plot_data$band.UCL, na.rm = TRUE) + 2
	}

	# Colour letters for Dunnett mode
	if (is_dunnett) {
	  plot_data$letter_colour <- ifelse(
	    plot_data$letter == "a", "#1b7837",
	    ifelse(plot_data$letter == "b", "#2166ac", "#c0392b")
	  )
	} else {
	  plot_data$letter_colour <- "black"
	}

	g <- ggplot(plot_data, aes(x = treatment, y = mean_pct, fill = treatment)) +
	  geom_col(alpha = 0.85, width = 0.7) +
	  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.25, linewidth = 0.8)

	if (!is_dunnett && isTRUE(input$show_group_bands) &&
	    all(c("band.LCL", "band.UCL") %in% names(plot_data))) {
	  g <- g + geom_errorbar(
	    aes(ymin = band.LCL, ymax = band.UCL),
	    width = 0.42, linewidth = 0.5, linetype = 2, colour = "black"
	  )
	}

	# Default subtitle
	if (is_dunnett) {
	  ctrl <- input$control_variant
	  sub_default <- sprintf(
	    "Dunnett vs %s: a = jako kontrola, b = lepší, c = horší. %d%% CI.",
	    ctrl, ci_pct
	  )
	} else if (isTRUE(input$show_group_bands)) {
	  sub_default <- sprintf("Silné úsečky = %d%% CI, tenké přerušované = simultánní grouping bands, písmenka = CLD.", ci_pct)
	} else {
	  sub_default <- sprintf("Modelové odhady (emmeans) ± %d%% CI. Písmenka = statistické skupiny.", ci_pct)
	}

	g +
	  geom_text(aes(label = letter, y = label_y),
		    size = 5, fontface = "bold", vjust = 0, colour = plot_data$letter_colour) +
	    geom_text(aes(label = sprintf("%.1f%%", mean_pct)),
		      vjust = -0.5, size = 3.5, colour = "grey30") +
	    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.18))) +
	    labs(
	      x = lbl("lbl_bar_x", NULL),
	      y = lbl("lbl_bar_y", "Vzcházivost (%)"),
	      title = lbl("lbl_bar_title", sprintf("Celková vzcházivost v posledním datu (%s)", last_time)),
	      subtitle = lbl("lbl_bar_subtitle", sub_default)
	    ) +
	    theme_minimal(base_size = 14) +
	    theme(
	      legend.position = "none",
	      plot.title = element_text(face = "bold", size = 18),
	      plot.subtitle = element_text(colour = "grey50", size = 12),
	      axis.text.x = element_text(size = 12, face = "bold"),
	      panel.grid.major.x = element_blank()
	    )
	})
  
	output$main_plot <- renderPlot({
	  req(rv$emm_plot)
	  df <- rv$emm_plot
	  band_df <- rv$emm_group_band
	  ci_pct <- round((1 - input$alpha) * 100)

	  df$treatment <- factor(df$treatment, levels = rv$treatment_order)

	  if (!is.null(band_df)) {
	    join_keys <- intersect(c("treatment", "time"), names(df))
	    band_keep <- band_df[, intersect(c("treatment", "time", "band.LCL", "band.UCL"), names(band_df)), drop = FALSE]
	    df <- dplyr::left_join(df, band_keep, by = join_keys)
	  }

	  if (!"time" %in% names(df)) {
	    g <- ggplot(df, aes(x = treatment, y = prob, fill = treatment)) +
	      geom_col(width = 0.7, alpha = 0.85) +
	      geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.25, linewidth = 0.8)

	    if (isTRUE(input$show_group_bands) &&
		all(c("band.LCL", "band.UCL") %in% names(df))) {
	      g <- g + geom_errorbar(
		aes(ymin = band.LCL, ymax = band.UCL),
		width = 0.42, linewidth = 0.5, linetype = 2, colour = "black"
	      )
	    }

	    g +
	      scale_y_continuous(labels = scales::percent_format()) +
	      labs(
		x = lbl("lbl_emm_x", NULL),
		y = lbl("lbl_emm_y", "Pravděpodobnost klíčení"),
		title = lbl("lbl_emm_title", "Odhady pravděpodobností klíčení dle varianty"),
		subtitle = lbl(
		  "lbl_emm_subtitle",
		  if (isTRUE(input$show_group_bands)) {
		    sprintf("Silné úsečky = %d%% CI, tenké přerušované = simultánní grouping bands", ci_pct)
		  } else {
		    sprintf("Chybové úsečky = %d%% interval spolehlivosti", ci_pct)
		  }
		)
	      ) +
	      theme_minimal(base_size = 14) +
	      theme(
		legend.position = "none",
		axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
		plot.title = element_text(face = "bold", size = 16),
		plot.subtitle = element_text(colour = "grey50"),
		panel.grid.major.x = element_blank()
	      )
	  } else {
	    dodge <- position_dodge(width = 0.8)

	    g <- ggplot(df, aes(x = treatment, y = prob, fill = time)) +
	      geom_col(position = dodge, width = 0.7, alpha = 0.85) +
	      geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
			    position = dodge, width = 0.25, linewidth = 0.8)

	    if (isTRUE(input$show_group_bands) &&
		all(c("band.LCL", "band.UCL") %in% names(df))) {
	      g <- g + geom_errorbar(
		aes(ymin = band.LCL, ymax = band.UCL),
		position = dodge, width = 0.42, linewidth = 0.5, linetype = 2, colour = "black"
	      )
	    }

	    g +
	      scale_y_continuous(labels = scales::percent_format()) +
	      scale_fill_brewer(palette = "Set2") +
	      labs(
		x = lbl("lbl_emm_x", NULL),
		y = lbl("lbl_emm_y", "Podmíněná pravděpodobnost klíčení"),
		fill = lbl("lbl_emm_legend", "Datum"),
		title = lbl("lbl_emm_title", "Podmíněná pravděpodobnost klíčení"),
		subtitle = lbl(
		  "lbl_emm_subtitle",
		  if (isTRUE(input$show_group_bands)) {
		    sprintf("Odhady dle varianty a data: silné = %d%% CI, tenké přerušované = simultánní grouping bands", ci_pct)
		  } else {
		    sprintf("Odhady pravděpodobností klíčení dle varianty a data (%d%% CI)", ci_pct)
		  }
		)
	      ) +
	      theme_minimal(base_size = 14) +
	      theme(
		axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
		plot.title = element_text(face = "bold", size = 16),
		plot.subtitle = element_text(colour = "grey50"),
		legend.position = "top",
		panel.grid.major.x = element_blank()
	      )
	  }
	})

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

      if (!is.null(rv$cum_group_band)) {
        band_keep <- rv$cum_group_band[, intersect(c("treatment", "time", "band.LCL", "band.UCL"), names(rv$cum_group_band)), drop = FALSE]
        timeline <- dplyr::left_join(timeline, band_keep, by = c("treatment", "time"))
        if ("band.LCL" %in% names(timeline)) timeline$band.LCL <- timeline$band.LCL * 100
        if ("band.UCL" %in% names(timeline)) timeline$band.UCL <- timeline$band.UCL * 100
      }

      sub_default <- if (isTRUE(input$show_group_bands)) {
        sprintf("Modelové odhady (emmeans) ± %d%% CI + simultánní grouping bands — orientační (korelovaná data)", ci_pct)
      } else {
        sprintf("Modelové odhady (emmeans) ± %d%% CI — orientační (korelovaná data)", ci_pct)
      }
    }

    timeline$treatment <- factor(timeline$treatment, levels = rv$treatment_order)
    dodge <- position_dodge(width = 0.8)

    g <- ggplot(timeline, aes(x = treatment, y = mean_pct, fill = time)) +
      geom_col(position = dodge, width = 0.7, alpha = 0.85) +
      geom_errorbar(
        aes(ymin = ci_lo, ymax = ci_hi),
        position = dodge, width = 0.25, linewidth = 0.8
      )

    if (isTRUE(input$show_group_bands) &&
        all(c("band.LCL", "band.UCL") %in% names(timeline))) {
      g <- g + geom_errorbar(
        aes(ymin = band.LCL, ymax = band.UCL),
        position = dodge, width = 0.42, linewidth = 0.5,
        linetype = 2, colour = "black"
      )
    }

    g +
      scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0))) +
      scale_fill_brewer(palette = "Set2") +
      labs(
        x = lbl("lbl_time_x", NULL),
        y = lbl("lbl_time_y", "Kumulativní klíčivost (%)"),
        fill = lbl("lbl_time_legend", "Datum"),
        title = lbl("lbl_time_title", "Postup klíčivosti dle varianty"),
        subtitle = lbl("lbl_time_subtitle", sub_default)
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(colour = "grey50", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        legend.position = "top",
        panel.grid.major.x = element_blank()
      )
  })
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

  output$final_sig_heatmap <- renderPlot({
    req(rv$final_posthoc)
    make_sig_heatmap_plot(rv$final_posthoc, "Jednoduchý test — matice signifikance", sprintf("Kontrasty po korekci %s; alpha = %.3f", input$p_method, input$alpha))
  })

  output$final_contrast_plot <- renderPlot({
    req(rv$final_posthoc)
    make_contrast_plot(rv$final_posthoc, "Jednoduchý test — forest plot kontrastů", sprintf("Adjustované %d%% intervaly kontrastů", round((1 - input$alpha) * 100)))
  })

  output$cum_sig_heatmap <- renderPlot({
    req(rv$cum_posthoc)
    make_sig_heatmap_plot(rv$cum_posthoc, "Kumulativní model — matice signifikance", sprintf("Po datech; kontrasty po korekci %s", input$p_method))
  })

  output$cum_contrast_plot <- renderPlot({
    req(rv$cum_posthoc)
    make_contrast_plot(rv$cum_posthoc, "Kumulativní model — forest plot kontrastů", sprintf("Po datech; adjustované %d%% intervaly kontrastů", round((1 - input$alpha) * 100)))
  })

  output$cond_sig_heatmap <- renderPlot({
    req(rv$posthoc)
    make_sig_heatmap_plot(rv$posthoc, "Podmíněný model — matice signifikance", sprintf("Kontrasty po korekci %s", input$p_method))
  })

  output$cond_contrast_plot <- renderPlot({
    req(rv$posthoc)
    make_contrast_plot(rv$posthoc, "Podmíněný model — forest plot kontrastů", sprintf("Adjustované %d%% intervaly kontrastů", round((1 - input$alpha) * 100)))
  })

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

  output$cum_emm_results <- renderDT({
    req(rv$cum_emm_marginal)
    df <- rv$cum_emm_marginal; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

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

  output$dunnett_table <- renderDT({
    req(rv$dunnett_posthoc)
    df <- rv$dunnett_posthoc; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$dunnett_cum_table <- renderDT({
    req(rv$dunnett_cum_posthoc)
    df <- rv$dunnett_cum_posthoc; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$dunnett_final_table <- renderDT({
    req(rv$dunnett_final_posthoc)
    df <- rv$dunnett_final_posthoc; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$dunnett_final_results <- renderDT({
    req(rv$dunnett_final_posthoc)
    df <- rv$dunnett_final_posthoc; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$dunnett_cum_results <- renderDT({
    req(rv$dunnett_cum_posthoc)
    df <- rv$dunnett_cum_posthoc; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$dunnett_cond_results <- renderDT({
    req(rv$dunnett_posthoc)
    df <- rv$dunnett_posthoc; nums <- names(df)[sapply(df, is.numeric)]
    datatable(df, options = list(pageLength = 30, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = nums, digits = 4)
  })

  output$dunnett_final_contrast_plot <- renderPlot({
    req(rv$dunnett_final_posthoc)
    ctrl <- input$control_variant
    make_contrast_plot(
      rv$dunnett_final_posthoc,
      sprintf("Jednoduchy test - Dunnett (vs %s)", ctrl),
      sprintf("Dunnettova korekce; %d%% intervaly kontrastu", round((1 - input$alpha) * 100))
    )
  })

  output$dunnett_cum_contrast_plot <- renderPlot({
    req(rv$dunnett_cum_posthoc)
    ctrl <- input$control_variant
    make_contrast_plot(
      rv$dunnett_cum_posthoc,
      sprintf("Kumulativni model - Dunnett (vs %s)", ctrl),
      sprintf("Dunnettova korekce; %d%% intervaly kontrastu (orientacni)", round((1 - input$alpha) * 100))
    )
  })

  output$dunnett_cond_contrast_plot <- renderPlot({
    req(rv$dunnett_posthoc)
    ctrl <- input$control_variant
    make_contrast_plot(
      rv$dunnett_posthoc,
      sprintf("Podmineny model - Dunnett (vs %s)", ctrl),
      sprintf("Dunnettova korekce; %d%% intervaly kontrastu", round((1 - input$alpha) * 100))
    )
  })

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

  output$log_text <- renderPrint({ cat(rv$log) })
}

shinyApp(ui, server)
