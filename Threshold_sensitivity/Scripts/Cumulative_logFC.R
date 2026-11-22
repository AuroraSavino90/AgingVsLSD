# ============================================================
# THRESHOLD-SENSITIVITY COPY of Scripts/Cumulative_logFC.R
# Changes vs the original: (1) gene selection threshold set by
# SELECTION_MODE for ALL datasets (P_ARRAY / P_SEQ columns);
# (2) outputs written to OUT (Threshold_sensitivity/Results_<mode>);
# (3) see sensitivity_setup.R. Run via Threshold_sensitivity/run_all.R
# ============================================================
# ============================================================
# CRB PANELS E-F
# Cumulative Reversal Burden for aging and dementia
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
source("Threshold_sensitivity/Scripts/sensitivity_setup.R")  # SELECTION_MODE, OUT, P_ARRAY/P_SEQ, shims


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

# written by save_mean_logFC_selected_genes.R
# SELECTION_MODE must match the one used there ("pvalue" or "fdr")
selection_mode <- tolower(Sys.getenv("SELECTION_MODE", "pvalue"))
stopifnot(selection_mode %in% c("pvalue", "fdr"))
out_suffix <- if (selection_mode == "fdr") "_FDR" else "_pvalue"

aging_file <-
  paste0(OUT, "/logFC_vs_aging_signatures", out_suffix, ".xlsx")

dementia_file <-
  paste0(OUT, "/logFC_vs_dementia_signatures", out_suffix, ".xlsx")


# ------------------------------------------------------------
# 2. Read detailed summaries
# ------------------------------------------------------------

aging <-
  read_excel(
    aging_file,
    sheet = "Detailed_summary"
  )

dementia <-
  read_excel(
    dementia_file,
    sheet = "Detailed_summary"
  )


# ------------------------------------------------------------
# 3. Required columns
# ------------------------------------------------------------

required_columns <-
  c(
    "DE",
    "intervention",
    "type",
    "direction",
    "reference_dataset",
    "reference_region",
    "reference_condition",
    "n_same_direction",
    "mean_logFC_same_direction",
    "n_opposite_direction",
    "mean_logFC_opposite_direction"
  )


missing_aging <-
  setdiff(
    required_columns,
    colnames(aging)
  )

missing_dementia <-
  setdiff(
    required_columns,
    colnames(dementia)
  )


if (length(missing_aging) > 0) {
  
  stop(
    paste(
      "Missing columns in aging table:",
      paste(
        missing_aging,
        collapse = ", "
      )
    )
  )
}


if (length(missing_dementia) > 0) {
  
  stop(
    paste(
      "Missing columns in dementia table:",
      paste(
        missing_dementia,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 4. Clean character columns
#
# Helps avoid encoding problems in PDF output.
# ------------------------------------------------------------

clean_text_columns <- function(dat) {
  
  dat %>%
    
    mutate(
      
      DE =
        enc2utf8(
          as.character(DE)
        ),
      
      intervention =
        enc2utf8(
          as.character(intervention)
        ),
      
      type =
        enc2utf8(
          as.character(type)
        ),
      
      direction =
        enc2utf8(
          as.character(direction)
        ),
      
      reference_dataset =
        enc2utf8(
          as.character(reference_dataset)
        ),
      
      reference_region =
        enc2utf8(
          as.character(reference_region)
        ),
      
      reference_condition =
        enc2utf8(
          as.character(reference_condition)
        )
    )
}


aging <-
  clean_text_columns(
    aging
  )

dementia <-
  clean_text_columns(
    dementia
  )


# ------------------------------------------------------------
# 5. Calculate CRB for each reference comparison
#
# concordant burden =
# n concordant genes *
# abs(mean logFC concordant genes)
#
# opposite burden =
# n opposite genes *
# abs(mean logFC opposite genes)
#
# CRB =
# opposite burden - concordant burden
# ------------------------------------------------------------

calculate_crb <- function(dat) {
  
  dat %>%
    
    mutate(
      
      concordant_burden =
        n_same_direction *
        abs(
          mean_logFC_same_direction
        ),
      
      opposite_burden =
        n_opposite_direction *
        abs(
          mean_logFC_opposite_direction
        ),
      
      CRB =
        opposite_burden -
        concordant_burden
    )
}


aging <-
  calculate_crb(
    aging
  )

dementia <-
  calculate_crb(
    dementia
  )


# ------------------------------------------------------------
# 6. Summarise CRB across reference datasets
#
# Produces one row per intervention:
#
# CRB_UP
# CRB_DN
# CRB_total
#
# Median values are also retained for inspection,
# but plots use mean CRB.
# ------------------------------------------------------------

summarise_crb <- function(dat) {
  
  # ----------------------------------------------------------
  # Mean CRB
  # ----------------------------------------------------------
  
  mean_tab <-
    dat %>%
    
    group_by(
      DE,
      intervention,
      type,
      direction
    ) %>%
    
    summarise(
      
      mean_CRB =
        mean(
          CRB,
          na.rm = TRUE
        ),
      
      .groups =
        "drop"
    ) %>%
    
    pivot_wider(
      
      id_cols =
        c(
          DE,
          intervention,
          type
        ),
      
      names_from =
        direction,
      
      values_from =
        mean_CRB,
      
      names_prefix =
        "CRB_"
    )
  
  
  # ----------------------------------------------------------
  # Median CRB
  # ----------------------------------------------------------
  
  median_tab <-
    dat %>%
    
    group_by(
      DE,
      intervention,
      type,
      direction
    ) %>%
    
    summarise(
      
      median_CRB =
        median(
          CRB,
          na.rm = TRUE
        ),
      
      .groups =
        "drop"
    ) %>%
    
    pivot_wider(
      
      id_cols =
        c(
          DE,
          intervention,
          type
        ),
      
      names_from =
        direction,
      
      values_from =
        median_CRB,
      
      names_prefix =
        "median_CRB_"
    )
  
  
  # ----------------------------------------------------------
  # Number of reference comparisons
  # ----------------------------------------------------------
  
  n_tests <-
    dat %>%
    
    group_by(
      DE,
      intervention,
      type
    ) %>%
    
    summarise(
      
      n_reference_tests =
        n_distinct(
          interaction(
            reference_dataset,
            reference_region,
            reference_condition,
            drop = TRUE
          )
        ),
      
      .groups =
        "drop"
    )
  
  
  # ----------------------------------------------------------
  # Merge summaries
  # ----------------------------------------------------------
  
  out <-
    mean_tab %>%
    
    left_join(
      median_tab,
      by =
        c(
          "DE",
          "intervention",
          "type"
        )
    ) %>%
    
    left_join(
      n_tests,
      by =
        c(
          "DE",
          "intervention",
          "type"
        )
    )
  
  
  # ----------------------------------------------------------
  # Safety checks for missing direction columns
  # ----------------------------------------------------------
  
  if (!"CRB_UP" %in% colnames(out)) {
    
    out$CRB_UP <-
      NA_real_
  }
  
  if (!"CRB_DN" %in% colnames(out)) {
    
    out$CRB_DN <-
      NA_real_
  }
  
  if (!"median_CRB_UP" %in% colnames(out)) {
    
    out$median_CRB_UP <-
      NA_real_
  }
  
  if (!"median_CRB_DN" %in% colnames(out)) {
    
    out$median_CRB_DN <-
      NA_real_
  }
  
  
  # ----------------------------------------------------------
  # Total CRB
  # ----------------------------------------------------------
  
  out <-
    out %>%
    
    mutate(
      
      CRB_total =
        CRB_UP +
        CRB_DN,
      
      median_CRB_total =
        median_CRB_UP +
        median_CRB_DN
    )
  
  
  return(
    out
  )
}


aging_CRB <-
  summarise_crb(
    aging
  )

dementia_CRB <-
  summarise_crb(
    dementia
  )


# ------------------------------------------------------------
# 7. Rename comparator classes
#
# Same biological classes as other plots:
#
# reverting comparator  = purple
# concordant comparator = orange
# psychoplastogen       = black
# ------------------------------------------------------------

aging_CRB <-
  aging_CRB %>%
  
  mutate(
    
    type =
      recode(
        type,
        
        "Positive Control" =
          "Aging-Reverting Comparator",
        
        "Negative Control" =
          "Aging-Concordant Comparator",
        
        .default =
          type
      )
  )


dementia_CRB <-
  dementia_CRB %>%
  
  mutate(
    
    type =
      recode(
        type,
        
        "Positive Control" =
          "Dementia-Reverting Comparator",
        
        "Negative Control" =
          "Dementia-Concordant Comparator",
        
        .default =
          type
      )
  )


# ------------------------------------------------------------
# 8. Remove incomplete rows from plot data
# ------------------------------------------------------------

aging_plot <-
  aging_CRB %>%
  
  filter(
    !is.na(CRB_UP),
    !is.na(CRB_DN),
    is.finite(CRB_UP),
    is.finite(CRB_DN)
  )


dementia_plot <-
  dementia_CRB %>%
  
  filter(
    !is.na(CRB_UP),
    !is.na(CRB_DN),
    is.finite(CRB_UP),
    is.finite(CRB_DN)
  )


# ------------------------------------------------------------
# 9. Print LSD values
# ------------------------------------------------------------

cat(
  "\n============================\n"
)

cat(
  "AGING - LSD\n"
)

cat(
  "============================\n"
)


print(
  
  aging_CRB %>%
    
    filter(
      grepl(
        "LSD",
        intervention,
        ignore.case = TRUE
      )
    )
)


cat(
  "\n============================\n"
)

cat(
  "DEMENTIA - LSD\n"
)

cat(
  "============================\n"
)


print(
  
  dementia_CRB %>%
    
    filter(
      grepl(
        "LSD",
        intervention,
        ignore.case = TRUE
      )
    )
)


# ------------------------------------------------------------
# 10. Check missing coordinates
# ------------------------------------------------------------

cat(
  "\nAGING rows missing CRB_UP or CRB_DN:\n"
)

print(
  
  aging_CRB %>%
    
    filter(
      is.na(CRB_UP) |
        is.na(CRB_DN)
    )
)


cat(
  "\nDEMENTIA rows missing CRB_UP or CRB_DN:\n"
)

print(
  
  dementia_CRB %>%
    
    filter(
      is.na(CRB_UP) |
        is.na(CRB_DN)
    )
)


# ============================================================
# 11. PANEL E — AGING
#
# ============================================================

set.seed(
  123
)


p_aging <-
  
  ggplot(
    
    aging_plot,
    
    aes(
      
      x =
        CRB_UP,
      
      y =
        CRB_DN,
      
      label =
        intervention,
      
      colour =
        type
    )
  ) +
  
  geom_point(
    size = 2.2
  ) +
  
  geom_text_repel(
    
    size =
      3.2,
    
    box.padding =
      0.7,
    
    point.padding =
      0.5,
    
    force =
      4,
    
    force_pull =
      0.2,
    
    max.overlaps =
      Inf,
    
    max.time =
      20,
    
    max.iter =
      100000,
    
    min.segment.length =
      0,
    
    segment.size =
      0.25,
    
    segment.colour =
      "grey55",
    
    seed =
      123
  ) +
  
  theme_classic() +
  
  geom_vline(
    
    xintercept =
      0,
    
    colour =
      "red",
    
    linetype =
      "dashed"
  ) +
  
  geom_hline(
    
    yintercept =
      0,
    
    colour =
      "red",
    
    linetype =
      "dashed"
  ) +
  
  scale_colour_manual(
    
    values =
      c(
        
        "Aging-Reverting Comparator" =
          "#8E24AA",
        
        "Aging-Concordant Comparator" =
          "#F57C00",
        
        "Psychoplastogen" =
          "black"
      )
  ) +
  
  scale_x_continuous(
    
    expand =
      expansion(
        mult =
          c(
            0.15,
            0.30
          )
      )
  ) +
  
  scale_y_continuous(
    
    expand =
      expansion(
        mult =
          c(
            0.15,
            0.25
          )
      )
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  labs(
    
    x =
      "Cumulative reversal burden - up-regulated genes",
    
    y =
      "Cumulative reversal burden - down-regulated genes",
    
    colour =
      NULL
  ) +
  
  theme(
    
    legend.position =
      "right",
    
    plot.margin =
      margin(
        t = 20,
        r = 60,
        b = 20,
        l = 30
      )
  )


print(
  p_aging
)


# ------------------------------------------------------------
# Save aging panel
# ------------------------------------------------------------

ggsave(
  
  filename =
    file.path(OUT, paste0("CRB_panel_E_aging", out_suffix, ".pdf")),
  
  plot =
    p_aging,
  
  width =
    8.5,
  
  height =
    7,
  
  device =
    cairo_pdf
)


ggsave(
  
  filename =
    file.path(OUT, paste0("CRB_panel_E_aging", out_suffix, ".png")),
  
  plot =
    p_aging,
  
  width =
    8.5,
  
  height =
    7,
  
  dpi =
    600
)


# ============================================================
# 12. PANEL F — DEMENTIA
# ============================================================

set.seed(
  123
)


p_dementia <-
  
  ggplot(
    
    dementia_plot,
    
    aes(
      
      x =
        CRB_UP,
      
      y =
        CRB_DN,
      
      label =
        intervention,
      
      colour =
        type
    )
  ) +
  
  geom_point(
    size = 2.2
  ) +
  
  geom_text_repel(
    
    size =
      3.2,
    
    box.padding =
      0.7,
    
    point.padding =
      0.5,
    
    force =
      4,
    
    force_pull =
      0.2,
    
    max.overlaps =
      Inf,
    
    max.time =
      20,
    
    max.iter =
      100000,
    
    min.segment.length =
      0,
    
    segment.size =
      0.25,
    
    segment.colour =
      "grey55",
    
    seed =
      123
  ) +
  
  theme_classic() +
  
  geom_vline(
    
    xintercept =
      0,
    
    colour =
      "red",
    
    linetype =
      "dashed"
  ) +
  
  geom_hline(
    
    yintercept =
      0,
    
    colour =
      "red",
    
    linetype =
      "dashed"
  ) +
  
  scale_colour_manual(
    
    values =
      c(
        
        "Dementia-Reverting Comparator" =
          "#8E24AA",
        
        "Dementia-Concordant Comparator" =
          "#F57C00",
        
        "Psychoplastogen" =
          "black"
      )
  ) +
  
  scale_x_continuous(
    
    expand =
      expansion(
        mult =
          c(
            0.15,
            0.30
          )
      )
  ) +
  
  scale_y_continuous(
    
    expand =
      expansion(
        mult =
          c(
            0.15,
            0.25
          )
      )
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  labs(
    
    x =
      "Cumulative reversal burden - up-regulated genes",
    
    y =
      "Cumulative reversal burden - down-regulated genes",
    
    colour =
      NULL
  ) +
  
  theme(
    
    legend.position =
      "right",
    
    plot.margin =
      margin(
        t = 20,
        r = 60,
        b = 20,
        l = 30
      )
  )


print(
  p_dementia
)


# ------------------------------------------------------------
# Save dementia panel
# ------------------------------------------------------------

ggsave(
  
  filename =
    file.path(OUT, paste0("CRB_panel_F_dementia", out_suffix, ".pdf")),
  
  plot =
    p_dementia,
  
  width =
    8.5,
  
  height =
    7,
  
  device =
    cairo_pdf
)


ggsave(
  
  filename =
    file.path(OUT, paste0("CRB_panel_F_dementia", out_suffix, ".png")),
  
  plot =
    p_dementia,
  
  width =
    8.5,
  
  height =
    7,
  
  dpi =
    600
)


# ============================================================
# 13. Save CRB values used in plots
# ============================================================

wb <-
  createWorkbook()


# ------------------------------------------------------------
# Aging summary
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Aging CRB"
)

writeData(
  wb,
  sheet =
    "Aging CRB",
  x =
    aging_CRB
)


# ------------------------------------------------------------
# Dementia summary
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Dementia CRB"
)

writeData(
  wb,
  sheet =
    "Dementia CRB",
  x =
    dementia_CRB
)


# ------------------------------------------------------------
# Aging reference-level CRB
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Aging reference CRB"
)

writeData(
  wb,
  sheet =
    "Aging reference CRB",
  x =
    aging
)


# ------------------------------------------------------------
# Dementia reference-level CRB
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Dementia reference CRB"
)

writeData(
  wb,
  sheet =
    "Dementia reference CRB",
  x =
    dementia
)


# ------------------------------------------------------------
# Save workbook
# ------------------------------------------------------------

saveWorkbook(
  
  wb,
  
  file =
    file.path(OUT, paste0("CRB_values_panels_E_F", out_suffix, ".xlsx")),
  
  overwrite =
    TRUE
)


# ============================================================
# 14. Final summary
# ============================================================

cat(
  "\n====================================\n"
)

cat(
  "CRB analysis completed\n"
)

cat(
  "====================================\n"
)


cat(
  "\nAging interventions plotted:",
  nrow(
    aging_plot
  ),
  "\n"
)


cat(
  "Dementia interventions plotted:",
  nrow(
    dementia_plot
  ),
  "\n"
)


cat(
  "\nGenerated files:\n"
)

cat(
  "CRB_panel_E_aging.pdf\n"
)

cat(
  "CRB_panel_E_aging.png\n"
)

cat(
  "CRB_panel_F_dementia.pdf\n"
)

cat(
  "CRB_panel_F_dementia.png\n"
)

cat(
  "CRB_values_panels_E_F.xlsx\n"
)