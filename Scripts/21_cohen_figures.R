# ============================================================
# Dementia Cohen / reversal figures 
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "ggrepel",
  "stringr", "forcats", "scales"
)

to_install <- packages[!packages %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(forcats)
library(scales)

# ------------------------------------------------------------
# Input / output
# ------------------------------------------------------------

input_xlsx <- "Results/dementia_cohen_all_tests.xlsx"   # written by 12_dementia_cohen.R

out_dir <- "figures_dementia_cohen"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw <- read_excel(input_xlsx, sheet = 2)
colnames(raw) <- trimws(colnames(raw))

# ------------------------------------------------------------
# Robust column detection
# ------------------------------------------------------------

find_col <- function(candidates, cols = colnames(raw), required = TRUE) {
  hit <- candidates[candidates %in% cols]
  if (length(hit) == 0) {
    if (required) {
      stop(
        "Could not find any of these columns: ",
        paste(candidates, collapse = ", ")
      )
    } else {
      return(NA_character_)
    }
  }
  hit[1]
}

intervention_dataset_col   <- find_col(c("intervention_dataset"))
intervention_condition_col <- find_col(c("intervention_condition", "intervention", "dataset"))
intervention_type_col      <- find_col(c("intervention_type", "type"), required = FALSE)

dementia_dataset_col <- find_col(
  c("dementia_dataset", "dementa_dataset", "dementia_dataset ", "dataset_id")
)

dementia_condition_col <- find_col(
  c("dementia_condition", "condition", "disease_condition")
)

direction_col <- find_col(c("direction"))
cohen_col     <- find_col(c("cohen", "Cohen", "cohen_d", "cohen_disease"))
pvalue_col    <- find_col(c("pvalue", "p_value", "pval"), required = FALSE)

p_onetail_rev_col <- find_col(
  c("p_onetail_rev", "p_one_tail_rev", "p_rev"),
  required = FALSE
)

p_onetail_mimick_col <- find_col(
  c("p_onetail_mimick", "p_onetail_mimic", "p_one_tail_mimick", "p_mimick"),
  required = FALSE
)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

first_non_na <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

safe_name <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+$", "")
}

neglog10_safe <- function(p) {
  p <- as.numeric(p)
  p <- pmax(p, 1e-300)
  -log10(p)
}

zero_line_col <- "red"

type_cols <- c(
  "Psychoplastogen"  = "black",
  "Positive Control" = "#8E24AA",
  "Negative Control" = "#F57C00"
)

theme_dementia <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid       = element_blank(),
      
      axis.line        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks       = element_line(colour = "black", linewidth = 0.4),
      axis.text        = element_text(colour = "black"),
      axis.title       = element_text(colour = "black"),
      
      plot.title       = element_text(
        colour = "black",
        face = "bold",
        hjust = 0.5,
        size = base_size
      ),
      
      legend.position   = "right",
      legend.title      = element_blank(),
      legend.text       = element_text(colour = "black"),
      legend.background = element_blank(),
      legend.key        = element_blank()
    )
}

# ------------------------------------------------------------
# Prepare long table
# ------------------------------------------------------------

dementia_long_raw <- raw %>%
  transmute(
    intervention_dataset = as.character(.data[[intervention_dataset_col]]),
    intervention_condition = as.character(.data[[intervention_condition_col]]),
    
    intervention_type = if (!is.na(intervention_type_col)) {
      as.character(.data[[intervention_type_col]])
    } else {
      NA_character_
    },
    
    dementia_dataset = as.character(.data[[dementia_dataset_col]]),
    dementia_condition = as.character(.data[[dementia_condition_col]]),
    
    direction = toupper(as.character(.data[[direction_col]])),
    cohen = as.numeric(.data[[cohen_col]]),
    
    pvalue = if (!is.na(pvalue_col)) {
      as.numeric(.data[[pvalue_col]])
    } else {
      NA_real_
    },
    
    p_onetail_rev = if (!is.na(p_onetail_rev_col)) {
      as.numeric(.data[[p_onetail_rev_col]])
    } else {
      NA_real_
    },
    
    p_onetail_mimick = if (!is.na(p_onetail_mimick_col)) {
      as.numeric(.data[[p_onetail_mimick_col]])
    } else {
      NA_real_
    }
  ) %>%
  mutate(
    direction = recode(direction, "DOWN" = "DN")
  ) %>%
  filter(
    direction %in% c("UP", "DN"),
    !is.na(cohen),
    !is.na(intervention_dataset),
    !is.na(intervention_condition),
    !is.na(dementia_dataset),
    !is.na(dementia_condition)
  )

# Collapse possible duplicated rows while keeping each dementia dataset,
# dementia condition, intervention condition, and intervention dataset distinct.
dementia_long <- dementia_long_raw %>%
  group_by(
    intervention_dataset,
    intervention_condition,
    dementia_dataset,
    dementia_condition,
    direction
  ) %>%
  summarise(
    cohen = mean(cohen, na.rm = TRUE),
    pvalue = suppressWarnings(min(pvalue, na.rm = TRUE)),
    p_onetail_rev = suppressWarnings(min(p_onetail_rev, na.rm = TRUE)),
    p_onetail_mimick = suppressWarnings(min(p_onetail_mimick, na.rm = TRUE)),
    intervention_type = first_non_na(intervention_type),
    .groups = "drop"
  ) %>%
  mutate(
    pvalue = ifelse(is.infinite(pvalue), NA_real_, pvalue),
    p_onetail_rev = ifelse(is.infinite(p_onetail_rev), NA_real_, p_onetail_rev),
    p_onetail_mimick = ifelse(is.infinite(p_onetail_mimick), NA_real_, p_onetail_mimick)
  )


entity_info <- dementia_long %>%
  distinct(
    intervention_condition,
    intervention_dataset,
    intervention_type
  ) %>%
  group_by(intervention_condition) %>%
  mutate(
    n_intervention_datasets = n_distinct(intervention_dataset)
  ) %>%
  ungroup() %>%
  mutate(
    intervention_entity = paste(
      intervention_condition,
      intervention_dataset,
      sep = " | "
    ),
    
    display_label = ifelse(
      n_intervention_datasets > 1,
      paste0(intervention_condition, "\n", intervention_dataset),
      intervention_condition
    ),
    
    intervention_group = case_when(
      intervention_type %in% c(
        "Psychoplastogen",
        "Positive Control",
        "Negative Control"
      ) ~ intervention_type,
      
      str_detect(intervention_condition, regex("Exercise|EE", ignore_case = TRUE)) ~
        "Positive Control",
      
      str_detect(intervention_condition, regex("Alcohol|Rat aging", ignore_case = TRUE)) ~
        "Negative Control",
      
      TRUE ~ "Psychoplastogen"
    ),
    
    intervention_group = factor(
      intervention_group,
      levels = c("Psychoplastogen", "Positive Control", "Negative Control")
    )
  ) %>%
  select(
    intervention_condition,
    intervention_dataset,
    intervention_entity,
    display_label,
    intervention_group
  )

dementia_long <- dementia_long %>%
  left_join(
    entity_info,
    by = c("intervention_condition", "intervention_dataset")
  ) %>%
  mutate(
    dementia_label = paste0(dementia_dataset, "\n", dementia_condition),
    dementia_entity = paste(dementia_dataset, dementia_condition, sep = " | ")
  )

dementia_wide <- dementia_long %>%
  select(
    intervention_condition,
    intervention_dataset,
    intervention_entity,
    display_label,
    intervention_group,
    dementia_dataset,
    dementia_condition,
    dementia_entity,
    dementia_label,
    direction,
    cohen
  ) %>%
  pivot_wider(
    names_from = direction,
    values_from = cohen,
    names_prefix = "cohen_"
  ) %>%
  filter(!is.na(cohen_UP), !is.na(cohen_DN))

all_entities <- entity_info %>%
  arrange(intervention_condition, intervention_dataset) %>%
  pull(intervention_entity)

# ============================================================
# 1. Scatterplot per intervention + intervention_dataset
# x = Cohen UP, y = Cohen DN
# ============================================================

plot_scatter_one_entity <- function(entity_id) {
  
  dat <- dementia_wide %>%
    filter(intervention_entity == entity_id)
  
  plot_title <- unique(dat$display_label)
  plot_title <- paste(plot_title, collapse = " / ")
  
  ggplot(
    dat,
    aes(x = cohen_UP, y = cohen_DN)
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      colour = zero_line_col,
      alpha = 0.45
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      colour = zero_line_col,
      alpha = 0.45
    ) +
    geom_point(
      aes(colour = intervention_group),
      size = 2.5
    ) +
    geom_text_repel(
      aes(label = dementia_label, colour = intervention_group),
      size = 2.8,
      max.overlaps = Inf,
      segment.size = 0.25,
      box.padding = 0.35,
      point.padding = 0.2,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = type_cols,
      breaks = c("Psychoplastogen", "Positive Control", "Negative Control")
    ) +
    labs(
      title = paste0(plot_title, ": Cohen UP vs Cohen DN"),
      x = "Cohen's d, UP with intervention",
      y = "Cohen's d, DN with intervention"
    ) +
    theme_dementia()
}

scatter_pdf <- file.path(
  out_dir,
  "scatter_Cohen_UP_vs_DN_by_intervention_dataset.pdf"
)

pdf(scatter_pdf, width = 7, height = 6)

for (entity_id in all_entities) {
  
  p <- plot_scatter_one_entity(entity_id)
  print(p)
  
  dat_entity <- entity_info %>%
    filter(intervention_entity == entity_id) %>%
    slice(1)
  
  out_name <- paste0(
    "scatter_Cohen_",
    safe_name(dat_entity$intervention_condition),
    "__",
    safe_name(dat_entity$intervention_dataset),
    ".png"
  )
  
  ggsave(
    filename = file.path(out_dir, out_name),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
}

dev.off()

# ============================================================
# 2. Boxplot Cohen UP
#    Boxes separated by intervention_condition + intervention_dataset
#    Ordered by mean Cohen UP
# ============================================================

cohen_up_df <- dementia_long %>%
  filter(direction == "UP")

up_order <- cohen_up_df %>%
  group_by(intervention_entity, display_label) %>%
  summarise(
    mean_cohen = mean(cohen, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_cohen)

cohen_up_df <- cohen_up_df %>%
  mutate(
    intervention_entity_ord = factor(
      intervention_entity,
      levels = up_order$intervention_entity
    )
  )

label_map_up <- setNames(
  up_order$display_label,
  up_order$intervention_entity
)

p_box_up <- ggplot(
  cohen_up_df,
  aes(
    x = intervention_entity_ord,
    y = cohen,
    colour = intervention_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    colour = zero_line_col,
    alpha = 0.45
  ) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    fill = "white"
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.4,
    alpha = 0.7
  ) +
  scale_colour_manual(
    values = type_cols,
    breaks = c("Psychoplastogen", "Positive Control", "Negative Control")
  ) +
  scale_x_discrete(labels = label_map_up) +
  labs(
    title = "Distribution of Cohen's d for UP genes across intervention datasets",
    x = "Intervention / intervention dataset ordered by mean Cohen's d",
    y = "Cohen's d, UP with intervention"
  ) +
  theme_dementia() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

p_box_up

ggsave(
  filename = file.path(out_dir, "boxplot_Cohen_UP_by_intervention_dataset_ordered_by_mean.png"),
  plot = p_box_up,
  width = 14,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "boxplot_Cohen_UP_by_intervention_dataset_ordered_by_mean.pdf"),
  plot = p_box_up,
  width = 14,
  height = 6
)

# ============================================================
# 3. Boxplot Cohen DN
#    Boxes separated by intervention_condition + intervention_dataset
#    Ordered by mean Cohen DN
# ============================================================

cohen_dn_df <- dementia_long %>%
  filter(direction == "DN")

dn_order <- cohen_dn_df %>%
  group_by(intervention_entity, display_label) %>%
  summarise(
    mean_cohen = mean(cohen, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_cohen)

cohen_dn_df <- cohen_dn_df %>%
  mutate(
    intervention_entity_ord = factor(
      intervention_entity,
      levels = dn_order$intervention_entity
    )
  )

label_map_dn <- setNames(
  dn_order$display_label,
  dn_order$intervention_entity
)

p_box_dn <- ggplot(
  cohen_dn_df,
  aes(
    x = intervention_entity_ord,
    y = cohen,
    colour = intervention_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    colour = zero_line_col,
    alpha = 0.45
  ) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    fill = "white"
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.4,
    alpha = 0.7
  ) +
  scale_colour_manual(
    values = type_cols,
    breaks = c("Psychoplastogen", "Positive Control", "Negative Control")
  ) +
  scale_x_discrete(labels = label_map_dn) +
  labs(
    title = "Distribution of Cohen's d for DN genes across intervention datasets",
    x = "Intervention / intervention dataset ordered by mean Cohen's d",
    y = "Cohen's d, DN with intervention"
  ) +
  theme_dementia() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

p_box_dn

ggsave(
  filename = file.path(out_dir, "boxplot_Cohen_DN_by_intervention_dataset_ordered_by_mean.png"),
  plot = p_box_dn,
  width = 14,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "boxplot_Cohen_DN_by_intervention_dataset_ordered_by_mean.pdf"),
  plot = p_box_dn,
  width = 14,
  height = 6
)

# ============================================================
# 4. Combined boxplot UP / DN
#    Boxes separated by intervention_condition + intervention_dataset
#    Ordered by global mean Cohen
# ============================================================

combined_order <- dementia_long %>%
  group_by(intervention_entity, display_label) %>%
  summarise(
    mean_cohen_global = mean(cohen, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_cohen_global)

cohen_combined_df <- dementia_long %>%
  mutate(
    intervention_entity_ord = factor(
      intervention_entity,
      levels = combined_order$intervention_entity
    ),
    direction = factor(direction, levels = c("UP", "DN"))
  )

label_map_combined <- setNames(
  combined_order$display_label,
  combined_order$intervention_entity
)

p_box_combined <- ggplot(
  cohen_combined_df,
  aes(
    x = intervention_entity_ord,
    y = cohen,
    fill = direction
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    colour = zero_line_col,
    alpha = 0.45
  ) +
  geom_boxplot(
    width = 0.7,
    position = position_dodge2(width = 0.75, preserve = "single"),
    outlier.size = 1.2,
    colour = "black"
  ) +
  scale_fill_manual(
    values = c("UP" = "red", "DN" = "blue"),
    name = "Direction"
  ) +
  scale_x_discrete(labels = label_map_combined) +
  labs(
    title = "Distribution of Cohen's d for UP and DN genes across intervention datasets",
    x = "Intervention / intervention dataset ordered by global mean Cohen's d",
    y = "Cohen's d"
  ) +
  theme_dementia() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top"
  )

p_box_combined

ggsave(
  filename = file.path(out_dir, "boxplot_Cohen_UP_DN_by_intervention_dataset_ordered_by_global_mean.png"),
  plot = p_box_combined,
  width = 16,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "boxplot_Cohen_UP_DN_by_intervention_dataset_ordered_by_global_mean.pdf"),
  plot = p_box_combined,
  width = 16,
  height = 6
)

# ============================================================
# 5. Scatterplot of mean Cohen
#    One point per intervention_condition + intervention_dataset
#    Legend aggregated by intervention_group
# ============================================================

mean_cohen <- dementia_wide %>%
  group_by(
    intervention_condition,
    intervention_dataset,
    intervention_entity,
    display_label,
    intervention_group
  ) %>%
  summarise(
    mean_cohen_UP = mean(cohen_UP, na.rm = TRUE),
    mean_cohen_DN = mean(cohen_DN, na.rm = TRUE),
    n_dementia_tests = n(),
    .groups = "drop"
  )

p_mean_scatter <- ggplot(
  mean_cohen,
  aes(
    x = mean_cohen_UP,
    y = mean_cohen_DN,
    colour = intervention_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = zero_line_col,
    linewidth = 0.35,
    alpha = 0.45
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = zero_line_col,
    linewidth = 0.35,
    alpha = 0.45
  ) +
  geom_point(
    size = 2.8
  ) +
  geom_text_repel(
    aes(label = display_label),
    size = 3,
    max.overlaps = Inf,
    segment.size = 0.25,
    box.padding = 0.35,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = type_cols,
    breaks = c("Psychoplastogen", "Positive Control", "Negative Control")
  ) +
  labs(
    title = "Mean Cohen's d per intervention dataset",
    x = "Mean Cohen's d, UP with intervention",
    y = "Mean Cohen's d, DN with intervention"
  ) +
  theme_dementia(base_size = 11)

p_mean_scatter

ggsave(
  filename = file.path(out_dir, "scatter_mean_Cohen_UP_vs_DN_by_intervention_dataset.png"),
  plot = p_mean_scatter,
  width = 7,
  height = 5.5,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "scatter_mean_Cohen_UP_vs_DN_by_intervention_dataset.pdf"),
  plot = p_mean_scatter,
  width = 7,
  height = 5.5
)

# ============================================================
# 6. Dementia-reversal p-ratio scatterplot
#    One point per intervention_condition + intervention_dataset
#
# DRP = -log10(p_reversal) - [-log10(p_mimicking)]
#     = log10(p_mimicking / p_reversal)
#
# Positive DRP = more dementia-reversal-like evidence
# Negative DRP = more dementia-mimicking-like evidence
# ============================================================

has_pratio <- all(
  !is.na(c(p_onetail_rev_col, p_onetail_mimick_col))
)

if (has_pratio) {
  
  drp_long <- dementia_long %>%
    filter(
      !is.na(p_onetail_rev),
      !is.na(p_onetail_mimick),
      p_onetail_rev > 0,
      p_onetail_mimick > 0
    ) %>%
    mutate(
      p_onetail_rev = pmax(p_onetail_rev, 1e-300),
      p_onetail_mimick = pmax(p_onetail_mimick, 1e-300),
      DRP = -log10(p_onetail_rev) - (-log10(p_onetail_mimick))
    )
  
  drp_entity <- drp_long %>%
    group_by(
      intervention_condition,
      intervention_dataset,
      intervention_entity,
      display_label,
      intervention_group,
      direction
    ) %>%
    summarise(
      DRP = sum(DRP, na.rm = TRUE),
      mean_DRP = mean(DRP, na.rm = TRUE),
      n_dementia_tests = n(),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = direction,
      values_from = DRP,
      names_prefix = "DRP_"
    ) %>%
    filter(!is.na(DRP_UP), !is.na(DRP_DN))
  
  p_drp <- ggplot(
    drp_entity,
    aes(
      x = DRP_UP,
      y = DRP_DN,
      colour = intervention_group
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      colour = zero_line_col,
      alpha = 0.45
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      colour = zero_line_col,
      alpha = 0.45
    ) +
    geom_point(size = 2.8) +
    geom_text_repel(
      aes(label = display_label),
      size = 3,
      max.overlaps = Inf,
      segment.size = 0.25,
      box.padding = 0.35,
      point.padding = 0.25,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = type_cols,
      breaks = c("Psychoplastogen", "Positive Control", "Negative Control")
    ) +
    labs(
      title = "Dementia-reversal p ratio per intervention dataset",
      x = "DRP, UP with intervention",
      y = "DRP, DN with intervention"
    ) +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = "dementia-reversal",
      hjust = 1.05,
      vjust = -0.8,
      size = 3
    ) +
    annotate(
      "text",
      x = -Inf,
      y = -Inf,
      label = "dementia-mimicking",
      hjust = -0.05,
      vjust = -0.8,
      size = 3
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = "dementia-reversal",
      hjust = -0.05,
      vjust = 1.3,
      angle = 90,
      size = 3
    ) +
    annotate(
      "text",
      x = -Inf,
      y = -Inf,
      label = "dementia-mimicking",
      hjust = -0.05,
      vjust = -0.3,
      angle = 90,
      size = 3
    ) +
    coord_cartesian(clip = "off") +
    theme_dementia()
  
  p_drp
  
  ggsave(
    filename = file.path(out_dir, "scatter_DRP_UP_vs_DN_by_intervention_dataset.png"),
    plot = p_drp,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(out_dir, "scatter_DRP_UP_vs_DN_by_intervention_dataset.pdf"),
    plot = p_drp,
    width = 7,
    height = 6
  )
  
  write.csv(
    drp_entity,
    file = file.path(out_dir, "DRP_by_intervention_dataset.csv"),
    row.names = FALSE
  )
}


# ============================================================
# Save processed tables
# ============================================================

write.csv(
  dementia_long,
  file = file.path(out_dir, "dementia_long_processed_by_intervention_dataset.csv"),
  row.names = FALSE
)

write.csv(
  dementia_wide,
  file = file.path(out_dir, "dementia_wide_Cohen_UP_DN_by_intervention_dataset.csv"),
  row.names = FALSE
)

write.csv(
  mean_cohen,
  file = file.path(out_dir, "mean_Cohen_by_intervention_dataset.csv"),
  row.names = FALSE
)

write.csv(
  entity_info,
  file = file.path(out_dir, "intervention_dataset_entity_info.csv"),
  row.names = FALSE
)

message("Done. Figures saved in: ", normalizePath(out_dir))