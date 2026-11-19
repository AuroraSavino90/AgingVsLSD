# ============================================================
# NES / APR figures
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

input_xlsx <- "Results/all_GSEA_aging.xlsx"   # written by 09_age_reversal_otherpsychedelics.R
out_dir <- "figures_NES"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw <- read_excel(input_xlsx) %>%
  rename(
    intervention_dataset = dataset,
    p_onetail_mimick     = p_two2one_meta,
    p_onetail_rev        = p_two2one_rev
  )

required_cols <- c(
  "aging_dataset",
  "intervention_dataset",
  "intervention",
  "direction",
  "NES"
)

missing_cols <- setdiff(required_cols, colnames(raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

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

zero_line_col <- "red"

type_cols <- c(
  "Psychoplastogen"  = "black",
  "Positive Control" = "purple",
  "Negative Control" = "orange"
)

theme_age_reversal <- function(base_size = 11) {
  theme(
    panel.background = element_blank(),
    plot.background  = element_rect(fill = "white", colour = NA),
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
# Prepare data
# ------------------------------------------------------------

gsea_long_raw <- raw %>%
  transmute(
    aging_dataset        = as.character(aging_dataset),
    intervention_dataset = as.character(intervention_dataset),
    intervention         = as.character(intervention),
    direction            = toupper(as.character(direction)),
    NES                  = as.numeric(NES),
    
    type_raw = if ("type" %in% colnames(raw)) {
      as.character(type)
    } else {
      NA_character_
    },
    
    p_onetail_mimick = if ("p_onetail_mimick" %in% colnames(raw)) {
      as.numeric(p_onetail_mimick)
    } else {
      NA_real_
    },
    
    p_onetail_rev = if ("p_onetail_rev" %in% colnames(raw)) {
      as.numeric(p_onetail_rev)
    } else {
      NA_real_
    }
  ) %>%
  mutate(
    direction = recode(direction, "DOWN" = "DN")
  ) %>%
  filter(
    direction %in% c("UP", "DN"),
    !is.na(NES),
    !is.na(aging_dataset),
    !is.na(intervention),
    !is.na(intervention_dataset)
  )

# Collapse possible duplicated rows, but keep intervention_dataset distinct
gsea_long <- gsea_long_raw %>%
  group_by(
    aging_dataset,
    intervention,
    intervention_dataset,
    direction
  ) %>%
  summarise(
    NES = mean(NES, na.rm = TRUE),
    type_raw = first_non_na(type_raw),
    p_onetail_mimick = suppressWarnings(min(p_onetail_mimick, na.rm = TRUE)),
    p_onetail_rev    = suppressWarnings(min(p_onetail_rev, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    p_onetail_mimick = ifelse(is.infinite(p_onetail_mimick), NA_real_, p_onetail_mimick),
    p_onetail_rev    = ifelse(is.infinite(p_onetail_rev), NA_real_, p_onetail_rev)
  )


entity_info <- gsea_long %>%
  distinct(intervention, intervention_dataset, type_raw) %>%
  group_by(intervention) %>%
  mutate(n_intervention_datasets = n_distinct(intervention_dataset)) %>%
  ungroup() %>%
  mutate(
    intervention_entity = paste(intervention, intervention_dataset, sep = " | "),
    
    display_label = ifelse(
      n_intervention_datasets > 1,
      paste0(intervention, "\n", intervention_dataset),
      intervention
    ),
    
    intervention_group = case_when(
      !is.na(type_raw) &
        type_raw %in% c("Psychoplastogen", "Positive Control", "Negative Control") ~ type_raw,
      
      intervention %in% c("EE (Old)", "EE (Young)", "Exercise") ~ "Positive Control",
      intervention %in% c("Alcohol", "Rat aging") ~ "Negative Control",
      
      TRUE ~ "Psychoplastogen"
    ),
    
    intervention_group = factor(
      intervention_group,
      levels = c("Psychoplastogen", "Positive Control", "Negative Control")
    )
  ) %>%
  select(
    intervention,
    intervention_dataset,
    intervention_entity,
    display_label,
    intervention_group
  )

gsea_long <- gsea_long %>%
  left_join(
    entity_info,
    by = c("intervention", "intervention_dataset")
  )

gsea_wide <- gsea_long %>%
  select(
    aging_dataset,
    intervention,
    intervention_dataset,
    intervention_entity,
    display_label,
    intervention_group,
    direction,
    NES
  ) %>%
  pivot_wider(
    names_from = direction,
    values_from = NES,
    names_prefix = "NES_"
  ) %>%
  filter(!is.na(NES_UP), !is.na(NES_DN))

all_entities <- entity_info %>%
  arrange(intervention, intervention_dataset) %>%
  pull(intervention_entity)

# ============================================================
# 1. Scatterplot per intervention + intervention_dataset
# x = NES UP, y = NES DN
# ============================================================

plot_scatter_one_entity <- function(entity_id) {
  
  dat <- gsea_wide %>%
    filter(intervention_entity == entity_id)
  
  plot_title <- unique(dat$display_label)
  plot_title <- paste(plot_title, collapse = " / ")
  
  ggplot(
    dat,
    aes(x = NES_UP, y = NES_DN)
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
      aes(label = aging_dataset, colour = intervention_group),
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
      title = paste0(plot_title, ": NES UP vs NES DOWN"),
      x = "NES UP",
      y = "NES DOWN"
    ) +
    theme_age_reversal()
}

scatter_pdf <- file.path(
  out_dir,
  "scatter_NES_UP_vs_DOWN_by_intervention_dataset.pdf"
)

pdf(scatter_pdf, width = 7, height = 6)

for (entity_id in all_entities) {
  
  p <- plot_scatter_one_entity(entity_id)
  print(p)
  
  dat_entity <- entity_info %>%
    filter(intervention_entity == entity_id) %>%
    slice(1)
  
  out_name <- paste0(
    "scatter_",
    safe_name(dat_entity$intervention),
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
# 2. Boxplot NES UP
#    Box separated by intervention + intervention_dataset
#    Ordered by mean NES UP
# ============================================================

gsea_up <- gsea_long %>%
  filter(direction == "UP")

up_order <- gsea_up %>%
  group_by(intervention_entity, display_label) %>%
  summarise(
    mean_NES = mean(NES, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_NES)

gsea_up <- gsea_up %>%
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
  gsea_up,
  aes(
    x = intervention_entity_ord,
    y = NES,
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
    title = "Distribution of NES UP across intervention datasets",
    x = "Intervention / intervention dataset ordered by mean NES UP",
    y = "NES UP"
  ) +
  theme_age_reversal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

p_box_up

ggsave(
  filename = file.path(out_dir, "boxplot_NES_UP_by_intervention_dataset_ordered_by_mean.png"),
  plot = p_box_up,
  width = 14,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "boxplot_NES_UP_by_intervention_dataset_ordered_by_mean.pdf"),
  plot = p_box_up,
  width = 14,
  height = 6
)

# ============================================================
# 3. Boxplot NES DOWN
#    Box separated by intervention + intervention_dataset
#    Ordered by mean NES DOWN
# ============================================================

gsea_dn <- gsea_long %>%
  filter(direction == "DN")

dn_order <- gsea_dn %>%
  group_by(intervention_entity, display_label) %>%
  summarise(
    mean_NES = mean(NES, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_NES)

gsea_dn <- gsea_dn %>%
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
  gsea_dn,
  aes(
    x = intervention_entity_ord,
    y = NES,
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
    title = "Distribution of NES DOWN across intervention datasets",
    x = "Intervention / intervention dataset ordered by mean NES DOWN",
    y = "NES DOWN"
  ) +
  theme_age_reversal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

p_box_dn

ggsave(
  filename = file.path(out_dir, "boxplot_NES_DOWN_by_intervention_dataset_ordered_by_mean.png"),
  plot = p_box_dn,
  width = 14,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "boxplot_NES_DOWN_by_intervention_dataset_ordered_by_mean.pdf"),
  plot = p_box_dn,
  width = 14,
  height = 6
)

# ============================================================
# 4. Combined boxplot UP / DOWN
#    Box separated by intervention + intervention_dataset
#    Ordered by global mean NES
# ============================================================

combined_order <- gsea_long %>%
  group_by(intervention_entity, display_label) %>%
  summarise(
    mean_NES_global = mean(NES, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_NES_global)

gsea_combined <- gsea_long %>%
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
  gsea_combined,
  aes(
    x = intervention_entity_ord,
    y = NES,
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
    title = "Distribution of NES UP and DOWN across intervention datasets",
    x = "Intervention / intervention dataset ordered by global mean NES",
    y = "NES"
  ) +
  theme_age_reversal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top"
  )

p_box_combined

ggsave(
  filename = file.path(out_dir, "boxplot_NES_UP_DOWN_by_intervention_dataset_ordered_by_global_mean.png"),
  plot = p_box_combined,
  width = 16,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "boxplot_NES_UP_DOWN_by_intervention_dataset_ordered_by_global_mean.pdf"),
  plot = p_box_combined,
  width = 16,
  height = 6
)

# ============================================================
# 5. Scatterplot of mean NES
#    One point per intervention + intervention_dataset
#    Legend aggregated by intervention_group
# ============================================================

mean_nes <- gsea_wide %>%
  group_by(
    intervention,
    intervention_dataset,
    intervention_entity,
    display_label,
    intervention_group
  ) %>%
  summarise(
    mean_NES_UP = mean(NES_UP, na.rm = TRUE),
    mean_NES_DN = mean(NES_DN, na.rm = TRUE),
    n_aging_datasets = n(),
    .groups = "drop"
  )

p_mean_scatter <- ggplot(
  mean_nes,
  aes(
    x = mean_NES_UP,
    y = mean_NES_DN,
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
    title = "Mean NES per intervention dataset",
    x = "Mean NES UP",
    y = "Mean NES DOWN"
  ) +
  theme_age_reversal(base_size = 11)

p_mean_scatter

ggsave(
  filename = file.path(out_dir, "scatter_mean_NES_UP_vs_DOWN_by_intervention_dataset.png"),
  plot = p_mean_scatter,
  width = 7,
  height = 5.5,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "scatter_mean_NES_UP_vs_DOWN_by_intervention_dataset.pdf"),
  plot = p_mean_scatter,
  width = 7,
  height = 5.5
)

# ============================================================
# 6. Optional APR-like scatterplot
#    One point per intervention + intervention_dataset
#
# APR = -log10(p_reversal) - [-log10(p_mimicking)]
#     = log10(p_mimicking / p_reversal)
#
# Positive APR = more reversal-like
# Negative APR = more mimicking-like
# ============================================================

if (all(c("p_onetail_mimick", "p_onetail_rev") %in% colnames(raw))) {
  
  apr_long <- gsea_long %>%
    filter(
      !is.na(p_onetail_mimick),
      !is.na(p_onetail_rev),
      p_onetail_mimick > 0,
      p_onetail_rev > 0
    ) %>%
    mutate(
      p_onetail_mimick = pmax(p_onetail_mimick, 1e-300),
      p_onetail_rev    = pmax(p_onetail_rev, 1e-300),
      APR = -log10(p_onetail_rev) - (-log10(p_onetail_mimick))
    )
  
  apr_entity <- apr_long %>%
    group_by(
      intervention,
      intervention_dataset,
      intervention_entity,
      display_label,
      intervention_group,
      direction
    ) %>%
    summarise(
      APR = sum(APR, na.rm = TRUE),
      n_aging_datasets = n(),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = direction,
      values_from = APR,
      names_prefix = "APR_"
    ) %>%
    filter(!is.na(APR_UP), !is.na(APR_DN))
  
  p_apr <- ggplot(
    apr_entity,
    aes(
      x = APR_UP,
      y = APR_DN,
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
      title = "Aging-reversal p ratio per intervention dataset",
      x = "APR, up-regulated genes",
      y = "APR, down-regulated genes"
    ) +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = "age-reversal",
      hjust = 1.05,
      vjust = -0.8,
      size = 3
    ) +
    annotate(
      "text",
      x = -Inf,
      y = -Inf,
      label = "age-mimicking",
      hjust = -0.05,
      vjust = -0.8,
      size = 3
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = "age-reversal",
      hjust = -0.05,
      vjust = 1.3,
      angle = 90,
      size = 3
    ) +
    annotate(
      "text",
      x = -Inf,
      y = -Inf,
      label = "age-mimicking",
      hjust = -0.05,
      vjust = -0.3,
      angle = 90,
      size = 3
    ) +
    coord_cartesian(clip = "off") +
    theme_age_reversal()
  
  p_apr
  
  ggsave(
    filename = file.path(out_dir, "scatter_APR_UP_vs_DOWN_by_intervention_dataset.png"),
    plot = p_apr,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(out_dir, "scatter_APR_UP_vs_DOWN_by_intervention_dataset.pdf"),
    plot = p_apr,
    width = 7,
    height = 6
  )
  
  write.csv(
    apr_entity,
    file = file.path(out_dir, "APR_by_intervention_dataset.csv"),
    row.names = FALSE
  )
}

# ============================================================
# Save processed tables
# ============================================================

write.csv(
  gsea_long,
  file = file.path(out_dir, "GSEA_long_processed_by_intervention_dataset.csv"),
  row.names = FALSE
)

write.csv(
  gsea_wide,
  file = file.path(out_dir, "GSEA_wide_NES_UP_DN_by_intervention_dataset.csv"),
  row.names = FALSE
)

write.csv(
  mean_nes,
  file = file.path(out_dir, "mean_NES_by_intervention_dataset.csv"),
  row.names = FALSE
)

write.csv(
  entity_info,
  file = file.path(out_dir, "intervention_dataset_entity_info.csv"),
  row.names = FALSE
)

message("Done. Figures saved in: ", normalizePath(out_dir))