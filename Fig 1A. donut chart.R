library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)

df <- Rejection@meta.data %>%
  as.data.frame() %>%
  dplyr::filter(
    !is.na(label),
    !is.na(group)
  ) %>%
  dplyr::mutate(
    label = as.character(label),
    group = as.character(group)
  ) %>%
  dplyr::count(
    label,
    group,
    name = "n"
  )

df <- df %>%
  dplyr::mutate(
    group = dplyr::recode(
      as.character(group),
      "Active_AMR" = "Active AMR",
      "CA_AMR" = "Chronic active AMR",
      "CA_AMR_borderline" = "Chronic active AMR with borderline",
      "ATN" = "ATN",
      "Borderline" = "Borderline for acute TCMR",
      "TCMR" = "TCMR"
    )
  )

group_order <- c(
  "ATN",
  "Borderline for acute TCMR",
  "TCMR",
  "Active AMR",
  "Chronic active AMR",
  "Chronic active AMR with borderline"
)

group_df <- df %>%
  group_by(group) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  mutate(
    group = factor(group, levels = group_order)
  ) %>%
  arrange(group) %>%
  mutate(
    end = cumsum(total),
    start = end - total,
    mid = (start + end) / 2
  )

label_df <- df %>%
  dplyr::left_join(
    group_df %>%
      dplyr::select(group, group_start = start),
    by = "group"
  ) %>%
  dplyr::mutate(
    group = factor(group, levels = group_order)
  ) %>%
  dplyr::arrange(group, label) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(
    end = min(group_start) + cumsum(n),
    start = end - n,
    mid = (start + end) / 2
  ) %>%
  dplyr::ungroup()

group_colors <- c(
  "ATN" = "chartreuse4",
  "Borderline for acute TCMR" = "skyblue",
  "TCMR" = "dodgerblue",
  "Active AMR" = "#E1B0FF",
  "Chronic active AMR" = "#DA6FD5",
  "Chronic active AMR with borderline" = "#801B80"
)

p <- ggplot() +
  geom_rect(
    data = label_df,
    aes(
      xmin = 0.8,
      xmax = 1.8,
      ymin = start,
      ymax = end,
      fill = group
    ),
    color = "white",
    linewidth = 0.5
  ) +
  geom_rect(
    data = group_df,
    aes(
      xmin = 1.8,
      xmax = 2.8,
      ymin = start,
      ymax = end,
      fill = group
    ),
    color = "white",
    linewidth = 0.8
  ) +
  geom_text(
    data = label_df,
    aes(
      x = 1.3,
      y = mid,
      label = gsub("_", "-", as.character(label))
    ),
    size = 3.5
  ) +
  geom_text(
    data = group_df,
    aes(
      x = 2.3,
      y = mid,
      label = dplyr::if_else(
        as.character(group) == "Chronic active AMR with borderline",
        "Chronic active AMR\nwith borderline",
        as.character(group)
      )
    ),
    size = 4,
    fontface = "bold"
  ) +  
  annotate(
    "text",
    x = 0,
    y = 0,
    label = "Xenium 5k",
    size = 4.5,
    fontface = "bold",
    color = "black"
  ) +  
  scale_fill_manual(
    values = group_colors,
    drop = FALSE
  ) +
  coord_polar(
    theta = "y",
    start = -pi/2
  ) +
  xlim(0, 3) +
  theme_void() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    legend.position = "right"
  ) 
print(p)
