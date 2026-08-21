library(readxl)
library(dplyr)

samples <- read_excel("terrapin_pond_survey_full.xlsx", sheet = "Samples")
str(samples$Terrapin_count)
samples <- samples %>%
  mutate(
    Terrapin_count = as.numeric(Terrapin_count)
  )
samples <- samples %>%
  mutate(
    Terrapins_seen = trimws(Terrapins_seen) 
  )
# extracting edna sampling visit
ves_day <- samples %>%
  filter(Visit_number == 1.0) %>%  
  group_by(Pond_ID) %>%
  summarise(
    ves_count = sum(Terrapin_count, na.rm = TRUE),
    ves_detected = if_else(any(Terrapins_seen == "Yes"), 1L, 0L),
    .groups = "drop"
  )

#edna
# reading barcodes

library(readr)
library(stringr)
library(purrr)

barcode_files <- list.files("C:/Users/nisar/OneDrive - Queen Mary, University of London/MSC Project/Sequencing Database/My data/barcode_files", pattern = "_COI_counts.tabular", full.names = TRUE)
metadata <- read_excel("terrapin_pond_survey_full.xlsx", sheet = "metadata")
extract_barcode <- function(fname) {
  basename(fname) |> str_replace("_COI_counts.tabular$", "")
}

all_counts <- map_dfr(barcode_files, ~ {
  bc <- read_tsv(.x, col_names = c("sseqid", "read_count"))
  bc$barcode <- extract_barcode(.x)
  bc
})
str(all_counts)
str(metadata)
metadata <- metadata %>%
  rename(barcode = Barcode) 
# adding pond and matrix info
all_counts <- all_counts %>%
  left_join(metadata, by = "barcode")  # metadata: barcode, Pond_ID, matrix

# terrapin genera

terrapin_counts <- all_counts %>%
  mutate(
    species_label = str_replace(sseqid, "_COI.*$", ""),
    Genus = str_split_fixed(species_label, "_", 2)[, 1]
  ) %>%
  filter(Genus %in% c("Trachemys", "Graptemys", "Malaclemys", "Pseudemys", 
                      "Chelydra", "Sternotherus", "Emys", "Mauremys", 
                      "Pelomedusa", "Apalone", "Macrochelys"))  
# detection

# water eDNA
water_terrapin <- terrapin_counts %>%
  filter(Matrix == "Water") %>%
  group_by(Pond_ID) %>%
  summarise(
    water_reads = sum(read_count),
    water_detected = if_else(water_reads > 0, 1L, 0L),
    .groups = "drop"
  )

# sediment eDNA (already combined by barcode)
sediment_terrapin <- terrapin_counts %>%
  filter(Matrix == "Sediment") %>%
  group_by(Pond_ID) %>%
  summarise(
    sediment_reads = sum(read_count),
    sediment_detected = if_else(sediment_reads > 0, 1L, 0L),
    .groups = "drop"
  )



detection_compare <- ves_day %>%
  left_join(water_terrapin, by = "Pond_ID") %>%
  left_join(sediment_terrapin, by = "Pond_ID")

# does edna detect terrapins where ves did not?

summary_compare <- detection_compare %>%
  summarise(
    ponds_total = n(),
    ponds_ves = sum(ves_detected),
    ponds_water_edna = sum(water_detected, na.rm = TRUE),
    ponds_sed_edna = sum(sediment_detected, na.rm = TRUE),
    ponds_both_water = sum(ves_detected == 1 & water_detected == 1, na.rm = TRUE),
    ponds_both_sed = sum(ves_detected == 1 & sediment_detected == 1, na.rm = TRUE),
    ponds_ves_only = sum(ves_detected == 1 & (water_detected == 0 & sediment_detected == 0), na.rm = TRUE),
    ponds_edna_only = sum((water_detected == 1 | sediment_detected == 1) & ves_detected == 0, na.rm = TRUE)
  )

#visualise

library(ggplot2)

ggplot(detection_compare, aes(x = Pond_ID, y = water_reads, fill = factor(ves_detected))) +
  geom_col() +
  labs(fill = "VES detected", y = "Water eDNA terrapin reads")

ggplot(detection_compare, aes(x = Pond_ID, y = sediment_reads, fill = factor(ves_detected))) +
  geom_col() +
  labs(fill = "VES detected", y = "Sediment eDNA terrapin reads")


# trying diff thresholds

thresholds <- c(1, 5, 10, 20, 50)

for (t in thresholds) {
  cat("\nThreshold:", t, "\n")
  print(
    detection_compare %>%
      mutate(
        water_call = if_else(water_reads >= t, 1L, 0L),
        sediment_call = if_else(sediment_reads >= t, 1L, 0L)
      ) %>%
      summarise(
        ponds_ves = sum(ves_detected),
        ponds_water_edna = sum(water_call),
        ponds_sed_edna = sum(sediment_call),
        ponds_both_water = sum(ves_detected == 1 & water_call == 1),
        ponds_both_sed = sum(ves_detected == 1 & sediment_call == 1),
        ponds_ves_only = sum(ves_detected == 1 & (water_call == 0 & sediment_call == 0)),
        ponds_edna_only = sum((water_call == 1 | sediment_call == 1) & ves_detected == 0)
      )
  )
}

library(dplyr)
library(purrr)

thresholds <- c(1, 5, 10, 20, 50)

threshold_summary <- map_dfr(thresholds, function(t) {
  detection_compare %>%
    mutate(
      water_call = if_else(water_reads >= t, 1L, 0L),
      sediment_call = if_else(sediment_reads >= t, 1L, 0L)
    ) %>%
    summarise(
      threshold = t,
      ponds_total = n(),
      ponds_ves = sum(ves_detected),
      ponds_water_edna = sum(water_call),
      ponds_sed_edna = sum(sediment_call),
      ponds_both_water = sum(ves_detected == 1 & water_call == 1),
      ponds_both_sed = sum(ves_detected == 1 & sediment_call == 1),
      ponds_ves_only = sum(ves_detected == 1 & (water_call == 0 & sediment_call == 0)),
      ponds_edna_only = sum((water_call == 1 | sediment_call == 1) & ves_detected == 0)
    )
})

threshold_summary

summary(detection_compare$water_reads)
summary(detection_compare$sediment_reads)


#cleaning the blast database.

terrapin_genus_summary <- terrapin_counts %>%
  group_by(Pond_ID, Matrix, Genus) %>%
  summarise(
    total_reads = sum(read_count),
    .groups = "drop"
  )

top_genera <- terrapin_genus_summary %>%
  group_by(Pond_ID) %>%
  slice_max(order_by = total_reads, n = 2, with_ties = FALSE) %>%
  ungroup()
#interpretation and summary

library(dplyr)
library(tidyr)

genus_presence <- terrapin_genus_summary %>%
  mutate(
    detected = if_else(total_reads > 10, 1L, 0L)  
  ) %>%
  group_by(Pond_ID, Genus) %>%        # combine water + sediment
  summarise(
    detected_any = if_else(sum(detected) > 0, 1L, 0L),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Genus,
    values_from = detected_any,
    values_fill = 0
  )
#compareing with ves
compare_genus_ves <- ves_day %>%
  left_join(genus_presence, by = "Pond_ID")



#building a ref summaary to mitigate bias

library(readxl)
library(dplyr)
library(stringr)

ref_db <- read_excel("terrapin_pond_survey_full.xlsx", sheet = "Reference database")

# Extract Genus from the Species name 
ref_db <- ref_db %>%
  mutate(
    Genus = word(Species, 1)
  )
ref_db <- ref_db %>%
  mutate(
    `Total Number of sequences` = as.numeric(`Total Number of sequences`)
  )
# Summarise number of reference sequences per genus
ref_summary <- ref_db %>%
  filter(`Taxon Groups` == "Terrapins") %>%  # focus on turtle genera
  group_by(Genus) %>%
  summarise(
    n_refs = sum(`Total Number of sequences`, na.rm = TRUE),
    .groups = "drop"
  )

# terrapin_genus_summary: Pond_ID, matrix, Genus, total_reads
genus_norm <- terrapin_genus_summary %>%
  left_join(ref_summary, by = "Genus") %>%
  mutate(
    n_refs = if_else(is.na(n_refs), 1, n_refs),  # fallback in case any genus is missing
    reads_per_ref = total_reads / n_refs
  )
genus_relative_norm <- genus_norm %>%
  group_by(Pond_ID, Matrix) %>%
  mutate(
    total_norm_reads = sum(reads_per_ref),
    rel_abundance_norm = reads_per_ref / total_norm_reads
  ) %>%
  ungroup()

genus_relative_norm %>% filter(Pond_ID %in% c("HG", "HGP1"))


library(tidyr)

genus_presence_norm <- genus_relative_norm %>%
  mutate(
    detected = if_else(rel_abundance_norm >= 0.05, 1L, 0L)
  ) %>%
  group_by(Pond_ID, Genus) %>%       # combine water + sediment
  summarise(
    detected_any = if_else(sum(detected) > 0, 1L, 0L),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Genus,
    values_from = detected_any,
    values_fill = 0
  )

compare_genus_ves_norm <- ves_day %>%
  left_join(genus_presence_norm, by = "Pond_ID")

compare_genus_ves_norm

final_table <- compare_genus_ves_norm %>%
  select(Pond_ID, ves_detected, ves_count,
         Trachemys, Graptemys, Pseudemys, Sternotherus, Emys)  # adjust to genera of interest


trachemys_summary <- compare_genus_ves_norm %>%
  summarise(
    ponds_total = n(),
    ponds_both_trachemys = sum(ves_detected == 1 & Trachemys == 1, na.rm = TRUE),
    ponds_edna_only_trachemys = sum(ves_detected == 0 & Trachemys == 1, na.rm = TRUE),
    ponds_ves_only_trachemys = sum(ves_detected == 1 & Trachemys == 0, na.rm = TRUE)
  )

graptemys_summary <- compare_genus_ves_norm %>%
  summarise(
    ponds_total = n(),
    ponds_both_graptemys = sum(ves_detected == 1 & Graptemys == 1, na.rm = TRUE),
    ponds_edna_only_graptemys = sum(ves_detected == 0 & Graptemys == 1, na.rm = TRUE),
    ponds_ves_only_graptemys = sum(ves_detected == 1 & Graptemys == 0, na.rm = TRUE)
  )


#making graphs
#graph 1: ves vs edna

library(ggplot2)
library(dplyr)

detection_compare <- detection_compare %>%
  mutate(
    water_reads = as.numeric(water_reads),
    sediment_reads = as.numeric(sediment_reads),
    total_edna_reads = water_reads + sediment_reads
  )

ggplot(detection_compare, aes(x = Pond_ID, y = total_edna_reads)) +
  geom_col(fill = "steelblue") +
  labs(
    x = "Pond",
    y = "Total terrapin eDNA reads",
    title = "Total eDNA reads per pond"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(detection_compare, aes(x = ves_count, y = total_edna_reads, label = Pond_ID)) +
  geom_point() +
  geom_text(nudge_y = 0.05 * max(detection_compare$total_edna_reads, na.rm = TRUE), size = 3) +
  labs(
    x = "Terrapins counted by VES",
    y = "Total terrapin eDNA reads (water + sediment)",
    title = "Relationship between VES terrapin counts and eDNA reads per pond"
  )

library(dplyr)
library(ggplot2)

detection_compare <- detection_compare %>%
  mutate(
    water_reads = as.numeric(water_reads),
    sediment_reads = as.numeric(sediment_reads),
    total_edna_reads = water_reads + sediment_reads
  )

# 2. Plot VES counts vs total eDNA reads with colours, jitter, and labels
ggplot(detection_compare,
       aes(x = ves_count,
           y = total_edna_reads,
           label = Pond_ID,
           colour = Pond_ID)) +
  geom_point(position = position_jitter(width = 0.1, height = 0),
             size = 3) +
  geom_text(position = position_jitter(width = 0.1, height = 0),
            nudge_y = 0.03 * max(detection_compare$total_edna_reads, na.rm = TRUE),
            size = 3,
            show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  labs(
    x = "Terrapins counted by visual encounter survey (VES)",
    y = "Total terrapin eDNA reads (water + sediment)",
    title = "Relationship between VES terrapin counts and total eDNA reads per pond"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )
#water edna and sediment edna vs ves
ggplot(detection_compare, aes(x = Pond_ID)) +
  geom_col(aes(y = water_reads), fill = "steelblue") +
  geom_point(aes(y = ves_count * 1000), colour = "red", size = 3) +
  labs(
    x = "Pond",
    y = "Water eDNA terrapin reads (bars) and VES count × 1000 (red points)",
    title = "Water eDNA vs VES per pond"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(detection_compare, aes(x = Pond_ID)) +
  geom_col(aes(y = sediment_reads), fill = "steelblue") +
  geom_point(aes(y = ves_count * 1000), colour = "red", size = 3) +
  labs(
    x = "Pond",
    y = "Sediment eDNA terrapin reads (bars) and VES count × 1000 (red points)",
    title = "Sediment eDNA vs VES per pond"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#genus heatmap
# use genus_relative_norm: Pond_ID, matrix, Genus, rel_abundance_norm

ggplot(genus_relative_norm, aes(x = Pond_ID, y = Genus, fill = rel_abundance_norm)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Normalised relative\nabundance") +
  labs(
    x = "Pond",
    y = "Genus",
    title = "Normalised terrapin genus abundances per pond (eDNA)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#by matrix
ggplot(genus_relative_norm, aes(x = Pond_ID, y = Genus, fill = rel_abundance_norm)) +
  geom_tile() +
  facet_wrap(~ Matrix) +
  scale_fill_viridis_c(name = "Normalised relative\nabundance") +
  labs(
    x = "Pond",
    y = "Genus",
    title = "Normalised terrapin genus abundances in water vs sediment"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#presence absence dot chart
library(dplyr)
library(tidyr)
library(ggplot2)


compare_genus_ves_manual <- compare_genus_ves_norm %>%
  mutate(
    # Manual VES genus presence based on your field notes
    Trachemys_VES = case_when(
      Pond_ID == "HG"    ~ 1L,
      Pond_ID == "HGP1"  ~ 1L,
      TRUE               ~ 0L
    ),
    Graptemys_VES = case_when(
      Pond_ID == "HG"    ~ 1L,
      Pond_ID == "HGP1"  ~ 0L,
      TRUE               ~ 0L
    )
  )
compare_long <- compare_genus_ves_manual %>%
  select(Pond_ID,
         Trachemys_VES, Graptemys_VES,   # VES genus presence
         Trachemys, Graptemys) %>%       # eDNA genus presence (0/1)
  pivot_longer(
    cols = c(Trachemys, Graptemys),
    names_to = "Genus",
    values_to = "edna_present"
  ) %>%
  mutate(
    # Match VES genus presence to the genus column
    ves_present_raw = if_else(
      Genus == "Trachemys", Trachemys_VES, Graptemys_VES
    ),
    ves_present = ves_present_raw == 1,
    edna_present = edna_present == 1,
    ves_present = factor(ves_present, levels = c(FALSE, TRUE),
                         labels = c("VES absent", "VES present")),
    edna_present = factor(edna_present, levels = c(FALSE, TRUE),
                          labels = c("eDNA absent", "eDNA present"))
  )

ggplot(compare_long, aes(x = Pond_ID, y = Genus)) +
  geom_point(aes(shape = ves_present, colour = edna_present), size = 5) +
  scale_shape_manual(values = c(1, 19), name = "VES") +
  scale_colour_manual(values = c("grey70", "darkgreen"), name = "eDNA") +
  labs(
    x = "Pond",
    y = "Genus",
    title = "Presence by VES and eDNA for Trachemys and Graptemys"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

#barplots
genus_relative_norm %>%
  filter(Pond_ID %in% c("HG", "HGP1")) %>%
  ggplot(aes(x = Genus, y = rel_abundance_norm, fill = Pond_ID)) +
  geom_col(position = "dodge") +
  labs(
    x = "Genus",
    y = "Normalised relative abundance",
    title = "Comparison of terrapin genera between HG and HGP1 (eDNA)"
  )


# relationship of ves detections and eDNA dectations by genus

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

ves_genus <- tibble(
  Pond_ID = c("HG", "HGP1", "HP2", "VP"),
  Trachemys_ves_count = c(6, 3, 0, 0),
  Graptemys_ves_count = c(1, 0, 0, 0)
)

ves_genus_long <- ves_genus %>%
  pivot_longer(
    cols = c(Trachemys_ves_count, Graptemys_ves_count),
    names_to = "Genus",
    values_to = "ves_genus_count"
  ) %>%
  mutate(
    Genus = recode(
      Genus,
      Trachemys_ves_count = "Trachemys",
      Graptemys_ves_count = "Graptemys"
    )
  )

edna_genus_total_norm <- genus_norm %>%
  filter(Genus %in% c("Trachemys", "Graptemys")) %>%
  group_by(Pond_ID, Genus) %>%
  summarise(
    normalised_edna_signal = sum(reads_per_ref, na.rm = TRUE),
    .groups = "drop"
  )
figure2_data <- ves_genus_long %>%
  left_join(edna_genus_total_norm, by = c("Pond_ID", "Genus")) %>%
  mutate(
    normalised_edna_signal = replace_na(normalised_edna_signal, 0),
    log_signal = normalised_edna_signal + 1
  )
figure2_data
ggplot(
  figure2_data,
  aes(
    x = ves_genus_count,
    y = log_signal,
    colour = Pond_ID,
    label = Pond_ID
  )
) +
  geom_point(size = 3.8) +
  geom_text_repel(
    size = 3.5,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  facet_wrap(~ Genus, scales = "free_y") +
  scale_y_log10() +
  labs(
    x = "Individuals observed by VES",
    y = "Normalised eDNA signal: reads per reference (water + sediment, log10 scale)",
    colour = "Pond",
    title = "Relationship between genus-level VES counts and normalised eDNA signal",
    subtitle = "eDNA reads were normalised by the number of COI reference sequences per genus"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )


#numeric version of the heatmap

library(dplyr)
library(ggplot2)
library(forcats)

bar_data_matrix <- genus_relative_norm %>%
  filter(rel_abundance_norm > 0) %>%
  mutate(
    Pond_ID = factor(Pond_ID, levels = c("HG", "HGP1", "HP2", "VP"))
  )

ggplot(
  bar_data_matrix,
  aes(
    x = Genus,
    y = rel_abundance_norm,
    fill = Genus
  )
) +
  geom_col(show.legend = FALSE) +
  geom_text(
    aes(label = round(rel_abundance_norm, 3)),
    hjust = -0.1,
    size = 2.8
  ) +
  coord_flip(clip = "off") +
  facet_grid(Matrix ~ Pond_ID, scales = "free_y") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    x = "Genus",
    y = "Normalised relative abundance",
    title = "Normalised terrapin genus abundances by pond and eDNA matrix",
    subtitle = "Bar values are identical to the matrix-separated heatmap"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(face = "italic")
  )
