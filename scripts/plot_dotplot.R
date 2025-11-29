# Figure 2E - Bibio_marci vs Ptychoptera_contaminata dotplot 
################################################################################
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('gsheet'))
suppressPackageStartupMessages(library('argparse'))

################################################################################
# root <- "/Users/ab66/Documents/sanger_work/diptera/diptera-ALGs/"
#root <- paste0(getwd(), "/")
################################################################################
# color palette
# pal <- c("M1" = "#1573afff", "M2" = "#e59d38ff", "M3" = "#f0e354ff", 
#          "M4" = "#169e73ff", "M5" = "#60b5e1ff", "M6" = "black", "unassigned" = "grey")
source('scripts/20250620_colour_pal.R')
pal <- c(pal, "unassigned" = "grey")
################################################################################


parser <- ArgumentParser()
parser$add_argument("-s1", 
    help="The name of the species 1 (x-axis)", default = "")
parser$add_argument("-s2", 
    help="The name of the species 2 (y-axis)", default = "")
parser$add_argument("-o", "--output",  
    dest="o", help="Base of the output name (.png will be attached)")
parser$add_argument("-a", "--alg-set",  
    help="File with BUSCO to ALG assignments", default = "data/diptera.no_plecia.mindist.m165_n1_n2.tsv")

args <- parser$parse_args()

target_species <- c(args$s1, args$s2) #c("Bibio_marci", "Ptychoptera_contaminata")
all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?gid=1940964825#gid=1940964825", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)

target_species_data <- all_genome_data %>% filter(species %in% target_species) %>%
  select(species, accession, chromosome, chromsome_size_b)

target_species_data$accession <- sub("\\.[0-9]", "", target_species_data$accession)

# load BUSCO data
busco_files_list <- paste0("data/busco_tables/", unique(target_species_data$species), ".syngraph.buscos.tsv")
buscos <- bind_rows(lapply(busco_files_list, read.table))
colnames(buscos) <- c("marker", "chromosome", "start", "end")

# read BUSCO to ALG file; "data/diptera.no_plecia.mindist.m165_n1_n2.tsv"
algs <- read.table(args$a, header = TRUE, sep = "\t", col.names = c("marker", "ALG"))

################################################################################
# merge everything with the main table
target_species_data <- left_join(target_species_data, buscos, by = "chromosome")
target_species_data <- left_join(target_species_data, algs, by = "marker")

# remove chromosomes that have no markers
target_species_data <- target_species_data %>% filter(!is.na(target_species_data$marker))
# mark buscos that are not in ALGs as unassigned
target_species_data$ALG[is.na(target_species_data$ALG)] <- "unassigned"

################################################################################
# convert coordinates into linear format
df_lin <- target_species_data %>% arrange(species, chromosome, start)
df_lin$linear_start <- NA
df_lin$linear_end <- NA

# for-loop with species + chromosome aware offset
chr_offset <- 0
current_species <- ""

for (i in 1:nrow(df_lin)) {
  
  if (i == 1 || df_lin$species[i] != df_lin$species[i - 1]) {
    # New species → reset offset
    chr_offset <- 0
    current_species <- df_lin$species[i]
    
  } else if (df_lin$chromosome[i] != df_lin$chromosome[i - 1]) {
    # Same species but new chromosome → increase offset
    chr_offset <- chr_offset + as.integer(df_lin$chromsome_size_b[i - 1])
  }
  
  # Assign linearized positions
  df_lin$linear_start[i] <- chr_offset + df_lin$start[i]
  df_lin$linear_end[i]   <- chr_offset + df_lin$end[i]
}

# remove duplicated markers
duplicates <- df_lin |>
  dplyr::summarise(n = dplyr::n(), .by = c(marker, species)) |>
  dplyr::filter(n > 1L)

# pivot wider including the requested columns
alg_df <- df_lin %>%
  select(marker, ALG) %>%
  distinct()

df_wide <- df_lin %>%
  filter(!marker %in% duplicates$marker) %>%
  select(marker, species, chromosome, linear_start, linear_end) %>%
  pivot_wider(
    id_cols = marker,
    names_from = species,
    values_from = c(chromosome, linear_start, linear_end),
    names_sep = "_") %>%
  
  # merge back ALGs
  left_join(alg_df, by = "marker") %>%
  
  # Optional: move stevens next to marker
  relocate(ALG, .after = marker)

################################################################################
# plotting

sp_x <- target_species[1]
sp_y <- target_species[2]

# Filter chrom info for species X and Y
chr_info_x <- df_lin %>%
  filter(species == sp_x) %>%
  arrange(chromosome) %>%
  distinct(chromosome, chromsome_size_b) %>%
  mutate(cum_end = cumsum(chromsome_size_b))

chr_info_x$order <- 1:length(chr_info_x$chromosome)

chr_info_y <- df_lin %>%
  filter(species == sp_y) %>%
  arrange(chromosome) %>%
  distinct(chromosome, chromsome_size_b) %>%
  mutate(cum_end = cumsum(chromsome_size_b))

chr_info_y$order <- 1:length(chr_info_y$chromosome)

# Chromosome label midpoints
chr_labels_x <- df_lin %>%
  filter(species == sp_x) %>%
  group_by(chromosome) %>%
  summarise(mid = mean(c(min(linear_start), max(linear_end))))

chr_labels_y <- df_lin %>%
  filter(species == sp_y) %>%
  group_by(chromosome) %>%
  summarise(mid = mean(c(min(linear_start), max(linear_end))))

# Build column names dynamically
start_x <- paste0("linear_start_", sp_x)
end_x   <- paste0("linear_end_", sp_x)
start_y <- paste0("linear_start_", sp_y)
end_y   <- paste0("linear_end_", sp_y)

# Buffer for rect thickness
buffer <- 1e5

# Build plot
p <- ggplot(df_wide) + 
  #geom_rect(aes_string(
  #  xmin = paste0(start_x, " - buffer"),
  #  xmax = paste0(end_x,   " + buffer"),
  #  ymin = paste0(start_y, " - buffer"),
  #  ymax = paste0(end_y,   " + buffer"),
  #  fill = "ALG",
  #  colour = "ALG"
  #), size = 1.5) +
  geom_point(aes_string(
    x = paste0(start_x, " - buffer"), 
    y = paste0(start_y, " - buffer"), 
    fill = "ALG", colour = "ALG")) +
  geom_vline(data = chr_info_x, aes(xintercept = cum_end), color = "grey80", linewidth = 0.2) +
  geom_hline(data = chr_info_y, aes(yintercept = cum_end), color = "grey80", linewidth = 0.2) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(
    expand = c(0,0),
    labels = function(x) x / 1e6,
    name = paste0("\n", gsub("_", " ", sp_x), " (Mb)"),
    sec.axis = sec_axis(~., breaks = chr_labels_x$mid, labels = chr_labels_x$chromosome)
  ) +
  scale_y_continuous(
    expand = c(0,0),
    labels = function(y) y / 1e6,
    name = paste0(gsub("_", " ", sp_y), " (Mb)\n"),
    sec.axis = sec_axis(~., breaks = chr_labels_y$mid, labels = chr_labels_y$chromosome)
  ) +
  theme_bw()

  p <- p + theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x.top = element_text(angle=90, size=14),
    axis.text.x.bottom = element_text(size=14),
    axis.title.x.bottom = element_text(size=18),
    axis.text.y.right = element_text(size=14),
    axis.text.y.left = element_text(size=14),
    axis.title.y.left = element_text(size=18),
    legend.position = "none"
  )

# ggsave("figures/two_species_dotplot.svg", plot = p, dpi = 600, height = 10.25, width = 11.58)
ggsave(paste0(args$o, ".png"), plot = p, dpi = 600, height = 10.25, width = 11.58)
