#Code obtained and modified from https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html

getwd()
setwd()

library("readr")
library("vegan")
library("splitstackshape")
library("dplyr")
library("taxa")
library("ggplot2")
library("agricolae")
library("metacoder")
library("phyloseq")

#Importing the necessary data tables 
otu_data <- read_tsv("feature-table.tsv")
tax_data <- read_tsv("taxonomy.tsv")
#separating the taxa data 
tax_data$Taxon <- concat.split(tax_data, split.col = 'Taxon', sep = ";", logical = TRUE)
tax_data <- as.data.frame(as.matrix(tax_data))
#Removing the second column
tax_data <- tax_data[-c(2)]
#Renaming the columns *Make sure the OTU_ID matches the OTU_ID from the OTU Table
names(tax_data) <- c("OTU_ID", "Taxonomy", "Kingdom", "Phylum",       "Class",           "Order",          "Family",         "Genus",            "Species")

#Joining the OTU Table and Taxonomy Tables based on OTU_ID
tax_data$OTU_ID <- as.character(tax_data$OTU_ID)
otu_data$OTU_ID <- as.character(otu_data$OTU_ID)
otu_data <- left_join(otu_data, tax_data,
                      by = c("OTU_ID" = "OTU_ID"))
tail(colnames(otu_data), n = 10)

#Importing sample data
sample_data <- read_tsv("Manifest.txt",
                        col_types = "cccccc")
obj <- parse_tax_data(otu_data,
                      class_cols = "Taxonomy",
                      class_sep = ";")
print(obj)
print(obj$data$tax_data)

#Converting to taxmap format
head(taxon_names(obj))

obj <- parse_tax_data(otu_data,
                      class_cols = "Taxonomy",
                      class_sep  = ";",
                      class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                      class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
head(taxon_names(obj))
obj$data$class_data
head(taxon_ranks(obj))
obj$data$class_data <- NULL
names(obj$data) <- "otu_counts"
print(obj)
head(taxon_names(obj))

obj <- taxa::filter_taxa(obj, taxon_names != "")

head(taxon_names(obj))
head(all_names(obj), 20)


#Heat Trees
obj <- taxa::filter_taxa(obj, taxon_names  == "Bacteria", subtaxa = TRUE)
print(obj)
obj$data$otu_counts <- obj$data$otu_counts[c("taxon_id", "OTU_ID", sample_data$sampleID)]

heat_tree(obj)
heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = n_obs)

obj$data$tax_abund <- calc_taxon_abund(obj, "otu_counts",
                                       cols = sample_data$sampleID,
                                       groups = sample_data$`symptom-type`)

obj %>%
  taxa::filter_taxa(taxon_ranks == "g", supertaxa = TRUE) %>%
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs,
            output_file = "plotnaive.pdf")

set.seed(2)
obj %>%
  taxa::filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>%
  taxa::filter_taxa(taxon_ranks =="g", supertaxa = TRUE) %>%
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = naive,
            initial_layout = "re", layout = "da",
            title = "Naive Read Depth",
            node_color_axis_label = "Sum of Naive Reads",
            node_size_axis_label =  "Number of OTUs",
            output_file = "plotnaive.pdf")

#Differential Heat Trees
print(obj$data$tax_abund)


obj$data$rel_abd <- calc_obs_props(obj, "otu_counts", other_cols = T)
obj$data$tax_abund <- calc_taxon_abund(obj, "rel_abd")
obj$data$diff_table <- compare_groups(obj, data = "tax_abund",
                                      cols = sample_data$sampleID,
                                      groups  = sample_data$location)
print(obj$data$diff_table)
obj <- mutate_obs(obj, "diff_table",
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))
obj$data$diff_tables$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0

set.seed(1)
heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label  = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportion",
                 layout  = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 output_file = "differential_heat_tree_location.pdf")

