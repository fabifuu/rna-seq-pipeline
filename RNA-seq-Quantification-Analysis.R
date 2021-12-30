# Install BioconductoR
if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install(version = "3.14")

# Install DESeq2 from BioconductoR
BiocManager::install("DESeq2")
# Install rhdf5 from BioconductoR
BiocManager::install("rhdf5")
# Install devtools
install.packages("devtools")
# Install sleuth
devtools::install_github("mschilli87/sleuth@drop-h5write-default-import")

# Set a working directory
setwd("/home/fabifuu/Bioinformatics/RNA-seq")

# Load sleuth and specify kallisto output dir
library("sleuth")
sample_id <- dir(file.path("./kallisto_quant"))
sample_id

# Make a list of paths to the kallisto results
kal_dir <- file.path("./kallisto_quant", sample_id)
kal_dir

# Load an auxiliary table that describes the experimental design
# You can create beforehand using Libreoffice, or manually type in R
# IMPORTANT!
# Sample header is "sample"
# Condition or treatment group header is "condition"
# Replication header is "reps"
exp_des <- read.csv(file.path("./experimental_design.csv"), header = TRUE)
exp_des

# Append kallisto dir with experimental design
exp_des <- dplyr::mutate(exp_des, path = kal_dir)
exp_des

# Construct a sleuth object
# Note: it is advisable to run sleuth from cli to take advantage of multiple core
slob <- sleuth_prep(exp_des, ~ condition, extra_bootstrap_summary = TRUE)

# Fitting a full model
slob <- sleuth_fit(slob)

# Fitting a reduced model
slob <- sleuth_fit(slob, ~1, 'reduced')

# Sleuth likelihood ratio test (LRT)
slob_lrt <- sleuth_lrt(slob, 'reduced', 'full')

# Examine a model
models(slob)
models(slob_lrt)

# Examine the results of slob
# The table shown below displays the top 20 significant genes
# with a Benjamini-Hochberg multiple testing corrected q-value <= 0.05
sleuth_table <- sleuth_results(slob_lrt, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) #only q<=0.05
head(sleuth_significant, 20) #show the first 20

# Export sleuth significant to csv file
write.csv(sleuth_significant, file = "./sleuth_significant.csv")

# Box plot for certain gene (est_results with the most significant q-value)
# The top results is gene ENST00000257570.9
plot_bootstrap(slob_lrt, "ENST00000257570.9", units = "est_counts", color_by = "condition")
plot_bootstrap(slob_lrt, "ENST00000371941.3", units = "est_counts", color_by = "condition")

# Including the gene name
# Install BiomaRt
BiocManager::install("biomaRt")
library("biomaRt")

# Collect the gene names
# Since we use ENSEMBL cDNA Homo sapiens (see "ENST- prefixes),
# we use ENSEMBL data set (URL may varies, consult Biostar)

# Download raw gene name
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'https://m.ensembl.org')

# Filter gene name (attributes)
gename <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                        "external_gene_name"), mart = mart)

# Modify header gene name into sleuth-specified header name
gename <- dplyr::rename(gename, target_id = ensembl_transcript_id,
                        ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# Adding gene name into sleuth object
# slob_gen is sleuth object that have gene name in the output.
# slob construction and model fitting is started all over again.
slob_gen <- sleuth_prep(exp_des, ~ condition, extra_bootstrap_summary = TRUE, target_mapping = gename)
slob_gen <- sleuth_fit(slob_gen, ~ condition, 'full')
slob_gen <- sleuth_fit(slob_gen, ~ 1, 'reduced')
slob_gen <- sleuth_lrt(slob_gen, 'reduced', 'full')

# Examine the results of slob_gen
sleuth_table_gen <- sleuth_results(slob_gen, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_gen <- dplyr::filter(sleuth_table_gen, qval <= 0.05) #only q<=0.05
head(sleuth_significant_gen, 20) #show the first 20

# Export slob_gen significant to csv file
write.csv(sleuth_significant_gen, file = "./sleuth_significant_gen.csv")

# Box plot with certain gene name
