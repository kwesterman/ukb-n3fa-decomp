# Perform initial processing of vQTL and main-effect summary stats
# (formatting, INFO filter, subsetting) and create a basic QQ plot


library(tidyverse)
library(data.table)


args <- commandArgs(trailingOnly=T)
filepath <- args[1]


### Define functions

calc_lambda <- function(x, p=0.5){
  # Calculate genomic inflation lambda value
  x = x[!is.na(x)]
  x.quantile <- quantile(x, p)
  round(qchisq(1 - x.quantile, 1) / qchisq(p, 1), 2)
}

make_qq <- function(data, pval_col, main=""){
  # Make a quantile-quantile plot
  data <- filter(data, data[[pval_col]] > 0)  # In case extremely low p-values are stored as zero
  
  # Process p-values
  y_vals <- sort(-log10(data[[pval_col]]))
  x_vals <- -log10(rev(ppoints(length(y_vals))))  # ppoints generates a uniform probability distribution
  
  # Trim points at higher p-values (credit to RaMWAS package for code snippet)
  levels = as.integer((x_vals - x_vals[1]) / (tail(x_vals, 1) - x_vals[1]) * 2000)
  keep = c(TRUE, diff(levels) != 0)
  levels = as.integer((y_vals - y_vals[1])/(tail(y_vals, 1) - y_vals[1]) * 2000)
  keep = keep | c(TRUE, diff(levels) != 0)
  keep = which(keep)
  
  par(ps=18)
  plot(x=x_vals[keep], y=y_vals[keep], 
       xlab=expression(-log[10](italic(p)) * " (Expected)"), 
       ylab=expression(-log[10](italic(p)) * " (Observed)"),
       main=main, cex=0.8, 
       cex.lab=0.8, cex.main=0.9, 
       pch=16, ylim=c(0, ceiling(max(y_vals))))
  abline(0, 1, lty=2)
  legend(x='topleft', y='topleft',
         bquote(lambda == .(calc_lambda(data[[pval_col]]))), 
         cex=0.9, bty="n")
}


### Read in summary stats and subset to columns of interest

ss_cols <- c(
  CHR="CHR", SNP="RSID", POS="POS",
  EA="Effect_Allele", NEA="Non_Effect_Allele", AF="AF",
  N="N_Samples", P_int="robust_P_Value_Interaction", P_joint="robust_P_Value_Joint"
)

high_qual_variants <- read_tsv("../data/processed/ukb_rsIDs_maf0.005_info0.5.txt", 
                               col_names=F, col_types="c")[[1]]

ss_df <- fread(filepath, stringsAsFactors=F, data.table=F) %>%
  select(all_of(ss_cols), matches("^Beta"), matches("^Var_Beta")) %>%
  mutate(P_int = as.numeric(P_int),
	 P_joint = as.numeric(P_joint)) %>%
  filter(SNP %in% high_qual_variants)
if ("AF" %in% names(ss_df)) ss_df <- filter(ss_df, pmin(AF, 1 - AF) > 0.01)


### Write particular subsets for downstream analysis

ss_df %>%
  filter(P_int < 0.05) %>%
  write_tsv(paste0(filepath, "_nom"))

ss_df %>%
  filter(P_int < 5e-8) %>%
  write_tsv(paste0(filepath, "_gw"))

magma_filepath <- gsub("_merged", "_magmaInput.tsv", filepath)
ss_df %>%
  select(SNP, N, P=P_int) %>%
  write_tsv(magma_filepath)

fuma_filepath <- gsub("_merged", "_fumaInput.tsv", filepath)
ss_df %>%
  select(SNP, CHR, POS, EA, NEA, AF, N, P=P_int) %>%
  write_tsv(fuma_filepath)
system(paste0("gzip -f ", fuma_filepath))


### Create Q-Q plot
  
qq_dir <- paste0(dirname(filepath), "/qq_plots/")
system(paste0("mkdir -p ", qq_dir))
write(calc_lambda(ss_df$P), paste0(qq_dir, gsub("_merged|\\.tbl", "_lambda", basename(filepath))))
plot_filepath <- paste0(qq_dir, gsub("_merged|\\.tbl", "_QQ.pdf", basename(filepath)))
pdf(file=plot_filepath)
make_qq(ss_df, "P_int")
dev.off()
