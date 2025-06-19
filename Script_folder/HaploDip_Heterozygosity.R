library(dplyr)
library(vcfR)
setwd("/home/paulos/PhD/Haplo-Dip_Model/")
# Import vcf
vcf <- read.vcfR("Fake_data/Caenea_FAKE_2contigs_2pops_5indvs.vcf")
# get only the genotypes from vcf file
gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = F)
head(gt_matrix)
# get positions and contigs (for sliding windows)
contig_vector <- vcf@fix[, "CHROM"]
positions <- as.numeric(vcf@fix[, "POS"])

# Get pop file
PopFile <- read.csv("Fake_data/Caenea_PopFile_Fake.txt", sep="\t", header= F)
colnames(PopFile) <- c("ID", "Pop")
head(PopFile)
# check if names from Popfile match names from vcf
colnames(gt_matrix)==PopFile$ID


# because it runs within window sizes it sums up the number of samples, males, females and genotypes of
## all loci within the window size (instead of doing an average)
## This means that if in population A, there are 1 male and two females a window (x1) that has two loci
## loci x1.1 has genotypes A, AA and AA and loci x1.2 has genotypes A, NA and Aa, 
## the number of samples in the window x1 is 5, the number of males is 2, the number of females is 3, 
## and the genotype count is A = 2, AA = 2 and Aa = 1. This way, NAs are accounted for.
compute_allele.freqs_SW <- function(geno.data, pop.file, contigs, positions, window.size, step.size) {
  results_list <- list()
  
  pops <- unique(pop.file$Pop)
  
  # for each poppulation do
  for (pop in pops) {
    # Get the sample names for the populations
    pops_samples <- pop.file$ID[pop.file$Pop == pop]
    gt_matrix_pop <- geno.data[, pops_samples, drop = FALSE]
    
    # Ensure everything is ordered the same
    df <- data.frame(contig = contigs, position = positions, stringsAsFactors = FALSE)
    df$order <- seq_along(positions)
    
    # Split by contig
    contig_names <- unique(df$contig)
    # for each contig do
    for (contig in contig_names) {
      contig_df <- df[df$contig == contig, ]
      contig_df <- contig_df[order(contig_df$position), ]
      # subset matrix for specific contig
      gt_sub <- gt_matrix_pop[contig_df$order, , drop = FALSE]
      pos_sub <- contig_df$position
      # get minimum and maximum postion of the contig
      min_pos <- min(pos_sub)
      max_pos <- max(pos_sub)
      
      # skip windows that are too small
      if ((max_pos - min_pos + 1) < window.size) next
      
      # window start and end positions
      windows_starts <- seq(min_pos, max_pos - window.size + 1, by = step.size)
      windows_ends <- windows_starts + window.size - 1
      
      # for each window do...
      for (i in seq_along(windows_starts)) {
        start_pos <- windows_starts[i]
        end_pos <- windows_ends[i]
        # get the loci within the window interval
        loci_in_window <- which(pos_sub >= start_pos & pos_sub <= end_pos)
        # if no loci within the window skip
        if (length(loci_in_window) < 1) next
        # subset gt_sub matrix for each loci within the window
        gt_window <- gt_sub[loci_in_window, , drop = FALSE]
        
        # Fun part
        # Lets calculate summary stats
        #gt_window <- t(as.matrix(geno.data))
        
        # 1st step: count each genotype in the window
        AA.f <- sum(gt_window %in% "0/0", na.rm= T)
        Aa.f <- sum(gt_window %in% c("0/1", "1/0"), na.rm= T)
        aa.f <- sum(gt_window %in% "1/1", na.rm= T)
        A.m <- sum(gt_window %in% "0", na.rm= T)
        a.m <- sum(gt_window %in% "1", na.rm= T)
        
        # 2nd step: Calculate number of females and males and total sample number in the window
        N.F <- AA.f + aa.f + Aa.f
        N.M <- A.m + a.m
        Total_Samples <- N.F + N.M
        
        # 3rd step: Calculate allele frequencies in the window
        F.A <- (2*AA.f + Aa.f + A.m) / (N.F * 2 + N.M)
        F.a <- (2*aa.f + Aa.f + a.m) / (N.F * 2 + N.M)
        
        # 4th step: Calculate expected genotype frequencies (assuming an equal sex ratio) in the window
        Exp.AA <- (F.A * F.A) * 0.5
        Exp.aa <- (F.a * F.a) * 0.5
        Exp.Aa <- 2 * (F.A * F.a) * 0.5
        Exp.A <- F.A  * 0.5
        Exp.a <- F.a  * 0.5
        
        # save values here
        results_window <- data.frame(
          Pop = pop, # which population
          Contig = contig, # which contig
          Window_starts = start_pos, # starting position of window
          Window_ends = end_pos, # ending position of window 
          N_samples = Total_Samples, # total number of samples in the window
          N_females = N.F, # total number of females in the window
          N_males = N.M, # total number of males in the window
          N_AA = AA.f, # total number of AA genotypes in the window
          N_Aa = Aa.f, # total number of Aa genotypes in the window
          N_aa = aa.f, # total number of aa genotypes in the window
          N_A = A.m, # total number of A genotypes in the window
          N_a = a.m, # total number of a genotypes in the window
          Freq.Ref = F.A, # frequency of the reference allele in the window
          Freq.Alt = F.a, # frequency of the alternative allele in the window
          Obs.Hom = (AA.f + aa.f) / Total_Samples, # observed frequency of homozygous in the window
          Obs.Het = Aa.f / Total_Samples, # observed frequency of heterozygous in the window
          Obs.M.Ref = A.m / Total_Samples, # observed frequency of males with reference allele in the window
          Exp.Hom = (Exp.AA + Exp.aa), # expected frequency of the homozygous in the window
          Exp.Het = Exp.Aa, # expected frequency of the heterozygous in the window
          Exp.M.Ref = Exp.A # expected frequency of males with reference allele in the window
        )
        
        results_list[[length(results_list) +1]] <- results_window
      }
    }
  }
  final_output <- do.call(rbind, results_list)
  return(final_output)
}

df <- compute_allele.freqs_SW(geno.data = gt_matrix, pop.file = PopFile, 
                              contigs = contig_vector, positions = positions,
                              window.size = 1000, step.size = 500)
head(df)


# This function summarises the output of the compute_allele.freqs_SW()
## the population average of each summary statistics is done by summing the values of that statistic of
## all windows and divided by the product of the number of windows and the summ of samples of all windows:
## Avg_He.i = summ(He.i) / (N.i x number of windows of i) where i refers to the population id
## It seems more logic to me than just do an average, since doing just and average of He per population
## would mean the average He over all windows not considering the number of samples.  

summary.He <- function(He_window_df) {
  # save here
  summary_results <- list()
  # unique populations
  pops <- unique(He_window_df$Pop)
  
  # for each population do
  for (pop in pops) {
    # split data.frame by pop
    He_df_pop <- split(He_window_df, He_window_df$Pop)
    
    N <- sum(He_df_pop[[pop]]$N_samples) # summ number of samples across all windows
    Males <- sum(He_df_pop[[pop]]$N_males) # summ number of males across all windows
    Windows <- nrow(He_df_pop[[pop]]) # summ number of windows
    
    
   
    results <- data_frame(
      Pop = pop, # which population
      N_Windows = Windows, # number of windows
      Avg_N.Samples = N / Windows, # average number of samples per window
      Avg_N.Males = Males / Windows, # average number of males per window
      Avg_Obs.Hom = sum(He_df_pop[[pop]]$Obs.Hom) / (N * Windows), # average number of observed homozygotes per window
      Avg_Obs.Het = sum(He_df_pop[[pop]]$Obs.Het) / (N * Windows), # average number of observed heterozygous per window
      Avg_Obs.M.Ref = sum(He_df_pop[[pop]]$Obs.M.Ref) / (N * Windows), # average number of observed reference allele males per window
      Avg_Exp.Hom = sum(He_df_pop[[pop]]$Exp.Hom) / (N * Windows), # average number of expected homozygotes per window
      Avg_Exp.Het = sum(He_df_pop[[pop]]$Exp.Het) / (N * Windows), # average number of expected heterozygous per window
      Avg_Exp.M.Ref = sum(He_df_pop[[pop]]$Exp.M.Ref) / (N * Windows) # average number of expected reference allele males per window
    )
    summary_results[[length(summary_results) +1]] <- results
  }
  summary_results_df <- as.data.frame(do.call(rbind, summary_results))
  return(summary_results_df)
}

summary.He(He_window_df = df)
