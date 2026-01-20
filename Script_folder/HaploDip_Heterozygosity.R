library(data.table)
library(dplyr)
library(vcfR)
setwd("/home/paulos/PhD/Haplo-Dip_Model/")
# Import vcf
vcf <- read.vcfR("SLiM/n500_2Sex.vcf")
# get only the gen# get only the genotypes from vcf file
gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = F)
head(gt_matrix)

# get positions and contigs (for sliding windows)
contig_vector <- vcf@fix[, "CHROM"]
positions <- as.numeric(vcf@fix[, "POS"])

remove(vcf)

# Get pop file
PopFile <- read.csv("SLiM/PopFile.csv", 
                    sep="\t", header= T)
head(PopFile)

# only two columns, one with indv names and other with populations names
#PopFile <- PopFile[-c(3:ncol(PopFile))]
colnames(PopFile) <- c("ID", "Pop")
PopFile <- PopFile[c(1:20),]

head(PopFile)
# check if names from Popfile match names from vcf
colnames(gt_matrix)==PopFile$ID


compute_allele.freqs_W <- function(geno.data, pop.file, contigs, positions, window.size, dip_freq) {
  
  all_results <- list()
  pops <- unique(pop.file$Pop)
  
  # for each population
  for (pop in pops) {
    
    cat("Processing population:", pop, "\n")
    
    # Get the sample names for the populations
    pops_samples <- pop.file$ID[pop.file$Pop == pop]
    gt_matrix_pop <- geno.data[, pops_samples, drop = FALSE]
    
    df <- data.table(
      contig = contigs, 
      pos = positions, 
      gt_matrix_pop
    )
    
    df$order <- seq_len(nrow(df)) # Preserve original row order
    
    pop_results <- list()
    
    # for each contig do
    for(contig_name in unique(df$contig)) {
      cat("  Contig:", contig_name, "\n")
      
      contig_data <- df[contig == contig_name]
      
      # Sort by genomic position
      setorder(contig_data, pos)
      
      # Define window start and end
      contig_start <- min(contig_data$pos)
      contig_end <- max(contig_data$pos)
      window_starts <- seq(contig_start, contig_end, by = window.size)
      
      results_list <- vector("list", length(window_starts))
      idx <- 1
      
      # for each window
      for (start_pos in window_starts) {
        end_pos <- start_pos + window.size - 1
        # Get rows in this window
        window_rows <- contig_data[pos >= start_pos & pos <= end_pos]
        
        if (nrow(window_rows) == 0) next  # Skip empty windows
        
        # Extract only genotype columns (excluding contig, pos, order)
        gt_window <- as.matrix(window_rows[, setdiff(names(window_rows),
                                                     c("contig", "pos", "order")), with = FALSE])
        
        
        # Flatten the genotype matrix for counting
        gt_flat <- as.vector(gt_window)
        
        # Fun part
        # Lets calculate summary stats
        
        # 1st step: count each genotype in the window
        AA.f <- sum(gt_flat %in% c("0/0", "0|0"), na.rm= T)
        Aa.f <- sum(gt_flat %in% c("0/1", "1/0", "0|1", "1|0"), na.rm= T)
        aa.f <- sum(gt_flat %in% c("1/1", "1|1"), na.rm= T)
        A.m <- sum(gt_flat %in% "0", na.rm= T)
        a.m <- sum(gt_flat %in% "1", na.rm= T)
        
        # 2nd step: Calculate the females, males and total number of sites in the window
        N.F <- AA.f + aa.f + Aa.f
        N.M <- A.m + a.m
        Total_Samples <- N.F + N.M
        
        # 3rd step: Calculate allele frequencies in the window
        F.A <- (2*AA.f + Aa.f + A.m) / (N.F * 2 + N.M)
        F.a <- (2*aa.f + Aa.f + a.m) / (N.F * 2 + N.M)
        
        # 4th step: Calculate expected genotype frequencies in the window
        # set haploid frequency from diploid (0.5 is highly recomended) 
        hap.freq <- 1 - dip_freq

        Exp.AA <- (F.A * F.A) * dip_freq
        Exp.aa <- (F.a * F.a) * dip_freq
        Exp.Aa <- 2 * (F.A * F.a) * dip_freq
        Exp.A <- F.A  * hap.freq
        Exp.a <- F.a  * hap.freq
        
        # save values here
        results_window <- data.table(
          Pop = pop, # which population
          Contig = contig_name, # which contig
          Window_starts = start_pos, # starting position of window
          Window_ends = end_pos, # ending position of window 
          N_sites = Total_Samples, # total number of samples in the window
          N_F = N.F, # Total number of female samples in the window
          N_M = N.M, # Total number of male samples in the window
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
          Exp.M.Ref = Exp.A, # expected frequency of males with reference allele in the window
          Fis = ifelse(Exp.Aa == 0, NA, 1 - ((Aa.f / Total_Samples) / Exp.Aa)) # Fis in the window
        )
        results_list[[idx]] <- results_window
        idx <- idx + 1
        
        # Clear window-specific objects
        rm(gt_window, gt_flat, window_rows)
        #gc(verbose = FALSE)
      }
      contig_results <- rbindlist(results_list[!sapply(results_list, is.null)])
      pop_results[[length(pop_results) + 1]] <- contig_results
      
      # Clear contig data
      rm(contig_data)
      #gc(verbose = FALSE)
    }
    all_results[[pop]] <- rbindlist(pop_results)
    # Clear population data
    rm(gt_matrix_pop, df, results_list)
    gc(verbose = FALSE)
  }
  
  final_output <- rbindlist(all_results)
  return(final_output)
}


df <- compute_allele.freqs_W(geno.data = gt_matrix, 
                             pop.file = PopFile, 
                             contigs = contig_vector, 
                             positions = positions, 
                             window.size = 10000,
                             dip_freq = 0.5)
head(df)
write.csv(df ,"Data/My_Data/complete_tablePseudoHap_mono_miss1.csv")

# weighted means and standard deviations
library(matrixStats)
summary_df <- df %>% group_by(Pop) %>% summarise(
  # weighted mean and sd of expected heterozygosity
  wMean.Exp.Het = weighted.mean(Exp.Het, N_sites, na.rm = T),
  wSD.Exp.Het = weightedSd(Exp.Het, N_sites, na.rm =T),
  # weighted mean and sd of observed heterozygosity
  wMean.Obs.Het = weighted.mean(Obs.Het, N_sites, na.rm = T),
  wSD.Obs.Het = weightedSd(Obs.Het, N_sites, na.rm =T),
  # weighted mean and sd of expected male ref. allele genotype
  wMean.Exp.M.Ref = weighted.mean(Exp.M.Ref, N_sites, na.rm = T),
  wSD.Exp.M.Ref = weightedSd(Exp.M.Ref, N_sites, na.rm =T),
  # weighted mean and sd of observed male ref. allele genotype
  wMean.Obs.M.Ref = weighted.mean(Obs.M.Ref, N_sites, na.rm = T),
  wSD.Obs.M.Ref = weightedSd(Obs.M.Ref, N_sites, na.rm =T),
)

summary_df %>% summarise(Mean = mean(wMean.Exp.Het), SD = sd(wMean.Exp.Het))

write.csv(summary_df, "Data/My_Data/summary_table_PseudoHap_mono_miss1.csv")






