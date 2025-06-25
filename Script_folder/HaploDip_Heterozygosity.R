library(data.table)
library(dplyr)
library(vcfR)
setwd("/home/paulos/PhD/Haplo-Dip_Model/Fake_data/")
# Import vcf
vcf <- read.vcfR("Caenea_FAKE_2contigs_2pops_5indvs.vcf")
# get only the gen# get only the genotypes from vcf file
gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = F)
head(gt_matrix)
# get positions and contigs (for sliding windows)
contig_vector <- vcf@fix[, "CHROM"]
positions <- as.numeric(vcf@fix[, "POS"])

# Get pop file
PopFile <- read.csv("Caenea_PopFile_Fake.txt", sep="\t", header= F, row.names = NULL)
head(PopFile)

# only two columns, one with indv names and other with populations names
PopFile <- PopFile[-c(3:ncol(PopFile))]
colnames(PopFile) <- c("ID", "Pop")

head(PopFile)
# check if names from Popfile match names from vcf
colnames(gt_matrix)==PopFile$ID


compute_allele.freqs_W <- function(geno.data, pop.file, contigs, positions, window.size) {
  
  all_results <- list()
  pops <- unique(pop.file$Pop)
  
  # for each poppulation
  for (pop in pops) {
    # Get the sample names for the populations
    pops_samples <- pop.file$ID[pop.file$Pop == pop]
    gt_matrix_pop <- geno.data[, pops_samples, drop = FALSE]
    
    df <- data.table(
      contig = contigs, 
      pos = positions, 
      gt_matrix
    )
    
    df$order <- seq_len(nrow(df)) # Preserve original row order
    
    results_list <- list()
    
    # for each contig do
    for(contig_name in unique(df$contig)) {
      contig_data <- df[contig == contig_name]
      # Sort by genomic position
      setorder(contig_data, pos)
      # Define window start and end
      contig_start <- min(contig_data$pos)
      contig_end <- max(contig_data$pos)
      window_starts <- seq(contig_start, contig_end, by = window.size)
      
      # for each window
      for (start_pos in window_starts) {
        end_pos <- start_pos + window.size - 1
        # Get rows in this window
        window_rows <- contig_data[pos >= start_pos & pos <= end_pos]
        
        if (nrow(window_rows) == 0) next  # Skip empty windows
        
        # Extract only genotype columns (excluding contig, pos, order)
        gt_window <- as.matrix(window_rows[, -(1:3), with = FALSE])
        
        # Flatten the genotype matrix for counting
        gt_flat <- as.vector(gt_window)
        
        # Fun part
        # Lets calculate summary stats
        #gt_window <- t(as.matrix(geno.data))
        
        # 1st step: count each genotype in the window
        AA.f <- sum(gt_flat %in% "0/0", na.rm= T)
        Aa.f <- sum(gt_flat %in% c("0/1", "1/0"), na.rm= T)
        aa.f <- sum(gt_flat %in% "1/1", na.rm= T)
        A.m <- sum(gt_flat %in% "0", na.rm= T)
        a.m <- sum(gt_flat %in% "1", na.rm= T)
        
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
        results_window <- data.table(
          Pop = pop, # which population
          Contig = contig_name, # which contig
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
    all_results[[pop]] <- rbindlist(results_list)
  }
  # Combine across populations
  final_output <- rbindlist(all_results)
  return(final_output)
}





df <- compute_allele.freqs_W(geno.data = gt_matrix, 
                             pop.file = PopFile, 
                             contigs = contig_vector, 
                             positions = positions, 
                             window.size = 10000)
head(df)
write.csv(df , "Diversity_Stats/complete_table.csv")

df_summary <- as.data.frame(df %>% group_by(Pop) %>% summarise(Exp.Het_Mean= mean(Exp.Het, na.rm = T),
                                                               Exp.Het_SD= sd(Exp.Het, na.rm = T),
                                                               Obs.Het_Mean= mean(Obs.Het, na.rm = T),
                                                               Obs.Het_SD= sd(Obs.Het, na.rm = T),
                                                               Exp.M.Ref_Mean = mean(Exp.M.Ref, na.rm = T),
                                                               Exp.M.Ref_SD = sd(Exp.M.Ref, na.rm = T),
                                                               Obs.M.Ref_Mean = mean(Obs.M.Ref, na.rm =T),
                                                               Obs.M.Ref_SD = sd(Obs.M.Ref, na.rm =T)))
df_summary

write.csv(df_summary, "Diversity_Stats/summary_table.csv")


df_summary %>% summarise(Exp.Het.mean = mean(Exp.Het_Mean), Exp.Het.sd = sd(Exp.Het_Mean),
                         Exp.Het.min = min(Exp.Het_Mean), Exp.Het.max = max(Exp.Het_Mean),
                         Obs.Het.mean = mean(Obs.Het_Mean), Obs.Het.sd = sd(Obs.Het_Mean),
                         Obs.Het.min = min(Obs.Het_Mean), Obs.Het.max = max(Obs.Het_Mean),
                         Exp.M.Ref.mean = mean(Exp.M.Ref_Mean), Exp.M.Ref.sd = sd(Exp.M.Ref_Mean),
                         Exp.M.Ref.min = min(Exp.M.Ref_Mean), Exp.M.Ref.max = max(Exp.M.Ref_Mean),
                         Obs.M.Ref.mean = mean(Obs.M.Ref_Mean), Obs.M.Ref.sd = sd(Obs.M.Ref_Mean),
                         Obs.M.Ref.min = min(Obs.M.Ref_Mean), Obs.M.Ref.max = max(Obs.M.Ref_Mean))

He.adj <- df[,c(1,6:7,19)]

He.adj <- as.data.frame(He.adj %>% mutate(Indv.1 = 2*N_females + N_males,
                                      Indv.2 = (2*N_females + N_males)-1 ))
He.adj <- He.adj %>%  mutate(uHe = Exp.Het * (Indv.1/Indv.2) )

He.adj %>% group_by(Pop) %>% summarise(Mean.uHe = mean(uHe, na.rm = T), SD.uHe = sd(uHe, na.rm = T),
                                       Min.uHe = min(uHe, na.rm = T), Max.uHe = max(uHe, na.rm = T),)




