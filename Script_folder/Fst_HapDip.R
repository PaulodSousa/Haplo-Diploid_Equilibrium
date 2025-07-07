################################### F-statistics
library(dplyr)
library(data.table)
library(vcfR)

## Fst

setwd("/home/paulos/PhD/WGS/Caenea/")
# Import vcf
vcf <- read.vcfR("Pruning/Caenea_99_Filtered_Pruned.vcf")
# get only the gen# get only the genotypes from vcf file
gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = F)
head(gt_matrix)

# get positions and contigs (for sliding windows)
contig_vector <- vcf@fix[, "CHROM"]
positions <- as.numeric(vcf@fix[, "POS"])

remove(vcf)

# Get pop file
PopFile <- read.csv("Caenea_PopFile99.txt", sep="\t", header= F)
head(PopFile)

# only two columns, one with indv names and other with populations names
PopFile <- PopFile[-c(3:ncol(PopFile))]
colnames(PopFile) <- c("ID", "Pop")

head(PopFile)
# check if names from Popfile match names from vcf
colnames(gt_matrix)==PopFile$ID

allele.freq.WS <- function(geno.data, pop.file, contigs, positions, window.size) {
  all_results <- list()
  
  # 1st section:  identified all possible population pairs
  pop.list <- unique(pop.file$Pop)
  
  # for each population
  for (pop in pop.list) {
    
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
        AA.f <- sum(gt_flat %in% "0/0", na.rm= T)
        Aa.f <- sum(gt_flat %in% c("0/1", "1/0"), na.rm= T)
        aa.f <- sum(gt_flat %in% "1/1", na.rm= T)
        A.m <- sum(gt_flat %in% "0", na.rm= T)
        a.m <- sum(gt_flat %in% "1", na.rm= T)
        
        # 2nd step: Calculate the females, males and total number of sites in the window
        N.F <- AA.f + aa.f + Aa.f
        N.M <- A.m + a.m
        Total_Samples <- N.F + N.M
        
        # 3rd step: Calculate allele frequencies in the window
        F.A <- (2*AA.f + Aa.f + A.m) / (N.F * 2 + N.M)
        F.a <- (2*aa.f + Aa.f + a.m) / (N.F * 2 + N.M)
        
        results_window <- data.table(
          Pop = pop, # which population
          Contig = contig_name, # which contig
          Window_starts = start_pos, # starting position of window
          Window_ends = end_pos, # ending position of window 
          N_sites = Total_Samples, # total number of samples in the window
          Freq.A = F.A # frequency of reference allele on window
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


alle.freq.df <- allele.freq.WS(geno.data = gt_matrix, pop.file = PopFile,
                               contigs = contig_vector, positions = positions,
                               window.size = 10000)


alle.freq.df

pairwise.fst <- function(allele.freq.table) {
  results <- list()
  
  # 1st section:  identified all possible population pairs
  pop.list <- unique(allele.freq.table$Pop)
  # generate all unique pairwise combinations as a matrix (2 rows, n columns)
  pairs <- combn(pop.list, 2)
  
  # convert to data.table with two columns
  pop.pairs <- data.table(
    Pop1 = pairs[1, ],
    Pop2 = pairs[2, ]
  )
  
  
  # for each population pair
  for (pair in 1:nrow(pop.pairs)) {
    cat("Processing population pair:", pair, "\n")
    # Get the population pair
    pops <- c(pop.pairs[pair, 1], pop.pairs[pair, 2])
    
    # subset data for desired populations
    allele.freqs_pop <- allele.freq.table %>% filter(Pop %in% c(pops))
    
    allele.freqs_pop$window_lims <- paste0(allele.freqs_pop$Window_starts, " - ", allele.freqs_pop$Window_ends)
    # calculate mean and population variance for the allele frequency, by window
    variables <- allele.freqs_pop %>% group_by(Contig, window_lims) %>% 
      summarise(Mean.p = mean(Freq.A), 
                Var.p = mean((Freq.A - mean(Freq.A))^2))
    
    # calculate the fst 
    Fst.by.window <- variables %>% group_by(Contig, window_lims) %>%
      summarise(Fst = ifelse(Mean.p == 0 | Mean.p == 1, 0, Var.p / (Mean.p * (1 - Mean.p))),
                .groups = "drop"
                )
    
    Fst.by.window$Pop_pair <- paste0(pops[1], " - ", pops[2])
    
    results[[pair]] <- Fst.by.window
    
    # Remove unnecessary objects and clean R environment
    rm(pops, allele.freqs_pop, variables, Fst.by.window)
    gc(verbose = FALSE)
  }
  
  final_output <- rbindlist(results)
  return(final_output)
}

fst <- pairwise.fst(alle.freq.df)

fst.summary <- fst %>% group_by(Pop_pair) %>% summarise(Mean = mean(Fst, na.rm= T),
                                                        SD = sd(Fst, na.rm= T))

