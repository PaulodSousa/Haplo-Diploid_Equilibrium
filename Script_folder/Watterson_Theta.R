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

remove(vcf)

# Get pop file
PopFile <- read.csv("Caenea_PopFile_Fake.txt", sep="\t", header= F)
head(PopFile)

# only three columns, one with indv names, other with populations names and other with sex of sample
colnames(PopFile) <- c("ID", "Pop", "Sex")

head(PopFile)
# check if names from Popfile match names from vcf
colnames(gt_matrix)==PopFile$ID

is.polymorphic <- function(site_genotype) {
  # Remove NAs and keep only valid genotypes
  valid_genotypes <- c("0/0", "0/1", "1/0", "1/1", "0", "1")
  site_genotype.clean <- site_genotype[site_genotype %in% valid_genotypes]
  
  unique_gts <- unique(site_genotype.clean)
  
  # If only one genotype of one sex is present -> not polymorphic
  if (length(unique_gts) == 1) {
    poly <- FALSE
    
    # If exactly two genotypes and both have just the  "A" or "a" alleles 
  } else if (length(unique_gts) == 2 &&
             (all(unique_gts %in% c("0/0", "0")) || all(unique_gts %in% c("1/1", "1")))) {
    poly <- FALSE
    
  } else {
    poly <- TRUE
  }
  
  return(poly)
}

compute_WTheta_W <- function(geno.data, pop.file, contigs, positions, window.size) {
  
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
        
        
        
        # Fun part
        # Apply to each site (row) of the window to get logical vector
        polymorphic_sites <- apply(gt_window, 1, is.polymorphic)
        
        # Count segregating sites
        S <- sum(polymorphic_sites)
        
        # Count number of samples (females and males)
        pop.table <- pop.file %>% group_by(Pop, Sex) %>% summarise(N_Samples = n())
        
        this.pop <- pop.table %>% filter(Pop == pop)
        
        n_f <- 2*this.pop$N[this.pop$Sex == "f"]
        n_m <- this.pop$N[this.pop$Sex == "m"]
        
        n_total <- n_f + n_m
        
        # Harmonic number
        a_n <- sum(1 / (1:(n_total - 1)))
        
        # Watterson's theta
        theta_watterson <- S / a_n
        
        # number of valid sites
        gt_flat <- as.vector(gt_window)
        gt_keep <- gt_flat[gt_flat %in% c("0/0", "0/1", "1/0", "1/1", "1", "0")]
        
        num_sites <- length(gt_keep)
        
        # Per-site theta
        theta_per_site <- theta_watterson / num_sites
           
        results_window <- data.table(
          Pop = pop, # which population
          Contig = contig_name, # which contig
          Window_starts = start_pos, # starting position of window
          Window_ends = end_pos, # ending position of window 
          N_sites = num_sites, # number of valid sites in the window
          Theta = theta_per_site, # Per-site theta of the window
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
      gc(verbose = FALSE)
    }
    all_results[[pop]] <- rbindlist(pop_results)
    # Clear population data
    rm(gt_matrix_pop, df, results_list)
    gc(verbose = FALSE)
  }
  
  final_output <- rbindlist(all_results)
  return(final_output)
}

df <- compute_WTheta_W(geno.data = gt_matrix, pop.file = PopFile, contigs = contig_vector,
                       positions = positions, window.size = 10000)
head(df)
