#' Compute per-window reference allele frequencies by sex
#'
#' Iterates over each population defined in \code{pop.file}, splits the
#' genotype data by contig, and slides a fixed-size window along each contig
#' to compute the reference allele frequency separately for diploid individuals
#' (females), haploid individuals (males), and both sexes combined. Sex is
#' inferred from ploidy: diploid genotypes (\code{"0/0"}, \code{"0/1"},
#' \code{"1/1"}) are assumed to belong to females and haploid genotypes
#' (\code{"0"}, \code{"1"}) to males.
#'
#' @param geno.data A character matrix of genotype strings with dimensions
#'   \code{n_sites x n_individuals}, as returned by [vcf2GT()].
#' @param pop.file A \code{data.frame} or \code{data.table} with at least two
#'   columns: \code{ID} (individual identifiers matching the column names of
#'   \code{geno.data}) and \code{Pop} (population labels).
#' @param contigs A character vector of length \code{n_sites} containing the
#'   contig (chromosome) name for each variant site, as returned by [vcf2GT()].
#' @param positions A numeric vector of length \code{n_sites} containing the
#'   physical position (bp) of each variant site, as returned by [vcf2GT()].
#' @param window.size A single positive integer giving the size of each
#'   sliding window in base pairs.
#'
#' @return A [data.table::data.table] with one row per population-contig-window
#'   combination and the following columns:
#'   \describe{
#'     \item{Pop}{Population label.}
#'     \item{Contig}{Contig (chromosome) name.}
#'     \item{Window_starts}{Genomic coordinate (bp) of the first position in
#'       the window.}
#'     \item{Window_ends}{Genomic coordinate (bp) of the last position in the
#'       window (\code{Window_starts + window.size - 1}).}
#'     \item{N_sites}{Total number of called genotype entries (diploid +
#'       haploid) within the window.}
#'     \item{Females.freq}{Reference allele frequency computed from diploid
#'       genotypes only: \code{(2*N_AA + N_Aa) / (2*N_dip)}.}
#'     \item{Males.freq}{Reference allele frequency computed from haploid
#'       genotypes only: \code{N_A / N_hap}.}
#'     \item{Total.freq}{Reference allele frequency computed from both sexes
#'       combined: \code{(2*N_AA + N_Aa + N_A) / (2*N_dip + N_hap)}.}
#'   }
#'
#' @examples
#' # Assuming vcf2GT() has already been run:
#' # result   <- vcf2GT("path/to/input.vcf")
#' # gt       <- result$gt_matrix
#' # contigs  <- result$contig_vector
#' # pos      <- result$positions
#' #
#' # pop.file <- data.frame(ID  = colnames(gt),
#' #                         Pop = c("PopA","PopA","PopB","PopB"))
#' #
#' # sex_ref <- compute.Female.Male.allele.W(geno.data   = gt,
#' #                                         pop.file    = pop.file,
#' #                                         contigs     = contigs,
#' #                                         positions   = pos,
#' #                                         window.size = 10000)
#'
#' @seealso [summarize_sex_ref()] for computing weighted genome-wide summary
#'   statistics from the output.
#'
#' @export compute.Female.Male.allele.W
compute.Female.Male.allele.W <- function(geno.data, pop.file, contigs, positions, window.size) {
  
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
        
        # 3rd step: Calculate ref. allele frequency for each sex and total in each window
        F.A <- (2*AA.f + Aa.f) / (N.F * 2) # females
        M.A <- (A.m) / (N.M) # males
        T.A <- (2*AA.f + Aa.f + A.m) / (N.F * 2 + N.M) # both sexes (total)
        
        # save values here
        results_window <- data.table(
          Pop = pop, # which population
          Contig = contig_name, # which contig
          Window_starts = start_pos, # starting position of window
          Window_ends = end_pos, # ending position of window 
          N_sites = Total_Samples, # total number of samples in the window
          Females.freq = F.A, # frequency on females
          Males.freq = M.A, # frequency on males
          Total.freq = T.A # frequency on both sexes
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


#' Summarize per-sex reference allele frequencies per population
#'
#' Computes the site-count-weighted mean and standard deviation of the
#' reference allele frequency for females, males, and both sexes combined,
#' across all windows for each population. Uses the per-window table produced
#' by [compute.Female.Male.allele.W()]. Weighting by \code{N_sites} ensures
#' that windows with more called genotypes contribute more to each estimate.
#'
#' @param allele_table A [data.table::data.table] produced by
#'   [compute.Female.Male.allele.W()], containing at minimum the columns
#'   \code{Pop}, \code{N_sites}, \code{Females.freq}, \code{Males.freq}, and
#'   \code{Total.freq}.
#'
#' @return A [tibble::tibble] with one row per population and the following
#'   columns:
#'   \describe{
#'     \item{Pop}{Population label.}
#'     \item{wMean.F.Ref}{Weighted mean of the female reference allele
#'       frequency across all windows.}
#'     \item{wSD.F.Ref}{Weighted standard deviation of the female reference
#'       allele frequency across all windows.}
#'     \item{wMean.M.Ref}{Weighted mean of the male reference allele frequency
#'       across all windows.}
#'     \item{wSD.M.Ref}{Weighted standard deviation of the male reference
#'       allele frequency across all windows.}
#'     \item{wMean.T.Ref}{Weighted mean of the combined reference allele
#'       frequency across all windows.}
#'     \item{wSD.T.Ref}{Weighted standard deviation of the combined reference
#'       allele frequency across all windows.}
#'   }
#'
#' @examples
#' # Assuming compute.Female.Male.allele.W() has already been run:
#' # sex_ref <- compute.Female.Male.allele.W(...)
#' # summary <- summarize_sex_ref(sex_ref)
#'
#' @seealso [compute.Female.Male.allele.W()] for computing the input
#'   per-window table.
#'
#' @export
summarize_sex_ref <- function(allele_table) {
  summary_df <- allele_table %>% group_by(Pop) %>% summarise(
  # females
  wMean.F.Ref = matrixStats::weighted.mean(Females.freq, N_sites, na.rm = TRUE),
  wSD.F.Ref = matrixStats::weightedSd(Females.freq, N_sites, na.rm =TRUE),
  # males
  wMean.M.Ref = matrixStats::weighted.mean(Males.freq, N_sites, na.rm = TRUE),
  wSD.M.Ref = matrixStats::weightedSd(Males.freq, N_sites, na.rm =TRUE),
  # both sexes
  wMean.T.Ref = matrixStats::weighted.mean(Total.freq, N_sites, na.rm = TRUE),
  wSD.T.Ref = matrixStats::weightedSd(Total.freq, N_sites, na.rm =TRUE)
  )
  return(summary_df)
}
