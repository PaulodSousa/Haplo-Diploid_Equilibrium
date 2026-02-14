#' Compute per-window Nei's H (gene diversity)
#'
#' Iterates over each population defined in \code{pop.file}, splits the
#' genotype data by contig, and slides a fixed-size window along each contig
#' to compute Nei's H (probability of sampling two different alleles) within
#' that window. Nei's H is calculated as \code{2pq}, where \code{p} and
#' \code{q} are the reference and alternative allele frequencies respectively.
#' Both diploid genotypes (\code{"0/0"}, \code{"0/1"}, \code{"1/1"}) and
#' haploid genotypes (\code{"0"}, \code{"1"}) are recognised when computing
#' allele frequencies. Despite the different ploidies, allele frequencies
#' should be the same between sexes, which means that Nei's H 
#' is agnostic to ploidy. 
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
#'     \item{Neis_H}{Nei's H (gene diversity) for the window, computed as
#'       \code{2 * Freq.Ref * Freq.Alt}.}
#'   }
#'
#' @examples
#' vcf_path <- system.file("extdata",
#'                         "example.vcf",
#'                         package = "HaploDiploidEquilibrium")
#'
#' result <- vcf2GT(vcf_path)
#' gt       <- result$gt_matrix
#' contigs  <- result$contig_vector
#' pos      <- result$positions
#' 
#' pop.file <- data.frame(ID  = colnames(gt),
#'                        Pop = c("PopA","PopA","PopB","PopB","PopB"))
#' 
#' hs <- compute_Hs_W(geno.data   = gt,
#'                    pop.file    = pop.file,
#'                    contigs     = contigs,
#'                    positions   = pos,
#'                    window.size = 10000)
#'
#' @seealso [summarize_NeisH()] for computing weighted genome-wide summary
#'   statistics from the output.
#'
#' @export
compute_Hs_W <- function(geno.data, pop.file, contigs, positions, window.size, verbose=TRUE) {
  
  all_results <- list()
  pops <- unique(pop.file$Pop)
  
  # for each population
  for (pop in pops) {
    
    if(verbose){
        cat("Processing population:", pop, "\n")
    }
    
    # Get the sample names for the populations
    pops_samples <- pop.file$ID[pop.file$Pop == pop]
    gt_matrix_pop <- geno.data[, pops_samples, drop = FALSE]
    
    df <- data.table::data.table(
      contig = contigs, 
      pos = positions, 
      gt_matrix_pop
    )
    
    df$order <- seq_len(nrow(df)) # Preserve original row order
    
    pop_results <- list()
    
    # for each contig do
    for(contig_name in unique(df$contig)) {
      if(verbose){
          cat("  Contig:", contig_name, "\n")
      }
      
      contig_data <- subset(df, contig == contig_name)
      
      # Sort by genomic position
      data.table::setorder(contig_data, pos)
      
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
        window_rows <- subset(contig_data, pos >= start_pos & pos <= end_pos)
        
        if (nrow(window_rows) == 0) next  # Skip empty windows
        
        # Extract only genotype columns (excluding contig, pos, order)
        gt_window <- as.matrix(window_rows[, setdiff(names(window_rows),
                                                     c("contig", "pos", "order")), with = FALSE])
        
        
        # Flatten the genotype matrix for counting
        gt_flat <- as.vector(gt_window)
        
        # Fun part
        # Lets calculate gene diveristy (Hs)
        
        # 1st step: count each genotype in the window
        AA.f <- sum(gt_flat %in% c("0/0", "0|0"), na.rm= T)
        Aa.f <- sum(gt_flat %in% c("0/1", "1/0", "1|0", "0|1"), na.rm= T)
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
        
        # 4th step: Calculate Hs in the window
        H <- 2* F.A * F.a
        
        # save values here
        results_window <- data.table::data.table(
          Pop = pop, # which population
          Contig = contig_name, # which contig
          Window_starts = start_pos, # starting position of window
          Window_ends = end_pos, # ending position of window 
          N_sites = Total_Samples, # total number of samples in the window
          Neis_H = H # Nei's H on the window
        )
        results_list[[idx]] <- results_window
        idx <- idx + 1
        
        # Clear window-specific objects
        rm(gt_window, gt_flat, window_rows)
        #gc(verbose = FALSE)
      }
      contig_results <- data.table::rbindlist(results_list[!sapply(results_list, is.null)])
      pop_results[[length(pop_results) + 1]] <- contig_results
      
      # Clear contig data
      rm(contig_data)
      gc(verbose = FALSE)
    }
    all_results[[pop]] <- data.table::rbindlist(pop_results)
    # Clear population data
    rm(gt_matrix_pop, df, results_list)
    gc(verbose = FALSE)
  }
  
  final_output <- data.table::rbindlist(all_results)
  return(final_output)
}


#' Summarize per-window Nei's H per population
#'
#' Computes the site-count-weighted mean and standard deviation of Nei's H
#' across all windows for each population, using the per-window table produced
#' by [compute_Hs_W()]. Weighting by \code{N_sites} ensures that windows with
#' more called genotypes contribute more to the estimate.
#'
#' @param neis_table A [data.table::data.table] produced by [compute_Hs_W()],
#'   containing at minimum the columns \code{Pop}, \code{N_sites}, and
#'   \code{Neis_H}.
#'
#' @return A [tibble::tibble] with one row per population and the following
#'   columns:
#'   \describe{
#'     \item{Pop}{Population label.}
#'     \item{wMean.Neis_H}{Weighted mean of Nei's H across all windows.}
#'     \item{wSD.Neis_H}{Weighted standard deviation of Nei's H across all
#'       windows.}
#'   }
#'
#' @examples
#' vcf_path <- system.file("extdata",
#'                         "example.vcf",
#'                         package = "HaploDiploidEquilibrium")
#'
#' result <- vcf2GT(vcf_path)
#' gt       <- result$gt_matrix
#' contigs  <- result$contig_vector
#' pos      <- result$positions
#' 
#' pop.file <- data.frame(ID  = colnames(gt),
#'                        Pop = c("PopA","PopA","PopB","PopB","PopB"))
#' 
#' hs <- compute_Hs_W(geno.data   = gt,
#'                    pop.file    = pop.file,
#'                    contigs     = contigs,
#'                    positions   = pos,
#'                    window.size = 10000)
#' summary <- summarize_NeisH(hs)
#'
#' @seealso [compute_Hs_W()] for computing the input per-window table.
#'
#' @export
summarize_NeisH <- function(neis_table) {
  summary_df <- neis_table |> dplyr::group_by(Pop) |> dplyr::summarise(
  # weighted mean and sd of genetic diversity
  wMean.Neis_H = stats::weighted.mean(Neis_H, N_sites, na.rm = TRUE),
  wSD.Neis_H = matrixStats::weightedSd(Neis_H, N_sites, na.rm = TRUE)
  )
  return(summary_df)
}
