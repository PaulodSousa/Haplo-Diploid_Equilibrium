#' Compute per-window reference allele frequencies across populations
#'
#' Iterates over each population defined in \code{pop.file}, splits the
#' genotype data by contig, and slides a fixed-size window along each contig to
#' compute the reference allele frequency within that window. Both diploid
#' genotypes (\code{"0/0"}, \code{"0/1"}, \code{"1/1"}) and haploid genotypes
#' (\code{"0"}, \code{"1"}) are recognised, making the calculation agnostic to
#' ploidy. The resulting per-window frequencies are the direct input expected
#' by [pairwise.fst()].
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
#'     \item{Freq.A}{Frequency of the reference allele in the window, computed
#'       as \code{(2*N_AA + N_Aa + N_A) / (2*N_dip + N_hap)}.}
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
#' # af <- allele.freq.WS(geno.data  = gt,
#' #                       pop.file   = pop.file,
#' #                       contigs    = contigs,
#' #                       positions  = pos,
#' #                       window.size = 10000)
#'
#' @seealso [pairwise.fst()] for computing Fst from the output of this
#'   function.
#'
#' @export
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


#' Compute pairwise Fst for all population pairs
#'
#' Takes the per-window allele frequency table produced by [allele.freq.WS()]
#' and computes Hudson's Fst for every unique pair of populations. Within each
#' contig-window, Fst is estimated as the ratio of the among-population
#' variance in reference allele frequency to its expected maximum under
#' panmixia:
#' \code{Fst = Var(p) / (mean(p) * (1 - mean(p)))}.
#' Windows in which the mean reference allele frequency is 0 or 1 (i.e.
#' monomorphic across the pair) are set to 0.
#'
#' @param allele.freq.table A [data.table::data.table] produced by
#'   [allele.freq.WS()], containing at minimum the columns \code{Pop},
#'   \code{Contig}, \code{Window_starts}, \code{Window_ends}, \code{N_sites},
#'   and \code{Freq.A}.
#'
#' @return A [data.table::data.table] with one row per population-pair-contig-
#'   window combination and the following columns:
#'   \describe{
#'     \item{Contig}{Contig (chromosome) name.}
#'     \item{window_lims}{A character string of the form
#'       \code{"Window_starts - Window_ends"} identifying the window.}
#'     \item{Sum.Sites}{Total number of called genotype entries summed across
#'       both populations in the window.}
#'     \item{Mean.p}{Mean reference allele frequency across the two
#'       populations in the window.}
#'     \item{Var.p}{Population variance of the reference allele frequency
#'       across the two populations in the window.}
#'     \item{Fst}{Estimated Fst for the window.}
#'     \item{Pop_pair}{A character string of the form \code{"Pop1 - Pop2"}
#'       identifying the population pair.}
#'   }
#'
#' @examples
#' # Assuming allele.freq.WS() has already been run:
#' # af  <- allele.freq.WS(...)
#' # fst <- pairwise.fst(af)
#'
#' @seealso [allele.freq.WS()] for computing the input allele frequency table,
#'   and [summarize_fst()] for computing weighted genome-wide summary
#'   statistics from the output.
#'
#' @export
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
    # Get the population pair
    pops <- c(pop.pairs[pair, 1], pop.pairs[pair, 2])
    cat("Processing population pair:", pops$Pop1, "-", pops$Pop2, "\n")
    
    # subset data for desired populations
    allele.freqs_pop <- allele.freq.table %>% filter(Pop %in% c(pops))
    
    allele.freqs_pop$window_lims <- paste0(allele.freqs_pop$Window_starts, " - ", allele.freqs_pop$Window_ends)
    # calculate A allele mean and variance by contig and window
    Fst.by.window <- allele.freqs_pop %>%
      group_by(Contig, window_lims) %>% # group by contig and window
      summarise(
        Sum.Sites = sum(N_sites), # number of sites across the same window of both populations
        Mean.p = mean(Freq.A), # mean frequency of A allele
        Var.p = mean((Freq.A - mean(Freq.A))^2), # # population variance in frequency of the A allele
        .groups = "drop"
      )
    # calculate Fst by contig and window
    Fst.by.window <- Fst.by.window %>%
      mutate(Fst = ifelse(Mean.p == 0 | Mean.p == 1, 0, Var.p / (Mean.p * (1 - Mean.p))))
    
    
    Fst.by.window$Pop_pair <- paste0(pops[1], " - ", pops[2])
    
    results[[pair]] <- Fst.by.window
    
    # Remove unnecessary objects and clean R environment
    rm(pops, allele.freqs_pop, Fst.by.window)
    gc(verbose = FALSE)
  }
  
  fst_table <- rbindlist(results)
  return(fst_table)
}


#' Summarize genome-wide Fst per population pair
#'
#' Computes the site-count-weighted mean and standard deviation of Fst across
#' all windows for each population pair, using the per-window Fst table
#' produced by [pairwise.fst()]. Weighting by \code{Sum.Sites} ensures that
#' windows with more called genotypes contribute more to the estimate.
#'
#' @param fst_table A [data.table::data.table] produced by [pairwise.fst()],
#'   containing at minimum the columns \code{Pop_pair}, \code{Fst}, and
#'   \code{Sum.Sites}.
#'
#' @return A [tibble::tibble] with one row per population pair and the
#'   following columns:
#'   \describe{
#'     \item{Pop_pair}{A character string identifying the population pair.}
#'     \item{wMean.Fst}{Weighted mean of Fst across all windows.}
#'     \item{wSD.Fst}{Weighted standard deviation of Fst across all windows.}
#'   }
#'
#' @examples
#' # Assuming pairwise.fst() has already been run:
#' # fst     <- pairwise.fst(af)
#' # summary <- summarize_fst(fst)
#'
#' @seealso [pairwise.fst()] for computing the input per-window Fst table.
#'
#' @export
summarize_fst <- function(fst_table) {
  summary_fst <- fst_table %>% group_by(Pop_pair) %>% summarise(
  # weighted mean and sd of genetic diversity
  wMean.Fst = matrixStats::weighted.mean(Fst, Sum.Sites, na.rm = TRUE),
  wSD.Fst = matrixStats::weightedSd(Fst, Sum.Sites, na.rm =TRUE),
  )  
}
