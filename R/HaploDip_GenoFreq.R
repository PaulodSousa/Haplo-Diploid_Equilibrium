#' Compute per-window genotype frequencies, allele frequencies, and Fis
#'
#' Iterates over each population defined in \code{pop.file}, splits the
#' genotype data by contig, and slides a fixed-size window along each contig
#' to compute observed and expected genotype frequencies, allele frequencies,
#' and the inbreeding coefficient (Fis). Expected genotype frequencies are
#' derived from the haplo-diploid equilibrium model, where the proportion of
#' diploid and haploid individuals in the population is controlled by
#' \code{dip_freq}. Sex is inferred from ploidy: diploid genotypes
#' (\code{"0/0"}, \code{"0/1"}, \code{"1/1"}) are assumed to belong to females
#' and haploid genotypes (\code{"0"}, \code{"1"}) to males. Fis is computed as
#' \code{1 - (Obs.Het / Exp.Het)} and is set to \code{NA} when expected
#' heterozygosity is zero.
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
#' @param dip_freq A single numeric value in the interval \code{(0, 1)}
#'   giving the expected proportion of diploid individuals in the population.
#'   The haploid proportion is set to \code{1 - dip_freq}. A value of
#'   \code{0.5} corresponds to an equal sex ratio and is recommended for
#'   standard haplo-diploid systems.
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
#'     \item{N_F}{Number of diploid (female) genotype entries in the window.}
#'     \item{N_M}{Number of haploid (male) genotype entries in the window.}
#'     \item{N_AA}{Count of homozygous reference (\code{AA}) genotypes.}
#'     \item{N_Aa}{Count of heterozygous (\code{Aa}) genotypes.}
#'     \item{N_aa}{Count of homozygous alternative (\code{aa}) genotypes.}
#'     \item{N_A}{Count of haploid reference (\code{A}) genotypes.}
#'     \item{N_a}{Count of haploid alternative (\code{a}) genotypes.}
#'     \item{Freq.Ref}{Overall reference allele frequency in the window.}
#'     \item{Freq.Alt}{Overall alternative allele frequency in the window.}
#'     \item{Obs.Hom}{Observed proportion of homozygous diploid genotypes
#'       (\code{AA + aa}) relative to total entries.}
#'     \item{Obs.Het}{Observed proportion of heterozygous diploid genotypes
#'       (\code{Aa}) relative to total entries.}
#'     \item{Obs.M.Ref}{Observed proportion of haploid reference genotypes
#'       (\code{A}) relative to total entries.}
#'     \item{Exp.Hom}{Expected frequency of homozygous diploid genotypes
#'       under the haplo-diploid equilibrium model.}
#'     \item{Exp.Het}{Expected frequency of heterozygous diploid genotypes
#'       under the haplo-diploid equilibrium model.}
#'     \item{Exp.M.Ref}{Expected frequency of haploid reference genotypes
#'       under the haplo-diploid equilibrium model.}
#'     \item{Fis}{Inbreeding coefficient for the window, or \code{NA} when
#'       expected heterozygosity is zero.}
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
#' # geno <- compute_allele.freqs_W(geno.data   = gt,
#' #                                pop.file    = pop.file,
#' #                                contigs     = contigs,
#' #                                positions   = pos,
#' #                                window.size = 10000,
#' #                                dip_freq    = 0.5)
#'
#' @seealso [summarize_geno()] for computing weighted genome-wide summary
#'   statistics from the output.
#'
#' @export
compute_allele.freqs_W <- function(geno.data, pop.file, contigs, positions, window.size, dip_freq) {
  
  all_results <- list()
  pops <- unique(pop.file$Pop)
  
  # for each population
  for (pop in pops) {
    
    cat("Processing population:", pop, "\n")
    
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
      cat("  Contig:", contig_name, "\n")
      
      contig_data <- df[contig == contig_name]
      
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
        results_window <- data.table::data.table(
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
      contig_results <- data.table::rbindlist(results_list[!sapply(results_list, is.null)])
      pop_results[[length(pop_results) + 1]] <- contig_results
      
      # Clear contig data
      rm(contig_data)
      #gc(verbose = FALSE)
    }
    all_results[[pop]] <- data.table::rbindlist(pop_results)
    # Clear population data
    rm(gt_matrix_pop, df, results_list)
    gc(verbose = FALSE)
  }
  
  final_output <- data.table::rbindlist(all_results)
  return(final_output)
}


#' Summarize per-window genotype frequencies per population
#'
#' Computes the site-count-weighted mean and standard deviation of observed
#' and expected heterozygosity, and of observed and expected haploid reference
#' allele frequency, across all windows for each population. Uses the
#' per-window table produced by [compute_allele.freqs_W()]. Weighting by
#' \code{N_sites} ensures that windows with more called genotypes contribute
#' more to each estimate.
#'
#' @param geno_table A [data.table::data.table] produced by
#'   [compute_allele.freqs_W()], containing at minimum the columns
#'   \code{Pop}, \code{N_sites}, \code{Exp.Het}, \code{Obs.Het},
#'   \code{Exp.M.Ref}, and \code{Obs.M.Ref}.
#'
#' @return A [tibble::tibble] with one row per population and the following
#'   columns:
#'   \describe{
#'     \item{Pop}{Population label.}
#'     \item{wMean.Exp.Het}{Weighted mean of expected heterozygosity across
#'       all windows.}
#'     \item{wSD.Exp.Het}{Weighted standard deviation of expected
#'       heterozygosity across all windows.}
#'     \item{wMean.Obs.Het}{Weighted mean of observed heterozygosity across
#'       all windows.}
#'     \item{wSD.Obs.Het}{Weighted standard deviation of observed
#'       heterozygosity across all windows.}
#'     \item{wMean.Exp.M.Ref}{Weighted mean of expected haploid reference
#'       allele frequency across all windows.}
#'     \item{wSD.Exp.M.Ref}{Weighted standard deviation of expected haploid
#'       reference allele frequency across all windows.}
#'     \item{wMean.Obs.M.Ref}{Weighted mean of observed haploid reference
#'       allele frequency across all windows.}
#'     \item{wSD.Obs.M.Ref}{Weighted standard deviation of observed haploid
#'       reference allele frequency across all windows.}
#'   }
#'
#' @examples
#' # Assuming compute_allele.freqs_W() has already been run:
#' # geno    <- compute_allele.freqs_W(...)
#' # summary <- summarize_geno(geno)
#'
#' @seealso [compute_allele.freqs_W()] for computing the input per-window
#'   table.
#'
#' @export
summarize_geno <- function(geno_table){
  summary_df <- geno_table |> dplyr::group_by(Pop) |> dplyr::summarise(
    # weighted mean and sd of expected heterozygosity
    wMean.Exp.Het = stats::weighted.mean(Exp.Het, N_sites, na.rm = TRUE),
    wSD.Exp.Het = matrixStats::weightedSd(Exp.Het, N_sites, na.rm = TRUE),
    # weighted mean and sd of observed heterozygosity
    wMean.Obs.Het = stats::weighted.mean(Obs.Het, N_sites, na.rm = TRUE),
    wSD.Obs.Het = matrixStats::weightedSd(Obs.Het, N_sites, na.rm = TRUE),
    # weighted mean and sd of expected male ref. allele genotype
    wMean.Exp.M.Ref = stats::weighted.mean(Exp.M.Ref, N_sites, na.rm = TRUE),
    wSD.Exp.M.Ref = matrixStats::weightedSd(Exp.M.Ref, N_sites, na.rm = TRUE),
    # weighted mean and sd of observed male ref. allele genotype
    wMean.Obs.M.Ref = stats::weighted.mean(Obs.M.Ref, N_sites, na.rm = TRUE),
    wSD.Obs.M.Ref = matrixStats::weightedSd(Obs.M.Ref, N_sites, na.rm = TRUE)
  )
  return(summary_df)
  # summary_df |> dplyr::summarise(Mean = mean(wMean.Exp.Het), SD = sd(wMean.Exp.Het))
}


