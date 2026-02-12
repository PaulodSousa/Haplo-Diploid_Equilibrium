#' Import a VCF file and extract genotype and positional data
#'
#' Reads a VCF file from disk and extracts three pieces of information: the
#' genotype calls (the \code{GT} field), the contig (chromosome) name for each
#' variant site, and its physical position. These are the three inputs required
#' by the windowed summary-statistic functions in this package.
#'
#' @param path_to_vcf A length-one character string giving the path to the
#'   input VCF file.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{contig_vector}{A character vector of length \code{n_sites}
#'       containing the contig (chromosome) name for each variant.}
#'     \item{positions}{A numeric vector of length \code{n_sites} containing
#'       the physical position (bp) of each variant.}
#'     \item{gt_matrix}{A character matrix of genotype strings with dimensions
#'       \code{n_sites x n_individuals} (e.g. \code{"0/0"}, \code{"0/1"},
#'       \code{"1"}). Row names are inherited from the VCF variant records and
#'       column names correspond to sample identifiers in the VCF header.}
#'   }
#'
#' @examples
#' vcf_path <- system.file("extdata",
#'                         "example.vcf",
#'                         package = "HaploDiploidEquilibrium")
#'
#' result <- vcf2GT(vcf_path)
#'
#' head(result$contig_vector)
#' head(result$positions)
#' head(result$gt_matrix)
#'
#' @seealso [vcfR::read.vcfR()] for full control over VCF parsing options.
#'
#' @export
vcf2GT <- function(path_to_vcf) {
  if (!length(path_to_vcf) || !nzchar(path_to_vcf) || !file.exists(path_to_vcf)) {
    stop("VCF file not found", call. = FALSE)
  }
  vcf <- vcfR::read.vcfR(path_to_vcf)
  contig_vector <- vcf@fix[, "CHROM"]
  positions <- as.numeric(vcf@fix[, "POS"])
  gt_matrix <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)
  return(list(contig_vector = contig_vector, positions = positions, gt_matrix = gt_matrix))
}

