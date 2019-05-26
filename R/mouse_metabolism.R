#' Mouse metabolism microarray data
#'
#' The purpose of the multi-tissue mouse metabolism transcriptomic data is to study how the gene expression changes with respect to the energy deficiency using mouse models. 
#' Very long-chain acyl-CoA dehydrogenase (VLCAD) deficiency was found to be associated with
#' energy metabolism disorder in children.
#' Two genotypes of the mouse model - wild type (VLCAD +/+) and VLCAD-deficient (VLCAD -/-) -
#' were studied for three types of tissues (brown fat, liver, heart) with 3 to 4 mice in each genotype group.
#' The sample size information is available in the table below.
#' A total of 6,883 genes are available in this example dataset.
#'
#' @format A list of data.frame with 6,883 genes (rows) and 3 - 4 mouse samples in each genotype group (columns).
#' \describe{
#'   \item{brown}{data for the brown fat tissue}
#'   \item{heart}{data for the heart tissue}
#'   \item{liver}{data for the liver tissue}
#' }
#' @source \url{https://projecteuclid.org/download/pdfview_1/euclid.aoas/1310562214}
#' @examples
#' data(data_mouseMetabolism)
"data_mouseMetabolism"


