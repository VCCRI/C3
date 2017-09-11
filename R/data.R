#' Human ENCODE data
#'
#' Data from the human ENCODE project, contains transcriptomic profiles of 144 different cell lines and/or tissues.
#'
#' @docType data
#'
#' @usage data(human.encode.data)
#'
#' @format An object of list containing 144 data frames of the gene expression profiles.
#'
#' @keywords datasets
#'
#' @references The ENCODE Project Consortium, Nature 2012 Sep 6;489(7414):57-74
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/22955616}{PubMed})
#'
#' @examples
#' summary(human.encode.data)
"human.encode.data"


#' Human ENCODE data descriptor
#'
#' Data from the human ENCODE project, contains experiment descriptors of the data sets.
#'
#' @docType data
#'
#' @usage data(human.encode.data.descriptor)
#'
#' @format An object of vector containing 144 experiment names.
#'
#' @keywords datasets
#'
#' @references The ENCODE Project Consortium, Nature 2012 Sep 6;489(7414):57-74
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/22955616}{PubMed})
#'
#' @examples
#' head(human.encode.data.descriptor)
"human.encode.data.descriptor"


#' Hawse human lens data
#'
#' Data from the Hawse et al human lens,
#' contains transcriptomic profiles data of lens epithelial cell (LEC) and lens fiber cell (LFC).
#' Each cell has 3 biological replicates.
#'
#' @docType data
#'
#' @usage data(hawse.human.lens.data)
#'
#' @format An object of matrix containing gene expression value of LEC and LFC.
#'
#' @keywords datasets
#'
#' @references Hawse et al, Mol Vis. 2005 Apr 18; 11: 274–283
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/15851978}{PubMed})
#'
#' @examples
#' head(hawse.human.lens.data)
"hawse.human.lens.data"


#' Hoang mouse lens data
#'
#' Data from the Hoang et al mouse lens,
#' contains transcriptomic profiles data of lens epithelial cell (LEC) and lens fiber cell (LFC).
#' Each cell contains average expression value of 3 biological replicates.
#'
#' @docType data
#'
#' @usage data(hoang.mouse.lens.data)
#'
#' @format An object of matrix containing average gene expression value of LEC and LFC.
#'
#' @keywords datasets
#'
#' @references Hoang et al, Mol Vis. 2014 Nov 4;20:1491-517
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/25489224}{PubMed})
#'
#' @examples
#' head(hoang.mouse.lens.data)
"hoang.mouse.lens.data"
