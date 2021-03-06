% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c3_functions.R
\name{c3_test}
\alias{c3_test}
\title{c3_test}
\usage{
c3_test(processed.queryData, marker.gene.compendium, min = 1, max = 5000,
  top.expressed.genes = NULL)
}
\arguments{
\item{processed.queryData}{A list containing the average gene expression values for the cell types and/or tissues with ensembl gene IDs and the species name, as returned by 'preprocess_querydata'.}

\item{marker.gene.compendium}{A list containing the data of marker gene compendium, as returned by 'make_marker_gene_compendium'.}

\item{min}{The minimum number of marker genes required in a marker gene set to be tested. Setting this higher may reduce false positives. Default is 0.}

\item{max}{The maximum number of genes allowed in a marker gene set to be tested. Default is '4000'.}

\item{top.expressed.genes}{The number of highest expressed genes from each cell/tissue that will be considered for marker gene status. Default is NULL, i.e., set the same number internally as the compendium.}
}
\value{
This function returns a list containing the XGSA p-values and the corresponding overlapping genes from the query data and the compendium. Each element of the list corresponds to one query sample vs. one compendium cell type.
}
\description{
This function performs the XGSA Fisher's exact tests between the marker genes of the query data and each cell type or tissue in the marker gene compendium.
}
\details{
This function performs the XGSA Fisher's exact tests between the marker genes of the query data and each cell type or tissue in the marker gene compendium.
}
\examples{
## Here we will use "human.encode.data" and "human.encode.data.descriptor" from C3 repository to make the compendium.
## These data sets are loaded automatically with the package.
## The gene expression data sets are in list format. Here every list element contains transcriptomic profile data of a cell or tissue type.
human.compendium<-make_gene_expression_compendium(gene.expression.data=human.encode.data, experiment.descriptor=human.encode.data.descriptor, expression.data.format="list")

## Then we will add "hawse.human.lens.data" to the compendium. This data is also loaded with the package.
## This data set has 2 different cell types - lens epithelial cell (LEC) and lens fiber cell (LFC). Each cell type has 3 biological replicates.
Hawse.data.descriptor<-c("Hawse_LEC","Hawse_LEC","Hawse_LEC","Hawse_LFC","Hawse_LFC","Hawse_LFC")
Hawse.human.compendium<-add_sample_into_compendium(compendium.data=human.compendium, sample.data=hawse.human.lens.data, species = "hsapiens", data.format = "matrix", geneID = "external_gene_name", experiment.descriptor = Hawse.data.descriptor)

## Next we will make the marker gene compendium from the Hawse.human.compendium
human.marker.gene.compendium<-make_marker_gene_compendium(compendium.data =Hawse.human.compendium, specific.cutoff = 0.05, top.expressed.genes = 500)

## Now we will pre-process the "hoang.mouse.lens.data" (query data). This data is also loaded with the package.
## This data set also contains 2 separate cell type: LEC and LFC. Each cell contains average expression value of 3 replicates.
Hoang.processed.query.data<-preprocess_querydata(cell.tissue.data = hoang.mouse.lens.data, species = "mmusculus", data.format = "matrix", geneID = "external_gene_name")

## Finally we will perform the test of the query data with the human.marker.gene.compendium
Hoang.data.test.result<-c3_test(processed.queryData = Hoang.processed.query.data, marker.gene.compendium = human.marker.gene.compendium)
head(sort(Hoang.data.test.result$Hoang_LEC$pvalue))
head(sort(Hoang.data.test.result$Hoang_LFC$pvalue))

}
