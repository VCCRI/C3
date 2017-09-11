#########################################################################################
#' @title C3: Cross-species compendium-based cell-type identification
#'
#' @description C3 is an R package for cross-species compendium-based cell-type identification. This package contains the neccessary R code to run C3 as described in "C3: An R package for compendium-based cross-species cell-type identification". C3 is implemented within a helpful framework to identify an unknown cell-type from a transcriptomic profile.
#'
#'
#' @author Md Humayun Kabir <humayun.mkabir@gmail.com>
#'
#' @docType package
#' @name c3-package
#' @aliases c3 C3
#'
#' @examples
#' ## Do a sample analysis using human compendium data and mouse lens query data.
#'
#'
NULL






##################################################################################
# Copyright Victor Chang Cardiac Research Institute 2016
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#





######################################################################################
# Users should also check the names of functions loaded by this script
#
# Sometimes you might get temporarily disconnected from the BioMart web server with an error message like:
#
# "Error in value[[3L]](cond) :
#  Request to BioMart web service failed. Verify if you are still connected to the internet.  Alternatively the BioMart web service is temporarily down."
#
# If this happens just try to re run the command that failed.
#
# C3 depends on three libraries - biomaRt, slam and xgsa
# To install these execute the following commands:
#
# source("http://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# install.packages("slam")
# install.packages('devtools')
# devtools::install_github('VCCRI/XGSA')
#
#
#
###############################
library(biomaRt)
library(slam)
library(xgsa)
#
#
#
#
################################################################
# Please note that it is recommended to pre-process both the compendium and query data in the same normalized format before utilizing these data with the C3 package,
# i.e., these data should be as consistent as possible in terms of pre-processing and log normalization.
# For micro-array data sets be aware that extra steps like background subtraction and quantile normalization are recommended.
# The C3 package does not perform any log normalizations.
#






###############################################################
# Use ENSEMBL to generate mapping  between ENSEMBL ids and gene symbols
# When Ensembl move servers (rarely) or it is down (less rare) these variable may need to be changed, for example to an archived version

#' @title c3_params
#' @description
#' When Ensembl move servers (rarely) or it is down (less rare) these variable may need to be changed, for example to an archived version
#' @rdname c3_params
#' @name c3_params
#' @export biomart.ID
#' @export host.ID
#'
biomart.ID <- "ENSEMBL_MART_ENSEMBL"
host.ID <- "www.ensembl.org"
#host.ID <- "asia.ensembl.org"
#host.ID <- "dec2015.archive.ensembl.org"






#########################################################################################
# For valid gene Ids of a species that have mapping to the ensembl_gene_id
# This function returns the valid gene ids of the species that have mapping to the ensembl_gene_id

#' @title find_valid_geneID
#'
#' @description
#' This helper function returns a data frame with two columns containing all valid gene IDs and their description that have mapping to Ensembl gene IDs for the given species. For supported species that have ensembl id mapping run 'supportedSpecies' (this command is from XGSA package), it will return all the names of the supported species. For each species there exists separate ensembl datasets.
#'
#' @rdname find_valid_geneID
#' @name find_valid_geneID
#'
#' @details This helper function returns a data frame with two columns containing all valid gene IDs and their description that have mapping to Ensembl gene IDs for the given species. For supported species that have ensembl id mapping run 'supportedSpecies' (this command is from XGSA package), it will return all the names of the supported species. For each species there exists separate ensembl dataset.
#' @return This helper function returns a data frame with two columns containing all valid gene IDs and their description that have mapping to Ensembl gene IDs for the given species
#'
#' @param species Species name in the form "hsapiens"
#'
#' @importFrom biomaRt useMart listAttributes
#'
#' @export
#'
#' @examples
#' human.valid.gene.IDs<-find_valid_geneID(species = "hsapiens")
#' head(human.valid.gene.IDs)

find_valid_geneID <- function(species = "hsapiens"){
  dataset.name <- paste(species, "_gene_ensembl", sep="")
  if(interactive()){
    ensembl = useMart("ensembl", dataset=dataset.name)
    return(listAttributes(ensembl))
  }
}






###############################################################################################
# Returns a data frame with two columns containing all Ensembl gene IDs and their corresponding gene IDs (eg. gene symbols) for the given species
# Input parameter is the species string like ('hsapiens')
# gene mapping for ensembl gene id with other gene id
# Here, geneID.forMapping is used for the existing gene ID (that has ensembl id mapping in the biomaRt database), such as "external_gene_name" for geneSymbol.
# For more details to get the valid gene ids for a species please see the find_valid_geneID() function of the C3 package.

#' @title make_ENSEMBL_symbol_map
#'
#' @description
#' This helper function returns a data frame with two columns containing all Ensembl gene IDs and their corresponding gene IDs for the given species
#'
#' @rdname make_ENSEMBL_symbol_map
#' @name make_ENSEMBL_symbol_map
#'
#' @details This helper function returns a data frame with two columns containing all Ensembl gene IDs and their corresponding gene IDs for the given species
#' @return This helper function returns a data frame with two columns containing all Ensembl gene IDs and their corresponding gene IDs for the given species
#'
#' @param species Species name in the form 'hsapiens'
#' @param ensemblId.toMapping The valid gene ID names to which the user cell data will be converted, default is "ensembl_gene_id"
#' @param geneID.forMapping The valid gene ID names of the user cell data, default is "external_gene_name"
#'
#' @importFrom biomaRt useMart getBM
#'
#' @export
#'
#' @examples
#' human.ensembl.symbol.map<-make_ENSEMBL_symbol_map('hsapiens', ensemblId.toMapping="ensembl_gene_id", geneID.forMapping="external_gene_name")
#' head(human.ensembl.symbol.map)
#'

make_ENSEMBL_symbol_map <- function(species, ensemblId.toMapping, geneID.forMapping){
  dataset.name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomart = biomart.ID, dataset = dataset.name, host = host.ID)
  ENSEMBL.symbol.map <- getBM(attributes=c(ensemblId.toMapping, geneID.forMapping), mart = ensembl)
  return(ENSEMBL.symbol.map)
}






############################################################################################################
# This function converts gene IDs to ensembl IDs
# 1st section is for marker genes - 2 options one for matrix and other for vector format
# for matrix- first create a list and then converting each column genes is stored in a sublist of the list
# 2nd section for cells/tissues gene expression profiles
# for example, used external_gene_name for geneID.forMapping and ensembl_gene_id for ensemblId.toMapping

#' @title convert_geneID_to_ensemblID
#'
#' @description
#' This function converts either a vector of gene IDs, a matrix of gene or a named cell / tissue expression profile or matrix to ensembl IDs
#'
#' @rdname convert_geneID_to_ensemblID
#' @name convert_geneID_to_ensemblID
#'
#' @details
#' TThis function converts a vector of gene IDs or a named cell/tissue expression profile to ensembl IDs
#'
#' @return This function returns either ensembl gene IDs for the corresponding gene IDs of marker genes or cell/tissue profiles with converted ensembl gene IDs
#'
#' @param user.cell.data Either valid gene IDs for marker genes or cell/tissue profile with valid gene IDs
#' @param species The species of the given gene IDs or cell/tissue profile, default is "hsapiens"
#' @param ensemblId.toMapping The valid gene ID names (that has ensembl id mapping in the biomaRt database) to which the user.cell.data will be converted, default is "ensembl_gene_id". For more details to get the valid gene ids for a species please see the find_valid_geneID() function of the C3 package.
#' @param geneID.forMapping The valid gene ID names (that has ensembl id mapping in the biomaRt database) of the user.cell.data, default is "external_gene_name". For more details to get the valid gene ids for a species please see the find_valid_geneID() function of the C3 package.
#' @param collapse.method Used when one ensembl_gene_id has more then one probes, usually two options are used, maximum probe value (i.e., "max") or average of probe values (i.e., "mean"), default is "max". Not used for marker gene option.
#' @param forMarkerGenes Boolean with either T for marker genes conversion or F for cell/tissue profile symbol conversion, default is F
#' @param markerGeneFormat Data format of marker genes either 'list' or 'vector', used only during marker genes conversion, default is "vector"
#'
#' @export
#'
#' @examples
#' marker.genes<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' convert_geneID_to_ensemblID(user.cell.data=marker.genes, species = "hsapiens", forMarkerGenes = T, markerGeneFormat = "vector")
#'

convert_geneID_to_ensemblID <- function(user.cell.data, species="hsapiens", ensemblId.toMapping="ensembl_gene_id", geneID.forMapping = "external_gene_name", collapse.method = "max", forMarkerGenes=F, markerGeneFormat="vector"){
  ensembl.symbol.map<-make_ENSEMBL_symbol_map(species, ensemblId.toMapping, geneID.forMapping)

  #for marker genes
  if(forMarkerGenes==T){
    if(markerGeneFormat=="list"){
      user.cell.data.ensemblID<-lapply(user.cell.data, function(x){
        ensembl.genes.exist.in.user.cell.data <- ensembl.symbol.map[ensembl.symbol.map[[2]] %in% x,]
        duplicate.ensembl.gene.list <- duplicated(ensembl.genes.exist.in.user.cell.data[[2]])
        ensembl.genes.exist.in.user.cell.data.without.duplication <- ensembl.genes.exist.in.user.cell.data[!duplicate.ensembl.gene.list,]
        return(ensembl.genes.exist.in.user.cell.data.without.duplication[[1]])
      })
      return(user.cell.data.ensemblID)
    }
    else if(markerGeneFormat=="vector"){
      ensembl.genes.exist.in.user.cell.data <- ensembl.symbol.map[ensembl.symbol.map[[2]] %in% user.cell.data,]
      duplicate.ensembl.gene.list <- duplicated(ensembl.genes.exist.in.user.cell.data[[2]])
      ensembl.genes.exist.in.user.cell.data.without.duplication <- ensembl.genes.exist.in.user.cell.data[!duplicate.ensembl.gene.list,]
      user.cell.data.ensemblID<- ensembl.genes.exist.in.user.cell.data.without.duplication[[1]]
      return(user.cell.data.ensemblID)
    }
    else{
      print("ERROR: Do not support other data format at this moment!!")
      return(NULL)
    }
  }
  #for cell/tissue gene expression data sets
  else{
    cell.tissue.names<-colnames(user.cell.data)
    list.gene.names<-rownames(user.cell.data)
    ensembl.genes.exist.in.user.cell.data <- ensembl.symbol.map[ensembl.symbol.map[[2]] %in% list.gene.names,]

    #duplicate probes / symbol IDS - multiple ENS genes match to the same probe or symbol - probably mostly NA or ""
    duplicate.ensembl.gene.list <- duplicated(ensembl.genes.exist.in.user.cell.data[[2]])
    ensembl.genes.exist.in.user.cell.data.without.duplication <- ensembl.genes.exist.in.user.cell.data[!duplicate.ensembl.gene.list,]

    idlist <- split(ensembl.genes.exist.in.user.cell.data.without.duplication[[2]], ensembl.genes.exist.in.user.cell.data.without.duplication[[1]])

    if(collapse.method == "mean"){
      #for mean
      collapsed <- lapply(idlist,function(X){
        return(colMeans(user.cell.data[unique(X),,drop=F]))
      })
    }
    else if(collapse.method == "max"){
      #for max
      collapsed <- lapply(idlist,function(X){
        return(user.cell.data[unique(X),,drop=F][which.max(rowMeans(user.cell.data[unique(X),,drop=F])),])
      })
    }
    else {
      print("ERROR: Other collapse methods not supported!!")
      return(NULL)
    }

    user.cell.data.clean <- do.call(rbind, collapsed)
    colnames(user.cell.data.clean)<-cell.tissue.names

    return(user.cell.data.clean)
  }
}






#################################################################################################
# calculate mean of the cells based on replicates

#' @title compute_mean
#'
#' @description
#' This function calculates average expression values of multiple replicate samples
#'
#' @rdname compute_mean
#' @name compute_mean
#'
#' @details
#' This function calculates average expression values of multiple replicate samples using rowMeans. If the data has no replicates (a single vector or a matrix with only one column) the original values are returned.
#'
#' @return This function returns a vector with the mean expression values.
#'
#' @param data.matrix A matrix with gene expression values, rows denote genes and columns denote replicate samples. Alternatively a vector of gene expression values.
#'
#' @export
#'
#' @examples
#' dd<-matrix(sample(1:10, 30, replace=TRUE), 10, 3)
#' rownames(dd)<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' compute_mean(dd)
#'

compute_mean <- function(data.matrix){
  if((is.vector(data.matrix)==TRUE) || (ncol(data.matrix)<2))
    data.matrix<-cbind(val1=data.matrix,vall2=data.matrix)
  return(rowMeans(data.matrix))
}






#################################################################################################
# generate matrix from a list (list can have more than one replicate for each cell/tissue types)

#' @title format_list_data
#'
#' @description
#' This function transforms a list of gene expression data (matrices or vectors) into a single matrix by averaging the replicates for each cell/tissue in the list.
#'
#' @rdname format_list_data
#' @name format_list_data
#'
#' @details
#' This function transforms a list of gene expression vectors or matrices with replicates into a single matrix by averaging the replicate values for each cell/tissue. The list elements should all have consistently named vectors, or rownames for matrices, representing the genes.
#'
#' @return This function returns a matrix with average gene expression values of the cells/tissues.
#'
#' @param list.data A list of cell type specific gene expression data. Each element of the list contains a named gene expression vector or matrix with replicates. The list names (denoting cell types) should be unique, otherwise the experiment.descriptor parameter value should be provided with unique names for each cell or tissue.
#' @param experiment.descriptor A vector correspoding to the elements of list.data that denotes the cell and/or tissue names of the data set. The names should be unique, defaults to the names of list.data.
#'
#' @export
#'
#' @examples
#' ll <- list()
#' ll[["cell1"]]<-matrix(sample(1:10, 30, replace=TRUE), 10, 3)
#' rownames(ll$cell1)<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' ll[["cell2"]]<-matrix(sample(1:10, 30, replace=TRUE), 10, 3)
#' rownames(ll$cell2)<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' ll[["cell3"]]<-matrix(sample(1:10, 30, replace=TRUE), 10, 3)
#' rownames(ll$cell1)<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' format_list_data(ll)
#'

format_list_data<-function(list.data, experiment.descriptor = names(list.data)){
  names(list.data)<-experiment.descriptor
  mat.data<-do.call(cbind,lapply(list.data,compute_mean))
  return(mat.data)
}






####################################################################################################
# format matrix data of gene expression profile
# generate a matrix from user input (make averages of the replicates)

#' @title format_matrix_data
#'
#' @description
#' This function re-formats gene expression matrix data by averaging replicate samples. The cell type identifiers should be the column names of the expression matrix, otherwise the experiment.descriptor parameter is provided to specify cell-types of the samples (i.e. from a separate metadata file).
#'
#' @rdname format_matrix_data
#' @name format_matrix_data
#'
#' @details
#' This function re-formats gene expression matrix data by averaging replicate samples. The cell type identifiers should be the column names of the expression matrix, otherwise the experiment.descriptor parameter is provided to specify cell-types of the samples (i.e. from a separate metadata file).
#'
#' @return This function returns a matrix with mean value of the replicates for each cell type.
#'
#' @param matrix.data A matrix with gene expression profiles where rows denote the genes and the columns denote the cells and/or tissues. Duplicate column names are expected in this case denoting replicate samples. All the replicate samples for a specific cell or tissue should have identical column names, otherwise the experiment.descriptor parameter should be used to identify replicate samples of a specific cell type or tissue.
#' @param experiment.descriptor A vector corresponding to the columns of matrix.data, containing the cell type or tissues of each sample. The names should be identical for a specific cell or tissue.  Defaults to the column names of the matrix.data.
#'
#' @export
#'
#' @examples
#' mm<-matrix(sample(1:10, 100, replace=TRUE),10,10)
#' rownames(mm)<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' colnames(mm)<-c("cell2", "cell3", "cell1", "cell2", "cell1", "cell3", "cell3", "cell2", "cell3", "cell1")
#' format_matrix_data(mm)
#'

format_matrix_data<-function(matrix.data, experiment.descriptor = colnames(matrix.data)){
  colnames(matrix.data)<-experiment.descriptor

  #sorting the matrix data sets according to the column names
  ordered.column.names<-colnames(matrix.data)[order(colnames(matrix.data))]
  matrix.data.ordered<-matrix.data[,order(colnames(matrix.data))]
  colnames(matrix.data.ordered)<-ordered.column.names

  repCellTissue=table(colnames(matrix.data.ordered))
  nCellTissue<-length(repCellTissue)
  nameCellTissue<-names(repCellTissue)


  #processing the matrix data to make a new matrix with number of rows same as before, but
  #number of columns equal to the number of cells/tissues by calculating the average value of the replicates
  if(ncol(matrix.data.ordered)==nCellTissue){
    f.mat<-matrix.data.ordered  #here each cell/tissue has only 1 replicate, with existing cell names
  }
  else{                         #here each cell/tissue has one or more replicates
    f.mat<-matrix(0, nrow=nrow(matrix.data.ordered), ncol=nCellTissue)
    r<-repCellTissue
    s<-1
    e<-r[1]
    for(i in 1:nCellTissue){
      mean.data<-compute_mean(matrix.data.ordered[,s:e])
      f.mat[,i]<-mean.data
      if(i != nCellTissue){
        s<-s+r[i]
        e<-e+r[i+1]
      }
    }
    rownames(f.mat)<-rownames(matrix.data.ordered)
    colnames(f.mat)<-nameCellTissue
  }

  return(f.mat)
}






#####################################################################################################
# For making compendium
# Make a compendium for gene expression data by averaging the replicate values, can be either in list or matrix format
# Compendium can be generated with any type of valid gene IDs. For more details to get the valid gene ids for a species please see the find_valid_geneID() function of the C3 package.
# Then can convert the gene IDs to ensembl_gene_id using function convert_geneID_to_ensemblID by setting appropriate parameters

#' @title make_gene_expression_compendium
#'
#' @description
#' This function makes a compendium (matrix) of gene expression data by transforming a list or matrix of gene expression profiles (with replicates) representing different cell types or tissues.
#'
#' @rdname make_gene_expression_compendium
#' @name make_gene_expression_compendium
#'
#' @details
#' This function makes a compendium (matrix) of gene expression data by transforming a list or matrix of gene expression profiles (with replicates) representing different cell types or tissues.
#' For both the input and output matrices, rows denote genes and columns denote cell/tissue types. The compendium can be generated with any type of valid gene IDs, and can then be converted to ensembl_gene_id using function convert_geneID_to_ensemblID with appropriate parameters.
#'
#' @return This function returns a list with two elements. The element 'counts' is the compendium matrix where rows denote genes and columns denote individual cell types or tissues. The element 'species' represents the compendium data species.
#'
#' @param gene.expression.data A list or matrix containing gene expression profiles. For list each element will be treated as an individual cell/tissue, therefore it should have a unique name and contain a vector of gene expression values, or a matrix with replicates. For a matrix input, the row names denote the genes and column names denote the cell type or tissue. Duplicate column names are expected in this case denoting replicate samples. All the replicate samples for a specific cell or tissue should have identical column names, otherwise the experiment.descriptor parameter should be used to identify replicate samples of a specific cell type or tissue.
#' @param species The species abbreviation of the compendium data (gene.expression.data). Default is "hsapiens".
#' @param experiment.descriptor A vector corresponding to the column names (matrix) or elements (list) of gene.expression.data, containing the cell type or tissues of each sample. The names should be identical for a specific cell or tissue.  Defaults to NULL.
#' @param expression.data.format The format of gene.expression.data, either 'list' or 'matrix'. Defaults to "list".
#'
#' @export
#'
#' @examples
#' exp.data<-matrix(sample(1:10, 100, replace=TRUE),10,10)
#' rownames(exp.data)<-c("ENSG00000100053", "ENSG00000109846", "ENSG00000244752", "ENSG00000138083", "ENSG00000181449", "ENSG00000160202", "ENSG00000007372", "ENSG00000117707", "ENSG00000101144", "ENSG00000277354")
#' colnames(exp.data)<-c("cell1", "cell1", "cell1", "cell2", "cell2", "cell2", "cell3", "cell3", "cell3", "cell3")
#' make_gene_expression_compendium(gene.expression.data=exp.data, expression.data.format="matrix")
#'

make_gene_expression_compendium<-function(gene.expression.data, species="hsapiens", experiment.descriptor=NULL, expression.data.format="list"){
  #make the compendium matrix according to the the data format
  if(expression.data.format=="list"){
    if(is.null(experiment.descriptor)) experiment.descriptor<-names(gene.expression.data)
    user.compendium.data<-format_list_data(gene.expression.data, experiment.descriptor)
  }
  else if(expression.data.format=="matrix"){
    if(is.null(experiment.descriptor)) experiment.descriptor<-colnames(gene.expression.data)
    user.compendium.data<-format_matrix_data(gene.expression.data, experiment.descriptor)
  }
  else{
    print("ERROR: could not support other data format at this moment!!")
    return(NULL)
  }

  #make a list 'compendium' with two elements - 'counts' as compendium matrix and 'species' as compendium data species
  compendium<-list()
  compendium[["avg.exp.matrix"]]<-user.compendium.data
  compendium[["species"]]<-species
  return(compendium)
}






#####################################################################################################################
# Add gene expression cell data to compendium
# 3 options for 3 formatted data - one for list data, next one  for matrix data, and final one for vector data
# for list data, assume that each sublist will contain each cell/tissue/development_stage gene expression data
# for matrix data, rows denotes the gene IDs and columns denote replicates of cells/tissues
# for vector data, elements will be expression values and the names of elements will be gene IDs

#' @title add_sample_into_compendium
#'
#' @description
#' This function adds cell/tissue gene expression data into a previously made compendium by making common genes for all cells/tissues.
#'
#' @rdname add_sample_into_compendium
#' @name add_sample_into_compendium
#'
#' @details
#' This function adds cell type or tissue specific gene expression data into an existing gene expression compendium and updates the compendium with a common gene universe.
#'
#' @return This function returns a  list by adding the new sample data into the gene expression compendium (i.e., into 'avg.exp.matrix').
#'
#' @param compendium.data A gene expression compendium (list), as returned by make_gene_expression_compendium.
#' @param sample.data A list, matrix or vector containing gene expression profiles. For a list, Each element of the list contains a named gene expression vector or matrix with replicates. The list names (denoting cell types) should be unique, otherwise the experiment.descriptor parameter value should be provided with unique names for each cell or tissue. For a matrix, rows denote the genes and the columns denote the cell types or tissues. Duplicate column names are expected in this case denoting replicate samples. All the replicate samples for a specific cell or tissue should have identical column names, otherwise the experiment.descriptor parameter should be used to identify replicate samples of a specific cell type or tissue. For a vector, it should be a contain expression values from a single cell type or tissue and the names of elements should be gene IDs.
#' @param species The species of sample.data, used to identify valid gene names. The species should be the same for the compendium. Default is "hsapiens"
#' @param data.format The format of sample.data, either 'list', 'matrix' or 'vector'. Default is "matrix".
#' @param geneID The code for the type of gene IDs used by sample.data, as used by the biomaRt database. To find the valid codes for gene IDs for a species, please see the find_valid_geneID() function of the C3 package. Default is "ensembl_gene_id".
#' @param experiment.descriptor A vector corresponding to the column names (matrix) or elements (list) of compendium.data, containing the cell type or tissues of each sample. The names should be identical for a specific cell or tissue.  Defaults to NULL.
#' @param collapse.method How to summarise values when one ensembl_gene_id has more then one value (multiple microarray probes or transcripts to one gene for example). Currently two options are implemented, 'max' or 'mean'. Default is "max".
#'
#' @export
#'
#' @examples
#' exp.data<-matrix(runif(100, min=0, max=10), 10, 10)
#' rownames(exp.data)<-c("ENSG00000100053", "ENSG00000109846", "ENSG00000244752", "ENSG00000138083", "ENSG00000181449", "ENSG00000160202", "ENSG00000007372", "ENSG00000117707", "ENSG00000101144", "ENSG00000277354")
#' colnames(exp.data)<-c("cell1", "cell2", "cell3", "cell4", "cell5", "cell6", "cell7", "cell8", "cell9", "cell10")
#' compendium.data<-make_gene_expression_compendium(gene.expression.data=exp.data, expression.data.format="matrix")
#' new.data<-matrix(sample(1:10, 100, replace=TRUE),10,10)
#' rownames(new.data)<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' colnames(new.data)<-c("cell11", "cell11", "cell11", "cell12", "cell12", "cell12", "cell13", "cell13", "cell13", "cell13")
#' add_sample_into_compendium(compendium.data = compendium.data, sample.data = new.data, geneID = "external_gene_name")
#'

add_sample_into_compendium<-function(compendium.data, sample.data, species = "hsapiens", data.format = "matrix", geneID = "ensembl_gene_id", experiment.descriptor = NULL, collapse.method = "max"){
  #checking for same species for compendium and sample data
  if(compendium.data$species!=species){
    print("ERROR: species must be same for the compendium and sample data!!")
    return(NULL)
  }

  #format the sample.data according to the type of format
  if(data.format=="matrix"){
    if(is.null(experiment.descriptor)) experiment.descriptor<-colnames(sample.data)
    cell.data.formatted<-format_matrix_data(sample.data, experiment.descriptor)
  }
  else if(data.format=="list"){
    if(is.null(experiment.descriptor)) experiment.descriptor<-names(sample.data)
    cell.data.formatted<-format_list_data(sample.data, experiment.descriptor)
  }
  else if(data.format=="vector"){
    cell.data.formatted<-data.frame("newly_added"=sample.data)
    if(!(is.null(experiment.descriptor))) colnames(cell.data.formatted)<-experiment.descriptor
  }
  else{
    print("ERROR: could not support other data format at this moment!!")
    return(NULL)
  }

  #conversion of gene IDs
  if(geneID=="ensembl_gene_id"){
    cell.data.ensemblID<-cell.data.formatted
  }
  else{
    cell.data.ensemblID<-convert_geneID_to_ensemblID(cell.data.formatted, species, geneID.forMapping = geneID, collapse.method = collapse.method)
  }

  #finding common genes and then taking only common genes data sets
  commonGenes<-intersect(rownames(cell.data.ensemblID),rownames(compendium.data$avg.exp.matrix))
  compendium.data$avg.exp.matrix<-compendium.data$avg.exp.matrix[commonGenes,]
  cell.data.merged <- cell.data.ensemblID[commonGenes,]

  #merging new data sets to the compendium
  compendium.data$avg.exp.matrix<-cbind(cell.data.merged, compendium.data$avg.exp.matrix)

  #checking for only one cell/tissue, i.e., a vector is generated after "cell.data.merged <- cell.data.ensemblID[commonGenes,]" for a single cell/tissue
  #then setting the names of the cell/tissue after combining the data sets
  if(is.vector(cell.data.merged))
    colnames(compendium.data$avg.exp.matrix)[1]<-colnames(cell.data.ensemblID)

  return(compendium.data)
}






####################################################################################################################
# generate marker gene compendium from the compendium data based on specific.cutoff value

#' @title make_marker_gene_compendium
#'
#' @description
#' This function transforms a gene expression compendium into a marker gene compendium, using certain threshold parameters.
#'
#' @rdname make_marker_gene_compendium
#' @name make_marker_gene_compendium
#'
#' @details
#' This function transforms a gene expression compendium into a marker gene compendium, using certain threshold parameters. The two parameters decide what is considered 'highly expressed' and how few cell types is considered 'unique' enough to be a marker gene.
#'
#' @return This function returns a  list with  cell/tissue marker gene matrix, compendium common genes, species and top number of expressed genes. In the marker gene compendium (marker.genes.matrix), rows denote genes and columns denote cell types / tissues, containing presence (1) or absence (0) of marker genes for each cell type / tissue.
#'
#' @param compendium.data A list consisiting of gene expression compendium (matrix) and compendium data species name, as returned by make_gene_expression_compendium.
#' @param specific.cutoff A percentage threshold representing the maximum number of cell types in which a gene can be highly expressed in order for it to be considered a marker gene. Default is 0.05.
#' @param top.expressed.genes The number of highest expressed genes from each cell/tissue that will be considered for marker gene status. Default is 500.
#'
#' @export
#'
#' @examples
#' ## Here we will use "human.encode.data" and "human.encode.data.descriptor" from C3 repository to make the compendium.
#' ## These data sets are loaded automatically with the package.
#' ## The gene expression data sets are in list format. Here every list element contains transcriptomic profile data of a cell or tissue type.
#' human.compendium<-make_gene_expression_compendium(gene.expression.data=human.encode.data, experiment.descriptor=human.encode.data.descriptor, expression.data.format="list")
#'
#' ## Finally we will make the marker gene compendium from the human.compendium
#' human.marker.gene.compendium<-make_marker_gene_compendium(compendium.data=human.compendium, specific.cutoff = 0.05, top.expressed.genes = 500)
#' summary(human.marker.gene.compendium)
#'

make_marker_gene_compendium<-function(compendium.data, specific.cutoff=0.05, top.expressed.genes=500){
  #compendium.data.log<-log2(compendium.data+1)

  #taking top number of expressed genes from each cell/tissue
  #top.expressed.genes.for.all.cells<-apply(compendium.data.log,2,function(x) order(x,decreasing=T)[1:top.expressed.genes])
  top.expressed.genes.for.all.cells<-apply(compendium.data$avg.exp.matrix,2,function(x) order(x,decreasing=T)[1:top.expressed.genes])


  #setting the value of top expressed genes in each cell/tissue as 1 and for other genes 0
  ExpMat<-matrix(0,nrow=dim(compendium.data$avg.exp.matrix)[1],ncol=dim(compendium.data$avg.exp.matrix)[2]);
  dimnames(ExpMat)<-dimnames(compendium.data$avg.exp.matrix);
  for(i in 1:length(top.expressed.genes.for.all.cells[1,])){
    ExpMat[top.expressed.genes.for.all.cells[,i],i]<-1;
  }

  #show the genes specificity among cells/tissues
  genes.specificity <-apply(ExpMat,1,sum);
  barplot(table(genes.specificity), main = "Genes specificity among cells", xlab = "Number of cells", ylab = "Number of expressed genes")

  #taking only cell/tissue specific genes. value 1 denotes the gene is specific in the respective cell/tissue and 0 not specific.
  indices.of.specific.genes <- which((genes.specificity > 0) & (genes.specificity < length(top.expressed.genes.for.all.cells[1,])*specific.cutoff))
  ExpMat.specific<-ExpMat[indices.of.specific.genes,]   ##This contains only cell/tissue specific genes

  #getting removed common genes
  rmGenes <- which((genes.specificity == 0) | (genes.specificity >= length(top.expressed.genes.for.all.cells[1,])*specific.cutoff))
  names.RM.Genes<-names(rmGenes)

  #making a list with specific matrix and common removal genes
  marker.gene.compendium<-list()
  marker.gene.compendium[["marker.genes.matrix"]]<-ExpMat.specific
  marker.gene.compendium[["compendium.common.genes"]]<-names.RM.Genes
  marker.gene.compendium[["species"]]<-compendium.data$species
  marker.gene.compendium[["top.expressed.genes"]]<-top.expressed.genes
  return(marker.gene.compendium)
}






#######################################################################################################################
# for adding marker genes to the marker gene compendium
# process the genes based on list or vector format
# for list - Id conversion is done at once, return a list where each element is a sublist that contains genes of each column of the matrix
# the loop runs number of matrix columns, each sublist is assigned in vector y, then find out unique existing genes and common genes in compendium
# then added it to the compendium using an array with value 1 that is previously created with default value 0
# for marker genes as ensembl id, no conversion needed and treated the original matrix for addition
# for vector - convert Id and then simply added to the compendium with data.format="vector" value

#' @title add_marker_genes_into_compendium
#'
#' @description
#' This function adds a set of marker genes to a marker gene compendium. The marker gene data can be a matrix with more than one cell type / tissue or a vector for a single cell type / tissue.
#'
#' @rdname add_marker_genes_into_compendium
#' @name add_marker_genes_into_compendium
#'
#' @details
#' This function adds a set of marker genes to a marker gene compendium. The marker gene data can be a list with more than one cell type / tissue or a vector for a single cell type / tissue. This function adds only those marker genes that are exist in the marker gene compendium.
#'
#' @return This function returns a list by adding the new sample data into the marker gene compendium ('marker.genes.matrix').
#'
#' @param marker.gene.compendium A list consisting of marker gene compendium, species, compendium common genes, as returned by make_marker_gene_compendium.
#' @param marker.genes A list or vector containing marker genes. For a list, each element denotes cell type or tissue and contains a vector of marker gene names. For a vector, it is the names of the marker genes of a single cell or tissue. The name of this cell type or tissue will need to be manually updated.
#' @param species The species abbreviation name of the marker genes.
#' @param data.format The data format of the marker.genes, either 'vector' or 'list'. Default is "vector"
#' @param geneID The code for the type of gene IDs used by sample.data, as used by the biomaRt database. To find the valid codes for gene IDs for a species, please see the find_valid_geneID() function of the C3 package. Default is "ensembl_gene_id".
#'
#' @export
#'
#' @examples
#' ## Here we will use "human.encode.data" and "human.encode.data.descriptor" from C3 repository to make the compendium.
#' ## These data sets are loaded automatically with the package.
#' ## The gene expression data sets are in list format. Here every list element contains transcriptomic profile data of a cell or tissue type.
#' human.compendium<-make_gene_expression_compendium(gene.expression.data=human.encode.data, experiment.descriptor=human.encode.data.descriptor, expression.data.format="list")
#'
#' ## Next we will make the marker gene compendium from the human.compendium
#' human.marker.gene.compendium<-make_marker_gene_compendium(compendium.data=human.compendium, specific.cutoff = 0.05, top.expressed.genes = 500)
#'
#' ## Finally we will add the marker genes of a cell type into the compendium
#' marker.genes<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "CRYBA1", "BFSP2", "CRYBB2", "BIRC6")
#' human.marker.gene.compendium.update<-add_marker_genes_into_compendium(marker.gene.compendium = human.marker.gene.compendium, marker.genes = marker.genes, data.format = "vector", geneID = "external_gene_name")
#'
#' ## Here, we will see that all marker genes are added as a 'newly added' cell type in the marker gene compendium except 'PROX1' as it does not exist in the universe of the marker gene compendium
#' names(which(human.marker.gene.compendium.update$marker.genes.matrix[,1]==1))
#'

add_marker_genes_into_compendium<-function(marker.gene.compendium, marker.genes, species = "hsapiens", data.format="vector", geneID="ensembl_gene_id"){
  #checking for same species for compendium and sample data
  if(marker.gene.compendium$species!=species){
    print("ERROR: species must be same for the compendium and marker genes data!!")
    return(NULL)
  }

  #add the marker.genes according to the type of format
  if(data.format=="list"){
    cNames<-names(marker.genes)

    #for gene ID conversion to ensembl ID
    if(geneID=="ensembl_gene_id")
      x<-marker.genes                 # here x is a list, same as marker.genes
    else
      x<-convert_geneID_to_ensemblID(marker.genes, species, geneID.forMapping = geneID, forMarkerGenes=T, markerGeneFormat=data.format)

    #for adding the marker genes list into the compendium
    #the loop runs for the number of cells/tissues in the list, each element (i.e., vector) of the list is assigned in vector y
    #then find out unique existing genes and common genes exists in the compendium
    #finally, added it to the compendium using an array with value 1 that is previously created with default value 0
    for(i in 1:length(marker.genes)){
      y<-x[[i]]                      # x is a list, each element is a vector consisting of the marker genes for a cell type
      ux<-unique(y[y != ""])         # find out unique existing genes
      ux<-intersect(ux, rownames(marker.gene.compendium$marker.genes.matrix))  #taking only existing ones in compendium
      a.marker<-array(0, nrow(marker.gene.compendium$marker.genes.matrix))
      names(a.marker)<-rownames(marker.gene.compendium$marker.genes.matrix)
      marker.gene.compendium$marker.genes.matrix<-cbind(a.marker,marker.gene.compendium$marker.genes.matrix)
      marker.gene.compendium$marker.genes.matrix[ux,1]<-1
      colnames(marker.gene.compendium$marker.genes.matrix)[1]<-cNames[i]
    }
    return(marker.gene.compendium)
  }
  else if(data.format=="vector"){
    #for gene ID conversion to ensembl ID
    if(geneID=="ensembl_gene_id")
      x<-marker.genes
    else
      x<-convert_geneID_to_ensemblID(marker.genes, species, geneID.forMapping = geneID, forMarkerGenes=T, markerGeneFormat=data.format)

    #adding the marker genes vector to the compendium
    #find out unique existing genes and common genes in compendium
    #then added it to the compendium using an array with value 1 that is previously created with default value 0
    ux<-unique(x[x != ""])
    ux<-intersect(ux, rownames(marker.gene.compendium$marker.genes.matrix))
    a.marker<-array(0, nrow(marker.gene.compendium$marker.genes.matrix))
    names(a.marker)<-rownames(marker.gene.compendium$marker.genes.matrix)
    marker.gene.compendium$marker.genes.matrix<-cbind(a.marker,marker.gene.compendium$marker.genes.matrix)
    marker.gene.compendium$marker.genes.matrix[ux,1]<-1
    colnames(marker.gene.compendium$marker.genes.matrix)[1]<-"newly_added"
    return(marker.gene.compendium)
  }
  else{
    print("ERROR: do not accept other format at this moment!!")
    return(NULL)
  }
}






###########################################################################################################
# For ID conversion and compute mean; return a list with elements - average expression value with ensembl id and species
# some ensembl associated gene ids -
#     external_gene_name -	Associated Gene Name
#     ensembl_gene_id -	Gene ID
#     ensembl_transcript_id -	Transcript ID
#     entrezgene -	EntrezGene ID

#' @title preprocess_querydata
#'
#' @description
#' This function preprocesses the query data for ID conversion and calculates average value of the replicates
#'
#' @rdname preprocess_querydata
#' @name preprocess_querydata
#'
#' @details
#' This function preprocesses the query data for ID conversion and calculates average value of the replicates
#'
#' @return This function returns a list with average value of the replicates with ensembl gene ID and species name
#'
#' @param cell.tissue.data A matrix, list or vector containing cell type / tissue specific gene expression data. For a list, Each element of the list contains a named gene expression vector or matrix with replicates. The list names (denoting cell types) should be unique, otherwise the experiment.descriptor parameter value should be provided with unique names for each cell or tissue. For a matrix, rows denote the genes and the columns denote the cell types or tissues. Duplicate column names are expected in this case denoting replicate samples. All the replicate samples for a specific cell or tissue should have identical column names, otherwise the experiment.descriptor parameter should be used to identify replicate samples of a specific cell type or tissue. For a vector, it should be a contain expression values from a single cell type or tissue and the names of elements should be gene IDs.
#' @param species The species abbreviation of the query data (cell.tissue.data). Default is "hsapiens".
#' @param data.format Format of cell.tissue.data, either 'list', 'matrix' or 'vector', Default is "matrix".
#' @param geneID The code for the type of gene IDs used by cell.tissue.data, as used by the biomaRt database. To find the valid codes for gene IDs for a species, please see the find_valid_geneID() function of the C3 package. Default is "ensembl_gene_id".
#' @param experiment.descriptor A vector corresponding to the column names (matrix) or elements (list) of cell.tissue.data, containing the cell type or tissues of each sample. The names should be identical for a specific cell or tissue.  Defaults to NULL.
#' @param collapse.method How to summarise values when one ensembl_gene_id has more then one value (multiple microarray probes or transcripts to one gene for example). Currently two options are implemented, 'max' or 'mean'. Default is "max".
#'
#' @export
#'
#' @examples
#' query.data<-matrix(sample(1:10, 100, replace=TRUE),10,10)
#' rownames(query.data)<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' colnames(query.data)<-c("cell1", "cell1", "cell1", "cell2", "cell2", "cell2", "cell3", "cell3", "cell3", "cell3")
#' preprocess_querydata(cell.tissue.data=query.data, geneID="external_gene_name")
#'

preprocess_querydata<-function(cell.tissue.data, species = "hsapiens", data.format = "matrix", geneID = "ensembl_gene_id", experiment.descriptor = NULL, collapse.method = "max"){
  if(data.format=="matrix"){
    if(is.null(experiment.descriptor)) experiment.descriptor<-colnames(cell.tissue.data)
    cell.data.formatted<-format_matrix_data(cell.tissue.data, experiment.descriptor)
  }
  else if(data.format=="list"){
    if(is.null(experiment.descriptor)) experiment.descriptor<-names(cell.tissue.data)
    cell.data.formatted<-format_list_data(cell.tissue.data, experiment.descriptor)
  }
  else if(data.format=="vector"){
    cell.data.formatted<-data.frame("newly_added"=cell.tissue.data)
    if(!(is.null(experiment.descriptor))) colnames(cell.data.formatted)<-experiment.descriptor
  }
  else{
    print("ERROR: could not support other data format at this moment!!")
    return(NULL)
  }

  #for gene id conversion
  if(geneID=="ensembl_gene_id"){
    cell.data.ensemblID<-cell.data.formatted
  }
  else{
    cell.data.ensemblID<-convert_geneID_to_ensemblID(user.cell.data = cell.data.formatted, species = species, geneID.forMapping = geneID, collapse.method = collapse.method)
  }

  #making a list with processed query data and the species name
  processed.queryData<-list()
  processed.queryData[["expressionData"]]<-cell.data.ensemblID
  processed.queryData[["species"]]<-species
  return(processed.queryData)
}






##################################################################################################################################
### Make user's query data ready for XGSA test
## takes input user's query data and species name

#' @title select_query_sample_specific_genes
#'
#' @description
#' This function converts the query data into a list of XGSA data sets containing the marker genes for each sample.
#'
#' @rdname select_query_sample_specific_genes
#' @name select_query_sample_specific_genes
#'
#' @details
#' This function converts the query data into a list of XGSA data sets containing the marker genes for each sample.
#'
#' @return This function returns a list containing XGSA data sets ready to test with the compendium data.
#'
#' @param processed.queryData A list containing the average gene expression values for the cell types and/or tissues with ensembl gene IDs and the species name, as returned by 'preprocess_querydata'.
#' @param comSpecies The species of the compendium which the query is being compared to. Default is "hsapiens".
#' @param topExpGenes The number of highest expressed genes from each cell/tissue that will be considered for marker gene status.
#' @param compendiumCommonGenes The list of commonly expressed (non-marker) genes as determined from the compendium, which will be removed from the list of highly expressed sample genes.
#'
#' @importFrom xgsa get_homology_matrix new_XGSA_dataset
#' @importFrom slam row_sums
#'
#' @export
#'
#' @examples
#' Used within C3 framework
#'

select_query_sample_specific_genes<-function(processed.queryData, comSpecies, topExpGenes, compendiumCommonGenes){
  specific.query.data<-list()

  ctnames<-colnames(processed.queryData$expressionData)
  qSpecies<-processed.queryData$species

  #process one cell/tissue data in one loop and then added to the 'specific.query.data' list
  for(i in 1:ncol(processed.queryData$expressionData)){
    expData<-processed.queryData$expressionData[,i]
    names(expData)<-rownames(processed.queryData$expressionData)

    #getting the top number of expressed genes
    topExpressed.expData<-expData[order(expData, decreasing = T)][1:topExpGenes]
    userData<-names(topExpressed.expData)

    #finding the orthologs for the compendium common removal genes for query species
    hkGenes<-compendiumCommonGenes
    if(qSpecies!=comSpecies){
      ## homology matrix with rows=qSpecies and cols=comSpecies
      homolog.mat<-get_homology_matrix(qSpecies, comSpecies)
      col.genes <- unique(as.character(hkGenes[hkGenes%in%colnames(homolog.mat)]))
      hkGenes <- rownames(homolog.mat)[which(row_sums(homolog.mat[,col.genes, drop=FALSE])>0)]
    }

    #removing the compendium common removal genes from the top expressed genes of query data
    userDataFiltered<-setdiff(userData, hkGenes)
    userUniverse<-names(expData)

    #making ready the query data for testing
    specific.query.data[[ctnames[i]]] <- new_XGSA_dataset(qSpecies, data = list(userGenes = userDataFiltered), type = 'genesetlist', name = 'userGenes', universe = unique(userUniverse))
  }

  return(specific.query.data)
}






##############################################################################################################################
# Making compendium data from the gene expression specific matrix
# takes the matrix as input and then build lists for each column of the matrix
# the list has the same number of elements as the matrix columns
# each element of the list contains the same number of genes of the matrix columns expressed genes, i.e., =1

#' @title convert_compendium_to_gene_sets
#'
#' @description
#' This function converts a marker gene matrix into marker gene gene-sets and gene universe.
#'
#' @rdname convert_compendium_to_gene_sets
#' @name convert_compendium_to_gene_sets
#'
#' @details
#' This function converts a marker gene matrix into marker gene gene-sets and gene universe.
#'
#' @return This function returns a list  with two sub-list consisting of marker gene gene-sets and the gene universe for each cell type or tissue, with elements of the sub-lists corresponding to columns of the input matrix.
#'
#' @param specific.compendium A marker gene matrix as returned by 'make_marker_gene_compendium' (i.e., marker.genes.matrix). That is, a matrix containing the presence (1) / absence (0) information for marker genes in each cell type / tissue. Rows denote genes, columns denote individual tissues or cell types.
#' @export
#'
#' @examples
#' specific.compendium<-matrix(sample(0:1, 100, replace=T), 10, 10)
#' rownames(specific.compendium)<-c("ENSG00000100053", "ENSG00000109846", "ENSG00000244752", "ENSG00000138083", "ENSG00000181449", "ENSG00000160202", "ENSG00000007372", "ENSG00000117707", "ENSG00000101144", "ENSG00000277354")
#' colnames(specific.compendium)<-c("cell1", "cell2", "cell3", "cell4", "cell5", "cell6", "cell7", "cell8", "cell9", "cell10")
#' convert_compendium_to_gene_sets(specific.compendium)
#'

convert_compendium_to_gene_sets<-function(specific.compendium){
  column.names<-colnames(specific.compendium)
  explist<-list()

  #takes the matrix as input and then build lists for each column of the matrix
  #the list has the same number of elements as the number of matrix columns
  for(i in 1:ncol(specific.compendium)){explist[[column.names[i]]]<-specific.compendium[,i]}

  compendium.universe<-list()
  #taking the names of the genes (both values for 1 and 0) and making the compendium universe for each cell/tissue
  for(i in 1:length(explist)){compendium.universe[[column.names[i]]]<-names(explist[[i]])}

  #making each element of the list contains only the genes that has only value 1
  for(i in 1:length(explist)){explist[[i]]<-which(explist[[i]]==1)}

  compendium.data<-list()
  #taking the names of the genes and making the compendium data for each cell/tissue
  for(i in 1:length(explist)){compendium.data[[column.names[i]]]<-names(explist[[i]])}

  compendium.gene.sets<-list()
  compendium.gene.sets[["data"]]<-compendium.data
  compendium.gene.sets[["universe"]]<-compendium.universe
  return(compendium.gene.sets)
}






################################################################################################################################
# perform test using the XGSA's  "run_XGSA_test" function

#' @title c3_test
#'
#' @description
#' This function performs the XGSA Fisher's exact tests between the marker genes of the query data and each cell type or tissue in the marker gene compendium.
#'
#' @rdname c3_test
#' @name c3_test
#'
#' @details
#' This function performs the XGSA Fisher's exact tests between the marker genes of the query data and each cell type or tissue in the marker gene compendium.
#'
#' @return This function returns a list containing the XGSA p-values and the corresponding overlapping genes from the query data and the compendium. Each element of the list corresponds to one query sample vs. one compendium cell type.
#'
#' @param processed.queryData A list containing the average gene expression values for the cell types and/or tissues with ensembl gene IDs and the species name, as returned by 'preprocess_querydata'.
#' @param marker.gene.compendium A list containing the data of marker gene compendium, as returned by 'make_marker_gene_compendium'.
#' @param min The minimum number of marker genes required in a marker gene set to be tested. Setting this higher may reduce false positives. Default is 0.
#' @param max The maximum number of genes allowed in a marker gene set to be tested. Default is '4000'.
#' @param top.expressed.genes The number of highest expressed genes from each cell/tissue that will be considered for marker gene status. Default is NULL, i.e., set the same number internally as the compendium.
#'
#' @importFrom xgsa new_XGSA_dataset run_XGSA_test
#'
#' @export
#'
#' @examples
#' ## Here we will use "human.encode.data" and "human.encode.data.descriptor" from C3 repository to make the compendium.
#' ## These data sets are loaded automatically with the package.
#' ## The gene expression data sets are in list format. Here every list element contains transcriptomic profile data of a cell or tissue type.
#' human.compendium<-make_gene_expression_compendium(gene.expression.data=human.encode.data, experiment.descriptor=human.encode.data.descriptor, expression.data.format="list")
#'
#' ## Then we will add "hawse.human.lens.data" to the compendium. This data is also loaded with the package.
#' ## This data set has 2 different cell types - lens epithelial cell (LEC) and lens fiber cell (LFC). Each cell type has 3 biological replicates.
#' Hawse.data.descriptor<-c("Hawse_LEC","Hawse_LEC","Hawse_LEC","Hawse_LFC","Hawse_LFC","Hawse_LFC")
#' Hawse.human.compendium<-add_sample_into_compendium(compendium.data=human.compendium, sample.data=hawse.human.lens.data, species = "hsapiens", data.format = "matrix", geneID = "external_gene_name", experiment.descriptor = Hawse.data.descriptor)
#'
#' ## Next we will make the marker gene compendium from the Hawse.human.compendium
#' human.marker.gene.compendium<-make_marker_gene_compendium(compendium.data =Hawse.human.compendium, specific.cutoff = 0.05, top.expressed.genes = 500)
#'
#' ## Now we will pre-process the "hoang.mouse.lens.data" (query data). This data is also loaded with the package.
#' ## This data set also contains 2 separate cell type: LEC and LFC. Each cell contains average expression value of 3 replicates.
#' Hoang.processed.query.data<-preprocess_querydata(cell.tissue.data = hoang.mouse.lens.data, species = "mmusculus", data.format = "matrix", geneID = "external_gene_name")
#'
#' ## Finally we will perform the test of the query data with the human.marker.gene.compendium
#' Hoang.data.test.result<-c3_test(processed.queryData = Hoang.processed.query.data, marker.gene.compendium = human.marker.gene.compendium)
#' head(sort(Hoang.data.test.result$Hoang_LEC$pvalue))
#' head(sort(Hoang.data.test.result$Hoang_LFC$pvalue))
#'

c3_test<-function(processed.queryData, marker.gene.compendium, min=1, max=5000, top.expressed.genes=NULL){
  comSpecies = marker.gene.compendium$species

  #make query data ready for testing
  if(is.null(top.expressed.genes)) top.expressed.genes = marker.gene.compendium$top.expressed.genes
  compendiumCommonGenes<-marker.gene.compendium$compendium.common.genes
  queryData.forTest <- select_query_sample_specific_genes(processed.queryData, comSpecies, top.expressed.genes, compendiumCommonGenes)

  #make the specific compendium ready for testing according to the XGSA format
  compendium.gene.sets<-convert_compendium_to_gene_sets(marker.gene.compendium$marker.genes.matrix)
  marker.gene.compendium.XGSA <- new_XGSA_dataset(comSpecies, data = compendium.gene.sets$data, type = 'genesetlist', name = "CompendiumData", universe = unique(unlist(compendium.gene.sets$universe)))
  compendiumData.forTest<-marker.gene.compendium.XGSA[1:5]


  individual.test.result<-list()
  ctnames<-colnames(processed.queryData$expressionData)

  #in each loop, one query cell/tissue type data is tested with every cell/tissue type of the compendium
  #the result for each query cell/tissue type is stored as a list element in the 'individual.test.result'
  for(i in 1:ncol(processed.queryData$expressionData)){
    #testing the query data with the compendium
    query.vs.compendium <- run_XGSA_test(queryData.forTest[[i]], compendiumData.forTest, min=min, max=max)

    #getting the p-values of the test result for every cell/tissue type of the compendium
    resulting.pvals <- lapply(query.vs.compendium, function(X){ X[["pvals"]] })
    #getting the overlap genes of the test result for every cell/tissue type of the compendium
    resulting.overlap.genes <- lapply(query.vs.compendium, function(X){ X[["genes"]] })

    #perform Benjamini Hochberg multiple hypothesis testing correction to the pvalues
    adjusted.pvals <- p.adjust(unlist(resulting.pvals), method = "BH")

    #make the names of the results interpretable for humans
    names(adjusted.pvals) <- unlist(lapply(strsplit(names(adjusted.pvals) ,"\\."), function(X){return(X[[2]])}))

    #making a list with adjusted p-values and overlap genes for every cell/tissue type of the compendium
    test.result<-list()
    test.result[["pvalue"]]<-adjusted.pvals
    test.result[["overlapGenes"]]<-resulting.overlap.genes$userGenes

    #finally store the test result as an individual element at the 'individual.test.result'
    individual.test.result[[ctnames[i]]]<-test.result
  }

  return(individual.test.result)
}






####################################################################################################################################
# display the overlapped pvalue by making -log10, reverse, sorted, i.e., show the lowest pvalue with the highest number and in top position

#' @title c3_display_result
#'
#' @description
#' This function displays the C3 results (gene set enrichment p-values for overlap between each query sample and the compendium data) in a barplot.
#'
#' @rdname c3_display_result
#' @name c3_display_result
#'
#' @details
#' This function displays the most significant C3 results (gene set enrichment p-values for overlap between each query sample and the compendium data) in a barplot. X axis is -log10 p-value, Y axis shows compendium cell types tested.
#'
#' @return NULL
#'
#' @param c3.result C3 / XGSA testing results, as returned by 'c3_test'.
#'
#' @export
#'
#' @examples
#' ## Here we will use "human.encode.data" and "human.encode.data.descriptor" from C3 repository to make the compendium.
#' ## These data sets are loaded automatically with the package.
#' ## The gene expression data sets are in list format. Here every list element contains transcriptomic profile data of a cell or tissue type.
#' human.compendium<-make_gene_expression_compendium(gene.expression.data=human.encode.data, experiment.descriptor=human.encode.data.descriptor, expression.data.format="list")
#'
#' ## Then we will add "hawse.human.lens.data" to the compendium. This data is also loaded with the package.
#' ## This data set has 2 different cell types - lens epithelial cell (LEC) and lens fiber cell (LFC). Each cell type has 3 biological replicates.
#' Hawse.data.descriptor<-c("Hawse_LEC","Hawse_LEC","Hawse_LEC","Hawse_LFC","Hawse_LFC","Hawse_LFC")
#' Hawse.human.compendium<-add_sample_into_compendium(compendium.data=human.compendium, sample.data=hawse.human.lens.data, species = "hsapiens", data.format = "matrix", geneID = "external_gene_name", experiment.descriptor = Hawse.data.descriptor)
#'
#' ## Next we will make the marker gene compendium from the Hawse.human.compendium
#' human.marker.gene.compendium<-make_marker_gene_compendium(compendium.data =Hawse.human.compendium, specific.cutoff = 0.05, top.expressed.genes = 500)
#'
#' ## Now we will pre-process the "hoang.mouse.lens.data" (query data). This data is also loaded with the package.
#' ## This data set also contains 2 separate cell type: LEC and LFC. Each cell contains average expression value of 3 replicates.
#' Hoang.processed.query.data<-preprocess_querydata(cell.tissue.data = hoang.mouse.lens.data, species = "mmusculus", data.format = "matrix", geneID = "external_gene_name")
#'
#' ## Next we will perform the test of the query data with the human.marker.gene.compendium
#' Hoang.data.test.result<-c3_test(processed.queryData = Hoang.processed.query.data, marker.gene.compendium = human.marker.gene.compendium)
#'
#' ## Finally we will display the test result in a bar plot where the top result shows the matched query data
#' par(mar=c(8,10,5,5))
#' c3_display_result(Hoang.data.test.result)
#'

c3_display_result<-function(c3.result){
  #in each loop, the test result of one query cell/tissue type is processed and showed in the bar plot
  for(i in 1:length(c3.result)){
    cell.tissue.names<-names(c3.result[i])
    result.pvalue<-c3.result[[i]]$pvalue
    #getting sorted p-values
    sorted.result.pvalue<-sort(result.pvalue)
    #getting significant sorted p-values
    significant.pvalue<-sorted.result.pvalue[sorted.result.pvalue < 0.05]

    #for significant test results
    if(length(significant.pvalue)>0){
      #for setting the title
      if(nchar(cell.tissue.names)>40){
        title<-paste("C3 result of:\n", cell.tissue.names, sep = " ")
      }
      else{
        title<-paste("C3 result of", cell.tissue.names, sep = " ")
      }
      #making negative log10 of the significant sorted p-values
      significant.pvalue.log <- -log10(significant.pvalue)
      #taking reverse of the negative log10 p-values and displayed in the bar plot
      barplot(rev(significant.pvalue.log), main = title, xlab = "-log10(p-value)", las=2, horiz=TRUE, cex.names = 0.7)
    }
    #for non-significant test results
    else{
      #for setting the title
      if(nchar(cell.tissue.names)>30){
        title<-paste("No significant test result..\n Top 10 C3 result of:\n", cell.tissue.names, sep = " ")
      }
      else{
        title<-paste("No significant test result..\n Top 10 C3 result of", cell.tissue.names, sep = " ")
      }
      #taking the top 10 p-values from the sorted resultant p-values
      top10.sorted.result.pvalue<-sorted.result.pvalue[1:10]
      #making negative log10 of the top 10 p-values
      top10.sorted.result.pvalue.log <- -log10(top10.sorted.result.pvalue)
      #taking reverse of the negative log10 p-values and displayed in the bar plot
      barplot(rev(top10.sorted.result.pvalue.log), main = title, xlab = "-log10(p-value)", las=2, horiz=TRUE, cex.names = 0.7)
    }
    # for console message for each cell/tissue type ploting
    consol.msg<-paste(cell.tissue.names, "-- result plotting done!!", sep = " ")
    print(consol.msg)
  }
}

