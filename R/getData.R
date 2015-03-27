#' This function formats your data for downstream analysis within the package
#' @param infoFile Data.frame of your sample info file.  Samples should be in rows with sample info in columns.
#' @param dataFile Data.frame with your data. Samples should be in rows and gene/cpgs in columns.
#' @return A list which contains your infoFile and Datafile for downstream analysis
#' @export


getData <- function(infoFile, dataFile) {
  #create list of both files
  dataBind <- list(infoFile, dataFile)
  #name list items
  names(dataBind)<- c("info", "data")
  return(dataBind)
}