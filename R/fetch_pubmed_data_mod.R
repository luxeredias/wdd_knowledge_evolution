#fetch_pubmed_data_mod is a modified function of the original function in
#the easyPubMed package. The modification allows for more than 500 documents
#to be retrieved from PubMed at a time

fetch_pubmed_data_2 <- function (pubmed_id_list, retstart = 0, retmax = 500, format = "xml", 
          encoding = "UTF8") 
{
  myIDlist <- pubmed_id_list
  if ((!is.list(myIDlist)) | is.na(myIDlist$WebEnv) | is.na(myIDlist$QueryKey) | 
      is.na(myIDlist$Count) | !is.integer(as.integer(retstart)) | 
      !is.integer(as.integer(retmax))) {
    message("There is an issue with the PubMed ID list you supplied. Please, call the function again and supply the result of a <get_pubmed_ids()> call as argument. Thank you.")
    return(NULL)
  }
  else {
    myWebEnv <- myIDlist$WebEnv
    myKey <- myIDlist$QueryKey
    myCount <- as.numeric(as.character(myIDlist$Count))
    myRetstart = as.integer(retstart)
    if (myRetstart < 0) {
      myRetstart = 0
    }
    myRetmax <- as.integer(retmax)
    if (myRetmax > 5000) {
      myRetmax = 5000
    }
    if (myRetmax < 1) {
      myRetmax = 1
    }
    if (format[1] %in% c("medline", "uilist", "abstract", 
                         "asn.1", "xml")) {
      myFormat <- format[1]
    }
    else {
      myFormat <- "xml"
    }
    typeMode <- switch(myFormat, asn.1 = c("null", "asn.1"), 
                       xml = c("null", "xml"), medline = c("medline", "text"), 
                       uilist = c("uilist", "text"), abstract = c("abstract", 
                                                                  "text"))
    efetch_url = paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?", 
                       "db=pubmed&WebEnv=", myWebEnv, "&query_key=", myKey, 
                       "&retstart=", myRetstart, "&retmax=", myRetmax, 
                       "&rettype=", typeMode[1], "&retmode=", typeMode[2], 
                       sep = "")
    api_key <- pubmed_id_list$APIkey
    if (!is.null(api_key)) {
      efetch_url <- paste(efetch_url, "&api_key=", api_key, 
                          sep = "")
    }
    out.data <- NULL
    try_num <- 1
    t_0 <- Sys.time()
    while (is.null(out.data)) {
      if (try_num > 1) 
        Sys.sleep(time = 2)
      t_1 <- Sys.time()
      if (as.numeric(difftime(t_1, t_0, units = "mins")) > 
          100) {
        message("Killing the request! Something is not working. Please, try again later")
        return(NULL)
      }
      out.data <- tryCatch({
        tmpConnect <- suppressWarnings(url(efetch_url, 
                                           open = "rb", encoding = "UTF8"))
        suppressWarnings(readLines(tmpConnect, warn = FALSE, 
                                   encoding = "UTF8"))
      }, error = function(e) {
        NULL
      }, finally = {
        try(suppressWarnings(close(tmpConnect)), silent = TRUE)
      })
      if (!is.null(out.data) && class(out.data) == "character" && 
          grepl("<ERROR>", substr(paste(utils::head(out.data, 
                                                    n = 100), collapse = ""), 1, 250))) {
        out.data <- NULL
      }
      try_num <- try_num + 1
    }
    if (is.null(out.data)) {
      message("Killing the request! Something is not working. Please, try again later")
      return(NULL)
    }
    if (encoding != "UTF8") 
      out.data <- base::iconv(out.data, from = "UTF8", 
                              to = encoding, sub = ".")
    if (format[1] == "xml") {
      out.data <- paste(out.data, collapse = "")
    }
    return(out.data)
  }
}
