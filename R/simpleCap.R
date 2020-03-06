#function to captilize all first letters of a string with multiple words
#by @RevoAndrie in the Rstudio development team
simpleCap <- function(x) {
  x <- tolower(as.character(x))
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}
