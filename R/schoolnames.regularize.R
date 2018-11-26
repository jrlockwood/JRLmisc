schoolnames.regularize <- function(x){
  ## x is a vector of character strings representing school names.
  ## result is a vector of character strings where we try make school
  ## names closer to regularized

  ## replace all punctuation with spaces, squeeze spaces, and lowercase everything
  a <- tolower(gsub('[[:space:]]+', ' ', gsub('[[:punct:]]', ' ', as.character(x))))
  ## split on spaces into list of character vectors
  a <- strsplit(a," ")
  ## drop all words with one or fewer characters
  a <- lapply(a, function(x){x[nchar(x) >= 2]})
  ## delete words that look like "school" and other useless words
  a <- lapply(a , function(x){x[ !(x %in% c("sc","sch","scho","schoo","school","schl","scl","an","the","of","cntr","center")) ]})
  ## make various expansions of words that look like something
  a <- lapply(a , function(x){x[x %in% c("es","el","ele","elem","eleme","elemen","element","elementa","elementar")] <- "elementary"; x})
  a <- lapply(a , function(x){x[x %in% c("ms","mid","midd","middl","midl","mdl")] <- "middle"; x})
  a <- lapply(a , function(x){x[x %in% c("hs","hi","hig")] <- "high"; x})
  a <- lapply(a , function(x){x[x %in% c("lear","learni","learnin","lrng")] <- "learning"; x})
  a <- lapply(a , function(x){x[x %in% c("aca","acad","acade","academ")] <- "academy"; x})
  a <- lapply(a , function(x){x[x %in% c("comm","commu","commun","communi","communit")] <- "community"; x})
  a <- lapply(a , function(x){x[x %in% c("educ","educa","educat","educati","educatio")] <- "education"; x})
  ## paste the strings back with underscore and return as vector
  return(sapply(a, paste, collapse="_"))
}
