soundex <- function(x){
    ## Taken from Splus manual.  not sure it works perfectly in R but it doesn't appear to be bad
    
    ## 1. extract the last word of surnames and translate to all upper case
    base <- gsub("[^A-Z]", "", toupper(gsub("^.*[ \t]","", gsub("[ \t]*$", "", x))))
    ## 2. encode the surnames (last word) using the soundex algorithm
    basecode <- gsub("[AEIOUY]", "", gsub("[R]+", "6", gsub("[MN]+", "5", gsub("[L]+", "4", gsub("[DT]+", "3", gsub("[CGJKQSXZ]+", "2", gsub("[BFPV]+", "1", gsub("[HW]", "", base))))))))
    ## 3. deal with the 1st letter and generate the final coding padded with 0
    sprintf("%4.4s", paste(substring(base, 1, 1), ifelse(regexpr("^[HWAEIOUY]", base) == 1, basecode, substring(basecode, 2)), "000", sep = ""))
}
