#' CLSPT function
#'
#' @description The main function in CSESA package
#'
#' @param in.file1 one of the input sequence file, the default value is NULL.
#' @param in.file2 the other one of the input sequence file, the default value is NULL.
#' @param out.file the file where the result will save. If not set the value, the result will show in the screen.
#'
#' @return the subset framedata which indicate the predictation of salmonella.
#'
#' @example
#'  CSESA("./testdata/1.fasta", "./testdata/2.fasta")
#'
#' @importFrom utils read.table
#' @export
#' 
CSESA <- function(in.file1 = NULL, in.file2 = NULL, out.file = NULL) {
    if (is.null(in.file1) && is.null(in.file2)) {
        # shiny
        print("No such files!")
        exit()
    }
    
    InitGlobal()
    
    csesa.result <- list()
    seq1 <- ReadInFile(in.file1)
    seq2 <- ReadInFile(in.file2)
    
    csesa.result$spacer1 <- GetAllNewSpacers(seq1)
    csesa.result$spacer2 <- GetAllNewSpacers(seq2)
    
    csesa.result$serotype <- FindSerotype(csesa.result$spacer1, csesa.result$spacer2)
    class(csesa) <- "CSESA"
    
    if (is.null(out.file)) {
        PrintClspt(csesa.result)
    }
    else {
        save(csesa.result, file = out.file)
    }
}

#' Get the new spacers from the molecular sequence and its reverse complement sequence
#'
#' @param molecular.seq the molecular sequence
#' @return The vector of the new spacers, which is extracted from the molecular sequence and its reverse complement sequence
#'
#' @note if there doesn't exist any new spacer, the function would return the whole mapping table.
#' @example
#' GetAllNewSpacers("GCGCCGGGAACACCAACGTCGGTTTATCCCCGCTGGCGCCGGGAACACAGGCGGACCGAAAAACCGTTTTCAGCCAACGTCGGTTTATCCCCGCTGGCGCCGGGAACACCAACGTCGGTTT")
#'
GetAllNewSpacers <- function(molecular.seq = NULL) {
    if (is.null(molecular.seq) || is.na(molecular.seq) || molecular.seq == "") 
        return (NA)
    molecular.seq.rev <- GetReverseComplement(molecular.seq)
    new.spacer <- GetNewSpacerCode(molecular.seq)
    new.spacer.rev <- GetNewSpacerCode(molecular.seq.rev)
    new.spacer.arr <- character()
    if (is.null(new.spacer) == FALSE && is.na(new.spacer) == FALSE)
        new.spacer.arr <- c(new.spacer.arr, new.spacer)
    if (is.null(new.spacer.rev) == FALSE && is.na(new.spacer.rev) == FALSE)
        new.spacer.arr <- c(new.spacer.arr, new.spacer.rev)
    if (length(new.spacer.arr) == 0) {
        return (NA)
    }
    else return (new.spacer.arr)
}


#' Find the serotype based on the analysis of the new spacers.
#'
#' @param csesa1 the new spacers of the first sequence
#' @param csesa2 the new spacers of the second sequence
#' 
#' @return the data frame which reprents the serotype
#'
FindSerotype <- function(csesa1 = NA, csesa2 = NA) {
    if (is.na(csesa1) == TRUE && is.na(csesa2) == TRUE) {
        print("Sorry. We did not find any corresponding serotype in the lib!")
        exit()
    }
    V1 = V2 = V3 = V4 = NULL
    if (is.na(csesa1) == TRUE || is.na(csesa2) == TRUE) {
        csesa <- csesa1
        if (is.na(csesa1))
            csesa <- csesa2
        serotype <- subset(mapping.table.global, is.element(V1, csesa) | is.element(V2, csesa), select = V3)
        if (nrow(serotype) > 1)
            serotype <- unique(serotype)
    }
    else {
        serotype <- subset(mapping.table.global, is.element(V1, csesa1) & is.element(V2, csesa2) | 
                   is.element(V1, csesa2) & is.element(V2, csesa1), select = c(V3, V4))
    }
    return (serotype)
}

#' Get the new spacer from the molecular sequence and map it to the code
#'
#' @param molecular.seq the molecular sequence
#' @return the string which is the new spacer code
#'
GetNewSpacerCode <- function(molecular.seq = NULL) {
    if (is.null(molecular.seq))
        return (NULL)
    new.spacer <- GetNewSpacer(molecular.seq)
    if (is.null(new.spacer))
        return (NULL)
    spacer.char <- as.character(subset(spacers.table.global, V3 == new.spacer, select = V1)[1, 1])
}


#' Get the new spacer from the molecular sequence
#'
#' @param molecular.seq the molecular sequence
#' @return the string which is the new spacer
#'
#' @example
#' GetNewSpacer("GCGCCGGGAACACCAACGTCGGTTTATCCCCGCTGGCGCCGGGAACACAGGCGGACCGAAAAACCGTTTTCAGCCAACGTCGGTTTATCCCCGCTGGCGCCGGGAACACCAACGTCGGTTT")
#'
#' @export
GetNewSpacer <- function(molecular.seq = NULL) {
    if (is.null(molecular.seq))
        return (NULL)
    InitGlobal()
    max.match <- "-"
    max.count <- -1
    for (x in dr.table.global$V3[-1]) {
        t <- gregexpr(pattern = x, text = molecular.seq)
        if (t[[1]][1] != -1 && length(t[[1]]) > max.count) {
            max.match = x
            max.count = length(t[[1]])
        }
    }
    if (max.count < 1)
        return (NULL)
    spacers <- strsplit(molecular.seq, max.match)[[1]]
    if (substr(molecular.seq, nchar(molecular.seq) - nchar(max.match) + 1, nchar(molecular.seq)) == max.match) {
        return (spacers[length(spacers)])
    }
    spacers[length(spacers) - 1]
}


#' Print the result which is got from CSESA
#'
#' @param  csesa the S3 object CSESA
#'
PrintClspt <- function(csesa) {
    if (is.null(csesa)) {
        print("The csesa object should be set!")
        exit()
    }
    
    print(paste("The new spacer in first sequence:", csesa$spacer1))
    print(paste("The new spacer in second sequence:", csesa$spacer2))
    
    if (is.na(csesa$spacer1) || is.na(csesa$spacer2)) {
        result <- ""
        if (is.atomic(csesa$serotype))
            result <- csesa$serotype
        else 
            result <- paste(csesa$serotype[, 1], collapse = "] [")
        result <- paste("The possible serotype: [", result, sep = "")
        print(paste(result, "]", sep = ""))
    }
    else {
        result <- paste(csesa$serotype[, 1], csesa$serotype[, 2])
        result <- paste("The possible serotype: [", paste(result, collapse = "] ["), sep = "")
        print(paste(result, "]", sep = ""))
    }
}


#' Read the three types of input file
#'
#' @param file.name the input file name
#' @return the string which represents the molecular sequence
#'
ReadInFile <- function(file.name) {
    if (is.null(file.name) || is.na(file.name))
        return (NA)
    if (file.exists(file.name) == FALSE) {
        print(paste0(file.name, " is not existed!"))
        exit()
    }
    data <- scan(file.name, what = "", quiet = TRUE)
    if (substring(data[1], 1, 1) == '>')
        data <- data[-1]
    toupper(paste(data, collapse = ""))
}


#' Return the reverse complement of the sequence
#'
#' @param x the input sequence
#' @return the reverse complement sequence string
#'
GetReverseComplement <- function(x) {
    a <- chartr("ATGC","TACG",x)
    paste(rev(substring(a, 1 : nchar(a), 1 : nchar(a))), collapse = "")
}


#' The exit function which help to exit from the script
#' 
exit <- function() {
    .Internal(.invokeRestart(list(NULL, NULL), NULL))
}


#' Setting the global varible reading from files
#' 
#' importFrom("utils", "read.table")
#' 
InitGlobal <- function() {
    # map.file <- "E:/code/CSESA/inst/extdata/mapping_tbl.txt"
    map.file <- system.file("extdata", "mapping_tbl.txt", package = "CSESA")
    # DR.file <- "E:/code/CSESA/inst/extdata/DR_tbl.txt"
    DR.file <- system.file("extdata", "DR_tbl.txt", package = "CSESA")
    # spacer.file <- "E:/code/CSESA/inst/extdata/spacer_tbl.txt"
    spacer.file <- system.file("extdata", "spacer_tbl.txt", package = "CSESA")

    if (exists("mapping.table.global") == FALSE)
        mapping.table.global <<- read.table(map.file, sep = "\t")
    if (exists("patterns.table.global") == FALSE)
        dr.table.global <<- read.table(DR.file)
    if (exists("spacers.table.global") == FALSE)
        spacers.table.global <<- read.table(spacer.file)

}
