#' CLSPT function
#'
#' @description The main function in CLSPT package
#'
#' @param in.file1 one of the input sequence file, the default value is NULL.
#' @param in.file2 the other one of the input sequence file, the default value is NULL.
#' @param out.file the file where the result will save. If not set the value, the result will show in the screen.
#'
#' @return the subset framedata which indicate the predictation of salmonella.
#'
#' @example
#'  CLSPT("./1.fasta", "./2.fasta")
#'
#' @importFrom utils read.table
#' @export
CLSPT <- function(in.file1 = NULL, in.file2 = NULL, out.file = NULL) {
    if (is.null(in.file1) && is.null(in.file2)) {
        # shiny
        stop("No such files!")
    }
    InitGlobal()
    # CLSPT("./testdata/1.fasta", "./testdata/2.fasta")
    file1.subset <- GetMapIterm(in.file1)
    file2.subset <- GetMapIterm(in.file2)
    data.result <- merge(file1.subset, file2.subset)
    if (is.null(out.file)) {
        print(data.result)
    }
    else {
        save(data.result, file = out.file)
    }
}


#' Return the mapping interval_seq_ID
#'
#' @param file.name the input file name
#' @return the subset framedata which include the input new spacer code of the sequence or the reverse complement sequence
#'
#' @example
#' GetMapIterm("./testdata/1.fasta")
#'
GetMapIterm <- function(file.name = NULL) {
    if (is.null(file.name))
        return (mapping.table.global)
    molecular.seq <- ReadInFile(file.name)
    molecular.seq.rev <- GetReverseComplement(molecular.seq)
    if (is.null(molecular.seq))
        return (mapping.table.global)
    new.spacer <- GetNewSpacerCode(molecular.seq)
    new.spacer.rev <- GetNewSpacerCode(molecular.seq.rev)
    new.spacer.arr <- character()
    if (is.null(new.spacer) == FALSE && is.na(new.spacer) == FALSE)
        new.spacer.arr <- c(new.spacer.arr, new.spacer)
    if (is.null(new.spacer.rev) == FALSE && is.na(new.spacer.rev) == FALSE)
        new.spacer.arr <- c(new.spacer.arr, new.spacer.rev)
    if (length(new.spacer.arr) == 0)
        return (mapping.table.global)
    subset(mapping.table.global, is.element(V1, new.spacer.arr) | is.element(V2, new.spacer.arr))
}


#' Get the new spacer from the molecular sequence and map it to the code
#'
#' @param molecular.seq the molecular sequence
#' @return the string which is the new spacer code\
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
#' GetNewSpacer(GCGCCGGGAACACCAACGTCGGTTTATCCCCGCTGGCGCCGGGAACACAGGCGGACCGAAAAACCGTTTTCAGCCAACGTCGGTTTATCCCCGCTGGCGCCGGGAACACCAACGTCGGTTT)
#'
#' @export
GetNewSpacer <- function(molecular.seq = NULL) {
    if (is.null(molecular.seq))
        return (NULL)
    molecular.seq <- toupper(molecular.seq)
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
    spacers[length(spacers) - 1]
}

#' Check if the new spacer exist in the sequence
#'
#' @param x the input sequence
#' @return TRUE or FALSE to represent the existance
#'
FindNewSpacer <- function(molecular.seq = NULL) {
    if (is.null(molecular.seq))
        return (FALSE)
    molecular.seq <- toupper(molecular.seq)
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
        return (FALSE)
    spacers <- strsplit(molecular.seq, max.match)[[1]]
    new.spacer <- spacers[length(spacers) - 1]
    return (is.element(new.spacer, spacers.table.global$V3))
}

#' Read the three types of input file
#'
#' @param file.name the input file name
#' @return the string which represents the molecular sequence
#'
ReadInFile <- function(file.name) {
    data <- readLines(file.name)
    if (substring(data[1], 1, 1) == '>')
        data <- data[-1]
    paste(data, collapse = "")
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


#' Setting the global varible reading from files
#'
InitGlobal <- function() {
    map.file <- "./inst/extdata/mapping_tbl.txt"
    # map.file <- "./mapping_tbl.txt" # in the package
    DR.file <- "./inst/extdata/dr_tbl.txt"
    # DR.file <- "./dr_tbl.txt" # in the package
    spacer.file <- "./inst/extdata/spacer_tbl.txt"
    # spacer.file <- "./spacer_tbl.txt" # in the package

    if (exists("mapping.table.global") == FALSE)
        mapping.table.global <<- read.table(map.file, sep = "\t")
    if (exists("patterns.table.global") == FALSE)
        dr.table.global <<- read.table(DR.file)
    if (exists("spacers.table.global") == FALSE)
        spacers.table.global <<- read.table(spacer.file)

}
