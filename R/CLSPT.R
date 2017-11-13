#' Mixed graph model estimate
#'
#' @description MixedGraph is the main function of the package, it uses the BRAIL to do the regression.
#'
CLSPT <- function(in.file1 = NULL, in.file2 = NULL, out.file = NULL) {
    if (is.null(in.file1) && is.null(in.file2)) {
        # shiny
    }
    file1.subset <- GetMapIterm(in.file1)
    file2.subset <- GetMapIterm(in.file2)
    merge(file1.subset, file2.subset)
}

GetMapIterm <- function(file.name = NULL) {
    map.file <- "./lib/mapping_tbl.txt"
    mapping.table <- read.table(map.file, sep = "\t")

    if (is.null(file.name))
        return (mapping.table)
    molecular.seq <- ReadInFile(file.name)
    molecular.seq.rev <- GetReverseComplement(molecular.seq)
    if (is.null(molecular.seq))
        return (mapping.table)
    new.spacer <- GetNewSpacerCode(molecular.seq)
    new.spacer.rev <- GetNewSpacerCode(molecular.seq.rev)
    new.spacer.arr <- character()
    if (is.null(new.spacer) == FALSE && is.na(new.spacer) == FALSE)
        new.spacer.arr <- c(new.spacer.arr, new.spacer)
    if (is.null(new.spacer.rev) == FALSE && is.na(new.spacer.rev) == FALSE)
        new.spacer.arr <- c(new.spacer.arr, new.spacer.rev)
    subset(mapping.table, is.element(V1, new.spacer.arr) | is.element(V2, new.spacer.arr))
}

GetNewSpacerCode <- function(molecular.seq = NULL) {
    if (is.null(molecular.seq))
        return (NULL)
    new.spacer <- GetNewSpacer(molecular.seq)
    if (is.null(new.spacer))
        return (NULL)
    spacer.file <- "./lib/spacer_tbl.txt"
    match.spacer <- read.table(spacer.file)
    spacer.char <- as.character(subset(match.spacer, V3 == new.spacer, select = V1)[1, 1])
}

GetNewSpacer <- function(molecular.seq = NULL) {
    if (is.null(molecular.seq))
        return(NULL)
    molecular.seq <- toupper(molecular.seq)
    DR.file <- "./lib/dr_tbl.txt"
    patterns <- read.table(DR.file)
    max.match <- "-"
    max.count <- -1
    for (x in patterns$V3[-1]) {
        t <- gregexpr(pattern = x, text = molecular.seq)
        if (t[[1]][1] != -1 && length(t[[1]]) > max.count) {
            max.match = x
            max.count = length(t[[1]])
        }
    }
    if (max.count < 1)
        return(NULL)
    spacers <- strsplit(molecular.seq, max.match)[[1]]
    spacers[length(spacers) - 1]
}

ReadInFile <- function(file.name) {
    data <- readLines(file.name)
    if (substring(data[1], 1, 1) == '>')
        data <- data[-1]
    paste(data, sep = "")
}

GetReverseComplement <- function(x) {
    rev(chartr("ATGC","TACG",x))
}
