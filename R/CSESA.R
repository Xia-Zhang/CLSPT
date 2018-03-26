#' CSESA (CRISPR-based Salmonella enterica Serotype Analyzer).
#'
#' @description The main function in CSESA package.
#'
#' @param in.file1 The first input file, the default value is NULL.
#' @param in.file2 The second input file (optional), the default value is NULL.
#' @param out.file Into which results will be saved if this value is set. Otherwise results will be displayed on the screen.
#' @param method The method to handle the input file(s), can only be "PCR" or "WGS".
#'
#' @return The subset framedata which indicate the predicted serotype of salmonella.
#' 
#' @note If you use the "WGS" method, please make sure you have install the BLAST software.
#'
#' @examples
#' \dontrun{
#'   CSESA("./testdata/1.fasta", "./testdata/2.fasta")
#' }
#' @importFrom utils read.table
#' @export
#' 
CSESA <- function(in.file1 = NULL, in.file2 = NULL, out.file = NULL, method = c("PCR", "WGS")) {
    tryCatch({
        if (is.null(in.file1) && is.null(in.file2)) {
            stop("No such files!")
        } 
        method <- toupper(method)
        method <- match.arg(method)
        
        if (method == "PCR") {
            seq1 <- ReadInFile(in.file1)
            seq2 <- ReadInFile(in.file2)
            PCR(seq1, seq2, out.file)
        }
        else {
            # read WGS
            # different, shouldn't union
            # WGS
        }
        
    }, error = function(e) {
        cat("ERROR :",conditionMessage(e),"\n")
    })
}


PCR <- function(seq1, seq2, out.file) {
    csesa.result <- list()
    
    csesa.result$spacer1 <- GetAllNewSpacers(seq1)
    csesa.result$spacer2 <- GetAllNewSpacers(seq2)
    
    csesa.result$serotype <- FindSerotype(csesa.result$spacer1, csesa.result$spacer2)
    class(csesa.result) <- "CSESA"
    
    if (is.null(out.file)) {
        PrintClspt(csesa.result)
    }
    else {
        save(csesa.result, file = out.file)
    }
}

WGS <- function(data) {
    path <- Sys.which("blastn")
    if (all(path == "")) {
        stop("Blast is not exist!")
    }
    
    # loading the database
    db <- ""
    blastcmd <- Sys.which("blastdbcmd")
    
    tmpwd <- tempdir()
    curwd <- getwd()
    tmp.prefix <- basename(tempfile(tmpdir = wd))
    on.exit({
        file.remove(Sys.glob(paste(tmp.prefix, "*")))
    })
    
    infile <- paste(temp_file, ".fasta")
    outfile <- paste(temp_file, ".out")
    
    writeXStringSet(data, infile, append=FALSE, format="fasta")
    system(paste(blastcmd, "-db", db, "-query", infile, "-out", outfile, "-outfmt 10", "-task blastn-short"), ignore.stdout = TRUE, ignore.stderr = FALSE)
    
    result.table <- read.table(outfile, sep=",", quote = "")
    colnames(result.table) <- c("QueryID",  "SubjectID", "Perc.Ident",
                                "Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
                                "S.start", "S.end", "E", "Bits" )
    
    # pick the useful sequence
}

#' Get the new spacers from the molecular sequence and its reverse complement.
#'
#' @param molecular.seq The molecular sequence.
#' @return The vector of the new spacers, which is extracted from the molecular sequence and its reverse complement.
#'
#' @note If there doesn't exist any new spacer, the function would return NA.
#'
GetAllNewSpacers <- function(molecular.seq = NULL) {
    if (is.null(molecular.seq) || is.na(molecular.seq) || molecular.seq == "") 
        return (NA)
    molecular.seq.rev <- GetReverseComplement(molecular.seq)
    
    # handle the Typhi special case
    typhi <- "ACGGCTATCCTTGTTGACGTGGGGAATACTGCTACACGCAAAAATTCCAGTCGTTGGCGCA"
    if (endsWith(molecular.seq, typhi) || endsWith(molecular.seq.rev, typhi))
        return c("EntB0var1")
    
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
#' @param csesa1 The new spacer of the first sequence.
#' @param csesa2 The new spacer of the second sequence.
#' 
#' @return The data frame which represents the serotype.
#'
FindSerotype <- function(csesa1 = NA, csesa2 = NA) {
    if (is.na(csesa1) == TRUE && is.na(csesa2) == TRUE) {
        stop("Sorry. We did not find any corresponding serotype in the lib!")
    }
    V1 = V2 = V3 = V4 = NULL
    mapping.table <- read.table(system.file("extdata", "mapping_tbl.txt", package = "CSESA"), sep = "\t")
    if (is.na(csesa1) == TRUE || is.na(csesa2) == TRUE) {
        csesa <- csesa1
        if (is.na(csesa1))
            csesa <- csesa2
        serotype <- subset(mapping.table, is.element(V1, csesa) | is.element(V2, csesa), select = V3)
        if (nrow(serotype) > 1)
            serotype <- unique(serotype)
    }
    else {
        serotype <- subset(mapping.table, is.element(V1, csesa1) & is.element(V2, csesa2) | 
                   is.element(V1, csesa2) & is.element(V2, csesa1), select = c(V3, V4))
    }
    return (serotype)
}

#' Get the new spacer from the molecular sequence and map it to the code.
#'
#' @param molecular.seq The molecular sequence.
#' @return The new spacer code as a string.
#'
GetNewSpacerCode <- function(molecular.seq = NULL) {
    if (is.null(molecular.seq))
        return (NULL)
    new.spacer <- GetNewSpacer(molecular.seq)
    if (is.null(new.spacer))
        return (NULL)
    V1 = V3 = NULL
    spacers.table <- read.table(system.file("extdata", "spacer_tbl.txt", package = "CSESA"))
    spacer.char <- as.character(subset(spacers.table, V3 == new.spacer, select = V1)[1, 1])
}


#' Get the new spacer from the molecular sequence.
#'
#' @param molecular.seq The molecular sequence.
#' @return The new spacer sequence as a string.
#'
#' @examples
#' GetNewSpacer("AGAGGCGGACCGAAAAACCGTTTTCAGCCAACGTAT")
#'
#' @export
GetNewSpacer <- function(molecular.seq = NULL) {
    if (is.null(molecular.seq))
        return (NULL)
    max.match <- "-"
    max.count <- -1
    dr.table <- read.table(system.file("extdata", "DR_tbl.txt", package = "CSESA"))
    for (x in dr.table$V3[-1]) {
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


#' Print the result of CSESA.
#'
#' @param  csesa The S3 object CSESA.
#'
PrintClspt <- function(csesa) {
    if (is.null(csesa)) {
        stop("The csesa object should be set!")
    }
    
    print(paste("The new spacer in the first sequence:", csesa$spacer1))
    print(paste("The new spacer in the second sequence:", csesa$spacer2))
    
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


#' Read the three types of input file.
#'
#' @param file.name The input file name.
#' @return The molecular sequence as a string.
#'
ReadInFile <- function(file.name) {
    if (is.null(file.name) || is.na(file.name))
        return (NA)
    if (file.exists(file.name) == FALSE) {
        stop(paste0(file.name, " is not existed!"))
    }
    data <- scan(file.name, what = "", quiet = TRUE)
    if (substring(data[1], 1, 1) == '>')
        data <- data[-1]
    toupper(paste(data, collapse = ""))
}


#' Return the reverse complement of the sequence.
#'
#' @param x The input sequence.
#' @return The reverse complement sequence as a string.
#'
GetReverseComplement <- function(x) {
    a <- chartr("ATGC","TACG",x)
    paste(rev(substring(a, 1 : nchar(a), 1 : nchar(a))), collapse = "")
}
