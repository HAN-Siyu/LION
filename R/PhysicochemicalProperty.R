#' Computation of the AAindex-Based Physicochemical Features of Protein Sequences
#' @description The function \code{computePhysChem_AAindex} computes the physicochemical features of
#' protein sequences with the help of AAindex provided by \code{\link[seqinr]{aaindex}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}).
#' @param seqPro protein sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr".
#' Each sequence should be a vector of single characters.
#' @param entry.index The entry number of AAindex. For example the entry number of the
#' Kyte & Doolittle Hydrophaty index is 151, and the entry number of Hopp-Woods Hydrophilicity value is 115. And
#' the corresponding \code{entry.index} is \code{c(151, 115)} for these two physicochemical indices. See examples and \code{\link[seqinr]{aaindex}} of
#' "seqinr" package.
#' @param Fourier.len postive integer specifying the Fourier series length that will be used as features.
#' The \code{Fourier.len} should be >= the length of the input sequence. Default: \code{10}.
#' @param parallel.cores an integer that indicates the number of cores for parallel computation. Default: \code{2}.
#' Set \code{parallel.cores = -1} to run with all the cores.
#'
#' @return This function returns a data frame.
#'
#'
#' @section References:
#' [1] Kawashima S, Kanehisa M.
#' AAindex: amino acid index database.
#' Nucleic Acids Res. 2000; 28:374
#'
#' [2] Charif D, Lobry JR.
#' SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis.
#' In: Structural approaches to sequence evolution: Molecules, networks, populations. 2007; 207-232
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom seqinr a
#' @importFrom utils data
#'
#' @seealso \code{\link{computePhysChem}}
#' @examples
#' data(demoPositiveSeq)
#' seqs <- demoPositiveSeq$Pro.positive
#' data(aaindex, package = "seqinr")
#'
#' # Check the aaindex list provided by package "seqinr":
#' # ?seqinr::aaindex
#'
#' # Supose that we need the Kyte & Doolittle Hydrophaty index,
#' # and we can find the entries with Kyte as author:
#'
#' which(sapply(aaindex, function(x) length(grep("Kyte", x$A)) != 0))
#'
#' # This returns entry number 151.
#' # And 151 is the entry.index used by function computePhysChem_AAindex.
#' # Users can also obtain the entry number by retrieving information in other
#' # fields such as x$D (Data description) or x$T (Title of the article).
#' # See documentation of seqinr::aaindex for detailed information.
#'
#' feature_aaindex_1 <- computePhysChem_AAindex(seqPro = seqs, entry.index = 151,
#'                                              Fourier.len = 10, parallel.cores = 2)
#'
#' # If you need to compute features with more than aaindex data:
#'
#' feature_aaindex_2 <- computePhysChem_AAindex(seqPro = seqs, entry.index = c(151, 68),
#'                                              Fourier.len = 10, parallel.cores = 2)
#'
#' @export

computePhysChem_AAindex <- function(seqPro, entry.index = NULL,
                                    Fourier.len = 10, parallel.cores = 2) {

        message("+ Initializing...  ", Sys.time())

        if (min(lengths(seqPro)) < Fourier.len) stop("The profile length of Fourier Series (Fourier.len) cannot larger than the minimum length of the input sequences!")

        data(aaindex, package = "seqinr", envir = environment())
        aaindex <- get("aaindex")

        if (is.null(entry.index)) stop("Please enter the entry index (entry.index) for AA Index!")

        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)

        cl <- parallel::makeCluster(parallel.cores)

        parallel::clusterExport(cl, varlist = c("Internal.FourierSeries",
                                                "Internal.convertSeq", "aaindex"), envir = environment())

        message("\n", "+ Calculating the physicochemical information of protein sequences...  ", Sys.time())

        output <- parallel::parSapply(cl, seqPro, Internal.featureAAindex, aaindex = aaindex,
                                      entry.index = entry.index, Fourier.len = Fourier.len)
        parallel::stopCluster(cl)

        message("\n", "+ Completed.  ", Sys.time(), "\n")

        output <- data.frame(t(output))
}

#' Computation of the Physicochemical Features of RNA or Protein Sequences
#' @description The function \code{computePhysChem} computes the physicochemical features of RNA
#' or protein sequences.
#' @param seqs sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of RNA/protein sequences.
#' RNA sequences will be converted into lower case letters, but
#' protein sequences will be converted into upper case letters.
#' Each sequence should be a vector of single characters.
#' @param seqType a string that specifies the nature of the sequence: \code{"RNA"} or \code{"Pro"} (protein).
#' If the input is DNA sequence and \code{seqType = "RNA"}, the DNA sequence will be converted to RNA sequence automatically.
#' Default: \code{"RNA"}.
#' @param Fourier.len postive integer specifying the Fourier series length that will be used as features.
#' The \code{Fourier.len} should be >= the length of the input sequence. Default: \code{10}.
#' @param physchemRNA strings specifying the physicochemical properties that are computed in RNA sequences. Ignored if \code{seqType = "Pro"}.
#' Options: \code{"hydrogenBonding"} for Hydrogen-bonding and \code{"vanderWaal"} for Van der Waal's interaction
#' Multiple elements can be selected at the same time. (Ref: [2])
#' @param physchemPro strings specifying the physicochemical properties that are computed in protein sequences. Ignored if \code{seqType = "RNA"}.
#' Options: \code{"polarity.Grantham"}, \code{"polarity.Zimmerman"}, \code{"bulkiness.Zimmerman"}, \code{"isoelectricPoint.Zimmerman"},
#' \code{"hphob.BullBreese"}, \code{"hphob.KyteDoolittle"}, \code{"hphob.Eisenberg"}, and
#' \code{"hphob.HoppWoods"}.
#' Multiple elements can be selected at the same time. See details below. (Ref: [3-9])
#' @param parallel.cores an integer that indicates the number of cores for parallel computation. Default: \code{2}.
#' Set \code{parallel.cores = -1} to run with all the cores.
#' @param as.list logical. The result will be returned as a list or data frame.
#'
#' @return This function returns a data frame if \code{as.list = FALSE} or returns a list if \code{as.list = TRUE}.
#'
#' @details
#' The default physicochemical properties are selected or derived from tool "catRAPID" (Ref: [10])
#' and "lncPro" (Ref: [11]).
#' In "catRAPID", \code{Fourier.len = 50}; in "lncPro", \code{Fourier.len} is set as \code{10}.
#'
#' \itemize{
#' \item The physicochemical properties of RNA
#'
#' \enumerate{
#' \item Hydrogen-bonding (\code{"hydrogenBonding"}) (Ref: [2])
#' \item Van der Waal's interaction (\code{"vanderWaal"}) (Ref: [2])
#' }
#'
#' \item The physicochemical properties of protein sequence
#'
#' \enumerate{
#' \item polarity \code{"polarity.Grantham"} (Ref: [3])
#' \item polarity \code{"polarity.Zimmerman"} (Ref: [4])
#' \item bulkiness \code{"bulkiness.Zimmerman"} Ref: [4]
#' \item isoelectric point \code{"isoelectricPoint.Zimmerman"} (Ref: [4])
#' \item hydropathicity \code{"hphob.BullBreese"} (Ref: [5])
#' \item hydropathicity \code{"hphob.KyteDoolittle"} (Ref: [6])
#' \item hydropathicity \code{"hphob.Eisenberg"} (Ref: [7])
#' \item hydropathicity \code{"hphob.HoppWoods"} (Ref: [8])
#' }}
#'
#' @section References:
#' [1] Han S, Liang Y, Li Y, \emph{et al}.
#' ncProR: an integrated R package for effective ncRNA-protein interaction prediction.
#' (\emph{Submitted})
#'
#' [2] Morozova N, Allers J, Myers J, \emph{et al}.
#' Protein-RNA interactions: exploring binding patterns with a three-dimensional superposition analysis of high resolution structures.
#' Bioinformatics 2006; 22:2746-52
#'
#' [3] Grantham R.
#' Amino acid difference formula to help explain protein evolution.
#' Science 1974; 185:862-4
#'
#' [4] Zimmerman JM, Eliezer N, Simha R.
#' The characterization of amino acid sequences in proteins by statistical methods.
#' J. Theor. Biol. 1968; 21:170-201
#'
#' [5] Bull HB, Breese K.
#' Surface tension of amino acid solutions: a hydrophobicity scale of the amino acid residues.
#' Arch. Biochem. Biophys. 1974; 161:665-670
#'
#' [6] Kyte J, Doolittle RF.
#' A simple method for displaying the hydropathic character of a protein.
#' J. Mol. Biol. 1982; 157:105-132
#'
#' [7] Eisenberg D, Schwarz E, Komaromy M, \emph{et al}.
#' Analysis of membrane and surface protein sequences with the hydrophobic moment plot.
#' J. Mol. Biol. 1984; 179:125-42
#'
#' [8] Hopp TP, Woods KR.
#' Prediction of protein antigenic determinants from amino acid sequences.
#' Proc. Natl. Acad. Sci. U. S. A. 1981; 78:3824-8
#'
#' [9] Kawashima S, Kanehisa M. AAindex: amino acid index database. Nucleic Acids Res. 2000; 28:374
#'
#' [10] Bellucci M, Agostini F, Masin M, \emph{et al}.
#' Predicting protein associations with long noncoding RNAs.
#' Nat. Methods 2011; 8:444-445
#'
#' [11] Lu Q, Ren S, Lu M, \emph{et al}.
#' Computational prediction of associations between long non-coding RNAs and proteins.
#' BMC Genomics 2013; 14:651
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom seqinr a
#' @importFrom utils data
#'
#' @seealso \code{\link{featurePhysChem}}
#' @examples
#' data(demoPositiveSeq)
#' seqsRNA <- demoPositiveSeq$RNA.positive
#' seqsPro <- demoPositiveSeq$Pro.positive
#'
#' physChemRNA <- computePhysChem(seqs = seqsRNA, seqType = "RNA",
#'                                Fourier.len = 10, as.list = FALSE)
#'
#' physChemPro <- computePhysChem(seqs = seqsPro, seqType = "Pro", Fourier.len = 8,
#'                                physchemPro = c("polarity.Grantham",
#'                                                "polarity.Zimmerman",
#'                                                "hphob.BullBreese",
#'                                                "hphob.KyteDoolittle",
#'                                                "hphob.Eisenberg",
#'                                                "hphob.HoppWoods"),
#'                                as.list = TRUE)
#'
#' @export

computePhysChem <- function(seqs, seqType = c("RNA", "Pro"), Fourier.len = 10,
                            physchemRNA = c("hydrogenBonding", "vanderWaal"),
                            physchemPro = c("polarity.Grantham", "polarity.Zimmerman",
                                            "bulkiness.Zimmerman", "isoelectricPoint.Zimmerman",
                                            "hphob.BullBreese", "hphob.KyteDoolittle",
                                            "hphob.Eisenberg", "hphob.HoppWoods"),
                            parallel.cores = 2, as.list = TRUE) {

        message("+ Initializing...  ", Sys.time())

        if (min(lengths(seqs)) < Fourier.len) stop("The profile length of Fourier Series (Fourier.len) cannot larger than the minimum length of the input sequences!")

        data(aaindex, package = "seqinr", envir = environment())
        aaindex <- get("aaindex")

        seqType <- match.arg(seqType)
        physchemRNA <- match.arg(physchemRNA, several.ok = TRUE)
        physchemPro <- match.arg(physchemPro, several.ok = TRUE)

        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)

        cl <- parallel::makeCluster(parallel.cores)

        parallel::clusterExport(cl, varlist = c("Internal.checkRNA", "Internal.FourierSeries", "Internal.convertSeq_customize", "Internal.convertSeq", "aaindex"), envir = environment())

        if (seqType == "RNA") message("\n", "+ Calculating the physicochemical information of RNA sequences...  ", Sys.time()) else message("\n", "+ Calculating the physicochemical information of protein sequences...  ", Sys.time())

        output <- parallel::parLapply(cl, seqs, Internal.computePhysChem, seqType = seqType,
                                      Fourier.len = Fourier.len, as.list = as.list, aaindex = aaindex,
                                      physchemRNA.mode = physchemRNA, physchemPro.mode = physchemPro)
        parallel::stopCluster(cl)

        if (!as.list) output <- as.data.frame(t(data.frame(output, check.names = FALSE)))

        message("\n", "+ Completed.  ", Sys.time(), "\n")

        output
}

#' Extraction of the Physicochemical Features of RNA and Protein Sequences
#' @description Basically a wrapper for \code{\link{computePhysChem}} function.
#' This function can extract physicochemical features of RNA and protein sequences at the same time
#' and format the results as the dataset that can be used to build classifier.
#' @param seqRNA RNA sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of RNA sequences.
#' RNA sequences will be converted into lower case letters.
#' Each sequence should be a vector of single characters.
#' @param seqPro protein sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of protein sequences.
#' Protein sequences will be converted into upper case letters.
#' Each sequence should be a vector of single characters.
#' @param label optional. A string that indicates the class of the samples such as
#' "Interact", "Non.Interact". Default: \code{NULL}
#' @param parallel.cores an integer that indicates the number of cores for parallel computation.
#' Default: \code{2}. Set \code{parallel.cores = -1} to run with all the cores.
#' @param ... arguments (\code{Fourier.len}, \code{physchemRNA}, and \code{physchemPro})
#' passed to \code{\link{computePhysChem}}. See examples below.
#' @return This function returns a data frame.
#'
#' @details see \code{\link{computePhysChem}}.
#' @section References:
#' [1] Han S, Liang Y, Li Y, \emph{et al}.
#' ncProR: an integrated R package for effective ncRNA-protein interaction prediction.
#' (\emph{Submitted})
#'
#' [2] Morozova N, Allers J, Myers J, \emph{et al}.
#' Protein-RNA interactions: exploring binding patterns with a three-dimensional superposition analysis of high resolution structures.
#' Bioinformatics 2006; 22:2746-52
#'
#' [3] Grantham R.
#' Amino acid difference formula to help explain protein evolution.
#' Science 1974; 185:862-4
#'
#' [4] Zimmerman JM, Eliezer N, Simha R.
#' The characterization of amino acid sequences in proteins by statistical methods.
#' J. Theor. Biol. 1968; 21:170-201
#'
#' [5] Bull HB, Breese K.
#' Surface tension of amino acid solutions: a hydrophobicity scale of the amino acid residues.
#' Arch. Biochem. Biophys. 1974; 161:665-670
#'
#' [6] Kyte J, Doolittle RF.
#' A simple method for displaying the hydropathic character of a protein.
#' J. Mol. Biol. 1982; 157:105-132
#'
#' [7] Eisenberg D, Schwarz E, Komaromy M, \emph{et al}.
#' Analysis of membrane and surface protein sequences with the hydrophobic moment plot.
#' J. Mol. Biol. 1984; 179:125-42
#'
#' [8] Hopp TP, Woods KR.
#' Prediction of protein antigenic determinants from amino acid sequences.
#' Proc. Natl. Acad. Sci. U. S. A. 1981; 78:3824-8
#'
#' [9] Kawashima S, Kanehisa M. AAindex: amino acid index database. Nucleic Acids Res. 2000; 28:374
#'
#' [10] Bellucci M, Agostini F, Masin M, \emph{et al}.
#' Predicting protein associations with long noncoding RNAs.
#' Nat. Methods 2011; 8:444-445
#'
#' [11] Lu Q, Ren S, Lu M, \emph{et al}.
#' Computational prediction of associations between long non-coding RNAs and proteins.
#' BMC Genomics 2013; 14:651
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom seqinr a
#' @importFrom utils data
#'
#' @seealso \code{\link{computePhysChem}}
#' @examples
#' data(demoPositiveSeq)
#' seqsRNA <- demoPositiveSeq$RNA.positive
#' seqsPro <- demoPositiveSeq$Pro.positive
#'
#' dataset1 <- featurePhysChem(seqRNA = seqsRNA, seqPro = seqsPro,
#'                             label = "Interact", Fourier.len = 10,
#'                             physchemRNA = c("hydrogenBonding", "vanderWaal"),
#'                             physchemPro = c("polarity.Grantham", "polarity.Zimmerman",
#'                                             "hphob.BullBreese", "hphob.KyteDoolittle",
#'                                             "hphob.Eisenberg", "hphob.HoppWoods"))
#'
#' dataset2 <- featurePhysChem(seqRNA = seqsRNA, seqPro = seqsPro, Fourier.len = 12)
#'
#' @export

featurePhysChem <- function(seqRNA, seqPro, label = NULL, parallel.cores = 2, ...) {

        if (length(seqRNA) != length(seqPro)) stop("The number of RNA sequences should match the number of protein sequences!")

        PhysChem.RNA <- computePhysChem(seqs = seqRNA, seqType = "RNA", as.list = FALSE, parallel.cores = parallel.cores, ...)
        PhysChem.Pro <- computePhysChem(seqs = seqPro, seqType = "Pro", as.list = FALSE, parallel.cores = parallel.cores, ...)

        sequenceName <- paste(row.names(PhysChem.RNA), row.names(PhysChem.Pro), sep = ".")
        features <- cbind(PhysChem.RNA, PhysChem.Pro, row.names = sequenceName)

        if (!is.null(label)) features <- data.frame(label = label, features)
        features
}
