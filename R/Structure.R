#' Run RNAsubopt in R
#' @description This function can run and capture the result of RNAsubopt (ViennaRNA Package) in R
#' by invoking the OS command.
#' The results can be formatted and used as the features of tool "catRAPID", "lncPro", etc.
#' ViennaRNA Package is required.
#'
#' @param seqs RNA sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of RNA/protein sequences.
#' RNA sequences will be converted into lower case letters.
#' Each sequence should be a vector of single characters.
#' @param structureRNA.num integer. Number of random samples of suboptimal structures. Default: \code{6}.
#' @param path.RNAsubopt string specifying the location of RNAsubopt program.
#' @param returnSum logical. If \code{FALSE}, the raw output of RNAsubopt (Dot-Bracket Notation of RNA structure) will be returned.
#' If \code{TRUE}, the results of RNAsubopt will be converted into numeric sequence:
#' In any RNA structure sequence (Dot-Bracket Notation), "." (Dot) will be replaced with 0;
#' "(" and ")" (Bracket) will be replaced with 1. And \emph{structureRNA.num} binary sequences will be added to obtain a
#' new numeric sequence.
#' @param verbose logical. Should the relevant information be printed during the calculation? (Only available on UNIX/Linux.)
#' @param parallel.cores an integer that indicates the number of cores for parallel computation. Default: \code{2}.
#' Set \code{parallel.cores = -1} to run with all the cores.
#'
#' @return This function returns a data frame that contains the raw outputs when \code{returnSum = FALSE}.
#' Or a list of numeric results of RNAsubopt when \code{returnSum = TRUE}.
#'
#' @details
#' This function depends on the program "RNAsubopt" of software "ViennaRNA".
#' (\url{http://www.tbi.univie.ac.at/RNA/index.html})
#'
#' Parameter \code{path.RNAsubopt} can be simply defined as \code{"RNAsubopt"} as
#' default when the OS is UNIX/Linux. However, for some OS, such as Windows, users
#' need to specify the \code{path.RNAsubopt} if the path of "RNAsubopt" haven't been
#' added in environment variables (e.g. \code{path.RNAsubopt = '"C:/Program Files/ViennaRNA/RNAsubopt.exe"'}).
#'
#' This function can print the related information when the OS is UNIX/Linux,
#' such as:
#'
#' \code{"25 of 100, length: 695 nt"},
#'
#' which means around 100 sequences are assigned to this node, and the program is
#' computing the 25th sequence. The length of this sequence is 695 nt.
#'
#' @section References:
#' Lorenz R, Bernhart SH, Honer zu Siederdissen C, \emph{et al}.
#' ViennaRNA Package 2.0.
#' Algorithms Mol. Biol. 2011; 6:26
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom parallel parLapply
#' @importFrom seqinr s2c
#'
#' @seealso \code{\link{runPredator}}
#' @examples
#' \dontrun{
#' data(demoPositiveSeq)
#' seqsRNA <- demoPositiveSeq$RNA.positive
#'
#' RNAsubopt.res <- runRNAsubopt(seqs = seqsRNA, path.RNAsubopt = "RNAsubopt",
#'                               returnSum = FALSE, parallel.cores = -1)
#' }
#' @export

runRNAsubopt <- function(seqs, structureRNA.num = 6, path.RNAsubopt = "RNAsubopt",
                         returnSum = FALSE, verbose = FALSE, parallel.cores = 2) {

        message("+ Initializing...  ", Sys.time())
        if (parallel.cores == 2) message("- Users can try to set parallel.cores = -1 to use all cores!")
        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)

        message("\n", "+ Calculating structural information of RNA sequences...  ", Sys.time())

        if (verbose) {
                cl <- parallel::makeCluster(parallel.cores, outfile = "")
        } else{
                cl <- parallel::makeCluster(parallel.cores)
        }

        seqValidate <- seqs[which(lengths(seqs) <= 4095)]
        if (length(seqValidate) < length(seqs)) {
                message("- Due to the limitation of RNAsubopt,")
                message("- sequences with length more than 4095 nt will be omitted.")
                message("- ", length(seqs) - length(seqValidate), " of ", length(seqs), " sequences have been removed.", "\n")
        }
        seqs <- seqValidate

        message("- Sequences Number: ", length(seqs))
        index <- 1
        seqs <- parallel::parLapply(cl, seqs, Internal.checkRNA)
        seq.string <- parallel::parSapply(cl, seqs, seqinr::getSequence, TRUE)
        info <- paste(round(length(seqs) / parallel.cores), ",", sep = "")

        parallel::clusterExport(cl, varlist = c("index"), envir = environment())

        RNAsubopt.seq <- parallel::parLapply(cl, seq.string, Internal.runRNAsubopt, info = info,
                                             structureRNA.num = structureRNA.num, verbose = verbose,
                                             path.RNAsubopt = path.RNAsubopt, outputNum = returnSum)
        parallel::stopCluster(cl)

        if (!returnSum) RNAsubopt.seq <- data.frame(RNAsubopt.seq, check.names = FALSE)

        message("\n", "+ Completed.  ", Sys.time(), "\n")

        RNAsubopt.seq
}

#' Run Predator in R
#' @description This function can run and capture the result of Predator in R
#' by invoking the OS command.
#' The results can be formatted and used as the features of tool "lncPro".
#' Predator is required (NOTE: Predator does not support 64-bit MS Windows OS. Linux OS is recommended).
#'
#' @param seqs RNA sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of RNA/protein sequences.
#' RNA sequences will be converted into lower case letters.
#' Each sequence should be a vector of single characters. (Non-AA letters will be ignored.)
#' @param path.Predator string specifying the location of Predator program.
#' @param path.stride string specifying the location of file "stride.dat" required by program Predator.
#' @param workDir string specifying the directory for temporary files.
#' The temp files will be deleted automatically when
#' the calculation is completed. The directory does not exist will be created automatically.
#' @param verbose logical. Should the relevant information be printed during the calculation? (Only available on UNIX/Linux.)
#' @param parallel.cores an integer that indicates the number of cores for parallel computation. Default: \code{2}.
#' Set \code{parallel.cores = -1} to run with all the cores.
#'
#' @return This function returns a data frame of the raw outputs of Predator.
#'
#' @details
#' This function depends on the program Predator.
#' Program Predator is only available on UNIX/Linux and 32-bit Windows OS.
#'
#' This function can print the related information when the OS is UNIX/Linux,
#' such as:
#'
#' \code{"25 of 100, length: 50 aa"},
#'
#' which means around 100 sequences are assigned to this node, and the program is
#' computing the 25th sequence. The length of this sequence is 50 aa.
#'
#' @section References:
#' Frishman D, Argos P.
#' Incorporation of non-local interactions in protein secondary structure prediction from the amino acid sequence.
#' Protein Eng. 1996; 9:133-42
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom parallel parLapply
#' @importFrom utils data
#' @importFrom seqinr c2s
#' @importFrom seqinr s2c
#' @importFrom seqinr write.fasta
#'
#' @seealso \code{\link{runRNAsubopt}}
#' @examples
#' \dontrun{
#' data(demoPositiveSeq)
#' seqsPro <- demoPositiveSeq$Pro.positive
#'
#' path.predator = "/mnt/external_drive_1/hansy/R_tmp/LncProInter/Predator/predator"
#' path.stride = "/mnt/external_drive_1/hansy/R_tmp/LncProInter/Predator/stride.dat"
#'
#' Predator.res <- runPredator(seqs = Pro.positive, path.Predator = path.Predator,
#'                             path.stride = path.stride, workDir = "tmp",
#'                             as.string = TRUE, parallel.cores = 4)
#' }
#' @export

runPredator <- function(seqs, path.Predator, path.stride, workDir = getwd(),
                        verbose = FALSE, parallel.cores = 2) {

        if (!file.exists(path.stride)) stop("The path of stride.dat is not correct! Please check parameter path.stride.")

        message("+ Initializing...  ", Sys.time())
        if (parallel.cores == 2) message("- Users can try to set parallel.cores = -1 to use all cores!")
        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        if (verbose) {
                cl <- parallel::makeCluster(parallel.cores, outfile = "")
        } else{
                cl <- parallel::makeCluster(parallel.cores)
        }

        message("\n", "+ Calculating structural information of protein sequences...  ", Sys.time())

        seqValidate <- seqs[which(lengths(seqs) >= 30)]
        if (length(seqValidate) < length(seqs)) {
                message("- Due to the limitation of predator,")
                message("- sequences with length less than 30 amino acids will be omitted.")
                message("- ", length(seqs) - length(seqValidate), " of ", length(seqs), " sequences have been removed.", "\n")
        }
        seqs <- seqValidate

        message("- Sequences Number: ", length(seqs))
        index <- 1

        info <- paste0(round(length(seqs) / parallel.cores), ",")

        parallel::clusterExport(cl, varlist = c("index", "Internal.randomName"), envir = environment())

        if (!dir.exists(workDir)) {
                dir.create(workDir, recursive = TRUE)
                createDir <- TRUE
        } else {
                createDir <- FALSE
        }

        predator.seq <- parallel::parSapply(cl, seqs, Internal.runPredator, info = info,
                                            workDir = workDir, path.Predator = path.Predator,
                                            path.stride = path.stride, as.string = TRUE,
                                            outputRaw = TRUE, verbose = verbose)
        parallel::stopCluster(cl)

        predator.seq <- data.frame(predator.seq, check.names = FALSE)

        if (createDir & length(dir(workDir, all.files = TRUE, recursive = TRUE)) == 0) unlink(workDir, recursive = TRUE)

        message("\n", "+ Completed.  ", Sys.time(), "\n")

        predator.seq
}

#' Computation of the Secondary Structural Features of RNA or Protein Sequences
#' @description The function \code{computeStructure} computes the secondary structural features of RNA
#' or protein sequences.
#' @param seqs sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of RNA/protein sequences.
#' RNA sequences will be converted into lower case letters, but
#' protein sequences will be converted into upper case letters, and non-AA letters will be ignored.
#' Each sequence should be a vector of single characters.
#' @param seqType a string that specifies the nature of the sequence: \code{"RNA"} or \code{"Pro"} (protein).
#' If the input is DNA sequence and \code{seqType = "RNA"}, the DNA sequence will be converted to RNA sequence automatically.
#' Default: \code{"RNA"}.
#' @param structureRNA.num integer. The number of random samples of suboptimal structures. Default: \code{6}.
#' @param structurePro strings specifying the secondary structural information that are extracted from protein sequences.
#' Ignored if \code{seqType = "RNA"}.
#' Options: \code{"ChouFasman"}, \code{"DeleageRoux"}, and \code{"Levitt"}. See details below.(Ref: [2-4])
#' Multiple elements can be selected at the same time.
#' @param Fourier.len postive integer specifying the Fourier series length that will be used as features.
#' The \code{Fourier.len} should be >= the length of the input sequence. Default: \code{10}.
#' @param workDir.Pro string specifying the directory for temporary files.
#' The temp files will be deleted automatically when
#' the calculation is completed.
#' @param as.list logical. The result will be returned as a list or data frame.
#' @param path.RNAsubopt string specifying the location of RNAsubopt program. (Ref: [5])
#' @param path.Predator string specifying the location of Predator program. (Ref: [6])
#' @param path.stride string specifying the location of file "stride.dat" required by program Predator.
#' @param verbose logical. Should the relevant information be printed during the calculation? (Only available on UNIX/Linux.)
#' @param parallel.cores an integer that indicates the number of cores for parallel computation. Default: \code{2}.
#' Set \code{parallel.cores = -1} to run with all the cores.
#'
#' @return This function returns a data frame if \code{as.list = FALSE} or returns a list if \code{as.list = TRUE}.
#'
#' @details
#' The secondary structures of RNA and protein are computed by RNAsubopt and Predator, respectively.
#' And the protein secondary features are encoded using three amino acid scales:
#' \enumerate{
#' \item Chou & Fasman conformational parameter (Ref: [2])
#' \item Deleage & Roux conformational parameter (Ref: [3])
#' \item Levitt normalised frequency (Ref: [4])
#' }
#' The feature encoding strategy is based on lncPro (Ref: [7]).
#'
#' @section References:
#' [1] Han S, Liang Y, Li Y, \emph{et al}.
#' ncProR: an integrated R package for effective ncRNA-protein interaction prediction.
#' (\emph{Submitted})
#'
#' [2] Chou PY, Fasman GD.
#' Prediction of the secondary structure of proteins from their amino acid sequence.
#' Adv. Enzymol. Relat. Areas Mol. Biol. 1978; 47:45-148
#'
#' [3] Deleage G, Roux B.
#' An algorithm for protein secondary structure prediction based on class prediction.
#' Protein Eng. Des. Sel. 1987; 1:289-294
#'
#' [4] Levitt M.
#' Conformational preferences of amino acids in globular proteins.
#' Biochemistry 1978; 17:4277-85
#'
#' [5] Frishman D, Argos P.
#' Incorporation of non-local interactions in protein secondary structure prediction from the amino acid sequence.
#' Protein Eng. 1996; 9:133-42
#'
#' [6] Lorenz R, Bernhart SH, Honer zu Siederdissen C, \emph{et al}.
#' ViennaRNA Package 2.0.
#' Algorithms Mol. Biol. 2011; 6:26
#'
#' [7] Lu Q, Ren S, Lu M, \emph{et al}.
#' Computational prediction of associations between long non-coding RNAs and proteins.
#' BMC Genomics 2013; 14:651
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom parallel parLapply
#' @importFrom seqinr getSequence
#' @importFrom seqinr write.fasta
#' @importFrom seqinr s2c
#' @importFrom seqinr a
#' @importFrom utils data
#'
#' @seealso \code{\link{runRNAsubopt}}, \code{\link{runPredator}}, \code{\link{featureStructure}}
#' @examples
#' \dontrun{
#' data(demoNegativeSeq)
#' seqsRNA <- demoNegativeSeq$RNA.negative
#' seqsPro <- demoNegativeSeq$Pro.negative
#'
#' StructureRNA <- computeStructure(seqsRNA, seqType = "RNA", structureRNA.num = 8,
#'                                  Fourier.len = 12, as.list = FALSE,
#'                                  path.RNAsubopt = "RNAsubopt", parallel.cores = 2)
#'
#' StructurePro <- computeStructure(seqsPro, seqType = "Pro",
#'                                  structurePro = c("ChouFasman", "DeleageRoux",
#'                                                   "Levitt"),
#'                                  Fourier.len = 10, workDir.Pro = getwd(),
#'                                  as.list = TRUE,
#'                                  path.Predator = "Predator/predator",
#'                                  path.stride = "Predator/stride.dat",
#'                                  parallel.cores = 2)
#' }
#' @export

computeStructure <- function(seqs, seqType = c("RNA", "Pro"), structureRNA.num = 6,
                             structurePro = c("ChouFasman", "DeleageRoux", "Levitt"),
                             Fourier.len = 10, workDir.Pro = getwd(), as.list = TRUE,
                             path.RNAsubopt = "RNAsubopt", path.Predator = "Predator/predator",
                             path.stride = "Predator/stride.dat", verbose = FALSE, parallel.cores = 2) {

        # path.stride <- normalizePath(path.stride)
        if (!file.exists(path.stride)) stop("The path of stride.dat is not correct! Please check parameter path.stride.")
        if (min(lengths(seqs)) < Fourier.len) stop("The profile length of Fourier Series (Fourier.len) cannot larger than the minimum length of the input sequences!")

        message("+ Initializing...  ", Sys.time())
        seqType <- match.arg(seqType)

        if (parallel.cores == 2) message("- Users can try to set parallel.cores = -1 to use all cores!")
        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        if (verbose) {
                cl <- parallel::makeCluster(parallel.cores, outfile = "")
        } else {
                cl <- parallel::makeCluster(parallel.cores)
        }

        if (seqType == "RNA") {
                message("\n", "+ Calculating structural information of RNA sequences...  ", Sys.time())
                seqValidate <- seqs[which(lengths(seqs) <= 4095)]
                if (length(seqValidate) < length(seqs)) {
                        message("- Due to the limitation of RNAsubopt,")
                        message("- sequences with length more than 4095 nt will be omitted.")
                        message("- ", length(seqs) - length(seqValidate), " of ", length(seqs), " sequences have been removed.", "\n")
                }
                seqs <- seqValidate

                message("- Sequences Number: ", length(seqs))
                index <- 1
                seqs <- parallel::parLapply(cl, seqs, Internal.checkRNA)
                seq.string <- parallel::parSapply(cl, seqs, seqinr::getSequence, TRUE)
                info <- paste(round(length(seqs) / parallel.cores), ",", sep = "")

                parallel::clusterExport(cl, varlist = c("index"), envir = environment())
                RNAsubopt.seq <- parallel::parLapply(cl, seq.string, Internal.runRNAsubopt, info = info,
                                                     structureRNA.num = structureRNA.num, verbose = verbose,
                                                     path.RNAsubopt = path.RNAsubopt, outputNum = TRUE)
                out.seq <- parallel::parLapply(cl, RNAsubopt.seq, Internal.FourierSeries, profile.length = Fourier.len)

                parallel::stopCluster(cl)

                if (!as.list) {
                        tmp.df <- t(as.data.frame(out.seq))
                        row.names(tmp.df) <- names(out.seq)
                        tmp.df <- as.data.frame(tmp.df)
                        names(tmp.df) <- paste0("RNAstructure", c(1:ncol(tmp.df)))
                        out.seq <- tmp.df
                }

        } else {
                data(aaindex, package = "seqinr", envir = environment())
                aaindex <- get("aaindex")

                message("\n", "+ Calculating structural information of protein sequences...")
                structurePro <- match.arg(structurePro, several.ok = TRUE)
                seqValidate <- seqs[which(lengths(seqs) >= 30)]
                if (length(seqValidate) < length(seqs)) {
                        message("- Due to the limitation of predator,")
                        message("- sequences with length less than 30 amino acids will be omitted.")
                        message("- ", length(seqs) - length(seqValidate), " of ", length(seqs), " sequences have been removed.", "\n")
                }
                seqs <- seqValidate

                message("- Sequences Number: ", length(seqs))
                index <- 1

                info <- paste0(round(length(seqs) / parallel.cores), ",")

                parallel::clusterExport(cl, varlist = c("index", "Internal.randomName", "aaindex", "Internal.convertSeq",
                                                        "Internal.convertSeq_customize", "Internal.FourierSeries"), envir = environment())

                if (!dir.exists(workDir.Pro)) {
                        dir.create(workDir.Pro, recursive = TRUE)
                        createDir <- TRUE
                } else {
                        createDir <- FALSE
                }

                out.seq <- parallel::parLapply(cl, seqs, Internal.convertPro_Struct,
                                               info = info, path.Predator = path.Predator,
                                               path.stride = path.stride, workDir = workDir.Pro,
                                               structurePro = structurePro, as.list = as.list,
                                               Fourier.len = Fourier.len, verbose = verbose, aaindex)
                parallel::stopCluster(cl)

                if (createDir & length(dir(workDir.Pro, all.files = TRUE, recursive = TRUE)) == 0) unlink(workDir.Pro, recursive = TRUE)

                if (!as.list) out.seq <- as.data.frame(t(data.frame(out.seq, check.names = FALSE)))
        }
        message("\n", "+ Completed.  ", Sys.time(), "\n")

        out.seq
}

#' Extraction of the Secondary Structural Features of RNA and Protein Sequences
#' @description Basically a wrapper for \code{\link{computeStructure}} function.
#' This function can extract secondary structure features of RNA and protein sequences at the same time
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
#' "Interact", "NonInteract". Default: \code{NULL}
#' @param parallel.cores an integer that indicates the number of cores for parallel computation.
#' Default: \code{2}. Set \code{parallel.cores = -1} to run with all the cores.
#' @param ... arguments (\code{structureRNA.num}, \code{structurePro}, \code{Fourier.len},
#' \code{workDir.Pro}, \code{path.RNAsubopt}, \code{path.Predator}, \code{path.stride}, and \code{verbose})
#' passed to function \code{\link{computeStructure}}. See example below.
#' @return This function returns a data frame.
#'
#' @details see \code{\link{computeStructure}}.
#' @section References:
#' [1] Han S, Liang Y, Li Y, \emph{et al}.
#' ncProR: an integrated R package for effective ncRNA-protein interaction prediction.
#' (\emph{Submitted})
#'
#' [2] Chou PY, Fasman GD.
#' Prediction of the secondary structure of proteins from their amino acid sequence.
#' Adv. Enzymol. Relat. Areas Mol. Biol. 1978; 47:45-148
#'
#' [3] Deleage G, Roux B.
#' An algorithm for protein secondary structure prediction based on class prediction.
#' Protein Eng. Des. Sel. 1987; 1:289-294
#'
#' [4] Levitt M.
#' Conformational preferences of amino acids in globular proteins.
#' Biochemistry 1978; 17:4277-85
#'
#' [5] Frishman D, Argos P.
#' Incorporation of non-local interactions in protein secondary structure prediction from the amino acid sequence.
#' Protein Eng. 1996; 9:133-42
#'
#' [6] Lorenz R, Bernhart SH, Honer zu Siederdissen C, \emph{et al}.
#' ViennaRNA Package 2.0.
#' Algorithms Mol. Biol. 2011; 6:26
#'
#' [7] Lu Q, Ren S, Lu M, \emph{et al}.
#' Computational prediction of associations between long non-coding RNAs and proteins.
#' BMC Genomics 2013; 14:651
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom parallel parLapply
#' @importFrom seqinr getSequence
#' @importFrom seqinr write.fasta
#' @importFrom seqinr s2c
#' @importFrom seqinr a
#' @importFrom utils data
#'
#' @seealso \code{\link{runRNAsubopt}}, \code{\link{runPredator}}, \code{\link{computeStructure}}
#' @examples
#' \dontrun{
#' data(demoNegativeSeq)
#' seqsRNA <- demoNegativeSeq$RNA.negative
#' seqsPro <- demoNegativeSeq$Pro.negative
#'
#' # Use your own paths of the program Predator and file "stride.dat". For example:
#'
#' path.stride = "/mnt/external_drive_1/hansy/R_tmp/LncProInter/Predator/stride.dat"
#' path.Predator = "/mnt/external_drive_1/hansy/R_tmp/LncProInter/Predator/predator"
#'
#' dataset <- featureStructure(seqRNA = seqsRNA, seqPro = seqsPro, label = "Non.Interact",
#'                             parallel.cores = 2, structureRNA.num = 6,
#'                             structurePro = c("ChouFasman", "DeleageRoux", "Levitt"),
#'                             Fourier.len = 10, workDir.Pro = "tmpDir",
#'                             path.RNAsubopt = "RNAsubopt", path.Predator = path.Predator,
#'                             path.stride = path.stride)
#' }
#' @export

featureStructure <- function(seqRNA, seqPro, label = NULL, parallel.cores = 2, ...) {

        if (length(seqRNA) != length(seqPro)) stop("The number of RNA sequences should match the number of protein sequences!")

        RNA.velidate.idx <- which(lengths(seqRNA) <= 4095)
        seqValidate <- seqRNA[RNA.velidate.idx]
        if (length(seqValidate) < length(seqRNA)) {
                message("- Due to the limitation of RNAsubopt,")
                message("- sequences with length more than 4095 nt will be omitted.")
                message("- ", length(seqRNA) - length(seqValidate), " of ", length(seqRNA), " pairs have been removed.", "\n")
        }
        seqRNA <- seqValidate
        seqPro <- seqPro[RNA.velidate.idx]

        Pro.velidate.idx <- which(lengths(seqPro) >= 30)
        seqValidate <- seqPro[Pro.velidate.idx]
        if (length(seqValidate) < length(seqPro)) {
                message("- Due to the limitation of predator,")
                message("- sequences with length less than 30 amino acids will be omitted.")
                message("- ", length(seqPro) - length(seqValidate), " of ", length(seqPro), " pairs have been removed.", "\n")
        }
        seqPro <- seqValidate
        seqRNA <- seqRNA[Pro.velidate.idx]

        Structure.RNA <- computeStructure(seqs = seqRNA, seqType = "RNA", as.list = FALSE, parallel.cores = parallel.cores, ...)
        Structure.Pro <- computeStructure(seqs = seqPro, seqType = "Pro", as.list = FALSE, parallel.cores = parallel.cores, ...)

        sequenceName <- paste(row.names(Structure.RNA), row.names(Structure.Pro), sep = ".")
        features <- cbind(Structure.RNA, Structure.Pro, row.names = sequenceName)

        if (!is.null(label)) features <- data.frame(label = label, features)
        features
}
