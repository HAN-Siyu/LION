#' Counting the Number of Motifs in RNA or Protein Sequences
#' @description Counts the number of motifs occurring in RNA/protein sequences. Motifs employed by tool "rpiCOOL"
#' can be selected. New motifs can also be defined.
#'
#' @param seqs sequences loaded by function \code{\link[seqinr]{read.fasta}} from \code{\link[seqinr]{seqinr-package}}. Or a list of RNA/protein sequences.
#' RNA sequences will be converted into lower case letters, but
#' protein sequences will be converted into upper case letters.
#' Each sequence should be a vector of single characters.
#' @param seqType a string that specifies the nature of the sequence: \code{"RNA"} or \code{"Pro"} (protein).
#' If the input is DNA sequence and \code{seqType = "RNA"}, the DNA sequence will be converted to RNA sequence automatically.
#' Default: \code{"RNA"}.
#' @param motifRNA strings specifying the motifs that are counted in RNA sequences. Ignored if \code{seqType = "Pro"}.
#' Options: \code{"rpiCOOL"}, \code{"selected5"},
#' \code{"Fox1"}, \code{"Nova"}, \code{"Slm2"}, \code{"Fusip1"}, \code{"PTB"}, \code{"ARE"}, \code{"hnRNPA1"},
#' \code{"PUM"}, \code{"U1A"}, \code{"HuD"}, \code{"QKI"}, \code{"U2B"}, \code{"SF1"}, \code{"HuR"}, \code{"YB1"},
#' \code{"AU"}, and \code{"UG"}. Multiple elements can be selected at the same time.
#' If \code{"rpiCOOL"}, all default motifs will be counted.
#' \code{"selected5"} indicates the total number of the occurrences of: PUM, Fox-1, U1A, Nova, and ARE which
#' are regarded as the five most over-presented binding motifs. See details below.
#' @param motifPro strings specifying the motifs that are counted in protein sequences. Ignored if \code{seqType = "RNA"}.
#' Options: \code{"rpiCOOL"}, \code{"E"}, \code{"H"}, \code{"K"}, \code{"R"}, \code{"H_R"},
#' \code{"EE"}, \code{"KK"}, \code{"HR_RH"}, \code{"RS_SR"},
#' \code{"RGG"}, and \code{"YGG"}. Multiple elements can be selected at the same time.
#' \code{"H_R"} indicates the total number of the occurrences of: H and R.
#' \code{"HR_RH"} indicates the total number of the occurrences of: HR and RH.
#' \code{"RS_SR"} indicates the total number of the occurrences of: RS and SR.
#' If \code{"rpiCOOL"}, the default motifs of rpiCOOL (\code{"E"}, \code{"K"}, \code{"H_R"},
#' \code{"EE"}, \code{"KK"}, \code{"RS_SR"}, \code{"RGG"}, and \code{"YGG"}) will be counted.
#' See details below.
#' @param newMotif list defining new motifs not listed above. New motifs are counted in RNA or protein sequences.
#' For example, \code{newMotif = list(hnRNPA1 = c("UAGGGU", "UAGGGA"), SF1 = "UACUAAC")}.
#' This parameter can be used together with parameter \code{motifRNA} or \code{motifPro}. Default: \code{NULL}.
#' @param newMotifOnly logical. If \code{TRUE}, only the new motifs defined in \code{newMotif} will be counted.
#' Default: \code{FALSE}.
#' @param parallel.cores an integer specifying the number of cores for parallel computation. Default: \code{2}.
#' Set \code{parallel.cores = -1} to run with all the cores. \code{parallel.cores} should be == -1 or >= 1.
#' @param cl parallel cores to be passed to this function.
#'
#' @return This function returns a data frame. Row names are sequences names, and column names are motif names.
#'
#' @details This function can count the motifs in RNA or protein sequences.
#'
#' The default motifs are selected or derived from tool "rpiCOOL" (Ref: [2]).
#'
#' \itemize{
#' \item Motifs of RNA
#'
#' \enumerate{
#' \item Fox1: UGCAUGU;
#' \item Nova: UCAUUUCAC, UCAUUUCAU, CCAUUUCAC, CCAUUUCAU;
#' \item Slm2: UAAAC, UAAAA, UAAUC, UAAUA;
#' \item Fusip1: AAAGA, AAAGG, AGAGA, AGAGG, CAAGA, CAAGG, CGAGA, CGAGG;
#' \item PTB: UUUUU, UUUCU, UCUUU, UCUCU;
#' \item ARE: UAUUUAUU;
#' \item hnRNPA1: UAGGGU, UAGGGA;
#' \item PUM: UGUAAAUA, UGUAGAUA, UGUAUAUA, UGUACAUA;
#' \item U1A: AUUGCAC;
#' \item HuD: UUAUUU;
#' \item QKI: AUUAAU, AUUAAC, ACUAAU, ACUAAC;
#' \item U2B: AUUGCAG;
#' \item SF1: UACUAAC;
#' \item HuR: UUUAUUU, UUUGUUU, UUUCUUU, UUUUUUU;
#' \item YB1: CCUGCG, UCUGCG;
#' \item AU: AU;
#' \item UG: UG.
#'
#' If \code{"rpiCOOL"}, all default motifs will be counted, and there is no need to input other default motifs.
#' \code{"selected5"} indicates the total number of the occurrences of: PUM, Fox-1, U1A, Nova, and ARE which
#' are regarded as the five most over-represented binding motifs.}
#'
#' \item Motifs of protein
#'
#' \enumerate{
#' \item E: E;
#' \item H: H;
#' \item K: K;
#' \item R: R;
#' \item EE: EE;
#' \item KK: KK;
#' \item HR (\code{"H_R"}): H, R;
#' \item HR (\code{"HR_RH"}): HR, RH;
#' \item RS (\code{"RS_SR"}): RS, SR;
#' \item RGG: RGG;
#' \item YGG: YGG.
#'
#' If \code{"rpiCOOL"}, default motifs of rpiCOOL (\code{"E"}, \code{"K"}, \code{"H_R"},
#' \code{"EE"}, \code{"KK"}, \code{"RS_SR"}, \code{"RGG"}, and \code{"YGG"}) will be counted.
#' } }
#'
#' There are some minor differences between this function and the extraction scheme of rpiCOOL.
#' In this function, motifs will be scanned directly.
#' As to the extraction scheme of rpiCOOL, some motifs (\code{"UG"}, \code{"AU"}, and \code{"H_R"})
#' are scanned in a 10 nt/aa sliding-window.
#'
#' New motif patterns are also supported. Users can pass new patterns to argument "newMotif" as a list.
#' Format:
#'
#' \code{newMotif = list(*motif_name* = c("*motif_pattern_1*", "*motif_pattern_2*"))}.
#'
#' For example: \code{newMotif = list(HR_RH = c("HR", "RH"), RGG = "RGG")}.
#' *HR_RH* is the name of this motif which contains two patterns: "HR" and "RH".
#'
#' @section References:
#' [1] Han S, Liang Y, Ma Q, \emph{et al}.
#' ncProR: an integrated R package for effective ncRNA-protein interaction prediction.
#' (\emph{Submitted})
#'
#' [2] Akbaripour-Elahabad M, Zahiri J, Rafeh R, \emph{et al}.
#' rpiCOOL: A tool for In Silico RNA-protein interaction detection using random forest.
#' J. Theor. Biol. 2016; 402:1-8
#'
#' [3] Pancaldi V, Bahler J.
#' In silico characterization and prediction of global protein-mRNA interactions in yeast.
#' Nucleic Acids Res. 2011; 39:5826-36
#'
#' [4] Castello A, Fischer B, Eichelbaum K, \emph{et al}.
#' Insights into RNA Biology from an Atlas of Mammalian mRNA-Binding Proteins.
#' Cell 2012; 149:1393-1406
#'
#' [5] Ray D, Kazan H, Cook KB, \emph{et al}.
#' A compendium of RNA-binding motifs for decoding gene regulation.
#' Nature 2013; 499:172-177
#'
#' [6] Jiang P, Singh M, Coller HA.
#' Computational assessment of the cooperativity between RNA binding proteins and MicroRNAs in Transcript Decay.
#' PLoS Comput. Biol. 2013; 9:e1003075
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom seqinr getSequence
#' @seealso \code{\link{featureMotifs}}
#' @examples
#' data(demoPositiveSeq)
#' seqsRNA <- demoPositiveSeq$RNA.positive
#' seqsPro <- demoPositiveSeq$Pro.positive
#'
#' motifRNA1 <- computeMotifs(seqsRNA, seqType = "RNA", motifRNA = "rpiCOOL",
#'                            parallel.cores = 2)
#'
#' motifRNA2 <- computeMotifs(seqsRNA, seqType = "RNA",
#'                            motifRNA = c("Fox1", "HuR", "ARE"), parallel.cores = 2)
#'
#' motifPro1 <- computeMotifs(seqsPro, seqType = "Pro",
#'                            motifPro = c("rpiCOOL", "HR_RH"), parallel.cores = 2)
#'
#' # Customized motifs are also supported and can be extracted with default motifs.
#' # Pass new motif patterns to "newMotif" argument as a list:
#'
#' motifPro2 <- computeMotifs(seqsPro, seqType = "Pro", motifPro = c("E", "K", "KK"),
#'                            newMotif = list(HR_RH = c("HR", "RH"), RGG = "RGG"),
#'                            parallel.cores = 2)
#'
#' # set "newMotifOnly = TRUE", if compute customized motifs only:
#'
#' motifPro3 <- computeMotifs(seqsPro, seqType = "Pro",
#'                            newMotif = list(HR_RH = c("HR", "RH"), RGG = "RGG"),
#'                            newMotifOnly = TRUE, parallel.cores = 2)
#' @export

computeMotifs <- function(seqs, seqType = c("RNA", "Pro"),
                          motifRNA = c("rpiCOOL", "Fox1", "Nova", "Slm2", "Fusip1", "PTB", "ARE", "hnRNPA1", "PUM",
                                       "U1A", "HuD", "QKI", "U2B", "SF1", "HuR", "YB1", "AU", "UG", "selected5"),
                          motifPro = c("rpiCOOL", "E", "H", "K", "R", "H_R", "EE", "KK", "HR_RH", "RS_SR",
                                       "RGG", "YGG"),
                          newMotif = NULL, newMotifOnly = FALSE, parallel.cores = 2, cl = NULL) {

        seqType <- match.arg(seqType)

        close_cl <- FALSE
        if (is.null(cl)) {
                message("+ Initializing...  ", Sys.time())
                message("- Creating cores...  ")
                parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
                cl <- parallel::makeCluster(parallel.cores)
                close_cl <- TRUE
        }

        if (seqType == "RNA") seqs <- parallel::parLapply(cl, seqs, Internal.checkRNA)

        seqs <- sapply(seqs, seqinr::getSequence, as.string = TRUE)

        motifPatterns <- NULL

        if (!is.null(newMotif)) {
                display_seqType <- ifelse(seqType == "RNA", "RNA sequences", "protein sequences")
                message("\n+ Formatting the new defined motifs of ", display_seqType, "...  ", Sys.time(),
                        "\n", "- ", paste(names(newMotif), collapse = ", "))
                motifPatterns <- lapply(newMotif, function(x) {
                        motifExpr <- paste(x, collapse = ")|(?=")
                        motifExpr <- paste0("(?=", motifExpr, ")")
                        motifExpr
                })
        }

        if (seqType == "RNA" & !newMotifOnly) {
                candidateMotif <- unique(match.arg(motifRNA, several.ok = TRUE))
                if ("rpiCOOL" %in% candidateMotif) {
                        candidateMotif <- c("Fox1", "Nova", "Slm2", "Fusip1", "PTB", "ARE", "hnRNPA1", "PUM", "U1A",
                                            "HuD", "QKI", "U2B", "SF1", "HuR", "YB1", "AU", "UG", "selected5")
                }
                message("\n", "+ Processing the default motifs of RNA sequences...  ", Sys.time(), "\n", "- ", paste(candidateMotif, collapse = ", "))

                Fox1 <- "(?=UGCAUGU)"
                Nova <- c("(?=UCAUUUCAC)|(?=UCAUUUCAU)|(?=CCAUUUCAC)|(?=CCAUUUCAU)")
                Slm2 <- c("(?=UAAAC)|(?=UAAAA)|(?=UAAUC)|(?=UAAUA)")
                Fusip1 <- c("(?=AAAGA)|(?=AAAGG)|(?=AGAGA)|(?=AGAGG)|(?=CAAGA)|(?=CAAGG)|(?=CGAGA)|(?=CGAGG)")
                PTB <- c("(?=UUUUU)|(?=UUUCU)|(?=UCUUU)|(?=UCUCU)")
                ARE <- "(?=UAUUUAUU)"
                hnRNPA1 <- c("(?=UAGGGU)|(?=UAGGGA)")
                PUM <- c("(?=UGUAAAUA)|(?=UGUAGAUA)|(?=UGUAUAUA)|(?=UGUACAUA)")
                U1A <- "(?=AUUGCAC)"
                HuD <- "(?=UUAUUU)"
                QKI <- c("(?=AUUAAU)|(?=AUUAAC)|(?=ACUAAU)|(?=ACUAAC)")
                U2B <- "(?=AUUGCAG)"
                SF1 <- "(?=UACUAAC)"
                HuR <- c("(?=UUUAUUU)|(?=UUUGUUU)|(?=UUUCUUU)|(?=UUUUUUU)")
                YB1 <- c("(?=CCUGCG)|(?=UCUGCG)")
                AU <- "(?=AU)"
                UG <- "(?=UG)"
                selected5 <- paste(Fox1, Nova, ARE, PUM, U1A, sep = "|")

                Patterns <- mget(candidateMotif)
                motifPatterns <- c(motifPatterns, Patterns)
        }

        if (seqType == "Pro" & !newMotifOnly) {
                candidateMotif <- match.arg(motifPro, several.ok = TRUE)
                if ("rpiCOOL" %in% candidateMotif) {
                        candidateMotif <- c(candidateMotif, "E", "K", "EE", "KK", "RS_SR", "RGG", "YGG", "H_R")
                        candidateMotif <- unique(candidateMotif)
                        candidateMotif <- candidateMotif[!candidateMotif %in% "rpiCOOL"]
                }

                message("\n", "+ Processing the default motifs of protein sequences...  ", Sys.time(), "\n", "- ", paste(candidateMotif, collapse = ", "))

                E <- "E"
                H <- "H"
                K <- "K"
                R <- "R"
                H_R <- c("(?=H)|(?=R)")
                EE <- "(?=EE)"
                KK <- "(?=KK)"
                HR_RH <- c("(?=HR)|(?=RH)")
                RS_SR <- c("(?=RS)|(?=SR)")
                RGG <- "(?=RGG)"
                YGG <- "(?=YGG)"

                Patterns <- mget(candidateMotif)
                motifPatterns <- c(motifPatterns, Patterns)
        }

        motifCounts <- parallel::parLapply(cl, seqs, Internal.computeMotifs, motifs = motifPatterns)
        if (close_cl) parallel::stopCluster(cl)

        motifCounts <- as.data.frame(t(data.frame(motifCounts, check.names = F)))
        formatNames <- paste("motif_", names(motifCounts), sep = "")
        names(motifCounts) <- formatNames

        # message("+ Completed.  ", Sys.time(), "\n")
        motifCounts
}

#' Extraction of the Motif Features of RNA and Protein Sequences
#' @description Basically a wrapper for \code{\link{computeMotifs}} function.
#' This function can count the motifs of RNA and protein sequences at the same time
#' and format the results as the dataset that can be used to build classifier.
#' @param seqRNA RNA sequences loaded by function \code{\link[seqinr]{read.fasta}} from \code{\link[seqinr]{seqinr-package}}. Or a list of RNA sequences.
#' RNA sequences will be converted into lower case letters.
#' @param seqPro protein sequences loaded by function \code{\link[seqinr]{read.fasta}} from \code{\link[seqinr]{seqinr-package}}. Or a list of protein sequences.
#' Protein sequences will be converted into upper case letters.
#' @param label optional. A string that indicates the class of the samples such as
#' "Interact", "Non.Interact". Default: \code{NULL}
#' @param featureMode a string that can be \code{"concatenate"} or \code{"combine"}.
#' If \code{"concatenate"}, the motif features of RNA and proteins will be simply concatenated.
#' If \code{"combine"}, the returned dataset will be formed by combining the motif features of RNA and proteins.
#' See details below. Default: \code{"concatenate"}.
#' @param newMotif.RNA a list specifying the motifs that are counted in RNA sequences. Default: \code{NULL}.
#' For example, \code{newMotif = list(hnRNPA1 = c("UAGGGU", "UAGGGA"), SF1 = "UACUAAC")}.
#' Can be used with parameter \code{motifRNA} (see parameter \code{...}) to count motifs in RNA sequences.
#' @param newMotif.Pro a list specifying the motifs that are counted in protein sequences. Default: \code{NULL}.
#' For example, \code{newMotif = list(YGG = "YGG", E = "E")}.
#' Can be used with parameter \code{motifPro} (see parameter \code{...}) to count motifs in protein sequences.
#' @param newMotifOnly.RNA logical. If \code{TRUE}, only the new motifs defined in \code{newMotif.RNA} will be counted.
#' Default: \code{FALSE}.
#' @param newMotifOnly.Pro logical. If \code{TRUE}, only the new motifs defined in \code{newMotif.Pro} will be counted.
#' Default: \code{FALSE}.
#' @param parallel.cores an integer that indicates the number of cores for parallel computation. Default: \code{2}.
#' Set \code{parallel.cores = -1} to run with all the cores. \code{parallel.cores} should be == -1 or >= 1.
#' @param cl parallel cores to be passed to this function.
#' @param ... argument \code{motifRNA} and \code{motifPro} to be passed to \code{\link{computeMotifs}}.
#' Used to compute the default motifs. See examples below.
#'
#' @return This function returns a data frame. Row names are the sequences names, and column names are the motif names.
#' The names of RNA and protein sequences are separated with ".",
#' i.e. row names format: "*RNASequenceName*.*proteinSequenceName*" (e.g. "YDL227C.YOR198C").
#' If \code{featureMode = "combine"}, the motif names of RNA and protein sequences are also separated with ".",
#' i.e. column names format: "motif_*RNAMotifName*.motif_*proteinMotifName*" (e.g. "motif_PUM.motif_EE").
#'
#' @details
#' If \code{featureMode = "concatenate"}, \emph{m} RNA motif features will be simply
#' concatenated with \emph{n} protein motif features, and the final result has \emph{m} + \emph{n} features.
#' If \code{featureMode = "combine"}, \emph{m} RNA motif features will be
#' combined with \emph{n} protein motif features, resulting in \emph{m} * \emph{n} possible combinations.
#'
#' \code{...} can be used to pass the default motif patterns of RNA and protein sequences.
#' See arguments \code{motifRNA} and \code{motifPro} in \code{\link{computeMotifs}}.
#'
#' @section References:
#' [1] Han S, Liang Y, Ma Q, \emph{et al}.
#' ncProR: an integrated R package for effective ncRNA-protein interaction prediction.
#' (\emph{Submitted})
#'
#' [2] Akbaripour-Elahabad M, Zahiri J, Rafeh R, \emph{et al}.
#' rpiCOOL: A tool for In Silico RNA-protein interaction detection using random forest.
#' J. Theor. Biol. 2016; 402:1-8
#'
#' [3] Pancaldi V, Bahler J.
#' In silico characterization and prediction of global protein-mRNA interactions in yeast.
#' Nucleic Acids Res. 2011; 39:5826-36
#'
#' [4] Castello A, Fischer B, Eichelbaum K, \emph{et al}.
#' Insights into RNA Biology from an Atlas of Mammalian mRNA-Binding Proteins.
#' Cell 2012; 149:1393-1406
#'
#' [5] Ray D, Kazan H, Cook KB, \emph{et al}.
#' A compendium of RNA-binding motifs for decoding gene regulation.
#' Nature 2013; 499:172-177
#'
#' [6] Jiang P, Singh M, Coller HA.
#' Computational assessment of the cooperativity between RNA binding proteins and MicroRNAs in Transcript Decay.
#' PLoS Comput. Biol. 2013; 9:e1003075
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom seqinr getSequence
#' @seealso \code{\link{computeMotifs}}
#' @examples
#' data(demoPositiveSeq)
#' seqsRNA <- demoPositiveSeq$RNA.positive
#' seqsPro <- demoPositiveSeq$Pro.positive
#'
#' dataset1 <- featureMotifs(seqRNA = seqsRNA, seqPro = seqsPro, featureMode = "conc",
#'                           newMotif.RNA = list(motif1 = c("cc", "cu")),
#'                           newMotif.Pro = list(motif2 = "KK"),
#'                           motifRNA = c("Fusip1", "AU", "UG"),
#'                           motifPro = c("E", "K", "HR_RH"))
#'
#' dataset2 <- featureMotifs(seqRNA = seqsRNA, seqPro = seqsPro, featureMode = "comb",
#'                           newMotif.RNA = list(motif1 = c("cc", "cu")),
#'                           newMotif.Pro = list(motif2 = c("R", "H")),
#'                           newMotifOnly.RNA = TRUE, newMotifOnly.Pro = FALSE)
#'
#' @export

featureMotifs <- function(seqRNA, seqPro, label = NULL, featureMode = c("concatenate", "combine"),
                          newMotif.RNA = NULL, newMotif.Pro = NULL, newMotifOnly.RNA = FALSE,
                          newMotifOnly.Pro = FALSE, parallel.cores = 2, cl = NULL, ...) {

        if (length(seqRNA) != length(seqPro)) stop("The number of RNA sequences should match the number of protein sequences!")

        featureMode <- match.arg(featureMode)

        close_cl <- FALSE
        if (is.null(cl)) {
                message("+ Initializing...  ", Sys.time())
                message("- Creating cores...  ")
                parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
                cl <- parallel::makeCluster(parallel.cores)
                close_cl <- TRUE
        }

        featureRNA <- computeMotifs(seqs = seqRNA, seqType = "RNA", newMotif = newMotif.RNA,
                                    newMotifOnly = newMotifOnly.RNA, cl = cl, ...)
        featurePro <- computeMotifs(seqs = seqPro, seqType = "Pro", newMotif = newMotif.Pro,
                                    newMotifOnly = newMotifOnly.Pro, cl = cl, ...)
        if (close_cl) parallel::stopCluster(cl)

        sequenceName <- paste(row.names(featureRNA), row.names(featurePro), sep = ".")

        if (featureMode == "combine") {
                featureName <- sapply(names(featureRNA), function(nameRNA) {
                        names <- paste(nameRNA, names(featurePro), sep = ".")
                })
                featureName <- as.character(featureName)
                featureValue <- mapply(Internal.combineMotifs, oneRNA = as.data.frame(t(featureRNA)),
                                       onePro = as.data.frame(t(featurePro)))
                features <- as.data.frame(t(featureValue), row.names = sequenceName)
                names(features) <- featureName
        } else {
                features <- cbind(featureRNA, featurePro, row.names = sequenceName)
        }

        if (!is.null(label)) features <- data.frame(label = label, features)
        features
}
