#' Predict RNA-Protein Interaction Using lncPro Method
#' @description This function can predict lncRNA/RNA-protein interactions using lncPro method.
#' Programs "RNAsubopt" from software "ViennaRNA Package" and "Predator" is required. Please also note that
#' "Predator" is only available on UNIX/Linux and 32-bit Windows OS.
#'
#' @param seqRNA RNA sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of RNA/protein sequences.
#' RNA sequences will be converted into lower case letters.
#' @param seqPro protein sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of protein sequences.
#' Protein sequences will be converted into upper case letters.
#' Each sequence should be a vector of single characters.
#' @param mode set \code{"original"} to use original lncPro algorithm, and \code{"retrained"} to call retrained model.
#' The retrained model is constructed with the same features as the original version, but random classifier is employed to
#' build classifier.
#' @param label string or \code{NULL}. Used to give label or just a note to the output result. Default: \code{NULL}.
#' @param path.RNAsubopt string specifying the location of "RNAsubopt" program.
#' @param path.Predator string specifying the location of "Predator" program.
#' @param path.stride string specifying the location of file "stride.dat" required by program Predator.
#' @param workDir.Pro string specifying the directory for temporary files.
#' The temp files will be deleted automatically when
#' the calculation is completed. If the directory does not exist, it will be created automatically.
#' @param parallel.cores an integer that indicates the number of cores for parallel computation.
#' Default: \code{2}. Set \code{parallel.cores = -1} to run with all the cores.
#'
#' @return This function returns a data frame that contains the predicted results.
#' @details
#' The method is proposed by lncPro. This function, \code{runlncPro}, has
#' improved and fixed the original code.
#'
#' \code{runlncPro} depends on the program "RNAsubopt" of software "ViennaRNA".
#' (\url{http://www.tbi.univie.ac.at/RNA/index.html})
#'
#' Parameter \code{path.RNAsubopt} can be simply defined as \code{"RNAsubopt"} as
#' default when the OS is UNIX/Linux. However, for some OS, such as Windows, users may
#' need to specify the \code{path.RNAsubopt} if the path of "RNAsubopt" haven't been
#' added in environment variables (e.g. \code{path.RNAsubopt = '"C:/Program Files/ViennaRNA/RNAsubopt.exe"'}).
#'
#' Program "Predator" is only available on UNIX/Linux and 32-bit Windows OS.
#'
#' @section References:
#' Lu Q, Ren S, Lu M, \emph{et al}.
#' Computational prediction of associations between long non-coding RNAs and proteins.
#' BMC Genomics 2013; 14:651
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom parallel mcmapply
#' @importFrom seqinr a
#' @importFrom seqinr s2c
#' @importFrom seqinr getSequence
#' @importFrom seqinr write.fasta
#' @importFrom utils data
#' @importFrom stats predict
#'
#' @examples
#' \dontrun{
#' data(demoPositiveSeq)
#' seqRNA <- demoPositiveSeq$RNA.positive
#' seqPro <- demoPositiveSeq$Pro.positive
#'
#' Res_lncPro <- run_lncPro(seqRNA = seqRNA, seqPro = seqPro,
#'                          path.RNAsubopt = "RNAsubopt",
#'                          path.Predator = "~/predator/predator",
#'                          path.stride = "~/predator/stride.dat",
#'                          workDir.Pro = "tmp", parallel.cores = 10)
#' }
#'
#' @export

run_lncPro <- function(seqRNA, seqPro, mode = c("original", "retrained"), label = NULL,
                       path.RNAsubopt = "RNAsubopt", path.Predator = "Predator/predator",
                       path.stride = "Predator/stride.dat", workDir.Pro = getwd(),
                       parallel.cores = 2) {

        mode <- match.arg(mode)
        if (length(seqRNA) != length(seqPro)) stop("The number of RNA sequences should match the number of protein sequences!")

        message("+ Initializing...  ", Sys.time())

        if (parallel.cores == 2) message("- Users can try to set parallel.cores = -1 to use all cores!")

        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        cl <- parallel::makeCluster(parallel.cores)

        seqRNA <- parallel::parSapply(cl, seqRNA, Internal.checkRNA)

        data(aaindex, package = "seqinr", envir = environment())
        aaindex <- get("aaindex")

        message("\n", "+ Extracting features from RNA sequences...  ", Sys.time())
        lncProFeatures.RNA <- Internal.lncPro_extractFeatures(cl = cl, seqs = seqRNA, seqType = "RNA",
                                                              path.RNAsubopt = path.RNAsubopt, aaindex = aaindex)

        message("\n", "+ Extracting features from protein sequences...  ", Sys.time())
        lncProFeatures.Pro <- Internal.lncPro_extractFeatures(cl = cl, seqs = seqPro, seqType = "Pro",
                                                              path.Predator = path.Predator,
                                                              path.stride = path.stride,
                                                              workDir.Pro = workDir.Pro, aaindex = aaindex)
        parallel::stopCluster(cl)

        if (length(names(seqRNA)) != 0) names.seqRNA <- names(seqRNA) else names.seqRNA <- "RNA.noName"
        if (length(names(seqPro)) != 0) names.seqPro <- names(seqPro) else names.seqPro <- "Protein.noName"

        lncProFeatures.RNA.df <- data.frame(lncProFeatures.RNA, stringsAsFactors = FALSE)
        lncProFeatures.Pro.df <- data.frame(lncProFeatures.Pro, stringsAsFactors = FALSE)

        if (mode == "original") {
                message("\n", "+ Predicting pairs using lncPro (original algorithm)...  ", Sys.time())
                scores <- t(parallel::mcmapply(Internal.lncPro_finalScore, lncProFeatures.RNA.df,
                                               lncProFeatures.Pro.df, mc.cores = parallel.cores))
                Prediction <- ifelse(scores[,1] > 50, "Interact", "Non.Interact")

                if (is.null(label)) {
                        lncPro.result <- data.frame(RNA_Name = names.seqRNA, Pro_Name = names.seqPro,
                                                    lncPro_original = Prediction,
                                                    scores, stringsAsFactors = FALSE)
                } else {
                        lncPro.result <- data.frame(RNA_Name = names.seqRNA, Pro_Name = names.seqPro,
                                                    label = label, lncPro_original = Prediction,
                                                    scores, stringsAsFactors = FALSE)
                }
        } else {
                data(mod_lncPro, envir = environment())
                mod_lncPro <- get("mod_lncPro")

                featureSet <- cbind(lncProFeatures.RNA.df, lncProFeatures.Pro.df)
                res <- stats::predict(mod_lncPro, featureSet, type = "response")
                prob <- stats::predict(mod_lncPro, featureSet, type = "prob")

                if (is.null(label)) {
                        lncPro.result <- data.frame(RNA_Name = names.seqRNA, Pro_Name = names.seqPro,
                                                    lncPro_retrain = res,
                                                    scores, stringsAsFactors = FALSE)
                } else {
                        lncPro.result <- data.frame(RNA_Name = names.seqRNA, Pro_Name = names.seqPro,
                                                    label = label, lncPro_original = Prediction,
                                                    scores, stringsAsFactors = FALSE)
                }

                message("\n", "+ Completed.  ", Sys.time())
                lncPro.result
        }
}

#' Predict RNA-Protein Interaction Using RPISeq Method
#' @description This function can predict lncRNA/RNA-protein interactions using RPISeq method.
#' Both the web-based original version and retrained model are available. Network is required to use
#' the original version.
#'
#' @param seqRNA RNA sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of RNA/protein sequences.
#' RNA sequences will be converted into lower case letters.
#' @param seqPro protein sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of protein sequences.
#' Protein sequences will be converted into upper case letters.
#' Each sequence should be a vector of single characters.
#' @param mode set \code{"web"} to employ web-based original version.
#' Use \code{"retrained"} to call retrained model.
#' @param label string or \code{NULL}. Used to give label or just a note to the output result. Default: \code{NULL}.
#' @param parallel.cores an integer that indicates the number of cores for parallel computation.
#' Default: \code{2}. Set \code{parallel.cores = -1} to run with all the cores.
#'
#' @return This function returns a data frame that contains the predicted results.
#'
#' @section References:
#' Muppirala UK, Honavar VG, Dobbs D.
#' Predicting RNA-protein interactions using only sequence information.
#' BMC Bioinformatics 2011; 12:489
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom seqinr count
#' @importFrom seqinr getSequence
#' @importFrom RCurl postForm
#' @importFrom stats predict
#'
#' @examples
#'
#' data(demoPositiveSeq)
#' seqRNA <- demoPositiveSeq$RNA.positive
#' seqPro <- demoPositiveSeq$Pro.positive
#'
#' Res_RPISeq_1 <- run_RPISeq(seqRNA = seqRNA, seqPro = seqPro, mode = "web",
#'                            parallel.cores = 2)
#'
#' Res_RPISeq_2 <- run_RPISeq(seqRNA = seqRNA, seqPro = seqPro, mode = "retrained",
#'                            label = "Interact", parallel.cores = 2)
#'
#'
#' @export

run_RPISeq <- function(seqRNA, seqPro, mode = c("web", "retrained"), label = NULL,
                       parallel.cores = 2) {

        mode <- match.arg(mode)
        if (length(seqRNA) != length(seqPro)) stop("The number of RNA sequences should match the number of protein sequences!")

        message("+ Initializing...  ", Sys.time())

        seqRNA <- sapply(seqRNA, Internal.checkRNA)

        if (length(names(seqRNA)) != 0) names.seqRNA <- names(seqRNA) else names.seqRNA <- "RNA.noName"
        if (length(names(seqPro)) != 0) names.seqPro <- names(seqPro) else names.seqPro <- "Protein.noName"

        if (mode == "web") {
                if (parallel.cores == 2) message("- Users can try to set parallel.cores = -1 to use all cores!")
                parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
                cl <- parallel::makeCluster(parallel.cores)

                ProSeq <- sapply(seqPro, seqinr::getSequence, as.string = T)
                RNASeq <- sapply(seqRNA, seqinr::getSequence, as.string = T)
                message("\n", "+ Predicting pairs using RPISeq (web server)...  ", Sys.time())

                parallel::clusterExport(cl, c("ProSeq", "RNASeq",
                                              "Internal.get_RPISeqWeb_res"), envir = environment())
                RPISeq <- parallel::parSapply(cl, 1:length(ProSeq), function(pairNum) {
                        RPISeq_res <- Internal.get_RPISeqWeb_res(ProSeq[[pairNum]], RNASeq[[pairNum]])
                        RPISeq_res
                })
                RPISeq <- t(RPISeq)
                RPISeq <- data.frame(RPISeq, stringsAsFactors = F)
                RFres  <- ifelse(RPISeq$RF_prob  > 0.5, "Interact", "Non.Interact")
                SVMres <- ifelse(RPISeq$SVM_prob > 0.5, "Interact", "Non.Interact")

                parallel::stopCluster(cl)

                if (is.null(label)) {
                        outRes <- cbind(RNA_Name = names.seqRNA, Pro_Name = names.seqPro, RF_res = RFres, SVM_res = SVMres, RPISeq)
                } else {
                        outRes <- cbind(RNA_Name = names.seqRNA, Pro_Name = names.seqPro, label = label, RF_res = RFres, SVM_res = SVMres, RPISeq)
                }

        } else {
                message("\n", "+ Predicting pairs using RPISeq (retrained model)...  ", Sys.time())

                data(mod_RPISeq, envir = environment())
                mod_RPISeq <- get("mod_RPISeq")

                featureSet <- featureFreq(seqRNA = seqRNA, seqPro = seqPro,
                                          featureMode = "conc", computePro = "RPISeq", k.Pro = 3,
                                          k.RNA = 4, normalize = "row",
                                          parallel.cores = parallel.cores)
                res <- stats::predict(mod_RPISeq, featureSet, type = "response")
                prob <- stats::predict(mod_RPISeq, featureSet, type = "prob")

                if (is.null(label)) {
                        outRes <- cbind(RNA_Name = names.seqRNA, Pro_Name = names.seqPro,
                                        RPISeq_res = as.character(res), RPISeq_prob = prob[,1])
                        outRes <- data.frame(outRes, stringsAsFactors = F)
                } else {
                        outRes <- cbind(RNA_Name = names.seqRNA, Pro_Name = names.seqPro,
                                        RPISeq_res = as.character(res), RPISeq_prob = prob[,1])
                        outRes <- data.frame(outRes, stringsAsFactors = F)
                }
        }
        row.names(outRes) <- NULL
        message("\n", "+ Completed.  ", Sys.time(), "\n")
        outRes
}

#' Predict RNA-Protein Interaction Using rpiCOOL Method
#' @description This function can predict lncRNA/RNA-protein interactions using retrained model of rpiCOOL method.
#' The codes of this function slightly differ from the rpiCOOL's script.
#'
#' @param seqRNA RNA sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of RNA/protein sequences.
#' RNA sequences will be converted into lower case letters.
#' @param seqPro protein sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of protein sequences.
#' Protein sequences will be converted into upper case letters.
#' Each sequence should be a vector of single characters.
#' @param label string or \code{NULL}. Used to give label or just a note to the output result. Default: \code{NULL}.
#' @param parallel.cores an integer that indicates the number of cores for parallel computation.
#' Default: \code{2}. Set \code{parallel.cores = -1} to run with all the cores.
#'
#' @return This function returns a data frame that contains the predicted results.
#'
#' @section References:
#' Akbaripour-Elahabad M, Zahiri J, Rafeh R, \emph{et al}.
#' rpiCOOL: A tool for In Silico RNA-protein interaction detection using random forest.
#' J. Theor. Biol. 2016; 402:1-8
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom seqinr count
#' @importFrom seqinr getSequence
#' @importFrom stats predict
#'
#' @examples
#'
#' data(demoPositiveSeq)
#' seqRNA <- demoPositiveSeq$RNA.positive
#' seqPro <- demoPositiveSeq$Pro.positive
#'
#' Res_rpiCOOL <- run_rpiCOOL(seqRNA = seqRNA, seqPro = seqPro,
#'                            label = "rpiCOOL_res", parallel.cores = 2)
#'
#'
#' @export

run_rpiCOOL <- function(seqRNA, seqPro, label = NULL, parallel.cores = 2) {

        if (length(seqRNA) != length(seqPro)) stop("The number of RNA sequences should match the number of protein sequences!")

        message("+ Initializing...  ", Sys.time())

        seqRNA <- sapply(seqRNA, Internal.checkRNA)

        if (length(names(seqRNA)) != 0) names.seqRNA <- names(seqRNA) else names.seqRNA <- "RNA.noName"
        if (length(names(seqPro)) != 0) names.seqPro <- names(seqPro) else names.seqPro <- "Protein.noName"

        message("\n", "+ Predicting pairs using rpiCOOL (retrained model)...  ", Sys.time())

        featureSetFreq <- featureFreq(seqRNA = seqRNA, seqPro = seqPro,
                                      featureMode = "conc", computePro = "rpiCOOL", k.Pro = 3,
                                      k.RNA = 4, normalize = "none", parallel.cores = parallel.cores)

        featureSetMotif <- featureMotifs(seqRNA = seqRNA, seqPro = seqPro, featureMode = "comb",
                                         motifRNA = "rpiCOOL", motifPro = "rpiCOOL",
                                         parallel.cores = parallel.cores)
        featureSet <- cbind(featureSetFreq, featureSetMotif)

        data(mod_rpiCOOL, envir = environment())
        mod_rpiCOOL <- get("mod_rpiCOOL")

        res <- stats::predict(mod_rpiCOOL, featureSet, type = "response")
        prob <- stats::predict(mod_rpiCOOL, featureSet, type = "prob")

        if (is.null(label)) {
                outRes <- cbind(RNA_Name = names.seqRNA, Pro_Name = names.seqPro,
                                rpiCOOL_res = as.character(res), rpiCOOL_prob = prob[,1])
                outRes <- data.frame(outRes, stringsAsFactors = F, row.names = NULL)
        } else {
                outRes <- cbind(RNA_Name = names.seqRNA, Pro_Name = names.seqPro,
                                rpiCOOL_res = as.character(res), rpiCOOL_prob = prob[,1])
                outRes <- data.frame(outRes, stringsAsFactors = F, row.names = NULL)
        }
        message("\n", "+ Completed.  ", Sys.time(), "\n")
        outRes
}

#' Predict RNA-Protein Interaction Using ncProR Method
#' @description This function can predict lncRNA/RNA-protein interactions using ncProR method.
#'
#' @param seqRNA RNA sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of RNA/protein sequences.
#' RNA sequences will be converted into lower case letters.
#' @param seqPro protein sequences loaded by function \code{\link[seqinr]{read.fasta}} of package
#' "seqinr" (\code{\link[seqinr]{seqinr-package}}). Or a list of protein sequences.
#' Protein sequences will be converted into upper case letters.
#' Each sequence should be a vector of single characters.
#' @param label string or \code{NULL}. Used to give label or just a note to the output result. Default: \code{NULL}.
#' @param parallel.cores an integer that indicates the number of cores for parallel computation.
#' Default: \code{2}. Set \code{parallel.cores = -1} to run with all the cores.
#'
#' @return This function returns a data frame that contains the predicted results.
#' @section References:
#' Han S, Liang Y, Li Y, \emph{et al}.
#' ncProR: an integrated R package for effective ncRNA-protein interaction prediction.
#' (\emph{Submitted})
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom parallel parSapply
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom seqinr count
#' @importFrom seqinr getSequence
#' @importFrom stats predict
#'
#' @examples
#'
#' data(demoPositiveSeq)
#' seqRNA <- demoPositiveSeq$RNA.positive
#' seqPro <- demoPositiveSeq$Pro.positive
#'
#' Res_ncProR <- run_ncProR(seqRNA = seqRNA, seqPro = seqPro, parallel.cores = 2)
#'
#'
#' @export

run_ncProR <- function(seqRNA, seqPro, label = NULL, parallel.cores = 2) {
        if (length(seqRNA) != length(seqPro)) stop("The number of RNA sequences should match the number of protein sequences!")

        message("+ Initializing...  ", Sys.time())

        seqRNA <- sapply(seqRNA, Internal.checkRNA)

        if (length(names(seqRNA)) != 0) names.seqRNA <- names(seqRNA) else names.seqRNA <- "RNA.noName"
        if (length(names(seqPro)) != 0) names.seqPro <- names(seqPro) else names.seqPro <- "Protein.noName"

        message("\n", "+ Predicting pairs using ncProR...  ", Sys.time())

        featureSetFreq <- featureFreq(seqRNA = seqRNA, seqPro = seqPro,
                                      featureMode = "conc", computePro = "DeNovo", k.Pro = 3,
                                      k.RNA = 4, normalize = "none", parallel.cores = parallel.cores)

        featureSetMotif <- featureMotifs(seqRNA = seqRNA, seqPro = seqPro, featureMode = "conc",
                                         parallel.cores = parallel.cores)

        featureSetPhysChem <- featurePhysChem(seqRNA = seqRNA, seqPro = seqPro,
                                           physchemRNA = c("hydrogenBonding", "vanderWaal"),
                                           physchemPro = c("bulkiness.Zimmerman", "isoelectricPoint.Zimmerman"),
                                           parallel.cores = parallel.cores)

        featureSet <- cbind(featureSetFreq, featureSetMotif, featureSetPhysChem)

        data(mod_ncProR, envir = environment())
        mod_ncProR <- get("mod_ncProR")

        res <- stats::predict(mod_ncProR, featureSet, type = "response")
        prob <- stats::predict(mod_ncProR, featureSet, type = "prob")

        if (is.null(label)) {
                outRes <- cbind(RNA_Name = names.seqRNA, Pro_Name = names.seqPro,
                                ncProR_res = as.character(res), ncProR_prob = prob[,1])
                outRes <- data.frame(outRes, stringsAsFactors = F, row.names = NULL)
        } else {
                outRes <- cbind(RNA_Name = names.seqRNA, Pro_Name = names.seqPro,
                                ncProR_res = as.character(res), ncProR_prob = prob[,1])
                outRes <- data.frame(outRes, stringsAsFactors = F, row.names = NULL)
        }
        message("\n", "+ Completed.  ", Sys.time(), "\n")
        outRes
}

