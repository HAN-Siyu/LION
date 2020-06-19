Internal.randomName <- function(digit = 10, prefix = NULL, suffix = NULL) {
        randomName <- paste0(prefix, paste0(sample(c(LETTERS, letters, 0:9),
                                                   size = digit, replace = TRUE),
                                            collapse = ""), suffix)
        randomName
}

Internal.combineMotifs <- function(oneRNA, onePro) {
        newValues <- sapply(oneRNA, function(element) {
                newValue <- element * onePro
        })

        result <- as.numeric(newValues)
        result
}

Internal.checkRNA <- function(oneSeq) {
        if (length(oneSeq) == 1) stop("Error: Please check the format of the input!")
        oneSeq <- tolower(unlist(oneSeq))
        oneSeq[oneSeq %in% "t"] <- "u"
        oneSeq
}

Internal.convertSeq <- function(oneSeq, aa.index, aaindex) {
        oneSeq <- oneSeq[oneSeq %in% LETTERS]
        tmp.seq <- factor(oneSeq)
        aa.label <- seqinr::a(names(aaindex[[aa.index]]$I))
        names(aa.label) <- aaindex[[aa.index]]$I
        levels(tmp.seq) <- as.list(aa.label)
        num.seq <- as.character(tmp.seq)
        num.seq[is.na(num.seq)] <- 0
        num.seq <- as.numeric(num.seq)
        num.seq
}

# AA.LETTERS = c("G", "V", "L", "F", "P", "I", "A", "M", "T", "S",
#                "Y", "N", "Q", "W", "H", "K", "R", "E", "D", "C")

Internal.convertSeq_customize <- function(oneSeq, indexVal) {
        oneSeq <- oneSeq[oneSeq %in% LETTERS]
        tmp.seq <- factor(oneSeq)
        levels(tmp.seq) <- as.list(indexVal)
        num.seq <- as.character(tmp.seq)
        num.seq[is.na(num.seq)] <- 0
        num.seq <- as.numeric(num.seq)
        num.seq
}

Internal.convertPro_Struct <- function(oneSeq, info, path.Predator, path.stride,
                                       workDir, structurePro.mode, as.list,
                                       Fourier.len, aaindex, verbose,
                                       AA.LETTERS = c("G", "V", "L", "F", "P", "I", "A", "M", "T", "S",
                                                      "Y", "N", "Q", "W", "H", "K", "R", "E", "D", "C")) {

        index <- get("index")

        if (verbose) {
                # if(is.null(info)) {
                #         outinfo <- paste(" ", index, "length:", length(oneSeq), "aa", "\n")
                # } else {
                #         outinfo <- paste(" ", index, "of", info, "length:", length(oneSeq), "aa", "\n")
                # }

                outinfo <- paste(" ", index, "of", info, "length:", length(oneSeq), "aa", "\n")
                cat(outinfo)
        }

        randomName <- Internal.randomName(digit = 15, "/", ".fa")
        path.tmpFile <- paste0(workDir, randomName)
        predator.command <- paste0(path.Predator, " -a ", path.tmpFile, " -b", path.stride)

        oneSeq <- oneSeq[oneSeq %in% AA.LETTERS]

        seqinr::write.fasta(oneSeq, names = "AA Seq", file.out = path.tmpFile, as.string = FALSE)
        output <- system(predator.command, intern = TRUE)

        # if(file.exists(path.tmpFile)) file.remove(path.tmpFile)
        index <<- index + 1

        SS.Seq <- gsub(".", NA, output, fixed = T)
        SS.Seq <- SS.Seq[!is.na(SS.Seq)]
        SS.Seq <- gsub(" ", "", SS.Seq, fixed = T)
        SS.Seq <- SS.Seq[sapply(SS.Seq, nchar) != 0]
        idx.1  <- which(sapply(SS.Seq, substr, 1, 1) == ">")
        SS.Seq <- SS.Seq[idx.1:length(SS.Seq)]
        SS.Seq <- SS.Seq[seq(3, length(SS.Seq), 2)]
        SS.Seq <- paste0(SS.Seq, collapse = "")
        SS.Seq <- seqinr::s2c(SS.Seq)

        # if(length(SS.Seq) != length(unlist(oneSeq))) stop("Extracting Secondary Structure of Protein Failed!")
        if (length(SS.Seq) == length(unlist(oneSeq))) {
                if (file.exists(path.tmpFile)) file.remove(path.tmpFile)
                # write.fasta(oneSeq, names = "Question Seq", file.out = paste0("QuestionSeq.", index, ".fa"), as.string = FALSE)
        } else {
                message("Warning: Structural Seq length: ", length(SS.Seq), ", AA Seq length: ", length(unlist(oneSeq)))
        }

        alphaHelix.idx <- which(SS.Seq == "H")
        betaTurn.idx   <- which(SS.Seq == "_")
        betaSheet.idx  <- which(SS.Seq == "E")

        if ("ChouFasman" %in% structurePro.mode) {
                alphaHelix <- Internal.convertSeq(oneSeq[alphaHelix.idx], aa.index = 38, aaindex = aaindex)
                betaTurn   <- Internal.convertSeq(oneSeq[betaTurn.idx],   aa.index = 37, aaindex = aaindex)
                betaSheet  <- Internal.convertSeq(oneSeq[betaSheet.idx],  aa.index = 39, aaindex = aaindex)

                ChFa.seq <- oneSeq
                ChFa.seq[alphaHelix.idx] <- alphaHelix
                ChFa.seq[betaTurn.idx]   <- betaTurn
                ChFa.seq[betaSheet.idx]  <- betaSheet
                ChFa.seq <- as.numeric(ChFa.seq)
                ChFa.num <- Internal.FourierSeries(ChFa.seq, profile.length = Fourier.len)
                ChouFasman <- ChFa.num
        }

        if ("Levitt" %in% structurePro.mode) {
                alphaHelix <- Internal.convertSeq(oneSeq[alphaHelix.idx], aa.index = 160, aaindex = aaindex)
                betaTurn   <- Internal.convertSeq(oneSeq[betaTurn.idx],   aa.index = 162, aaindex = aaindex)
                betaSheet  <- Internal.convertSeq(oneSeq[betaSheet.idx],  aa.index = 161, aaindex = aaindex)

                Levi.seq <- oneSeq
                Levi.seq[alphaHelix.idx] <- alphaHelix
                Levi.seq[betaTurn.idx]   <- betaTurn
                Levi.seq[betaSheet.idx]  <- betaSheet
                Levi.seq <- as.numeric(Levi.seq)
                Levi.num <- Internal.FourierSeries(Levi.seq, profile.length = Fourier.len)
                Levitt <- Levi.num
        }

        if ("DeleageRoux" %in% structurePro.mode) {
                data(aaindex, package = "seqinr", envir = environment())
                DeleageRoux.alphaHelix <- seqinr::a(names(aaindex[[1]]$I))
                DeleageRoux.betaTurn   <- DeleageRoux.alphaHelix
                DeleageRoux.betaSheet  <- DeleageRoux.alphaHelix

                names(DeleageRoux.alphaHelix) <- c(1.489, 1.224, 0.772, 0.924, 0.966,
                                                   1.164, 1.504, 0.510, 1.003, 1.003,
                                                   1.236, 1.172, 1.363, 1.195, 0.492,
                                                   0.739, 0.785, 1.090, 0.787, 0.990)
                names(DeleageRoux.betaTurn)   <- c(0.788, 0.912, 1.572, 1.197, 0.965,
                                                   0.997, 1.149, 1.860, 0.970, 0.240,
                                                   0.670, 1.302, 0.436, 0.624, 1.415,
                                                   1.316, 0.739, 0.546, 0.795, 0.387)
                names(DeleageRoux.betaSheet)  <- c(0.709, 0.920, 0.604, 0.541, 1.191,
                                                   0.840, 0.567, 0.657, 0.863, 1.799,
                                                   1.261, 0.721, 1.210, 1.393, 0.354,
                                                   0.928, 1.221, 1.306, 1.266, 1.965)

                alphaHelix <- Internal.convertSeq_customize(oneSeq[alphaHelix.idx], DeleageRoux.alphaHelix)
                betaTurn   <- Internal.convertSeq_customize(oneSeq[betaTurn.idx],   DeleageRoux.betaTurn)
                betaSheet  <- Internal.convertSeq_customize(oneSeq[betaSheet.idx],  DeleageRoux.betaSheet)

                DeRo.seq <- oneSeq
                DeRo.seq[alphaHelix.idx] <- alphaHelix
                DeRo.seq[betaTurn.idx]   <- betaTurn
                DeRo.seq[betaSheet.idx]  <- betaSheet
                DeRo.seq <- as.numeric(DeRo.seq)
                DeRo.num <- Internal.FourierSeries(DeRo.seq, profile.length = Fourier.len)
                DeleageRoux <- DeRo.num
        }

        if (as.list) {
                out.seq <- mget(structurePro.mode)
        } else {
                out.seq <- unlist(mget(structurePro.mode))
        }
}

Internal.FourierSeries <- function(Num.Seq, profile.length = 10){
        LenSeq <- length(Num.Seq)
        out.val <- c()
        i <- 1
        while (i <= profile.length) {
                tmp <- 0

                sapply(c(1:LenSeq), function(x, .i = i, .LenSeq = LenSeq){
                        tmp <<- tmp + (Num.Seq[[x]] * cos(pi * (x - 0.5) * (.i - 0.5) / .LenSeq))
                })

                tmp <- sqrt(2 / LenSeq) * tmp
                out.val <- c(out.val, tmp)
                i <- i + 1
        }
        out.val
}

Internal.computeFreq <- function(oneSeq, seqType = c("RNA", "Pro"),
                                 computePro = c("RPISeq", "DeNovo", "rpiCOOL"),
                                 k = 3) {

        if (seqType == "RNA") {
                oneSeq <- Internal.checkRNA(oneSeq)
                freq.raw <- seqinr::count(oneSeq, wordsize = k, freq = TRUE, alphabet = c("a", "c", "g", "u"))
        } else {
                oneSeq <- toupper(unlist(oneSeq))
                if (computePro == "RPISeq") {
                        oneSeq[oneSeq %in% c("G", "V")] <- "A"
                        oneSeq[oneSeq %in% c("L", "F", "P")] <- "I"
                        oneSeq[oneSeq %in% c("M", "T", "S")] <- "Y"
                        oneSeq[oneSeq %in% c("N", "Q", "W")] <- "H"
                        oneSeq[oneSeq %in% c("K")] <- "R"
                        oneSeq[oneSeq %in% c("E")] <- "D"
                        freq.raw <- seqinr::count(oneSeq, wordsize = k, freq = TRUE,
                                                  alphabet = c("A", "I", "Y", "H", "R", "D", "C"))
                } else if (computePro == "DeNovo") {
                        oneSeq[oneSeq %in% c("E")] <- "D"
                        oneSeq[oneSeq %in% c("R", "K")] <- "H"
                        oneSeq[oneSeq %in% c("G", "N", "Q", "S", "T", "Y")] <- "C"
                        oneSeq[oneSeq %in% c("F", "I", "L", "M", "P", "V", "W")] <- "A"
                        freq.raw <- seqinr::count(oneSeq, wordsize = k, freq = TRUE,
                                                  alphabet = c("A", "C", "D", "H"))
                } else if (computePro == "rpiCOOL") {
                        oneSeq[oneSeq %in% c("E")] <- "A"
                        oneSeq[oneSeq %in% c("L", "F", "M", "V")] <- "I"
                        oneSeq[oneSeq %in% c("D", "T", "S")] <- "N"
                        oneSeq[oneSeq %in% c("H", "Q", "K")] <- "R"
                        oneSeq[oneSeq %in% c("W")] <- "Y"
                        freq.raw <- seqinr::count(oneSeq, wordsize = k, freq = TRUE,
                                                  alphabet = c("A", "I", "N", "G", "R", "P", "Y", "C"))
                }
        }
        return(freq.raw)
}

Internal.computeMotifs <- function(oneSeq, motifs = NULL) {
        countMotifs <- sapply(motifs, function(motif) {
                infoMotif <- unlist(gregexpr(motif, oneSeq, ignore.case = TRUE, perl = TRUE))
                countMotif <- ifelse(infoMotif[1] == -1, 0, length(infoMotif))
                countMotif
        })
        countMotifs
}

Internal.featureAAindex <- function(oneSeq, aaindex, entry.index = entry.index, Fourier.len = Fourier.len) {
        output_seq <- lapply(entry.index, function(aa.index) {
                num_seq <- Internal.convertSeq(oneSeq, aa.index = aa.index, aaindex = aaindex)
                Fourier_seq <- Internal.FourierSeries(num_seq, profile.length = Fourier.len)
                names(Fourier_seq) <- paste(aaindex[[aa.index]]$H, 1:Fourier.len, sep = ".")
                Fourier_seq
        })
        aa_seq <- unlist(output_seq)
        aa_seq
}

Internal.computePhysChem <- function(oneSeq, seqType = c("RNA", "Pro"), Fourier.len = 10, as.list = TRUE,
                                     physchemRNA.mode = c("hydrogenBonding", "vanderWaal"), aaindex,
                                     physchemPro.mode = c("polarity.Grantham", "polarity.Zimmerman",
                                                          "bulkiness.Zimmerman", "isoelectricPoint.Zimmerman",
                                                          "hphob.BullBreese", "hphob.KyteDoolittle",
                                                          "hphob.Eisenberg", "hphob.HoppWoods")) {

        if (seqType == "RNA") {
                oneSeq <- Internal.checkRNA(oneSeq)
                Hyd.seq <- rep("0", length(oneSeq))
                Van.seq <- rep("0", length(oneSeq))

                a.idx <- oneSeq %in% "a"
                c.idx <- oneSeq %in% "c"
                g.idx <- oneSeq %in% "g"
                u.idx <- oneSeq %in% "u"
                out.seq <- c()

                if ("hydrogenBonding" %in% physchemRNA.mode) {
                        Hyd.seq[a.idx] <- "24 17 40 10"
                        Hyd.seq[c.idx] <- "49 21 26"
                        Hyd.seq[g.idx] <- "21 86 17 41 29"
                        Hyd.seq[u.idx] <- "24 17 22"
                        # Hyd.seq[!(a.idx | c.idx | g.idx | u.idx)] <- "0"
                        # Hyd.seq <- paste(Hyd.seq, collapse = " ")
                        Hyd.seq <- unlist(strsplit(Hyd.seq, split = " "))
                        Hyd.seq <- as.numeric(Hyd.seq)
                        hydrogenBonding <- Internal.FourierSeries(Hyd.seq, profile.length = Fourier.len)
                }

                if ("vanderWaal" %in% physchemRNA.mode) {
                        Van.seq[a.idx] <- "79 98 69 40 53 37 84 62 49 28"
                        Van.seq[c.idx] <- "14 44 98 42 30 50 39 19"
                        Van.seq[g.idx] <- "26 74 24 37 22 21 19 67 48 44 21"
                        Van.seq[u.idx] <- "25 42 74 53 43 67 44 24"
                        # Van.seq[!(a.idx | c.idx | g.idx | u.idx)] <- "0"
                        # Van.seq <- paste(Van.seq, collapse = " ")
                        Van.seq <- unlist(strsplit(Van.seq, split = " "))
                        Van.seq <- as.numeric(Van.seq)
                        vanderWaal <- Internal.FourierSeries(Van.seq, profile.length = Fourier.len)
                }
                if (as.list) {
                        out.seq <- mget(physchemRNA.mode)
                } else {
                        out.seq <- unlist(mget(physchemRNA.mode))
                }
        } else {
                if ("hphob.BullBreese" %in% physchemPro.mode) {
                        BullBresse.index <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                                              "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
                        names(BullBresse.index) <- c(0.61, 0.69, 0.89, 0.61, 0.36,
                                                     0.97, 0.51, 0.81, 0.69, -1.45,
                                                     -1.65, 0.46, -0.66, -1.52, -0.17,
                                                     0.42, 0.29, -1.20, -1.43, -0.75)
                        oneSeq <- toupper(unlist(oneSeq))
                        BuBr.seq <- Internal.convertSeq_customize(oneSeq, BullBresse.index)
                        hphob.BullBreese <- Internal.FourierSeries(BuBr.seq, profile.length = Fourier.len)
                }

                if ("polarity.Grantham" %in% physchemPro.mode) polarity.Grantham <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 111, aaindex), profile.length = Fourier.len)
                if ("polarity.Zimmerman" %in% physchemPro.mode) polarity.Zimmerman <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 400, aaindex), profile.length = Fourier.len)
                if ("bulkiness.Zimmerman" %in% physchemPro.mode) bulkiness.Zimmerman <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 399, aaindex), profile.length = Fourier.len)
                if ("isoelectricPoint.Zimmerman" %in% physchemPro.mode) isoelectricPoint.Zimmerman <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 401, aaindex), profile.length = Fourier.len)
                if ("hphob.KyteDoolittle" %in% physchemPro.mode) hphob.KyteDoolittle <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 151, aaindex), profile.length = Fourier.len)
                if ("hphob.Eisenberg" %in% physchemPro.mode) hphob.Eisenberg <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 68, aaindex), profile.length = Fourier.len)
                if ("hphob.HoppWoods" %in% physchemPro.mode) hphob.HoppWoods <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 115, aaindex), profile.length = Fourier.len)

                if (as.list) {
                        out.seq <- mget(physchemPro.mode)
                } else {
                        out.seq <- unlist(mget(physchemPro.mode))
                }
        }
        out.seq
}

Internal.runRNAsubopt <- function(oneSeq, info = NULL, structureRNA.num = 6,
                                  path.RNAsubopt, outputNum = FALSE, verbose = FALSE){
        index <- get("index")

        if (verbose) {
                outinfo <- paste(" ", index, "of", info, "length:", nchar(oneSeq), "nt", "\n")
                cat(outinfo)
        }

        RNAsubopt.cmd <- paste(path.RNAsubopt, "-p ", structureRNA.num)
        RNAsubopt.seq <- system(RNAsubopt.cmd, intern = TRUE, input = oneSeq)
        index <<- index + 1
        if (outputNum) {
                RNAsubopt.seq <- RNAsubopt.seq[(length(RNAsubopt.seq) - structureRNA.num + 1):length(RNAsubopt.seq)]
                RNAsubopt.num <- sapply(RNAsubopt.seq, function(x) {
                        x <- seqinr::s2c(x)
                        x[x %in% "."] <- "0"
                        x[x %in% c("(", ")")] <- "1"
                        x <- as.numeric(x)
                        return(x)
                })
                RNAsubopt.sum <- rowSums(RNAsubopt.num)
                return(RNAsubopt.sum)
        } else  return(RNAsubopt.seq)
}

Internal.runPredator <- function(oneSeq, info = NULL, path.Predator, path.stride, workDir = getwd(),
                                 as.string = FALSE, outputRaw = FALSE, verbose = FALSE,
                                 AA.LETTERS = c("G", "V", "L", "F", "P", "I", "A", "M", "T", "S",
                                                "Y", "N", "Q", "W", "H", "K", "R", "E", "D", "C")) {

        index <- get("index")

        if (verbose) {
                # if(is.null(info)) {
                #         outinfo <- paste(" ", index, "length:", length(oneSeq), "aa", "\n")
                # } else {
                #         outinfo <- paste(" ", index, "of", info, "length:", length(oneSeq), "aa", "\n")
                # }

                outinfo <- paste(" ", index, "of", info, "length:", length(oneSeq), "aa", "\n")
                cat(outinfo)
        }

        randomName <- Internal.randomName(digit = 15, "/", ".fa")
        path.tmpFile <- paste0(workDir, randomName)
        predator.command <- paste0(path.Predator, " -a ", path.tmpFile, " -b", path.stride)

        oneSeq <- oneSeq[oneSeq %in% AA.LETTERS]

        seqinr::write.fasta(oneSeq, names = "AA Seq", file.out = path.tmpFile, as.string = FALSE)
        output <- system(predator.command, intern = TRUE)

        # if(file.exists(path.tmpFile)) file.remove(path.tmpFile)
        index <<- index + 1

        SS.Seq <- gsub(".", NA, output, fixed = T)
        SS.Seq <- SS.Seq[!is.na(SS.Seq)]
        SS.Seq <- gsub(" ", "", SS.Seq, fixed = T)
        SS.Seq <- SS.Seq[sapply(SS.Seq, nchar) != 0]
        idx.1  <- which(sapply(SS.Seq, substr, 1, 1) == ">")
        SS.Seq <- SS.Seq[idx.1:length(SS.Seq)]
        SS.Seq <- SS.Seq[seq(3, length(SS.Seq), 2)]
        SS.Seq <- paste0(SS.Seq, collapse = "")

        # if(nchar(SS.Seq) != length(unlist(oneSeq))) stop("Extracting Secondary Structure of Protein Failed!")
        if (nchar(SS.Seq) == length(unlist(oneSeq))) {
                if (file.exists(path.tmpFile)) file.remove(path.tmpFile)
                # write.fasta(oneSeq, names = "Question Seq", file.out = paste0("QuestionSeq.", index, ".fa"), as.string = FALSE)
        } else {
                message("Warning: Structural Seq length: ", length(SS.Seq), ", AA Seq length: ", length(unlist(oneSeq)))
        }

        if (outputRaw) SS.Seq <- c(SS.Seq, seqinr::c2s(oneSeq))

        if (!as.string) SS.Seq <- seqinr::s2c(SS.Seq)

        SS.Seq
}

Internal.convertPro_Struct.lncPro <- function(oneSeq, info, path.Predator, path.stride,
                                              workDir, Fourier.len, aaindex,
                                              AA.LETTERS = c("G", "V", "L", "F", "P", "I", "A", "M", "T", "S",
                                                             "Y", "N", "Q", "W", "H", "K", "R", "E", "D", "C")) {

        index <- get("index")

        randomName <- Internal.randomName(digit = 15, "/", ".fa")
        path.tmpFile <- paste0(workDir, randomName)
        predator.command <- paste0(path.Predator, " -a ", path.tmpFile, " -b", path.stride)

        oneSeq <- oneSeq[oneSeq %in% AA.LETTERS]

        seqinr::write.fasta(oneSeq, names = "AA Seq", file.out = path.tmpFile, as.string = FALSE)
        output <- system(predator.command, intern = TRUE)

        # if(file.exists(path.tmpFile)) file.remove(path.tmpFile)
        index <<- index + 1

        SS.Seq <- gsub(".", NA, output, fixed = T)
        SS.Seq <- SS.Seq[!is.na(SS.Seq)]
        SS.Seq <- gsub(" ", "", SS.Seq, fixed = T)
        SS.Seq <- SS.Seq[sapply(SS.Seq, nchar) != 0]
        idx.1  <- which(sapply(SS.Seq, substr, 1, 1) == ">")
        SS.Seq <- SS.Seq[idx.1:length(SS.Seq)]
        SS.Seq <- SS.Seq[seq(3, length(SS.Seq), 2)]
        SS.Seq <- paste0(SS.Seq, collapse = "")
        SS.Seq <- seqinr::s2c(SS.Seq)

        if (length(SS.Seq) == length(unlist(oneSeq))) {
                if (file.exists(path.tmpFile)) file.remove(path.tmpFile)
        }

        alphaHelix.idx <- which(SS.Seq == "H")
        betaTurn.idx   <- which(SS.Seq == "_")
        betaSheet.idx  <- which(SS.Seq == "E")

        alphaHelix <- Internal.convertSeq(oneSeq[alphaHelix.idx], aa.index = 38, aaindex = aaindex)
        betaTurn   <- Internal.convertSeq(oneSeq[betaTurn.idx],   aa.index = 37, aaindex = aaindex)
        betaSheet  <- Internal.convertSeq(oneSeq[betaSheet.idx],  aa.index = 39, aaindex = aaindex)

        ChFa.seq <- oneSeq
        ChFa.seq[alphaHelix.idx] <- alphaHelix
        ChFa.seq[betaTurn.idx]   <- betaTurn
        ChFa.seq[betaSheet.idx]  <- betaSheet
        ChFa.seq <- as.numeric(ChFa.seq)

        ChFa.seq
}

Internal.computePhysChem <- function(oneSeq, seqType = c("RNA", "Pro"), Fourier.len = 10, as.list = TRUE,
                                     physchemRNA.mode = c("hydrogenBonding", "vanderWaal"), aaindex,
                                     physchemPro.mode = c("polarity.Grantham", "polarity.Zimmerman",
                                                          "bulkiness.Zimmerman", "isoelectricPoint.Zimmerman",
                                                          "hphob.BullBreese", "hphob.KyteDoolittle",
                                                          "hphob.Eisenberg", "hphob.HoppWoods")) {
        if (seqType == "RNA") {
                oneSeq <- Internal.checkRNA(oneSeq)
                Hyd.seq <- oneSeq
                Van.seq <- oneSeq

                a.idx <- oneSeq %in% "a"
                c.idx <- oneSeq %in% "c"
                g.idx <- oneSeq %in% "g"
                u.idx <- oneSeq %in% "u"
                out.seq <- c()

                if ("hydrogenBonding" %in% physchemRNA.mode) {
                        Hyd.seq[a.idx] <- "24 17 40 10"
                        Hyd.seq[c.idx] <- "49 21 26"
                        Hyd.seq[g.idx] <- "21 86 17 41 29"
                        Hyd.seq[u.idx] <- "24 17 22"
                        Hyd.seq[!(a.idx | c.idx | g.idx | u.idx)] <- "0"
                        Hyd.seq <- paste(Hyd.seq, collapse = " ")
                        Hyd.seq <- unlist(strsplit(Hyd.seq, split = " "))
                        Hyd.seq <- as.numeric(Hyd.seq)
                        hydrogenBonding <- Internal.FourierSeries(Hyd.seq, profile.length = Fourier.len)

                        # print("hydrogenBonding: ")
                        # print(hydrogenBonding)
                }

                if ("vanderWaal" %in% physchemRNA.mode) {
                        Van.seq[a.idx] <- "79 98 69 40 53 37 84 62 49 28"
                        Van.seq[c.idx] <- "14 44 98 42 30 50 39 19"
                        Van.seq[g.idx] <- "26 74 24 37 22 21 19 67 48 44 21"
                        Van.seq[u.idx] <- "25 42 74 53 43 67 44 24"
                        Van.seq[!(a.idx | c.idx | g.idx | u.idx)] <- "0"
                        Van.seq <- paste(Van.seq, collapse = " ")
                        Van.seq <- unlist(strsplit(Van.seq, split = " "))
                        Van.seq <- as.numeric(Van.seq)
                        vanderWaal <- Internal.FourierSeries(Van.seq, profile.length = Fourier.len)

                        # print("vanderWaal: ")
                        # print(vanderWaal)
                }
                if (as.list) {
                        out.seq <- mget(physchemRNA.mode)
                } else {
                        out.seq <- unlist(mget(physchemRNA.mode))
                }
        } else {
                if ("hphob.BullBreese" %in% physchemPro.mode) {
                        BullBresse.index <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                                              "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
                        names(BullBresse.index) <- c(0.61, 0.69, 0.89, 0.61, 0.36,
                                                     0.97, 0.51, 0.81, 0.69, -1.45,
                                                     -1.65, 0.46, -0.66, -1.52, -0.17,
                                                     0.42, 0.29, -1.20, -1.43, -0.75)
                        oneSeq <- toupper(unlist(oneSeq))
                        BuBr.seq <- Internal.convertSeq_customize(oneSeq, BullBresse.index)
                        hphob.BullBreese <- Internal.FourierSeries(BuBr.seq, profile.length = Fourier.len)
                }

                if ("polarity.Grantham" %in% physchemPro.mode) polarity.Grantham <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 111, aaindex), profile.length = Fourier.len)
                if ("polarity.Zimmerman" %in% physchemPro.mode) polarity.Zimmerman <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 400, aaindex), profile.length = Fourier.len)
                if ("bulkiness.Zimmerman" %in% physchemPro.mode) bulkiness.Zimmerman <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 399, aaindex), profile.length = Fourier.len)
                if ("isoelectricPoint.Zimmerman" %in% physchemPro.mode) isoelectricPoint.Zimmerman <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 401, aaindex), profile.length = Fourier.len)
                if ("hphob.KyteDoolittle" %in% physchemPro.mode) hphob.KyteDoolittle <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 151, aaindex), profile.length = Fourier.len)
                if ("hphob.Eisenberg" %in% physchemPro.mode) hphob.Eisenberg <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 68, aaindex), profile.length = Fourier.len)
                if ("hphob.HoppWoods" %in% physchemPro.mode) hphob.HoppWoods <- Internal.FourierSeries(Internal.convertSeq(oneSeq, aa.index = 115, aaindex), profile.length = Fourier.len)

                # print("polarity.Grantham: ")
                # print(polarity.Grantham)
                #
                # print("polarity.Zimmerman: ")
                # print(polarity.Zimmerman)
                #
                # print("hphob.KyteDoolittle: ")
                # print(hphob.KyteDoolittle)
                #
                # print("hphob.BullBreese: ")
                # print(hphob.BullBreese)

                if (as.list) {
                        out.seq <- mget(physchemPro.mode)
                } else {
                        out.seq <- unlist(mget(physchemPro.mode))
                }
        }
        out.seq
}

Internal.lncPro_SS <- function(cl, seqs, seqType = c("RNA", "AA"), aaindex = NULL,
                               path.RNAsubopt = NULL, path.Predator = NULL,
                               path.stride = NULL, workDir.Pro = NULL) {
        if (seqType == "RNA") {
                index <- 1
                seqs <- parallel::parLapply(cl, seqs, Internal.checkRNA)
                seq.string <- parallel::parSapply(cl, seqs, seqinr::getSequence, TRUE)

                parallel::clusterExport(cl, varlist = c("index"), envir = environment())

                RNAsubopt.seq <- parallel::parLapply(cl, seq.string, Internal.runRNAsubopt,
                                                     structureRNA.num = 6, verbose = FALSE,
                                                     path.RNAsubopt = path.RNAsubopt, outputNum = TRUE)

                out.seq <- parallel::parLapply(cl, RNAsubopt.seq, Internal.FourierSeries, profile.length = 10)


        } else {
                index <- 1

                parallel::clusterExport(cl, varlist = c("index", "Internal.randomName", "aaindex", "Internal.convertSeq",
                                                        "Internal.convertSeq_customize", "Internal.FourierSeries"),
                                        envir = environment())

                if (!dir.exists(workDir.Pro)) {
                        dir.create(workDir.Pro, recursive = TRUE)
                        createDir <- TRUE
                } else {
                        createDir <- FALSE
                }

                Predator.seq <- parallel::parLapply(cl, seqs, Internal.convertPro_Struct.lncPro,
                                                    path.Predator = path.Predator,
                                                    path.stride = path.stride, workDir = workDir.Pro,
                                                    Fourier.len = 10, aaindex = aaindex)

                out.seq <- parallel::parLapply(cl, Predator.seq, Internal.FourierSeries, profile.length = 10)

                if (createDir & length(dir(workDir.Pro, all.files = TRUE, recursive = TRUE)) == 0) unlink(workDir.Pro, recursive = TRUE)
        }
        out.seq
}

Internal.lncPro_extractFeatures <- function(cl, seqs, seqType = c("RNA", "Pro"), path.RNAsubopt = NULL,
                                            path.Predator = NULL, path.stride = NULL,
                                            workDir.Pro = NULL, aaindex = NULL) {

        if (seqType == "RNA") {

                # seqValidate <- seqs[which(lengths(seqs) <= 4095)]
                # if (length(seqValidate) < length(seqs)) {
                #         message("- Due to the limitation of RNAsubopt,")
                #         message("- sequences with length more than 4095 nt will be omitted.")
                #         message("- ", length(seqs) - length(seqValidate), " of ", length(seqs), " sequences have been removed.")
                # }
                # seqs <- seqValidate

                parallel::clusterExport(cl, varlist = c("Internal.checkRNA", "Internal.FourierSeries",
                                                        "Internal.convertSeq_customize", "Internal.convertSeq",
                                                        "aaindex"), envir = environment())

                message("- RNA Sequences Number: ", length(seqs))

                message("- Extracting structural features...")

                output.RNA_Structure <- Internal.lncPro_SS(cl = cl, seqs = seqs, seqType = "RNA", aaindex = aaindex,
                                                           path.RNAsubopt = path.RNAsubopt)

                message("- Extracting physicochemical features...")

                output.RNA_PhysChem <- parallel::parSapply(cl, seqs, Internal.computePhysChem, seqType = "RNA",
                                                           Fourier.len = 10, as.list = TRUE, aaindex = aaindex,
                                                           physchemRNA.mode = c("hydrogenBonding", "vanderWaal"))

                outputFeatures <- rbind(RNAStructure = output.RNA_Structure, output.RNA_PhysChem)

        } else {

                # if (!file.exists(path.stride)) stop("The path of stride.dat is not correct! Please check parameter path.stride.")

                # seqValidate <- seqs[which(lengths(seqs) >= 30)]
                # if (length(seqValidate) < length(seqs)) {
                #         message("- Due to the limitation of predator,")
                #         message("- sequences with length less than 30 amino acids will be omitted.")
                #         message("- ", length(seqs) - length(seqValidate), " of ", length(seqs), " sequences have been removed.")
                # }
                # seqs <- seqValidate

                message("- Protein Sequences Number: ", length(seqs))

                message("- Extracting structural features...")

                output.Pro_Structure <- Internal.lncPro_SS(cl = cl, seqs = seqs, seqType = "Pro", aaindex = aaindex,
                                                           path.Predator = path.Predator, path.stride = path.stride,
                                                           workDir.Pro = workDir.Pro)

                message("- Extracting physicochemical features...")

                output.Pro_PhysChem <- parallel::parSapply(cl, seqs, Internal.computePhysChem, seqType = "Pro",
                                                           Fourier.len = 10, as.list = TRUE, aaindex = aaindex,
                                                           physchemPro = c("polarity.Grantham", "polarity.Zimmerman",
                                                                           "hphob.KyteDoolittle", "hphob.BullBreese"))

                outputFeatures <- rbind(output.Pro_PhysChem, ProteinStructure = output.Pro_Structure)
                outputFeatures <- outputFeatures[c(5, 1:4),,drop = FALSE]
        }
        outputFeatures
}

Internal.lncPro_computeScore <- function(Data1, Data2, weight.w) {
        Data1 <- unlist(Data1)
        Data2 <- unlist(Data2)
        tmpScore <- sapply(Data1, function(x) {
                .score <- x * Data2
        })
        tmpScore <- matrix(data = as.numeric(t(tmpScore)), nrow = 1)

        outScore <- as.numeric(tmpScore %*% weight.w)
}

Internal.lncPro_transferScore <- function(Score, m1, m2, w) {
        c1 <- m1 %*% w
        c2 <- m2 %*% w
        cc <- as.numeric((c1 + c2) / 2)

        normalScore = (100 / pi) * atan(2 * (Score - cc) / (c1 - c2)) + 50
}

Internal.lncPro_finalScore <- function(lncProFeatures.RNA, lncProFeatures.Pro, showDetails = TRUE) {

        Score.1 <- Internal.lncPro_computeScore(Data1 = lncProFeatures.RNA[1],
                                                Data2 = lncProFeatures.Pro[1], weight.w = w_1) # SS -SS

        Score.2 <- Internal.lncPro_computeScore(Data1 = lncProFeatures.RNA[2],
                                                Data2 = lncProFeatures.Pro[2], weight.w = w_2) # Hyd-bonding - Grantham

        Score.3 <- Internal.lncPro_computeScore(Data1 = lncProFeatures.RNA[2],
                                                Data2 = lncProFeatures.Pro[3], weight.w = w_3) # Hyd-bonding - Zimmerman

        Score.4 <- Internal.lncPro_computeScore(Data1 = lncProFeatures.RNA[3],
                                                Data2 = lncProFeatures.Pro[4], weight.w = w_4) # Van.der.Waal - KyteDoolittle

        Score.5 <- Internal.lncPro_computeScore(Data1 = lncProFeatures.RNA[3],
                                                Data2 = lncProFeatures.Pro[5], weight.w = w_5) # Van.der.Waal - BullBreese


        finalScores <- sapply(c(1:5), function(x) {
                weight.m1 <- get(paste0("m_", x, "_1"))
                weight.m2 <- get(paste0("m_", x, "_2"))
                weight.w  <- get(paste0("w_", x))
                .Score    <- get(paste0("Score.", x))

                normalizedScore <- Internal.lncPro_transferScore(Score = .Score, m1 = weight.m1, m2 = weight.m2, w = weight.w)
        })

        finalScore <- mean(finalScores)

        if (showDetails) {
                finalScore <- c(finalScore = finalScore,
                                structure_structure = finalScores[1],
                                hydrogenBonding_polarity.Grantham = finalScores[2],
                                hydrogenBonding_polarity.Zimmerman = finalScores[3],
                                vanderWaal_hphob.KyteDoolittle = finalScores[4],
                                vanderWaal_hphob.BullBreese = finalScores[5])
        }
        finalScore
}


Internal.get_RPISeqWeb_res <- function(onePro, oneRNA) {

        resPage <- RCurl::postForm("http://pridb.gdcb.iastate.edu/RPISeq/results.php",
                                   .params = list(p_input = onePro, r_input = oneRNA),
                                   .opts = list(referer = "http://pridb.gdcb.iastate.edu/RPISeq/results.php"))

        part_1 <- strsplit(resPage, '<table border="0"><tr><td><b><u>Interaction probabilities</u></b></td></tr><tr /><tr /><tr><td>Prediction using RF classifier </td> <td></td><td>')
        part_2 <- strsplit(part_1[[1]][[2]], "</td></tr><tr><td>Prediction using SVM classifier </td><td></td><td>")
        RF_res <- as.numeric(part_2[[1]][[1]])
        part_3 <- strsplit(part_2[[1]][[2]], "</td></tr></table>")
        SVM_res <- as.numeric(part_3[[1]][[1]])
        res <- c(RF_prob = RF_res, SVM_prob = SVM_res)
        res
}

Internal.randomForest_CV <- function(datasets = list(), all_folds, label.col = 1,
                                     positive.class = NULL, ntree = 1500,
                                     folds.num = folds.num, cl, ...) {

        perf.res <- parallel::parSapply(cl, 1:folds.num, function(x, datasets,
                                                                  all_folds,
                                                                  positive.class,
                                                                  ntree) {
                trainSet <- c()
                testSet <- c()
                for (i in 1:length(datasets)) {
                        dataset_tmp <- datasets[[i]]
                        fold_tmp <- all_folds[[i]]
                        train_tmp <- dataset_tmp[fold_tmp[[x]],]
                        test_tmp <- dataset_tmp[-fold_tmp[[x]],]

                        trainSet <- rbind(trainSet, train_tmp)
                        testSet <- rbind(testSet, test_tmp)
                }

                RF_mod <- randomForest::randomForest(label ~ ., data = trainSet,
                                                     ntree = ntree, ...)
                res <- predict(RF_mod, testSet, type = "response")

                confusion.res <- caret::confusionMatrix(data.frame(res)[,1], testSet$label,
                                                        positive = positive.class,
                                                        mode = "everything")
                TP <- confusion.res$table[1]
                FN <- confusion.res$table[2]
                FP <- confusion.res$table[3]
                TN <- confusion.res$table[4]
                N  <- sum(confusion.res$table)
                S  <- (TP + FN) / N
                P  <- (TP + FP) / N
                MCC <- ((TP / N) - (S * P)) / sqrt(P * S * (1 - S) * (1 - P))
                # MCC <- ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

                performance.res <- data.frame(TP = TP, TN = TN, FP = FP, FN = FN,
                                              Sensitivity = confusion.res$byClass[1],
                                              Specificity = confusion.res$byClass[2],
                                              Accuracy    = confusion.res$overall[1],
                                              F.Measure   = confusion.res$byClass[7],
                                              MCC         = MCC,
                                              Kappa       = confusion.res$overall[2])

        }, datasets = datasets, all_folds = all_folds, positive.class = positive.class, ntree = ntree)

        Ave.res <- apply(perf.res, 1, as.numeric)
        Ave.res <- as.data.frame(t(Ave.res))
        Ave.res$Ave.Res <- rowMeans(Ave.res)
        Ave.res
}

Internal.randomForest_tune <- function(datasets = list(), label.col = 1,
                                       positive.class = NULL, folds.num = 10,
                                       ntree.range = c(200, 500, 1000, 1500, 2000),
                                       seed = 1, return.model = TRUE,
                                       parallel.cores = 2, ...) {

        for (i in 1:length(datasets)) {
                names(datasets[[i]])[[label.col]] <- "label"
                datasets[[i]]$label <- as.factor(datasets[[i]]$label)
        }

        all_folds <- lapply(datasets, function(x) {
                set.seed(seed)
                folds <- caret::createFolds(x$label, k = folds.num, returnTrain = TRUE)
        })

        if (is.null(positive.class)) {
                positive.class <- as.character(datasets[[1]]$label[[1]])
        }

        ntree.range <- sort(unique(ntree.range))
        perf_tune <- data.frame()

        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        parallel.cores <- ifelse(parallel.cores > folds.num, folds.num, parallel.cores)

        cl <- parallel::makeCluster(parallel.cores)

        for (ntree in ntree.range) {
                message("- ntree - ", ntree, "   ", Sys.time())

                ntree_res <- Internal.randomForest_CV(datasets = datasets, label.col = label.col,
                                                      positive.class = positive.class, ntree = ntree,
                                                      folds.num = folds.num, all_folds = all_folds,
                                                      cl = cl, ...)
                ntree_perf <- t(ntree_res[folds.num + 1])
                row.names(ntree_perf) <- paste0("ntree_", ntree)
                perf_tune <- rbind(perf_tune, ntree_perf)
                # print(ntree_perf)
                print(round(ntree_perf, digits = 4)[,-c(1:4)])
        }
        parallel::stopCluster(cl)

        output <- ntree.range[which.max(perf_tune$Accuracy)]
        message("- Optimal ntree: ", output)

        if (return.model) {
                message("\n+ Training the model...   ", Sys.time())
                trainSet <- do.call("rbind", datasets)
                output <- randomForest::randomForest(label ~ ., data = trainSet,
                                                     ntree = output, ...)
        } else {
                output = list(ntree = output, performance = perf_tune)
        }
        output
}

Internal.checkNa <- function(dataset) {
        # dataset[is.na(dataset)] <- 0
        na_idx <- sapply(1:ncol(dataset), function(x, raw_dataset) {
                idx <- all(is.na(raw_dataset[,x]))
        }, raw_dataset = dataset)
        dataset[,which(na_idx)] <- 0
        dataset
}

# usethis::use_data(m_1_1, m_1_2, m_2_1, m_2_2, m_3_1, m_3_2, m_4_1, m_4_2, m_5_1,
#                   m_5_2, w_1, w_2, w_3, w_4, w_5, internal = T,
#                   overwrite = T, version = 2)

# save(mod_lncPro, file = "mod_lncPro.RData", compress = "xz")
# save(mod_ncProR, file = "mod_ncProR.RData", compress = "xz")
# save(mod_RPISeq, file = "mod_RPISeq.RData", compress = "xz")
# save(mod_rpiCOOL, file = "mod_rpiCOOL.RData", compress = "xz")

# devtools::check(document = T, cleanup = FALSE, manual = T, cran = F, check_dir = "E:\\ncProR")
# devtools::document(roclets = c('rd', 'collate', 'namespace'))
# devtools::build_manual()
#
# if (!library("devtools", logical.return = T)) install.packages("devtools")
# devtools::install_github("HAN-Siyu/ncProR")
