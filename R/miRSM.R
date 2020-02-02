## Internal function cluster from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
cluster <- function(graph, method = "MCL", expansion = 2, inflation = 2,
    hcmethod = "average", directed = FALSE, outfile = NULL, ...) {

    method <- match.arg(method)
    if (method == "FN") {
        graph <- simplify(graph)
        fc <- fastgreedy.community(graph, merges = TRUE, modularity = TRUE)
        membership <- membership(fc)
        if (!is.null(V(graph)$name)) {
            names(membership) <- V(graph)$name
        }
        if (!is.null(outfile)) {
            cluster.save(cbind(names(membership), membership), 
                         outfile = outfile)
        } else {
            return(membership)
        }
    } else if (method == "LINKCOMM") {
        edgelist <- get.edgelist(graph)
        if (!is.null(E(graph)$weight)) {
            edgelist <- cbind(edgelist, E(graph)$weight)
        }
        lc <- getLinkCommunities(edgelist, plot = FALSE, directed = directed,
            hcmethod = hcmethod)
        if (!is.null(outfile)) {
            cluster.save(lc$nodeclusters, outfile = outfile)
        } else {
            return(lc$nodeclusters)
        }
    } else if (method == "MCL") {
        adj <- matrix(rep(0, length(V(graph))^2), nrow = length(V(graph)),
            ncol = length(V(graph)))
        for (i in seq_along(V(graph))) {
            neighbors <- neighbors(graph, v = V(graph)$name[i], mode = "all")
            j <- match(neighbors$name, V(graph)$name, nomatch = 0)
            adj[i, j] = 1
        }
        lc <- mcl(adj, addLoops = TRUE, expansion = expansion, 
            inflation = inflation, allow1 = TRUE, 
            max.iter = 100, ESM = FALSE)
        lc$name <- V(graph)$name
        lc$Cluster <- lc$Cluster

        if (!is.null(outfile)) {
            cluster.save(cbind(lc$name, lc$Cluster), outfile = outfile)
        } else {
            result <- lc$Cluster
            names(result) <- V(graph)$name
            return(result)
        }
    } else if (method == "MCODE") {
        compx <- mcode(graph, vwp = 0.9, haircut = TRUE, fluff = TRUE,
            fdt = 0.1)
        index <- which(!is.na(compx$score))
        membership <- rep(0, vcount(graph))
        for (i in seq_along(index)) {
            membership[compx$COMPLEX[[index[i]]]] <- i
        }
        if (!is.null(V(graph)$name))
            names(membership) <- V(graph)$name
        if (!is.null(outfile)) {
            cluster.save(cbind(names(membership), membership), 
                         outfile = outfile)
            invisible(NULL)
        } else {
            return(membership)
        }
    }
}

## Internal function cluster.save from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
#' @importFrom utils write.table
cluster.save <- function(membership, outfile) {
    wd <- dirname(outfile)
    wd <- ifelse(wd == ".", paste(wd, "/", sep = ""), wd)
    filename <- basename(outfile)
    if ((filename == "") || (grepl(":", filename))) {
        filename <- "membership.txt"
    } else if (grepl("\\.", filename)) {
        filename <- sub("\\.(?:.*)", ".txt", filename)
    }
    write.table(membership, file = paste(wd, filename, sep = "/"), 
                row.names = FALSE, col.names = c("node", "cluster"),
                quote = FALSE)
}

## Internal function mcode.vertex.weighting from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.vertex.weighting <- function(graph, neighbors) {
    stopifnot(is.igraph(graph))
    weight <- lapply(seq_len(vcount(graph)), function(i) {
        subg <- induced.subgraph(graph, neighbors[[i]])
        core <- graph.coreness(subg)
        k <- max(core)
        ### k-coreness
        kcore <- induced.subgraph(subg, which(core == k))
        if (vcount(kcore) > 1) {
            if (any(is.loop(kcore))) {
                k * ecount(kcore)/choose(vcount(kcore) + 1, 2)
            } else {
                k * ecount(kcore)/choose(vcount(kcore), 2)
            }
        } else {
            0
        }
    })

    return(unlist(weight))
}

## Internal function mcode.find.complex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.find.complex <- function(neighbors, neighbors.indx, vertex.weight,
    vwp, seed.vertex, seen) {

    res <- .C("complex", as.integer(neighbors), as.integer(neighbors.indx),
        as.single(vertex.weight), as.single(vwp), as.integer(seed.vertex),
        seen = as.integer(seen), COMPLEX = as.integer(rep(0, length(seen))),
        PACKAGE = "miRSM")

    return(list(seen = res$seen, COMPLEX = which(res$COMPLEX != 0)))
}

## Internal function mcode.find.complexex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.find.complexex <- function(graph, neighbors, vertex.weight, vwp) {

    seen <- rep(0, vcount(graph))

    neighbors <- lapply(neighbors, function(item) {
        item[-1]
    })
    neighbors.indx <- cumsum(unlist(lapply(neighbors, length)))

    neighbors.indx <- c(0, neighbors.indx)
    neighbors <- unlist(neighbors) - 1

    COMPLEX <- list()
    n <- 1
    w.order <- order(vertex.weight, decreasing = TRUE)
    for (i in w.order) {
        if (!(seen[i])) {
            res <- mcode.find.complex(neighbors, neighbors.indx, vertex.weight,
                vwp, i - 1, seen)
            if (length(res$COMPLEX) > 1) {
                COMPLEX[[n]] <- res$COMPLEX
                seen <- res$seen
                n <- n + 1
            }
        }
    }
    rm(neighbors)
    return(list(COMPLEX = COMPLEX, seen = seen))
}

## Internal function mcode.fluff.complex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.fluff.complex <- function(graph, vertex.weight, fdt = 0.8, complex.g,
    seen) {

    seq_complex.g <- seq_along(complex.g)
    for (i in seq_complex.g) {
        node.neighbor <- unlist(neighborhood(graph, 1, complex.g[i]))
        if (length(node.neighbor) > 1) {
            subg <- induced.subgraph(graph, node.neighbor)
            if (graph.density(subg, loops = FALSE) > fdt) {
                complex.g <- c(complex.g, node.neighbor)
            }
        }
    }

    return(unique(complex.g))
}

## Internal function mcode.post.process from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.post.process <- function(graph, vertex.weight, haircut, fluff, fdt = 0.8,
    set.complex.g, seen) {

    indx <- unlist(lapply(set.complex.g, function(complex.g) {
        if (length(complex.g) <= 2)
            0 else 1
    }))
    set.complex.g <- set.complex.g[indx != 0]
    set.complex.g <- lapply(set.complex.g, function(complex.g) {
        coreness <- graph.coreness(induced.subgraph(graph, complex.g))
        if (fluff) {
            complex.g <- mcode.fluff.complex(graph, vertex.weight, fdt,
                complex.g, seen)
            if (haircut) {
                ## coreness needs to be recalculated
                coreness <- graph.coreness(induced.subgraph(graph, complex.g))
                complex.g <- complex.g[coreness > 1]
            }
        } else if (haircut) {
            complex.g <- complex.g[coreness > 1]
        }
        return(complex.g)
    })
    set.complex.g <- set.complex.g[lapply(set.complex.g, length) > 2]
    return(set.complex.g)
}

## Internal function mcode from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode <- function(graph, vwp = 0.5, haircut = FALSE, fluff = FALSE, fdt = 0.8,
    loops = TRUE) {

    stopifnot(is.igraph(graph))
    if (vwp > 1 | vwp < 0) {
        stop("vwp must be between 0 and 1")
    }
    if (!loops) {
        graph <- simplify(graph, remove.multiple = FALSE, remove.loops = TRUE)
    }
    neighbors <- neighborhood(graph, 1)
    W <- mcode.vertex.weighting(graph, neighbors)
    res <- mcode.find.complexex(graph, neighbors = neighbors, vertex.weight = W,
        vwp = vwp)
    COMPLEX <- mcode.post.process(graph, vertex.weight = W, haircut = haircut,
        fluff = fluff, fdt = fdt, res$COMPLEX, res$seen)
    score <- unlist(lapply(COMPLEX, function(complex.g) {
        complex.g <- induced.subgraph(graph, complex.g)
        if (any(is.loop(complex.g)))
            score <- ecount(complex.g)/choose(vcount(complex.g) + 1, 2) *
                vcount(complex.g) 
        else score <- ecount(complex.g)/choose(vcount(complex.g), 2) *
            vcount(complex.g)
        return(score)
    }))
    order_score <- order(score, decreasing = TRUE)
    return(list(COMPLEX = COMPLEX[order_score], score = score[order_score]))
}

## Internal function CandModgenes for extracting candidate module genes
CandModgenes <- function(ceRExp, mRExp, Modulegenes, num.ModuleceRs = 2, 
    num.ModulemRs = 2){
  
    ceR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in%
        colnames(ceRExp))))
    mR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in%
        colnames(mRExp))))

    index <- which(ceR_Num >= num.ModuleceRs & mR_Num >= num.ModulemRs)
    CandidateModulegenes <- lapply(index, function(i) Modulegenes[[i]])
    CandidateModulegenes <- lapply(seq_along(index), function(i) GeneSet(CandidateModulegenes[[i]], 
        setName = paste("Module", i, sep=" ")))
    CandidateModulegenes <- GeneSetCollection(CandidateModulegenes)
    
    return(CandidateModulegenes)
}

#' Identification of co-expressed gene modules from matched ceRNA and mRNA
#' expression data using WGCNA package
#'
#' @title module_WGCNA
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp A SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param RsquaredCut Desired minimum scale free topology fitting index 
#' R^2 with interval [0 1].
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @importFrom WGCNA pickSoftThreshold
#' @importFrom WGCNA adjacency
#' @importFrom WGCNA TOMdist
#' @importFrom WGCNA standardColors
#' @importFrom flashClust flashClust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom GSEABase GeneSet
#' @importFrom GSEABase GeneSetCollection
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples 
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp[, seq_len(80)], 
#'     mRExp[, seq_len(80)])
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
#' @references Langfelder P, Horvath S. WGCNA: an R package for weighted 
#' correlation network analysis. BMC Bioinformatics. 2008, 9:559.#' 
module_WGCNA <- function(ceRExp, mRExp, RsquaredCut = 0.9, num.ModuleceRs = 2,
    num.ModulemRs = 2) {

    ExpData <- cbind(assay(ceRExp), assay(mRExp))

    Optimalpower <- pickSoftThreshold(ExpData, RsquaredCut = RsquaredCut)$powerEstimate
    adjacencymatrix <- adjacency(ExpData, power = Optimalpower)
    dissTOM <- TOMdist(adjacencymatrix)
    hierTOM <- flashClust(as.dist(dissTOM), method = "average")

    # The function cutreeDynamic colors each gene by the branches that
    # result from choosing a particular height cutoff.
    colorh <- cutreeDynamic(hierTOM, method = "tree") + 1
    StandColor <- c("grey", standardColors(n = NULL))
    colorh <- unlist(lapply(seq_len(length(colorh)), function(i) StandColor[colorh[i]]))
    colorlevels <- unique(colorh)
    colorlevels <- colorlevels[-which(colorlevels == "grey")]

    Modulegenes <- lapply(seq_len(length(colorlevels)), function(i) colnames(ExpData)[which(colorh ==
        colorlevels[i])])
    
    CandidateModulegenes <- CandModgenes(ceRExp, mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
        num.ModulemRs = num.ModulemRs)

    return(CandidateModulegenes)
}


#' Identification of gene modules from matched ceRNA and mRNA 
#' expression data using GFA package
#'
#' @title module_GFA
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp A SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param StrengthCut Desired minimum strength (absolute value of 
#' association with interval [0 1]) for each bicluster.
#' @param iter.max The total number of Gibbs sampling steps 
#' (default 1000).
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @importFrom GFA normalizeData
#' @importFrom GFA getDefaultOpts
#' @importFrom GFA gfa
#' @importFrom GSEABase GeneSet
#' @importFrom GSEABase GeneSetCollection
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_GFA <- module_GFA(ceRExp[seq_len(20), seq_len(15)],
#'     mRExp[seq_len(20), seq_len(15)], iter.max = 2600)
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
#' @references Bunte K, Lepp\'{a}aho E, Saarinen I, Kaski S. 
#' Sparse group factor analysis for biclustering of multiple data sources. Bioinformatics. 2016, 32(16):2457-63.
#' @references Lepp\'{a}aho E, Ammad-ud-din M, Kaski S. GFA: 
#' exploratory analysis of multiple data sources with group factor 
#' analysis. J Mach Learn Res. 2017, 18(39):1-5.
module_GFA <- function(ceRExp, mRExp, StrengthCut = 0.9, iter.max = 5000,
    num.ModuleceRs = 2, num.ModulemRs = 2) {

    ExpData <- list(assay(ceRExp), assay(mRExp))
    names(ExpData) = c("ceRNA expression", "mRNA expression")

    # Normalize the data - here we assume that every feature is equally
    # important
    norm <- normalizeData(ExpData, type = "scaleFeatures")

    # Get the model options to detect bicluster structure
    opts <- getDefaultOpts(bicluster = TRUE)

    # Check for sampling chain convergence
    opts$convergenceCheck <- TRUE
    opts$iter.max <- iter.max

    # Infer the model
    res <- gfa(norm$train, opts = opts)

    # Extract gene index of each bicluster, using stength cutoff (absolute
    # value of association)
    BCresnum <- lapply(seq_len(dim(res$W)[2]), function(i) which(abs(res$W[,
        i]) >= StrengthCut))

    # Extract genes of each bicluster
    Modulegenes <- lapply(seq_along(BCresnum), function(i) colnames(cbind(assay(ceRExp),
        assay(mRExp)))[BCresnum[[i]]])

    CandidateModulegenes <- CandModgenes(ceRExp, mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
        num.ModulemRs = num.ModulemRs)

    return(CandidateModulegenes)
}


#' Identification of gene modules from matched ceRNA and mRNA 
#' expression data using igraph package
#'
#' @title module_igraph
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp A SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param cor.method The method of calculating correlation selected, 
#' including 'pearson' (default), 'kendall', 'spearman'.
#' @param pos.p.value.cutoff The significant p-value cutoff of 
#' positive correlation.
#' @param cluster.method The clustering method selected in 
#' \pkg{igraph} package, including 'betweenness', 'greedy' (default), 
#' 'infomap', 'prop', 'eigen', 'louvain', 'walktrap'.
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @importFrom igraph graph_from_incidence_matrix
#' @importFrom igraph cluster_edge_betweenness
#' @importFrom igraph cluster_fast_greedy
#' @importFrom igraph cluster_infomap
#' @importFrom igraph cluster_label_prop
#' @importFrom igraph cluster_leading_eigen
#' @importFrom igraph cluster_louvain
#' @importFrom igraph cluster_walktrap
#' @importFrom GSEABase GeneSet
#' @importFrom GSEABase GeneSetCollection
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_igraph <- module_igraph(ceRExp[, seq_len(10)],
#'     mRExp[, seq_len(10)])
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
#' @references Csardi G, Nepusz T. The igraph software package for 
#' complex network research, InterJournal, Complex Systems. 2006:1695.
module_igraph <- function(ceRExp, mRExp, cor.method = "pearson", pos.p.value.cutoff = 0.01,
    cluster.method = "greedy", num.ModuleceRs = 2, num.ModulemRs = 2) {

    cor.binary <- cor_binary(ceRExp, mRExp, cor.method = cor.method, pos.p.value.cutoff = pos.p.value.cutoff)
    cor.binary.graph <- graph_from_incidence_matrix(cor.binary)

    if (cluster.method == "betweenness") {
        Modulegenes <- cluster_edge_betweenness(cor.binary.graph)
    } else if (cluster.method == "greedy") {
        Modulegenes <- cluster_fast_greedy(cor.binary.graph)
    } else if (cluster.method == "infomap") {
        Modulegenes <- cluster_infomap(cor.binary.graph)
    } else if (cluster.method == "prop") {
        Modulegenes <- cluster_label_prop(cor.binary.graph)
    } else if (cluster.method == "eigen") {
        Modulegenes <- cluster_leading_eigen(cor.binary.graph)
    } else if (cluster.method == "louvain") {
        Modulegenes <- cluster_louvain(cor.binary.graph)
    } else if (cluster.method == "walktrap") {
        Modulegenes <- cluster_walktrap(cor.binary.graph)
    }

    CandidateModulegenes <- CandModgenes(ceRExp, mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
        num.ModulemRs = num.ModulemRs)

    return(CandidateModulegenes)
}


#' Identification of gene modules from matched ceRNA and mRNA 
#' expression data using ProNet package
#'
#' @title module_ProNet
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp A SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param cor.method The method of calculating correlation selected, 
#' including 'pearson' (default), 'kendall', 'spearman'.
#' @param pos.p.value.cutoff The significant p-value cutoff of 
#' positive correlation
#' @param cluster.method The clustering method selected in 
#' \pkg{ProNet} package, including 'FN', 'MCL' (default), 
#' 'LINKCOMM', 'MCODE'.
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @import igraph
#' @importFrom Rcpp evalCpp
#' @importFrom MCL mcl
#' @importFrom linkcomm getLinkCommunities
#' @importFrom GSEABase GeneSet
#' @importFrom GSEABase GeneSetCollection
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_ProNet <- module_ProNet(ceRExp[, seq_len(10)],
#'     mRExp[, seq_len(10)])
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
#' @references Clauset A, Newman ME, Moore C. Finding community 
#' structure in very large networks. Phys Rev E Stat Nonlin Soft 
#' Matter Phys., 2004, 70(6 Pt 2):066111.
#' @references Enright AJ, Van Dongen S, Ouzounis CA. An efficient 
#' algorithm for large-scale detection of protein families. 
#' Nucleic Acids Res., 2002, 30(7):1575-84.
#' @references Kalinka AT, Tomancak P. linkcomm: an R package 
#' for the generation, visualization, and analysis of link 
#' communities in networks of arbitrary size and type. 
#' Bioinformatics, 2011, 27(14):2011-2.
#' @references Bader GD, Hogue CW. An automated method for 
#' finding molecular complexes in large protein interaction 
#' networks. BMC Bioinformatics, 2003, 4:2.
module_ProNet <- function(ceRExp, mRExp, cor.method = "pearson", pos.p.value.cutoff = 0.01,
    cluster.method = "MCL", num.ModuleceRs = 2, num.ModulemRs = 2) {

    cor.binary <- cor_binary(ceRExp, mRExp, cor.method = cor.method, pos.p.value.cutoff = pos.p.value.cutoff)
    cor.binary.graph <- graph_from_incidence_matrix(cor.binary)

    if (cluster.method == "FN" | cluster.method == "MCL") {
        network_Cluster <- cluster(cor.binary.graph, method = cluster.method)

        Modulegenes <- lapply(seq_len(max(network_Cluster)), function(i) rownames(as.matrix(network_Cluster))[which(network_Cluster ==
            i)])
    } else if (cluster.method == "LINKCOMM") {
        edgelist <- get.edgelist(cor.binary.graph)
        network_Cluster <- getLinkCommunities(edgelist)$nodeclusters
        Modulegenes <- lapply(seq_len(max(c(network_Cluster$cluster))),
            function(i) as.character(network_Cluster$node[which(c(network_Cluster$cluster) ==
                i)]))
    } else if (cluster.method == "MCODE") {
        network_Cluster <- cluster(cor.binary.graph, method = cluster.method) +
            1
        Modulegenes <- lapply(seq_len(max(network_Cluster)), function(i) rownames(as.matrix(network_Cluster))[which(network_Cluster ==
            i)])
    }

    CandidateModulegenes <- CandModgenes(ceRExp, mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
        num.ModulemRs = num.ModulemRs)

    return(CandidateModulegenes)
}


#' Identification of gene modules from matched ceRNA and mRNA 
#' expression data using NMF package
#'
#' @title module_NMF
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp A SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param NMF.algorithm Specification of the NMF algorithm, 
#' including 'brunet' (default), 'Frobenius', 'KL', 'lee', 'nsNMF', 
#' 'offset', 'siNMF', 'snmf/l', 'snmf/r'.
#' @param num.modules The number of modules to be identified.
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @importFrom NMF nmf
#' @importFrom NMF predict
#' @importFrom NMF nneg
#' @importFrom GSEABase GeneSet
#' @importFrom GSEABase GeneSetCollection
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples
#' data(BRCASampleData)
#' # Reimport NMF package to avoid conflicts with DelayedArray package
#' library(NMF)
#' modulegenes_NMF <- module_NMF(ceRExp[, seq_len(10)],
#'     mRExp[, seq_len(10)])
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
#' @references  Gaujoux R, Seoighe C. A flexible R package for 
#' nonnegative matrix factorization. BMC Bioinformatics. 2010, 11:367.
module_NMF <- function(ceRExp, mRExp, NMF.algorithm = "brunet", num.modules = 10,
    num.ModuleceRs = 2, num.ModulemRs = 2) {

    ExpData <- cbind(assay(ceRExp), assay(mRExp))

    # Run NMF algorithm with rank num.modules, negative values are transformed
    # into 0 if exist in expression data
    res <- nmf(nneg(ExpData), rank = num.modules, method = NMF.algorithm)

    # Predict column clusters
    Cluster.membership <- predict(res)

    # Extract genes of each cluster
    Modulegenes <- lapply(seq_len(num.modules), function(i) colnames(ExpData)[which(Cluster.membership ==
        i)])

    CandidateModulegenes <- CandModgenes(ceRExp, mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
        num.ModulemRs = num.ModulemRs)

    return(CandidateModulegenes)
}

#' Identification of gene modules from matched ceRNA and mRNA 
#' expression data using a series of clustering packages, 
#' including stats, flashClust, dbscan, subspace, mclust, SOMbrero and ppclust packages.
#' 
#' @title module_clust 
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp A SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param cluster.method Specification of the clustering method, 
#' including 'kmeans'(default), 'hclust', 'dbscan' , 'clique', 
#' 'gmm', 'som' and 'fcm'.
#' @param num.modules Parameter of the number of modules to be identified
#' for the 'kmeans', 'hclust', 'gmm' and 'fcm' methods. Parameter of the number
#' of intervals for the 'clique' method. For the 'dbscan' and 'som' methods,
#' no need to set the parameter.
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @importFrom stats kmeans
#' @importFrom stats dist
#' @importFrom stats cutree
#' @importFrom flashClust flashClust
#' @importFrom dbscan optics
#' @importFrom dbscan dbscan
#' @importFrom subspace CLIQUE
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom SOMbrero trainSOM
#' @importFrom ppclust fcm
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_clust <- module_clust(ceRExp[, seq_len(30)],
#'     mRExp[, seq_len(30)])
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
#' @references Forgy EW. Cluster analysis of multivariate 
#' data: efficiency vs interpretability of classifications. 
#' Biometrics, 1965, 21:768-769.
#' @references Hartigan JA, Wong MA. 
#' Algorithm AS 136: A K-means clustering algorithm. 
#' Applied Statistics, 1979, 28:100-108.
#' @references Lloyd SP. Least squares quantization in PCM. 
#' Technical Note, Bell Laboratories. Published in 1982 
#' in IEEE Transactions on Information Theory, 1982, 28:128-137.
#' @references MacQueen J. Some methods for classification 
#' and analysis of multivariate observations. 
#' In Proceedings of the Fifth Berkeley Symposium on 
#' Mathematical Statistics and Probability, 
#' eds L. M. Le Cam & J. Neyman, 1967, 1, pp.281-297. 
#' Berkeley, CA: University of California Press.
#' @references Langfelder P, Horvath S. Fast R Functions for 
#' Robust Correlations and Hierarchical Clustering. 
#' Journal of Statistical Software. 2012, 46(11):1-17.
#' @references Ester M, Kriegel HP, Sander J, Xu X. A density-based 
#' algorithm for discovering clusters in large spatial databases with 
#' noise, Proceedings of 2nd International Conference on Knowledge Discovery and
#' Data Mining (KDD-96), 1996, 96(34): 226-231.
#' @references Campello RJGB, Moulavi D, Sander J. 
#' Density-based clustering based on hierarchical density estimates,
#' Pacific-Asia conference on knowledge discovery and data mining. 
#' Springer, Berlin, Heidelberg, 2013: 160-172.
#' @references Agrawal R, Gehrke J, Gunopulos D, Raghavan P. 
#' Automatic subspace clustering of high dimensional data for 
#' data mining applications. In Proc. ACM SIGMOD, 1998.
#' @references Scrucca L, Fop M, Murphy TB, Raftery AE. 
#' mclust 5: clustering, classification and density estimation using 
#' Gaussian finite mixture models The R Journal 8/1, 2016, pp. 205-233.
#' @references Kohonen T. Self-Organizing Maps. 
#' Berlin/Heidelberg: Springer-Verlag, 3rd edition, 2001. 
#' @references Dunn JC. A fuzzy relative of the ISODATA process 
#' and its use in detecting compact well-separated clusters. Journal of Cybernetics, 
#' 1973, 3(3):32-57.
#' @references Bezdek JC. Cluster validity with fuzzy sets. Journal of Cybernetics, 1974, 3: 58-73.
#' @references Bezdek JC. Pattern recognition with fuzzy objective function 
#' algorithms. Plenum, NY, 1981. 
module_clust <- function(ceRExp, mRExp, cluster.method = "kmeans", num.modules = 10,
                           num.ModuleceRs = 2, num.ModulemRs = 2) {
  
  ExpData <- cbind(assay(ceRExp), assay(mRExp))
  
  if (cluster.method == "kmeans") {
    res <- kmeans(t(ExpData), centers = num.modules, iter.max = 100)
  } else if (cluster.method == "hclust") {
    diss <- dist(t(ExpData))
    hc <- flashClust(diss, method = "average")
    res <- cutree(hc, k = num.modules)
  } else if (cluster.method == "dbscan") {
    eps <- optics(t(ExpData))$eps
    res <- dbscan(t(ExpData), eps = eps)
  } else if (cluster.method == "clique") {
    res <- CLIQUE(t(ExpData), xi = num.modules)
  } else if (cluster.method == "gmm") {
    res <- Mclust(t(ExpData), G = num.modules)$classification
  } else if (cluster.method == "som") {
    res <- trainSOM(t(ExpData))$clustering
  } else if (cluster.method == "fcm") {
    res <- fcm(t(ExpData), centers = num.modules)$cluster
  } 
  
  # Extract genes of each cluster
  if (cluster.method == "kmeans") {
    Cluster.membership <- res$cluster
    Modulegenes <- lapply(seq_len(num.modules), function(i) 
      colnames(ExpData)[which(Cluster.membership == i)])
  }
  
  if (cluster.method == "hclust" | cluster.method == "gmm") {
    Cluster.membership <- res
    Modulegenes <- lapply(seq_len(num.modules), function(i) 
      names(res)[which(Cluster.membership == i)])
  }
  
  if (cluster.method == "dbscan" ) {
    Cluster.membership <- res$cluster
    Modulegenes <- lapply(seq_len(max(Cluster.membership)), function(i) 
      colnames(ExpData)[which(Cluster.membership == i)])
  }
  
  if (cluster.method == "clique") {
    Modulegenes <- lapply(seq(length(res)), function(i) 
      colnames(ExpData)[res[[i]]$objects])
  }
  
  if (cluster.method == "som" ) {
    Cluster.membership <- res
    Modulegenes <- lapply(seq_len(max(Cluster.membership)), function(i) 
      names(res)[which(Cluster.membership == i)])
  }
  
  if (cluster.method == "fcm" ) {
    Cluster.membership <- res
    Modulegenes <- lapply(seq_len(num.modules), function(i) 
      colnames(ExpData)[which(Cluster.membership == i)])
  }
                                                                                          
  
  CandidateModulegenes <- CandModgenes(ceRExp, mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
                                       num.ModulemRs = num.ModulemRs)
  
  return(CandidateModulegenes)
}


#' Identification of gene modules from matched ceRNA and mRNA 
#' expression data using a series of biclustering packages, 
#' including biclust, runibic, iBBiG, fabia, BicARE, isa2, s4vd, 
#' BiBitR and rqubic
#'
#' @title module_biclust
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp A SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param BCmethod Specification of the biclustering method, 
#' including 'BCBimax', 'BCCC', 'BCPlaid' (default), 'BCQuest', 
#' 'BCSpectral', 'BCXmotifs', 'BCUnibic', iBBiG', 'fabia', 'fabiap', 
#' 'fabias', 'mfsc', 'nmfdiv', 'nmfeu', 'nmfsc', 'FLOC', 'isa', 
#' 'BCs4vd', 'BCssvd', 'bibit' and 'quBicluster'.
#' @param num.modules The number of modules to be identified. For the 'BCPlaid',
#' 'BCSpectral', 'isa' and 'bibit' methods, no need to set the parameter. For the 
#' 'quBicluster' method, the parameter is used to set the number of biclusters 
#' that should be reported. 
#' @param num.ModuleceRs The minimum number of ceRNAs in each module.
#' @param num.ModulemRs The minimum number of mRNAs in each module.
#' @import SummarizedExperiment
#' @importFrom biclust biclust
#' @importFrom biclust binarize
#' @importFrom biclust discretize
#' @importFrom biclust BCBimax
#' @importFrom biclust BCCC
#' @importFrom biclust BCPlaid
#' @importFrom biclust BCQuest
#' @importFrom biclust BCSpectral
#' @importFrom biclust BCXmotifs
#' @importFrom biclust biclusternumber
#' @importFrom runibic BCUnibic
#' @importFrom iBBiG iBBiG
#' @importFrom fabia fabia
#' @importFrom fabia fabiap
#' @importFrom fabia fabias
#' @importFrom fabia mfsc
#' @importFrom fabia nmfdiv
#' @importFrom fabia nmfeu
#' @importFrom fabia nmfsc
#' @importFrom fabia extractBic
#' @importFrom BicARE FLOC
#' @importFrom BicARE bicluster
#' @importFrom isa2 isa
#' @importFrom isa2 isa.biclust
#' @importFrom s4vd BCs4vd
#' @importFrom s4vd BCssvd
#' @importFrom BiBitR bibit
#' @importFrom BiBitR MaxBC
#' @importFrom rqubic quBicluster
#' @importFrom rqubic quantileDiscretize
#' @importFrom rqubic generateSeeds
#' @importFrom Biobase ExpressionSet
#' @importFrom GSEABase GeneSet
#' @importFrom GSEABase GeneSetCollection
#' @export
#' @return GeneSetCollection object: a list of module genes.
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_biclust <- module_biclust(ceRExp[, seq_len(30)],
#'     mRExp[, seq_len(30)])
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
#' @references Preli\'{c} A, Bleuler S, Zimmermann P, Wille A, 
#' B\'{u}hlmann P, Gruissem W, Hennig L, Thiele L, Zitzler E. 
#' A systematic comparison and evaluation of biclustering methods 
#' for gene expression data. Bioinformatics. 2006, 22(9):1122-9.
#' @references Cheng Y, Church GM. Biclustering of expression data. 
#' Proc Int Conf Intell Syst Mol Biol. 2000, 8:93-103.
#' @references Turner H, Bailey T, Krzanowski W. Improved 
#' biclustering of microarray data demonstrated through systematic 
#' performance tests. Comput Stat Data Anal. 2003, 48(2): 235-254.
#' @references Murali TM, Kasif S. Extracting conserved gene 
#' expression motifs from gene expression data. 
#' Pac Symp Biocomput. 2003:77-88.
#' @references Kluger Y, Basri R, Chang JT, Gerstein M. 
#' Spectral biclustering of microarray data: coclustering genes 
#' and conditions. Genome Res. 2003, 13(4):703-16.
#' @references Wang Z, Li G, Robinson RW, Huang X. UniBic: Sequential 
#' row-based biclustering algorithm for analysis of gene expression 
#' data. Sci Rep. 2016, 6:23466.
#' @references Gusenleitner D, Howe EA, Bentink S, Quackenbush J, 
#' Culhane AC. iBBiG: iterative binary bi-clustering of gene sets. 
#' Bioinformatics. 2012, 28(19):2484-92.
#' @references Hochreiter S, Bodenhofer U, Heusel M, Mayr A, 
#' Mitterecker A, Kasim A, Khamiakova T, Van Sanden S, Lin D, 
#' Talloen W, Bijnens L, G\'{o}hlmann HW, Shkedy Z, Clevert DA. 
#' FABIA: factor analysis for bicluster acquisition. 
#' Bioinformatics. 2010, 26(12):1520-7.
#' @references Yang J, Wang H, Wang W, Yu, PS. An improved 
#' biclustering method for analyzing gene expression. 
#' Int J Artif Intell Tools. 2005, 14(5): 771-789.
#' @references Bergmann S, Ihmels J, Barkai N. Iterative 
#' signature algorithm for the analysis of large-scale gene 
#' expression data. Phys Rev E Stat Nonlin Soft Matter Phys. 
#' 2003, 67(3 Pt 1):031902.
#' @references Sill M, Kaiser S, Benner A, Kopp-Schneider A. 
#' Robust biclustering by sparse singular value decomposition 
#' incorporating stability selection. Bioinformatics. 2011, 
#' 27(15):2089-97.
#' @references Lee M, Shen H, Huang JZ, Marron JS. Biclustering 
#' via sparse singular value decomposition. Biometrics. 2010, 
#' 66(4):1087-95.
#' @references Rodriguez-Baena DS, Perez-Pulido AJ, Aguilar-Ruiz JS. 
#' A biclustering algorithm for extracting bit-patterns from 
#' binary datasets. Bioinformatics. 2011, 27(19):2738-45.
#' @references Li G, Ma Q, Tang H, Paterson AH, Xu Y. 
#' QUBIC: a qualitative biclustering algorithm for analyses of 
#' gene expression data. Nucleic Acids Res. 2009, 37(15):e101.
module_biclust <- function(ceRExp, mRExp, BCmethod = "fabia", num.modules = 10,
    num.ModuleceRs = 2, num.ModulemRs = 2) {

    ExpData <- cbind(assay(ceRExp), assay(mRExp))

    if (BCmethod == "BCBimax") {
        ExpData <- binarize(ExpData)
        BCres <- biclust(ExpData, method = BCBimax(), number = num.modules)
    } else if (BCmethod == "BCCC") {
        BCres <- biclust(ExpData, method = BCCC(), number = num.modules)
    } else if (BCmethod == "BCPlaid") {
        BCres <- biclust(ExpData, method = BCPlaid())
    } else if (BCmethod == "BCQuest") {
        BCres <- biclust(ExpData, method = BCQuest(), number = num.modules)
    } else if (BCmethod == "BCSpectral") {
        BCres <- biclust(ExpData, method = BCSpectral())
    } else if (BCmethod == "BCXmotifs") {
        ExpData <- discretize(ExpData)
        BCres <- biclust(ExpData, method = BCXmotifs(), number = num.modules)
    } else if (BCmethod == "BCUnibic") {
        BCres <- biclust(t(ExpData), method = BCUnibic(), nbic = num.modules)
    } else if (BCmethod == "iBBiG") {
        ExpData <- binarize(ExpData)
        BCres <- iBBiG(ExpData, nModules = num.modules)
    } else if (BCmethod == "fabia") {
        BCres <- fabia(t(ExpData), p = num.modules)
    } else if (BCmethod == "fabiap") {
        BCres <- fabiap(t(ExpData), p = num.modules)
    } else if (BCmethod == "fabias") {
        BCres <- fabias(t(ExpData), p = num.modules)
    } else if (BCmethod == "mfsc") {
        BCres <- mfsc(t(ExpData), p = num.modules)
    } else if (BCmethod == "nmfdiv") {
        BCres <- nmfdiv(t(ExpData), p = num.modules)
    } else if (BCmethod == "nmfeu") {
        BCres <- nmfeu(t(ExpData), p = num.modules)
    } else if (BCmethod == "nmfsc") {
        BCres <- nmfsc(t(ExpData), p = num.modules)
    } else if (BCmethod == "FLOC") {
        ExpData <- t(ExpData)
        BCres <- FLOC(ExpData, k = num.modules)
    } else if (BCmethod == "isa") {
        BCres <- isa(ExpData)
        BCres <- isa.biclust(BCres)
    } else if (BCmethod == "BCs4vd") {
        BCres <- biclust(ExpData, method = BCs4vd(), nbiclust = num.modules)
    } else if (BCmethod == "BCssvd") {
        BCres <- biclust(ExpData, method = BCssvd(), K = num.modules)
    } else if (BCmethod == "bibit") {
        ExpData <- binarize(ExpData)
        BCres <- bibit(ExpData)
    } else if (BCmethod == "quBicluster") {
        ExpDataSet <- ExpressionSet(assayData = ExpData)
        ExpData.discret <- quantileDiscretize(ExpDataSet)
        ExpData.seeds <- generateSeeds(ExpData.discret)
        BCres <- quBicluster(ExpData.seeds, ExpData.discret, report.no = num.modules)
    }

    # Extract genes of each bicluster
    if (BCmethod == "BCBimax" | BCmethod == "BCCC" | BCmethod == "BCPlaid" |
        BCmethod == "BCQuest" | BCmethod == "BCSpectral" | BCmethod ==
        "BCXmotifs" | BCmethod == "iBBiG" | BCmethod ==
        "isa" | BCmethod == "BCs4vd" | BCmethod == "BCssvd" | 
        BCmethod == "quBicluster") {
        BCresnum <- biclusternumber(BCres)
        Modulegenes <- lapply(seq_along(BCresnum), function(i) colnames(ExpData)
            [BCresnum[[i]]$Cols])
    }
    
    if (BCmethod == "BCUnibic") {
      BCresnum <- biclusternumber(BCres)
      Modulegenes <- lapply(seq_along(BCresnum), function(i) colnames(ExpData)
                            [BCresnum[[i]]$Rows])
    }
    
    if (BCmethod == "bibit") {
      BCresnum <- biclusternumber(BCres)
      BCresnum <- lapply( which( names(BCresnum) %in% 
          gsub("BC", "Bicluster", colnames(MaxBC(BCres, top = num.modules)$column)) ), 
          function(i) BCresnum[[i]])
      Modulegenes <- lapply(seq_along(BCresnum), function(i) colnames(ExpData)[BCresnum[[i]]$Cols])
    }

    if (BCmethod == "fabia" | BCmethod == "fabiap" | BCmethod == "fabias" |
        BCmethod == "mfsc" | BCmethod == "nmfdiv" | BCmethod == "nmfeu" |
        BCmethod == "nmfsc") {
        Modulegenes <- lapply(seq_len(num.modules), function(i) extractBic(BCres)$bic[i,
            ]$bixn)
    }

    if (BCmethod == "FLOC") {
        Modulegenes <- lapply(seq_len(num.modules), function(i) rownames(bicluster(BCres,
            i, graph = FALSE)))
    }

    CandidateModulegenes <- CandModgenes(ceRExp, mRExp, Modulegenes, num.ModuleceRs = num.ModuleceRs, 
        num.ModulemRs = num.ModulemRs)

    return(CandidateModulegenes)
}


#' Generation of positively correlated binary matrix between 
#' ceRNAs and mRNAs
#'
#' @title cor_binary
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp A SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param cor.method The method of calculating correlation selected, 
#' including 'pearson' (default), 'kendall', 'spearman'.
#' @param pos.p.value.cutoff The significant p-value cutoff of 
#' positive correlation.
#' @import SummarizedExperiment
#' @importFrom WGCNA cor
#' @importFrom WGCNA corPvalueFisher
#' @export
#' @return A binary matrix.
#'
#' @examples
#' data(BRCASampleData)
#' cor_binary_matrix <- cor_binary(ceRExp, mRExp)
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
#' @references Langfelder P, Horvath S. WGCNA: an R package for 
#' weighted correlation network analysis. BMC Bioinformatics. 
#' 2008, 9:559.
cor_binary <- function(ceRExp, mRExp, cor.method = "pearson", 
                       pos.p.value.cutoff = 0.01) {

    cor.r <- cor(assay(ceRExp), assay(mRExp), method = cor.method)
    cor.pvalue <- corPvalueFisher(cor.r, nSamples = dim(ceRExp)[1])

    index1 <- which(cor.r > 0)
    index2 <- c(which(cor.r <= 0), which(cor.r %in% NA))
    index3 <- which(cor.pvalue < pos.p.value.cutoff)
    index4 <- c(which(cor.pvalue >= pos.p.value.cutoff), which(cor.pvalue %in%
        NA))

    cor.r[index1] <- 1
    cor.r[index2] <- 0
    cor.pvalue[index3] <- 1
    cor.pvalue[index4] <- 0

    cor.binary <- cor.r * cor.pvalue

    return(cor.binary)
}

## Identify miRNA sponge modules using sensitivity canonical correlation (SCC) method
miRSM_SCC <- function(miRExp, ceRExp, mRExp, miRTarget, CandidateModulegenes,  
                      typex = "standard", typez = "standard", nperms = 100, num_shared_miRNAs = 3,
                      pvalue.cutoff = 0.05, CC.cutoff = 0.8, SCC.cutoff = 0.3) {
  
  miRNames <- colnames(miRExp)
  ceRNames <- colnames(ceRExp)
  mRNames <- colnames(mRExp)
  CandidateModulegenes <- geneIds(CandidateModulegenes)
  
  miRTarget <- assay(miRTarget)
  miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in%
                                                    miRNames), which(miRTarget[, 2] %in% c(ceRNames, mRNames))), ]
  Res <- c()
  
  
  for (i in seq_along(CandidateModulegenes)) {
    # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
    tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                              intersect(CandidateModulegenes[[i]], ceRNames)), 1])
    M1 <- length(tmp1)
    tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                              intersect(CandidateModulegenes[[i]], mRNames)), 1])
    M2 <- length(tmp2)
    tmp3 <- intersect(tmp1, tmp2)
    M3 <- length(tmp3)
    M4 <- length(miRNames)
    M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
    
    if (M3 >= 3) {
      
      # Canonical correlation between a group of ceRNAs and a group of mRNAs
      perm.out_ceR_mR <- CCA.permute(assay(ceRExp)[, which(ceRNames %in%
                                                             CandidateModulegenes[[i]])], assay(mRExp)[, which(mRNames %in%
                                                                                                                 CandidateModulegenes[[i]])], typex = typex, typez = typez,
                                     nperms = nperms)
      out_ceR_mR <- CCA(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                        assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])],
                        typex = typex, typez = typez, penaltyx = perm.out_ceR_mR$bestpenaltyx,
                        penaltyz = perm.out_ceR_mR$bestpenaltyz, v = perm.out_ceR_mR$v.init)
      M6 <- out_ceR_mR$cor
      
      # Canonical correlation between a group of miRNAs and a group of mRNAs
      perm.out_miR_mR <- CCA.permute(assay(miRExp)[, which(miRNames %in%
                                                             tmp3)], assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])],
                                     typex = typex, typez = typez, nperms = nperms)
      out_miR_mR <- CCA(assay(miRExp)[, which(miRNames %in% tmp3)],
                        assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])],
                        typex = typex, typez = typez, penaltyx = perm.out_miR_mR$bestpenaltyx,
                        penaltyz = perm.out_miR_mR$bestpenaltyz, v = perm.out_miR_mR$v.init)
      M7 <- out_miR_mR$cor
      
      # Canonical correlation between a group of miRNAs and a group of
      # ceRNAs
      perm.out_miR_ceR <- CCA.permute(assay(miRExp)[, which(miRNames %in%
                                                              tmp3)], assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                                      typex = typex, typez = typez, nperms = nperms)
      out_miR_ceR <- CCA(assay(miRExp)[, which(miRNames %in% tmp3)],
                         assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                         typex = typex, typez = typez, penaltyx = perm.out_miR_ceR$bestpenaltyx,
                         penaltyz = perm.out_miR_ceR$bestpenaltyz, v = perm.out_miR_ceR$v.init)
      M8 <- out_miR_ceR$cor
      
      # Calculate partial canonical correlation between a group of ceRNAs
      # and a group of mRNAs on condition a group of miRNAs
      M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
      
      # Calculate sensitivity canonical correlation between a group of
      # ceRNAs and a group of mRNAs on condition a group of miRNAs
      M10 <- M6 - M9
    } else {
      M6 <- NA
      M7 <- NA
      M8 <- NA
      M9 <- NA
      M10 <- NA
    }
    
    tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)
    Res <- rbind(Res, tmp)
    
  }
  colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                     "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                     "Canonical correlation of ceRNAs:mRNAs", "Canonical correlation of miRNAs:mRNAs",
                     "Canonical correlation of miRNAs:ceRNAs", "Partial canonical correlation of ceRNAs:mRNAs",
                     "Sensitivity canonical correlation of ceRNAs:mRNAs")
  index <- which(Res[, "#Shared miRNAs"] > num_shared_miRNAs &
                   Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                   Res[, "Canonical correlation of ceRNAs:mRNAs"] > CC.cutoff &
                   Res[, "Sensitivity canonical correlation of ceRNAs:mRNAs"] > SCC.cutoff)
  
  if (length(index) == 0) {
    Result <- "No miRNA sponge modules identified"
  } else {
    miRSM_genes <- lapply(index, function(i) CandidateModulegenes[[i]])
    names(miRSM_genes) <- paste("miRSM", seq_along(index), sep=" ")
    Res <- Res[index, ]
    if (length(index) > 1) {
      rownames(Res) <- paste("miRSM", seq_along(index), sep = " ")
    }
    Result <- list(Res, miRSM_genes)
    names(Result) <- c("Group competition of miRNA sponge modules", "miRNA sponge modules")
  }
  return(Result)
}

## Identify miRNA sponge modules using sensitivity distance correlation (SDC) method 
miRSM_SDC <- function(miRExp, ceRExp, mRExp, miRTarget, CandidateModulegenes,  
                      num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                      DC.cutoff = 0.8, SDC.cutoff = 0.3) {    
  
  miRNames <- colnames(miRExp)
  ceRNames <- colnames(ceRExp)
  mRNames <- colnames(mRExp)
  CandidateModulegenes <- geneIds(CandidateModulegenes)
  
  miRTarget <- assay(miRTarget)
  miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in%
                                                    miRNames), which(miRTarget[, 2] %in% c(ceRNames, mRNames))), ]
  Res <- c()
  
  
  for (i in seq_along(CandidateModulegenes)) {
    # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
    tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                              intersect(CandidateModulegenes[[i]], ceRNames)), 1])
    M1 <- length(tmp1)
    tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                              intersect(CandidateModulegenes[[i]], mRNames)), 1])
    M2 <- length(tmp2)
    tmp3 <- intersect(tmp1, tmp2)
    M3 <- length(tmp3)
    M4 <- length(miRNames)
    M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
    
    if (M3 >= 3) {        
      
      # Calculate distance correlation between a group of ceRNAs
      # and a group of mRNAs
      M6 <- dcor(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                 assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
      
      # Calculate partial distance correlation between a group of ceRNAs
      # and a group of mRNAs on condition a group of miRNAs
      M7 <- abs(pdcor(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                      assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])],
                      assay(miRExp)[, which(miRNames %in% tmp3)]))
      
      # Calculate sensitivity distance correlation between a group of
      # ceRNAs and a group of mRNAs on condition a group of miRNAs
      M8 <- M6 - M7
      
    } else {
      M6 <- NA
      M7 <- NA
      M8 <- NA
    }
    
    tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8)
    Res <- rbind(Res, tmp)
    
  }
  colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                     "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                     "Distance correlation of ceRNAs:mRNAs", "Partial distance correlation of ceRNAs:mRNAs",
                     "Sensitivity distance correlation of ceRNAs:mRNAs")
  index <- which(Res[, "#Shared miRNAs"] > num_shared_miRNAs &
                   Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                   Res[, "Distance correlation of ceRNAs:mRNAs"] > DC.cutoff &
                   Res[, "Sensitivity distance correlation of ceRNAs:mRNAs"] > SDC.cutoff)
  
  if (length(index) == 0) {
    Result <- "No miRNA sponge modules identified"
  } else {
    miRSM_genes <- lapply(index, function(i) CandidateModulegenes[[i]])
    names(miRSM_genes) <- paste("miRSM", seq_along(index), sep=" ")
    Res <- Res[index, ]
    if (length(index) > 1) {
      rownames(Res) <- paste("miRSM", seq_along(index), sep = " ")
    }
    Result <- list(Res, miRSM_genes)
    names(Result) <- c("Group competition of miRNA sponge modules", "miRNA sponge modules")
  }
  return(Result)
}

## Identify miRNA sponge modules using sensitivity RV coefficient (SRVC) method
miRSM_SRVC <- function(miRExp, ceRExp, mRExp, miRTarget, CandidateModulegenes,  
                       num_shared_miRNAs = 3, pvalue.cutoff = 0.05, RVC.cutoff = 0.8, 
                       SRVC.cutoff = 0.3, RV_method = "RV") {    
  
  miRNames <- colnames(miRExp)
  ceRNames <- colnames(ceRExp)
  mRNames <- colnames(mRExp)
  CandidateModulegenes <- geneIds(CandidateModulegenes)
  
  miRTarget <- assay(miRTarget)
  miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in%
                                                    miRNames), which(miRTarget[, 2] %in% c(ceRNames, mRNames))), ]
  Res <- c()
  
  
  for (i in seq_along(CandidateModulegenes)) {
    # Calculate significance of miRNAs shared by each ceRNAs:mRNAs
    tmp1 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                              intersect(CandidateModulegenes[[i]], ceRNames)), 1])
    M1 <- length(tmp1)
    tmp2 <- unique(miRTargetCandidate[which(miRTargetCandidate[, 2] %in%
                                              intersect(CandidateModulegenes[[i]], mRNames)), 1])
    M2 <- length(tmp2)
    tmp3 <- intersect(tmp1, tmp2)
    M3 <- length(tmp3)
    M4 <- length(miRNames)
    M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
    
    if (M3 >= 3 & RV_method == "RV") {
      
      # RV coefficient between a group of ceRNAs and a group of mRNAs       
      M6 <- RV(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
               assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
      
      # RV coefficient between a group of miRNAs and a group of mRNAs        
      M7 <- RV(assay(miRExp)[, which(miRNames %in% tmp3)],
               assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
      
      # RV coefficient between a group of miRNAs and a group of ceRNAs         
      M8 <- RV(assay(miRExp)[, which(miRNames %in% tmp3)],
               assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])])
      
      # Calculate partial RV coefficient between a group of ceRNAs
      # and a group of mRNAs on condition a group of miRNAs
      M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
      
      # Calculate sensitivity RV coefficient between a group of
      # ceRNAs and a group of mRNAs on condition a group of miRNAs
      M10 <- M6 - M9
      
    } else if (M3 >= 3 & RV_method == "RV2") {
      
      # RV coefficient between a group of ceRNAs and a group of mRNAs       
      M6 <- RV2(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
      
      # RV coefficient between a group of miRNAs and a group of mRNAs        
      M7 <- RV2(assay(miRExp)[, which(miRNames %in% tmp3)],
                assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
      
      # RV coefficient between a group of miRNAs and a group of ceRNAs         
      M8 <- RV2(assay(miRExp)[, which(miRNames %in% tmp3)],
                assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])])
      
      # Calculate partial RV coefficient between a group of ceRNAs
      # and a group of mRNAs on condition a group of miRNAs
      M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
      
      # Calculate sensitivity RV coefficient between a group of
      # ceRNAs and a group of mRNAs on condition a group of miRNAs
      M10 <- M6 - M9
      
    } else if (M3 >= 3 & RV_method == "RVadjMaye") {
      
      # RV coefficient between a group of ceRNAs and a group of mRNAs       
      M6 <- RVadjMaye(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                      assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
      
      # RV coefficient between a group of miRNAs and a group of mRNAs        
      M7 <- RVadjMaye(assay(miRExp)[, which(miRNames %in% tmp3)],
                      assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
      
      # RV coefficient between a group of miRNAs and a group of ceRNAs         
      M8 <- RVadjMaye(assay(miRExp)[, which(miRNames %in% tmp3)],
                      assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])])
      
      # Calculate partial RV coefficient between a group of ceRNAs
      # and a group of mRNAs on condition a group of miRNAs
      M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
      
      # Calculate sensitivity RV coefficient between a group of
      # ceRNAs and a group of mRNAs on condition a group of miRNAs
      M10 <- M6 - M9
      
    } else if (M3 >= 3 & RV_method == "RVadjGhaziri") {
      
      # RV coefficient between a group of ceRNAs and a group of mRNAs       
      M6 <- RVadjGhaziri(assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])],
                         assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
      
      # RV coefficient between a group of miRNAs and a group of mRNAs        
      M7 <- RVadjGhaziri(assay(miRExp)[, which(miRNames %in% tmp3)],
                         assay(mRExp)[, which(mRNames %in% CandidateModulegenes[[i]])])
      
      # RV coefficient between a group of miRNAs and a group of ceRNAs         
      M8 <- RVadjGhaziri(assay(miRExp)[, which(miRNames %in% tmp3)],
                         assay(ceRExp)[, which(ceRNames %in% CandidateModulegenes[[i]])])
      
      # Calculate partial RV coefficient between a group of ceRNAs
      # and a group of mRNAs on condition a group of miRNAs
      M9 <- (M6 - M7 * M8)/(sqrt(1 - M7^2) * sqrt(1 - M8^2))
      
      # Calculate sensitivity RV coefficient between a group of
      # ceRNAs and a group of mRNAs on condition a group of miRNAs
      M10 <- M6 - M9
      
    } else {
      M6 <- NA
      M7 <- NA
      M8 <- NA
      M9 <- NA
      M10 <- NA
    }
    
    tmp <- c(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)
    Res <- rbind(Res, tmp)
    
  }
  colnames(Res) <- c("#miRNAs regulating ceRNAs", "#miRNAs regulating mRNAs",
                     "#Shared miRNAs", "#Background miRNAs", "Sig. p.value of sharing miRNAs",
                     "RV coefficient of ceRNAs:mRNAs", "RV coefficient of miRNAs:mRNAs",
                     "RV coefficient of miRNAs:ceRNAs", "Partial RV coefficient of ceRNAs:mRNAs",
                     "Sensitivity RV coefficient of ceRNAs:mRNAs")
  index <- which(Res[, "#Shared miRNAs"] > num_shared_miRNAs &
                   Res[, "Sig. p.value of sharing miRNAs"] < pvalue.cutoff & 
                   Res[, "RV coefficient of ceRNAs:mRNAs"] > RVC.cutoff &
                   Res[, "Sensitivity RV coefficient of ceRNAs:mRNAs"] > SRVC.cutoff)
  
  if (length(index) == 0) {
    Result <- "No miRNA sponge modules identified"
  } else {
    miRSM_genes <- lapply(index, function(i) CandidateModulegenes[[i]])
    names(miRSM_genes) <- paste("miRSM", seq_along(index), sep=" ")
    Res <- Res[index, ]
    if (length(index) > 1) {
      rownames(Res) <- paste("miRSM", seq_along(index), sep = " ")
    }
    Result <- list(Res, miRSM_genes)
    names(Result) <- c("Group competition of miRNA sponge modules", "miRNA sponge modules")
  }
  return(Result)
}


#' Identify miRNA sponge modules using sensitivity canonical correlation (SCC), sensitivity distance correlation (SDC),
#' and sensitivity RV coefficient (SRVC) methods.
#'
#' @title miRSM
#' @param miRExp A SummarizedExperiment object. miRNA expression data: 
#' rows are samples and columns are miRNAs.
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp A SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param miRTarget A SummarizedExperiment object. Putative 
#' miRNA-target binding information.
#' @param CandidateModulegenes List object: a list of candidate 
#' miRNA sponge modules.
#' @param typex The columns of x unordered (type='standard') or 
#' ordered (type='ordered'). Only for the SCC method.
#' @param typez The columns of z unordered (type='standard') or 
#' ordered (type='ordered'). Only for the SCC method.
#' @param nperms The number of permutations. Only for the SCC method.
#' @param method The method selected to identify miRNA sponge 
#' modules, including 'SCC', 'SDC' and 'SRVC'.
#' @param num_shared_miRNAs The number of common miRNAs shared 
#' by a group of ceRNAs and mRNAs.
#' @param pvalue.cutoff The p-value cutoff of significant sharing 
#' of common miRNAs by a group of ceRNAs and mRNAs.
#' @param MC.cutoff The cutoff of matrix correlation (canonical correlation, 
#' distance correlation and RV coefficient).
#' @param SMC.cutoff The cutoff of sensitivity matrix correlation
#' (sensitivity canonical correlation, sensitivity distance correlation 
#' and sensitivity RV coefficient).
#' @param RV_method the method of calculating RV coefficients. Select
#' one of 'RV', 'RV2', 'RVadjMaye' and 'RVadjGhaziri' methods.
#' Only for the SRVC method.
#' @import SummarizedExperiment
#' @importFrom PMA CCA.permute
#' @importFrom PMA CCA
#' @importFrom energy dcor
#' @importFrom energy pdcor
#' @importFrom MatrixCorrelation RV
#' @importFrom MatrixCorrelation RV2
#' @importFrom MatrixCorrelation RVadjMaye
#' @importFrom MatrixCorrelation RVadjGhaziri
#' @importFrom stats phyper
#' @importFrom GSEABase geneIds
#' @export
#' @return List object: Sensitivity correlation,
#' and genes of miRNA sponge modules.
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_igraph <- module_igraph(ceRExp[, seq_len(10)], 
#'     mRExp[, seq_len(10)])
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_igraph_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_igraph, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
#' @references Witten DM, Tibshirani R, Hastie T. A penalized matrix 
#' decomposition, with applications to sparse principal components 
#' and canonical correlation analysis. Biostatistics. 
#' 2009, 10(3):515-34.
#' @references Szekely GJ, Rizzo ML. Partial distance 
#' correlation with methods for dissimilarities. Annals of Statistics. 
#' 2014, 42(6):2382-2412.
#' @references Szekely GJ, Rizzo ML, Bakirov NK. 
#' Measuring and Testing Dependence by Correlation of Distances, 
#' Annals of Statistics, 2007, 35(6):2769-2794.
#' @references Robert P, Escoufier Y. A unifying tool for 
#' linear multivariate statistical methods: the RV-Coefficient. 
#' Applied Statistics, 1976, 25(3):257-265.
#' @references Smilde AK, Kiers HA, Bijlsma S, Rubingh CM, 
#' van Erk MJ. Matrix correlations for high-dimensional 
#' data: the modified RV-coefficient. Bioinformatics, 
#' 2009, 25(3):401-405.
#' @references Maye CD, Lorent J, Horgan GW. 
#' Exploratory analysis of multiple omics datasets using 
#' the adjusted RV coefficient". Stat Appl Genet Mol Biol., 
#' 2011, 10, 14.
#' @references EIGhaziri A, Qannari EM. Measures 
#' of association between two datasets; Application to sensory data, 
#' Food Quality and Preference, 2015, 40(A):116-124.
miRSM <- function(miRExp, ceRExp, mRExp, miRTarget, CandidateModulegenes,
    typex = "standard", typez = "standard", nperms = 100, method = c("SCC",
    "SDC", "SRVC"), num_shared_miRNAs = 3, pvalue.cutoff = 0.05, MC.cutoff = 0.8,
    SMC.cutoff = 0.3, RV_method = c("RV", "RV2", "RVadjMaye", "RVadjGhaziri")) {

    if (method == "SCC") {
        Res <- miRSM_SCC(miRExp, ceRExp, mRExp, miRTarget, CandidateModulegenes,
            typex = "standard", typez = "standard", nperms = nperms, num_shared_miRNAs = num_shared_miRNAs,
            pvalue.cutoff = pvalue.cutoff, CC.cutoff = MC.cutoff, 
            SCC.cutoff = SMC.cutoff)
    } else if (method == "SDC") {
        Res <- miRSM_SDC(miRExp, ceRExp, mRExp, miRTarget, CandidateModulegenes,
            num_shared_miRNAs = num_shared_miRNAs,
            pvalue.cutoff = pvalue.cutoff, DC.cutoff = MC.cutoff, 
            SDC.cutoff = SMC.cutoff)
    } else if (method == "SRVC") {
        Res <- miRSM_SRVC(miRExp, ceRExp, mRExp, miRTarget, CandidateModulegenes,
            num_shared_miRNAs = num_shared_miRNAs,
            pvalue.cutoff = pvalue.cutoff, RVC.cutoff = MC.cutoff, 
            SRVC.cutoff = SMC.cutoff, RV_method = RV_method)
    }

    return(Res)
}


#' Functional analysis of miRNA sponge modules, including functional 
#' enrichment and disease enrichment analysis
#'
#' @title module_FA
#' @param Modulelist List object: a list of miRNA sponge modules.
#' @param GOont One of 'MF', 'BP', and 'CC' subontologies.
#' @param Diseaseont One of 'DO', and 'DOLite' subontologies.
#' @param KEGGorganism Organism, supported organism listed 
#' in http://www.genome.jp/kegg/catalog/org_list.html.
#' @param Reactomeorganism Organism, one of 'human', 'rat', '
#' mouse', 'celegans', 'yeast', 'zebrafish', 'fly'.
#' @param OrgDb OrgDb
#' @param padjustvaluecutoff A cutoff value of adjusted p-values.
#' @param padjustedmethod Adjusted method of p-values, can select 
#' one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 
#' 'fdr', 'none'.
#' @param Analysis.type The type of functional analysis selected, 
#' including 'FEA' (functional enrichment analysis) and 'DEA' 
#' (disease enrichment analysis).
#' @importFrom miRspongeR moduleFEA
#' @importFrom miRspongeR moduleDEA
#' @export
#' @return List object: a list of enrichment analysis results.
#'
#' @examples
#' \dontrun{
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM_WGCNA_SRVC_FEA <- module_FA(miRSM_WGCNA_SRVC_genes, Analysis.type = 'FEA')
#' miRSM_WGCNA_SRVC_DEA <- module_FA(miRSM_WGCNA_SRVC_genes, Analysis.type = 'DEA')
#' }
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
#' @references Zhang J, Liu L, Xu T, Xie Y, Zhao C, Li J, Le TD (2019). 
#' “miRspongeR: an R/Bioconductor package for the identification and analysis of 
#' miRNA sponge interaction networks and modules.” BMC Bioinformatics, 20, 235.
#' @references Yu G, Wang L, Han Y, He Q (2012). 
#' “clusterProfiler: an R package for comparing biological themes among gene clusters.” 
#' OMICS: A Journal of Integrative Biology, 16(5), 284-287.
module_FA <- function(Modulelist, GOont = "BP", Diseaseont = "DO", KEGGorganism = "hsa",
    Reactomeorganism = "human", OrgDb = "org.Hs.eg.db", padjustvaluecutoff = 0.05,
    padjustedmethod = "BH", Analysis.type = c("FEA", "DEA")) {
    if (Analysis.type == "FEA") {
        Res <- moduleFEA(Modulelist, ont = GOont, KEGGorganism = KEGGorganism,
            Reactomeorganism = Reactomeorganism, OrgDb = OrgDb, padjustvaluecutoff = padjustvaluecutoff,
            padjustedmethod = padjustedmethod)
    } else if (Analysis.type == "DEA") {
        Res <- moduleDEA(Modulelist, OrgDb = OrgDb, ont = Diseaseont, padjustvaluecutoff = padjustvaluecutoff,
            padjustedmethod = padjustedmethod)
    }

    return(Res)
}

#' Cancer enrichment analysis of miRNA sponge modules using hypergeometric distribution test
#'
#' @title module_CEA
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp A SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs.
#' @param Cancergenes A SummarizedExperiment object: a list of cancer genes given.
#' @param Modulelist List object: a list of the identified miRNA sponge modules. 
#' @import SummarizedExperiment
#' @importFrom stats phyper
#' @export
#' @return Cancer enrichment significance p-values of the identified miRNA sponge modules
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM.CEA.pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, 
#'                               miRSM_WGCNA_SRVC_genes)
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
#' @references Johnson NL, Kotz S, Kemp AW (1992) 
#' "Univariate Discrete Distributions", Second Edition. New York: Wiley.
module_CEA <- function(ceRExp, mRExp, Cancergenes, Modulelist) {
  
  ExpData <- cbind(assay(ceRExp), assay(mRExp))      
  
  B <- ncol(ExpData)
  N <- length(intersect(colnames(ExpData), as.matrix(assay(Cancergenes))))
  M <- unlist(lapply(seq_along(Modulelist), function(i) length(Modulelist[[i]])))
  x <- unlist(lapply(seq_along(Modulelist), function(i) 
      length(intersect(Modulelist[[i]], as.matrix(assay(Cancergenes))))))    
  p.value <- 1 - phyper(x - 1, N, B - N, M)
  
  names(p.value) <- names(Modulelist)
  return(p.value)
}

#' Validation of miRNA sponge interactions in each miRNA sponge module
#'
#' @title module_Validate
#' @param Modulelist List object: a list of the identified miRNA sponge modules. 
#' @param Groundtruth Matrix object: a list of experimentally validated miRNA sponge interactions.
#' @export
#' @return List object: a list of validated miRNA sponge interactions in each miRNA sponge module
#'
#' @examples
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' library(miRspongeR)
#' Groundtruthcsv <- system.file("extdata", "Groundtruth.csv", package="miRspongeR")
#' Groundtruth <- read.csv(Groundtruthcsv, header=TRUE, sep=",") 
#' miRSM.Validate <- module_Validate(miRSM_WGCNA_SRVC_genes, Groundtruth)
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
module_Validate <- function(Modulelist, Groundtruth) {
  
  validate_res <- lapply(seq(Modulelist), function(i) 
      Groundtruth[intersect(which(as.matrix(Groundtruth[, 1]) %in% 
      Modulelist[[i]]), which(as.matrix(Groundtruth[, 2]) %in% Modulelist[[i]])), ])
  names(validate_res) <- names(Modulelist)
  
  return(validate_res)
}


#' Co-expression analysis of each miRNA sponge module and its corresponding random miRNA sponge module
#' 
#' @title module_Coexpress
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp A SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs. 
#' @param Modulelist List object: a list of the identified miRNA sponge modules. 
#' @param resample The number of random miRNA sponge modules generated, and 1000 times in default.
#' @param method The method used to evaluate the co-expression level of each miRNA sponge module.
#' Users can select "mean" or "median" to calculate co-expression value of each miRNA sponge module
#' and its corresponding random miRNA sponge module. 
#' @import SummarizedExperiment
#' @importFrom WGCNA cor
#' @importFrom stats median
#' @importFrom stats na.omit
#' @export
#' @return List object: co-expression values of miRNA sponge modules and their corresponding random miRNA sponge modules.
#'
#' @examples 
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM_WGCNA_Coexpress <-  module_Coexpress(ceRExp, mRExp, 
#'                                            miRSM_WGCNA_SRVC_genes, 
#'                                            resample = 10, method = "mean")
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
module_Coexpress <- function(ceRExp, mRExp, Modulelist, resample = 1000, method = c("mean", "median")) {
  
  ceRExp <- assay(ceRExp)
  mRExp <- assay(mRExp)
  module_ceRExp <- lapply(seq_along(Modulelist), function(i) 
    ceRExp[, which(colnames(ceRExp) %in% Modulelist[[i]])])
  module_mRExp <- lapply(seq_along(Modulelist), function(i) 
    mRExp[, which(colnames(mRExp) %in% Modulelist[[i]])])
  
  if (method == "mean"){
    module_avg_cor <- unlist(lapply(seq_along(Modulelist), function(i) 
      mean(abs(cor(module_ceRExp[[i]], module_mRExp[[i]])))))
  } else if (method == "median"){
    module_avg_cor <- unlist(lapply(seq_along(Modulelist), function(i) 
      median(abs(cor(module_ceRExp[[i]], module_mRExp[[i]])))))
  }
  
  module_avg_cor_resample <- c()
  for (i in seq_along(Modulelist)){
    temp1 <- replicate(resample, sample(seq_len(ncol(ceRExp)), size = ncol(module_ceRExp[[i]])))
    temp2 <- replicate(resample, sample(seq_len(ncol(mRExp)), size = ncol(module_mRExp[[i]])))
    module_ceRExp_resample <- lapply(seq_len(resample), function(i) ceRExp[, temp1[, i]])
    module_mRExp_resample <- lapply(seq_len(resample), function(i) mRExp[, temp2[, i]])
    
      if (method == "mean"){
          module_avg_cor_resample[i] <- mean(unlist(lapply(seq_len(resample), function(i) 
            mean(na.omit(abs(cor(module_ceRExp_resample[[i]], module_mRExp_resample[[i]])))))))
      } else if (method == "median"){
          module_avg_cor_resample[i] <- median(unlist(lapply(seq_len(resample), function(i) 
            median(na.omit(abs(cor(module_ceRExp_resample[[i]], module_mRExp_resample[[i]])))))))
      }
  }
  module_coexpress <- list(module_avg_cor, module_avg_cor_resample)
  names(module_coexpress) <- c("Real miRNA sponge modules", "Random miRNA sponge modules")
  return(module_coexpress)
}

#' Extract common miRNAs of each miRNA sponge module
#' 
#' @title share_miRs
#' @param miRExp A SummarizedExperiment object. miRNA expression data: 
#' rows are samples and columns are miRNAs.
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp A SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs. 
#' @param miRTarget A SummarizedExperiment object. Putative 
#' miRNA-target binding information.
#' @param Modulelist List object: a list of the identified miRNA sponge modules.
#' @import SummarizedExperiment 
#' @export
#' @return List object: a list of common miRNAs of each miRNA sponge module.
#'
#' @examples 
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM_WGCNA_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, 
#'                                       miRTarget, miRSM_WGCNA_SRVC_genes)
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
share_miRs <- function(miRExp, ceRExp, mRExp, miRTarget, Modulelist){
  
  miRExp <- assay(miRExp)
  ceRExp <- assay(ceRExp)
  mRExp <- assay(mRExp)
  miRTarget <- assay(miRTarget)
  miRTarget <- as.matrix(miRTarget)            
  miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% colnames(miRExp)),
                        which(miRTarget[, 2] %in% c(colnames(ceRExp),colnames(mRExp)))),]
  
  Res <- list()
  for (i in seq_along(Modulelist)){
    tmp1 <- unique(miRTargetCandidate[which( miRTargetCandidate[, 2] %in% 
                  intersect(Modulelist[[i]], colnames(ceRExp)) ), 1])
    tmp2 <- unique(miRTargetCandidate[which( miRTargetCandidate[, 2] %in% 
                  intersect(Modulelist[[i]], colnames(mRExp)) ), 1])
    tmp3 <- intersect( tmp1, tmp2 )
    Res[[i]] <- tmp3
  }
  
  names(Res) <- names(Modulelist)
  return(Res)
}

#' miRNA distribution analysis of sharing miRNAs by the identified miRNA sponge modules
#' 
#' @title module_miRdistribute
#' @param share_miRs List object: a list of common miRNAs of each miRNA sponge module 
#' generated by share_miRs function. 
#' @export
#' @return Matrix object: miRNA distribution in each miRNA sponge module.
#'
#' @examples 
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM_WGCNA_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, 
#'                                       miRTarget, miRSM_WGCNA_SRVC_genes)
#' miRSM_WGCNA_miRdistribute <- module_miRdistribute(miRSM_WGCNA_share_miRs)
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
module_miRdistribute <- function(share_miRs) {
  
  miRs <- unique(unlist(share_miRs))
  res <- NULL
  interin <- NULL
  for (i in seq_along(miRs)) {
    for (j in seq_along(share_miRs)) {
      if (length(which(miRs[i] %in% share_miRs[[j]]) == 1)) {
        interin <- c(interin, names(share_miRs)[j])
      }
    }
    res1 <- paste(interin, collapse = ", ")        
    res2 <- length(interin)
    res <- rbind(res, c(miRs[i], res1, res2))
    interin <- NULL
  }
  colnames(res) <- c("miRNA", "Module ID", "Number of modules")
    
  return(res)
}

#' Extract miRNA-target interactions of each miRNA sponge module
#' 
#' @title module_miRtarget
#' @param share_miRs List object: a list of common miRNAs of each miRNA sponge module 
#' generated by share_miRs function. 
#' @param Modulelist List object: a list of the identified miRNA sponge modules.
#' @export
#' @return List object: miRNA-target interactions of each miRNA sponge module.
#'
#' @examples 
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM_WGCNA_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, 
#'                                       miRTarget, miRSM_WGCNA_SRVC_genes)
#' miRSM_WGCNA_miRtarget <- module_miRtarget(miRSM_WGCNA_share_miRs, 
#'                                           miRSM_WGCNA_SRVC_genes)
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
module_miRtarget <- function(share_miRs, Modulelist){
  
  res_int <- list()
  for (k in seq(share_miRs)){
      CommonmiRs <- share_miRs[[k]]
      targets <- Modulelist[[k]]
      len_CommonmiRs <- length(CommonmiRs)
      len_targets <- length(targets)
      res_interin <- matrix(NA, len_CommonmiRs*len_targets, 2)
      for (i in seq_len(len_CommonmiRs)){
          for (j in seq_len(len_targets)){
              res_interin[(i-1)*len_targets+j, 1] <- CommonmiRs[i]
              res_interin[(i-1)*len_targets+j, 2] <- targets[j]
          }
      }
  res_int[[k]] <- res_interin
  }
  names(res_int) <- names(Modulelist)
  return(res_int)
}

#' Extract miRNA sponge interactions of each miRNA sponge module
#'
#' @title module_miRsponge
#' @param ceRExp A SummarizedExperiment object. ceRNA expression data: 
#' rows are samples and columns are ceRNAs.
#' @param mRExp A SummarizedExperiment object. mRNA expression data: 
#' rows are samples and columns are mRNAs. 
#' @param Modulelist List object: a list of the identified miRNA sponge modules.
#' @import SummarizedExperiment
#' @export
#' @return List object: miRNA sponge interactions of each miRNA sponge module.
#'
#' @examples 
#' data(BRCASampleData)
#' modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#' # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#' miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#'                         modulegenes_WGCNA, method = "SRVC",
#'                         SMC.cutoff = 0.01, RV_method = "RV")
#' miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#' miRSM_WGCNA_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, 
#'                                       miRTarget, miRSM_WGCNA_SRVC_genes)
#' miRSM_WGCNA_miRsponge <- module_miRsponge(ceRExp, mRExp, 
#'                                          miRSM_WGCNA_SRVC_genes)
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng_Zhang3})
module_miRsponge<- function(ceRExp, mRExp,  Modulelist){
  
    res_int <- list()
    for (k in seq(Modulelist)){
        Modulegenes <- Modulelist[[k]]
        ceRNAs <- Modulegenes[which(Modulegenes %in% colnames(ceRExp))]
        mRNAs <- Modulegenes[which(Modulegenes %in% colnames(mRExp))]
        len_ceRNAs <- length(ceRNAs)
        len_mRNAs <- length(mRNAs)
        res_interin <- matrix(NA, len_ceRNAs*len_mRNAs, 2)
        for (i in seq_len(len_ceRNAs)){
            for (j in seq_len(len_mRNAs)){
                res_interin[(i-1)*len_mRNAs+j, 1] <- ceRNAs[i]
                res_interin[(i-1)*len_mRNAs+j, 2] <- mRNAs[j]
            }
        }
    res_int[[k]] <- res_interin     
    } 
    names(res_int) <- names(Modulelist)
  return(res_int)
}
