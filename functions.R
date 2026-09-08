################################################################################
## functions.R — helper functions used by 01_analysis.R and 02_figures.R.
##
## Copied verbatim from the project's working functions.R; only the helpers the
## two published scripts actually call are kept (interactive/exploratory
## helpers have been dropped).
##
##   enrich2vector.V2       pairwise association (Fisher OR, AMI/ARI, Jaccard)
##                          between two categorical vectors  -> Figure 4C
##   computeAssociationVals single cell-of the above contingency table
##   fast_group_agg         sparse-aware group-wise row/column aggregation
##                          -> collapsing Ensembl IDs to gene symbols
##   dl_JAFAR_geneSyms_Human  human TF gene symbols from the JASPAR2024 SQLite
##   renameModulePerTF      rename gene modules after their top-contributing TF
##
## Requires: Matrix, aricode, JASPAR2024, RSQLite.
################################################################################

enrich2vector.V2  = function(a,b,varnames=c(deparse(substitute(a)),deparse(substitute(b)))){
    if(length(a) != length(b) ) stop("a and b must have the same size")
    
    levelsA<-as.factor(a) |> levels()
    levelsB<-as.factor(b) |> levels()
    
    binA<-sapply(levelsA, function(x) x == a); colnames(binA)<-levelsA
    binB<-sapply(levelsB, function(x) x == b); colnames(binB)<-levelsB
    
    allComb <- expand.grid(levelsA,levelsB)
    
    resList<-vector("list",length(levelsA)*length(levelsB))
    names(resList)<-apply(allComb,1,paste0,collapse = "_")
    
    for(lvlA in levelsA){
        for(lvlB in levelsB){
            vectA<-binA[,lvlA]
            vectB<-binB[,lvlB]
            resList[[paste0(lvlA,"_",lvlB)]]<-computeAssociationVals(vectA,vectB)
        }
    }
    res<-do.call(rbind,resList) |> data.frame()
    res$padj <- p.adjust(res$pval,method="BH")
    res[,varnames[1]] <- allComb[,1]
    res[,varnames[2]] <- allComb[,2]
    
    res
    
}

computeAssociationVals<-function(a,b){
    nA<-sum(a)
    nB<-sum(b)
    universeSize<-length(a)
    observedSize <- sum(a & b)
    expectedSize <- (nA * nB)/universeSize
    unionSize <- sum(a | b)
    fisherRes <- fisher.test(matrix(c(observedSize, nA - observedSize, 
                                      nB - observedSize, universeSize - unionSize), 
                                    nrow = 2))
    # pval <- min(phyper(q = observedSize - 0.5, m = nA,
    #            n = universeSize - nA, k = nB,
    #            lower.tail = TRUE),
    #     phyper(q = observedSize - 0.5, m = nA,
    #            n = universeSize - nA, k = nB,
    #            lower.tail = FALSE))*2
    # 
    
    c(
        "observed" = observedSize,
        "expected" = expectedSize,
        "nA" = nA,
        "nB" = nB,
        "margin" = unionSize,
        "universe" = universeSize,
        "OEdeviation" = (observedSize-expectedSize)/sqrt(universeSize),
        "AMI" =   aricode::AMI(a,b),
        "ARI" =  aricode::ARI(a,b),
        "Jacard" = observedSize/unionSize,
        "pval" = fisherRes$p.value,
        "OR" = as.vector(fisherRes$estimate)
    )
    
}


fast_group_agg <- function(x, by, FUN = mean, byRow = NULL) {
    stopifnot(is.matrix(x) || inherits(x, "Matrix"))
    nr <- nrow(x); nc <- ncol(x)
    
    if (is.null(byRow)) {
        if (length(by) == nr) byRow <- TRUE
        else if (length(by) == nc) byRow <- FALSE
        else stop("`by` length must match nrow(x) or ncol(x).")
    }
    
    # Prepare grouping as an integer factor (cheap & consistent)
    f <- factor(by, exclude = NULL)                # keep NAs separate if present
    g <- as.integer(f)
    labs <- levels(f)
    k <- nlevels(f)
    cnt <- tabulate(g, nbins = k)
    
    is_sparse <- inherits(x, "dgCMatrix")
    
    if (byRow) {
        if (is_sparse && (identical(FUN, mean) || identical(FUN, base::mean) || identical(FUN, "mean") ||
                          identical(FUN, sum)  || identical(FUN, base::sum)  || identical(FUN, "sum"))) {
            # Sparse row-group sum via G %*% x
            G <- Matrix::sparseMatrix(i = g, j = seq_len(nr), x = 1L, dims = c(k, nr))
            out <- G %*% x                              # (k x nr) %*% (nr x nc) -> (k x nc)
            if (identical(FUN, mean) || identical(FUN, base::mean) || identical(FUN, "mean")) {
                out <- Matrix::Diagonal(x = 1 / pmax(cnt, 1L)) %*% out
            }
            rownames(out) <- labs
            return(as.matrix(out))
        } else if (!is_sparse && (identical(FUN, mean) || identical(FUN, base::mean) || identical(FUN, "mean"))) {
            out <- rowsum(x, group = g, reorder = FALSE)  # dense, fast C code
            out <- out / pmax(cnt, 1L)
            rownames(out) <- labs
            return(as.matrix(out))
        } else if (!is_sparse && (identical(FUN, sum) || identical(FUN, base::sum) || identical(FUN, "sum"))) {
            out <- rowsum(x, group = g, reorder = FALSE)
            rownames(out) <- labs
            return(as.matrix(out))
        } else {
            # Generic FUN (slower): apply per group without copying whole matrix
            idx <- split(seq_len(nr), g, drop = TRUE)
            out <- vapply(idx, function(ix) apply(x[ix, , drop = FALSE], 2, FUN), numeric(nc))
            out <- t(out)
            rownames(out) <- labs
            return(as.matrix(out))
        }
    } else { # group columns
        if (is_sparse && (identical(FUN, mean) || identical(FUN, base::mean) || identical(FUN, "mean") ||
                          identical(FUN, sum)  || identical(FUN, base::sum)  || identical(FUN, "sum"))) {
            # Sparse col-group sum via x %*% t(G)
            Gt <- Matrix::sparseMatrix(i = seq_len(nc), j = g, x = 1L, dims = c(nc, k)) # this is t(G)
            out <- x %*% Gt                           # (nr x nc) %*% (nc x k) -> (nr x k)
            if (identical(FUN, mean) || identical(FUN, base::mean) || identical(FUN, "mean")) {
                out <- Matrix::t(Matrix::t(out) / pmax(cnt, 1L))  # sweep-like without dense copy
            }
            colnames(out) <- labs
            return(as.matrix(out))
        } else if (!is_sparse && (identical(FUN, mean) || identical(FUN, base::mean) || identical(FUN, "mean"))) {
            # Dense columns: avoid data.frames; one transpose (fast BLAS) is usually OK
            out <- t(rowsum(t(x), g, reorder = FALSE))   # sum per column-group
            out <- sweep(out, 2, pmax(cnt, 1L), "/")
            colnames(out) <- labs
            return(as.matrix(out))
        } else if (!is_sparse && (identical(FUN, sum) || identical(FUN, base::sum) || identical(FUN, "sum"))) {
            out <- t(rowsum(t(x), g, reorder = FALSE))
            colnames(out) <- labs
            return(as.matrix(out))
        } else {
            idx <- split(seq_len(nc), g, drop = TRUE)
            out <- vapply(idx, function(ix) apply(x[, ix, drop = FALSE], 1, FUN), numeric(nr))
            colnames(out) <- labs
            return(as.matrix(out))
        }
    }
}


dl_JAFAR_geneSyms_Human <- function(){
    require(JASPAR2024)
    JASPAR2024<-JASPAR2024()
    JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))
    geneNames <- RSQLite::dbGetQuery(JASPARConnect, 
                                     'SELECT DISTINCT m.name FROM matrix AS m
                        JOIN matrix_species AS ms ON ms.id = m.id
                        WHERE ms.tax_id = 9606;')[,"NAME"]
    strsplit(geneNames,split = ":", fixed = T) |> unlist() |> unique()
}


renameModulePerTF<-function(geneAnnot,TFgynSyms,ModuleNameCol="Module",interestCol="Contribution"){
    TFgynSyms <- intersect(TFgynSyms, rownames(geneAnnot))
    if(!is.factor(geneAnnot[,ModuleNameCol])) geneAnnot[,ModuleNameCol] <- as.factor(geneAnnot[,ModuleNameCol])
    modules <- levels(geneAnnot[,ModuleNameCol])
    newNames<- c()
    for(module in modules){
        genes <- rownames(geneAnnot)[geneAnnot[,ModuleNameCol] == module]
        genes <- genes[!is.na(geneAnnot[genes,interestCol])]
        genesTF <- intersect(genes, TFgynSyms)
        if(length(genesTF)<1) genesTF <- genes
        best <- genesTF[which.max(geneAnnot[genesTF,interestCol])]
        newNames[module] <- paste0(best,".Mod")
    }
    geneAnnot$oldModuleName <- geneAnnot[,ModuleNameCol]
    levels(geneAnnot[,ModuleNameCol]) <- newNames
    
    data.frame(geneAnnot)
}


