
geneset_score <- function(...) UseMethod("geneset_score")

geneset_score.Seurat <- function(
	object, 
	geneset, 
	method = 'average', 
	slot = 'data', 
	highly_variable = TRUE,
	gene_upper = TRUE,
	seurat_nbin = 24,
	seurat_control = 24,
	... # other args
) {
	slot <- match.arg(slot, c('counts', 'data', 'scale.data'))
	method <- match.arg(method, c('average', 'GSVA', 'Seurat'))

	if (method == 'Seurat') {
		Seurat_available <- require(Seurat)
		if (!Seurat_available)
			stop("Seurat analysis requires packages Seurat!")

		# get gene set
	    if (is(geneset, "GeneSetCollection")) {
	        list_of_genes <- lapply(
	        	geneset@.Data, 
	        	function(x) x@geneIds
	        )
	        names(list_of_genes) <- lapply(
	        	geneset@.Data, 
	        	function(x) x@setName
	        )
	    } else if (is(geneset, "list")) {
	    	if (is.null(names(geneset)))
	    		names(geneset) <- paste('geneset', seq(geneset), sep = "-")
	    	list_of_genes <- geneset
	    } else {
	    	stop("Input geneset should be the class of list or GeneSetCollection!")
	    }

	    if (gene_upper) {
	    	message("Use upper gene name.")
	    	exp_Names <- toupper(rownames(object))
	    	list_of_overlap_genes <- lapply(
	    		list_of_genes, 
	    		function(x) rownames(object)[exp_Names %in% x]
	    	)
	    } else {
	    	list_of_overlap_genes <- lapply(
	    		list_of_genes, 
	    		function(x) intersect(x, rownames(object))
	    	)
	    }
	    print("start")
	    for (i in names(list_of_overlap_genes)) {
	    	features <- list_of_overlap_genes[[i]]
	    	tmp_name <- 'XXX_Features'
	    	print(i)
	    	tryCatch(expr = {
    			object <- Seurat::AddModuleScore(
					object = object,
					features = features,
					nbin = seurat_nbin,
					ctrl = seurat_control,
					name = tmp_name,
					...
				)
				object@meta.data[[i]] <- object@meta.data[[paste0(tmp_name, "1")]]
				object@meta.data[[paste0(tmp_name, "1")]] <- NULL
			}, error = function(...) 
				warning("Problem may be caused by seurat_control set when ", i, " !")
			)
	    }
	} else {
	    assay <- object@active.assay
	    cat("Use data from default assay ", assay, "\n", sep = "")
	    exp_X <- methods::slot(
	    	object@assays[[assay]], slot
	    )
	    if (highly_variable) {
	    	var_features <- methods::slot(
	    		object@assays[[assay]], 'var.features'
	    	)
	    	if (length(var_features) == 0)
	    		stop("No highly variable genes in assay ", assay, " !")
	    	exp_X <- exp_X[var_features, , drop = F]
	    }
	    if (gene_upper) {
	    	message("Use upper gene name.")
	    	rownames(exp_X) <- toupper(rownames(exp_X))
	    }
	    
	    geneset_score_mtx <- geneset_score(
	    	exp_mtx = exp_X, 
	    	geneset = geneset, 
	    	method = method, 
	    	...
	    ) # row-samples, col-pathway
	    
	    object@meta.data <- cbind(object@meta.data, geneset_score_mtx)
	}


    return(object)
}

geneset_score.default <- function(
    exp_mtx, 
    geneset, 
    method = 'average', 
    gsva_method = NULL, # "gsva", "ssgsea", "zscore", "plage"
    gsva_kcdf = NULL,
    gsva_min.sz = 1,
    gsva_max.sz = Inf,
    gene_row = NULL, 
    sample_col = NULL,
    ... # GSVA::gsva args
) {
    # exp_mtx
    method <- match.arg(method, c('average', 'GSVA'))
    
    # get gene set
    if (is(geneset, "GeneSetCollection")) {
        list_of_genes <- lapply(
        	geneset@.Data, 
        	function(x) x@geneIds
        )
        names(list_of_genes) <- lapply(
        	geneset@.Data, 
        	function(x) x@setName
        )
    } else if (is(geneset, "list")) {
    	if (is.null(names(geneset)))
    		names(geneset) <- paste('geneset', seq(geneset), sep = "-")
    	list_of_genes <- geneset
    } else {
    	stop("Input geneset should be the class of list or GeneSetCollection")
    }
    set_of_genes <- unique(unlist(list_of_genes))
    
    # main computation
    if (method == 'average') {

    	# # Method by apply function to compute average value
        # pb <- utils::txtProgressBar(style = 3)
        # iii <- 0
        # set_of_genes_overlap <- intersect(set_of_genes, rownames(exp_mtx))
        # if (length(set_of_genes_overlap) <= 1)
        # 	stop("No enough genes in expression matrix! ",
        # 		"Case conversion may be required.")
        # exp_mtx <- exp_mtx[set_of_genes_overlap, , drop = F]
        # result <- sapply(list_of_genes, function(x) {
        #     gene_exist <- set_of_genes_overlap %in% x
        #     res <- apply(
        #         X = exp_mtx[gene_exist, , drop = F], 
        #         MARGIN = 2, 
        #         FUN = function(y) sum(y) / length(x)
        #     )
        #     iii <<- iii + 1
        #     utils::setTxtProgressBar(pb, iii / length(list_of_genes))
        #     res
        # }, USE.NAMES = T)
        # # row-samples, col-pathway
        # colnames(result) <- names(list_of_genes)
        # close(pb)

        # Method by matrix multiplication (Linear Algebra)
        list_of_overlap_genes <- lapply(
    		list_of_genes, 
    		function(x) intersect(x, rownames(exp_mtx))
    	)
    
    	geneset_mtx <- matrix(
    		0, nrow = nrow(exp_mtx), ncol = length(list_of_overlap_genes), 
    		dimnames = list(rownames(exp_mtx), names(list_of_overlap_genes))
    	)
    	for (i in colnames(geneset_mtx)) {
    		allgenes <- list_of_genes[[i]]
    		genes <- list_of_overlap_genes[[i]]
    		if (length(genes) == 0)
    			next
    		geneset_mtx[genes, i] <- 1 / length(allgenes)
    	}

    	result <- t(geneset_mtx) %*% exp_mtx 
    	# row-pathway, col-samples
    	result <- t(as.matrix(result))

    } else if (method == 'GSVA') {

    	GSVA_available <- require(GSVA) & require(GSEABase)
    	if (!GSVA_available)
    		stop("GSVA analysis requires packages GSVA and GSEABase!")
    	if (is.null(gsva_method))
    		gsva_method <- "gsva"
    	else
    		gsva_method <- match.arg(gsva_method, c("gsva", "ssgsea", "zscore", "plage"))
    	if (is.null(gsva_kcdf))
    		gsva_kcdf <- "Gaussian"
    	else
    		gsva_kcdf <- match.arg(gsva_kcdf, c("Gaussian", "Poisson"))
	
	exp_mtx_tmp <- if (gsva_method == 'ssgsea') exp_mtx else as.matrix(exp_mtx)
    	result <- GSVA::gsva(
		    exp_mtx_tmp, 
		    list_of_genes, 
		    min.sz = gsva_min.sz, 
		    max.sz = gsva_max.sz, 
		    verbose = TRUE, 
		    method = gsva_method, 
		    kcdf = gsva_kcdf,
		    ...
	    )
    	# row-pathway, col-samples
    	result <- t(result)

    }

    return(result) # row-samples, col-pathway
}

