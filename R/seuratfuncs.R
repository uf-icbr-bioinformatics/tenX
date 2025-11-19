# ariva, 11/4/2019

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(celldex)

# Utils

rs_saveplot <- function(plot, data, tag) {
    filename = paste(Project(data), "-", tag, ".png", sep="")
    print(paste("Plotting", filename))
    png(file=filename, width=800, height=600)
    print(plot)
    dev.off()
}

# Pipeline steps

rs_load10Xrun <- function(directory, project="myproject", mincells=3, minfeatures=200) {
    data = Read10X(data.dir = directory)
    sobj = CreateSeuratObject(counts = data, project=project, min.cells=mincells, min.features=minfeatures)
    return(sobj)
}

rs_loadCsv <- function(filename, project="myproject", mincells=3, minfeatures=200) {
    data = read.csv(filename, comment.char="#", row.names=1)
    sobj = CreateSeuratObject(counts = data, project=project, min.cells=mincells, min.features=minfeatures)
    return(sobj)
}

rs_preProcess <- function(data, mtpattern="^MT-", percentmt=5, mingenes=200, maxgenes=5000) {
    data[["percent.mt"]] = PercentageFeatureSet(data, pattern=mtpattern)

    saveplot(VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),
             data, "violin")

    data = subset(data, subset = nFeature_RNA > mingenes & nFeature_RNA < maxgenes & percent.mt < 5)

    plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    saveplot(CombinePlots(plots = list(plot1, plot2)), data, "scatter")

    return(data)
}

rs_normalize <- function(data, method="LogNormalize", factor=10000) {
    # Normalized values are stored in data[["RNA"]]@data
    
    data = NormalizeData(data, normalization.method = method, scale.factor = factor)
    return(data)
}

rs_selectFeatures <- function(data, nfeatures=2000, top=10) {

    data = FindVariableFeatures(data, selection.method = "vst", nfeatures = nfeatures)

    # Identify the 10 most highly variable genes
    topgenes = head(VariableFeatures(data), top)

    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(data)
    plot2 <- LabelPoints(plot = plot1, points = topgenes, repel = TRUE)
    saveplot(CombinePlots(plots = list(plot1, plot2)), data, "varfeat")
    return(data)
}

rs_doPCA <- function(data) {
    # Results of scaling are stored in data[["RNA"]]@scale.data
    all.genes = rownames(data)
    data = ScaleData(data, features = all.genes)
    data = RunPCA(data, features = VariableFeatures(object = data))
    saveplot(DimPlot(data, reduction="pca"), data, "pca")
    return(data)
}

rs_doClustering <- function(data, resolution=0.5, dimensions=10) {
    data <- FindNeighbors(data, dims = 1:dimensions)
    data <- FindClusters(data, resolution = resolution)

    data <- RunUMAP(data, dims = 1:dimensions)
    # note that you can set `label = TRUE` or use the LabelClusters function to help label
    # individual clusters
    saveplot(DimPlot(data, reduction = "umap"), data, "umap")

    data <- RunTSNE(data)
    saveplot(DimPlot(data, reduction = "tsne"), data, "tsne")

    return(data)
}

rs_findMarkers <- function(data, minpct=0.25, logfc=0.25) {
    data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = minpct, logfc.threshold = logfc)
    data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
    return(data.markers)
}

rs_main <- function(directory, project) {
    data = rs_load10Xrun(directory, project=project)
    data = rs_preProcess(data)
    data = rs_normalize(data)
    data = rs_selectFeatures(data)
    data = rs_doPCA(data)
    data = rs_doClustering(data)
    saveRDS(data, file = paste(Project(data), "rds", sep="."))
    return(data)
}

rs_genes <- function(data) {
    markers = rs_findMarkers(data)
    
    return(markers)
}

rs_writeMarkers <- function(data.markers, nclusters, ntop) {
    for (i in 0:(nclusters-1)) {
        m = subset(data.markers, subset = cluster == i)
        print(kable(m[1:ntop,], caption = paste("Top ", ntop, " marker genes for cluster ", i)))
        cat("\n\n")
    }
}

# Integration and differential analysis

make_seurat_object <- function(name, mtpatt="^mt:", maxmt=5, mingenes=200, maxgenes=5000, condition=FALSE, parse=FALSE, path=FALSE) {
    if ( condition == FALSE ) { condition = name }
    if (parse) {
        so = make_seurat_object_parse(name, path)
    } else {
        message(paste("Loading dataset", name))
        so = CreateSeuratObject(counts = Read10X_h5(paste(name, "/outs/filtered_feature_bc_matrix.h5", sep="")), project=name)
    }
    so[["percent.mt"]] = PercentageFeatureSet(so, pattern=mtpatt)
    genes.no.mt = rownames(so)[grep(mtpatt, rownames(so), invert=TRUE)]
    so = subset(so, subset = percent.mt < maxmt & nFeature_RNA > mingenes & nFeature_RNA < maxgenes, features = genes.no.mt) %>%
        NormalizeData() %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
        AddMetaData(col.name = "Sample", metadata = name) %>%
        AddMetaData(col.name = "Condition", metadata = condition)
    return(so)
}

make_seurat_object_parse <- function(name, path) {
    if (path == FALSE) { path = paste(name, "DGE_filtered", sep='/') }
    message(paste("Loading Parse dataset", name))
    mat = ReadParseBio(path)
    cell_meta = read.csv(paste0(path, "/cell_metadata.csv"), row.names=1)
    so = CreateSeuratObject(counts = mat, project=name, names.field=0, meta.data=cell_meta)
    return(so)
}


#
#
#integrated_analysis v 20240715
#
#
# integrated_analysis <- function(smp1, smp2, plots=TRUE, resolution=0.25, logfc=0.25, mtpatt="^mt:",
#                                 maxmt=5, outdir="",
#                                 mingenes=500, maxgenes=5000, condition1=FALSE, condition2=FALSE, parse=FALSE) {
#
integrated_analysis <- function(params, outdir="", classify.cells=FALSE) {
    #
    log = c()

    smp1<- params$sample1
    smp2<- params$sample2
    condition1<- params$condition1
    condition2<- params$condition2
    if (condition1 == FALSE) { condition1 = smp1 }
    if (condition2 == FALSE) { condition2 = smp2 }
    label = paste(condition1, "_", condition2, sep="")

    solist = c()
    for (smp in unlist(strsplit(smp1, ","))) {
        so = make_seurat_object(
          smp, mtpatt=params$mtpatt, maxmt=params$maxmt, mingenes=params$mingenes,
          maxgenes=params$maxgenes, condition=condition1, parse=params$parse)
        solist = c(solist, so)
        log = c(log, paste(smp, "original number of cells:", nrow(so@meta.data)))
    }
    for (smp in unlist(strsplit(smp2, ","))) {
        so = make_seurat_object(
          smp, mtpatt=params$mtpatt, maxmt=params$maxmt,
          mingenes=params$mingenes, maxgenes=params$maxgenes, condition=condition2,
          parse=params$parse)
        solist = c(solist, so)
        log = c(log, paste(smp, "original number of cells:", nrow(so@meta.data)))
    }
    features <- SelectIntegrationFeatures(object.list = solist)
    anchors <- FindIntegrationAnchors(object.list = solist, anchor.features = features)
    combined <- IntegrateData(anchorset = anchors)

    log = c(log, paste("Integrated number of cells:", nrow(combined@meta.data)))
    
    DefaultAssay(combined) <- "integrated"

    # Run the standard workflow for visualization and clustering
    combined <- ScaleData(combined, verbose = FALSE) %>%
        RunPCA(npcs = 30, verbose = FALSE) %>%
        RunUMAP(reduction = "pca", dims = 1:params$dimensions) %>%
        FindNeighbors(reduction = "pca", dims = 1:params$dimensions) %>%
        FindClusters(resolution = params$resolution)

    log = c(log, paste("Combined number of cells:", nrow(combined@meta.data)))
    
    cluster.names = sort(unique(Idents(combined)))
    #
    # Save marker genes
    combined.markers <- rs_findMarkers(combined, minpct=0.25, logfc=params$logfc)
    markersfile = paste("markers", "csv", sep=".")
    write.table(combined.markers, file=paste(outdir, markersfile, sep='/'),
                sep='\t', quote=FALSE, col.names=NA, row.names=TRUE)

    # Save all genes
    all.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
    all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
    markersfile = paste(label, "markers", "csv", sep=".")
    write.table(all.markers, file=markersfile, sep='\t', quote=FALSE, col.names=NA, row.names=TRUE)
    #
    #
    cluster_counts = t(table(combined$orig.ident, Idents(combined)))
    #
    # Convert to percentages
    for (row in 1:nrow(cluster_counts)) {
        cluster_counts[row,] = 100.0 * cluster_counts[row,]/sum(cluster_counts[row,])
    }
    write.table(
      cluster_counts, file=paste(label, "counts", "csv", sep="."),
      quote=FALSE, col.names=TRUE, row.names=TRUE)
                                        #
    #
    out_list<- list(
      seurat=combined, markers=combined.markers, markersfile=markersfile,
        clusters=cluster.names)
    #
    # saveRDS(object=out_list, file=paste(outdir, '/inan.list1.', smp1, '.', smp2, '.rds', sep = ''))
                                        #
    
    if (params$celldex != FALSE) {
        out_list<- run_singler(params, out_list, by.cell=classify.cells, wanted.assays=c("integrated"))
    }
    names(out_list)[names(out_list)=='seurat']<- 'combined'
    
    #
    ##ERRATA 20240716v2
    saveRDS(object = out_list, file = paste(
                                   outdir, '/inan.list.', smp1, '.', smp2, '.rds', sep = ''))

    writeLines(log, "integr.log")
    #
    #####End of ERRATA 20240716v2
    #
    return(out_list)
    #
    #return(list(combined=combined, markers=combined.markers, markersfile=markersfile, clusternames=cluster.names))
}








diff_analysis <- function(integ, smp1, smp2, mincells=50, outdir="") {
    filenames = c()
    clustids = c()
    markers = c()
    label = paste(outdir, smp1, "_", smp2, sep="")
    DefaultAssay(integ) = "RNA"
    clusters = sort(unique(Idents(integ)))
    saveRDS(c(clusters,integ), "res.rds")
    # Add label consisting of condition and cluster number
    integ$clust.sample = paste("clust", Idents(integ), "_", integ$Condition, sep="")
    integ$celltype = Idents(integ)
    Idents(integ) = "clust.sample"

    for (cl in clusters) {
        id1 = paste("clust", cl, "_", smp1, sep="")
        id2 = paste("clust", cl, "_", smp2, sep="")

        if ( sum(integ$clust.sample == id1) > mincells & sum(integ$clust.sample == id2) > mincells) {
            message("*** Processing cluster ", cl)

            # Write all fold changes
            foldchanges = FoldChange(integ, ident.1 = id1, ident.2 = id2)
            outfile = paste(label, paste("clust", cl, sep=""), "fc", "csv", sep=".")
            write.table(foldchanges, file=outfile, sep='\t', quote=FALSE, col.names=NA, row.names=TRUE)

            response = FindMarkers(integ, ident.1 = id1, ident.2 = id2)
            outfile = paste(label, paste("clust", cl, sep=""), "csv", sep=".")
            filenames = c(filenames, outfile)
            clustids = c(clustids, cl)
            markers = c(markers, response)
            write.table(response, file=outfile, sep='\t', quote=FALSE, col.names=NA, row.names=TRUE)

        }
    }
    return(list(filenames, markers, clustids))
}

# Put everything together

make_params <- function(sample) {
    params = list()
    params$sample = sample
    params$mingenes = 500
    params$maxgenes = 5000
    params$maxmt = 10
    params$mtpatt = "^MT-"
    params$nfeatures = 2000
    params$method = "LogNormalize"
    params$ntop = 12
    params$dimensions = 10
    params$resolution = 0.5
    params$minpct = 0.25
    params$logfc = 0.25
    params$celldex = "FALSE"
    params$parse = FALSE
    return(params)
}

seurat_main <- function(params, save.markers=TRUE, save.rds=TRUE, classify.cells=FALSE) {
    result = list()
    
    # Load dataset, no mt filtering for now.
    message("/// Starting Seurat analysis.")
    data = make_seurat_object(
      params$sample, maxmt=100, mtpatt=params$mtpatt,
      mingenes=0, maxgenes=100000, parse=params$parse)
    AddMetaData(data, col.name = "percent.mt", metadata = PercentageFeatureSet(
      data, pattern=params$mtpatt))
    result$cells.before = dim(data)[2]
    message(paste("// Cells before filtering:", result$cells.before))
    result$violin = VlnPlot(data, features = c(
      "nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

    # Remove mt genes, filter by mingenes and maxgenes
    data = subset(
      data,
      subset = nFeature_RNA > params$mingenes &
        nFeature_RNA < params$maxgenes & percent.mt < params$maxmt)
    result$cells.after = dim(data)[2]
    message(paste("// Cells after filtering:", result$cells.after))
    
    # Normalization and feature selection
    message("/// Performing normalization.")
    data = NormalizeData(data, normalization.method=params$method, scale.factor=10000) %>%
        FindVariableFeatures(selection.method="vst", nfeatures=params$nfeatures)
    result$top.genes = head(VariableFeatures(data), params$ntop)

    # Scale data and run PCA
    data = ScaleData(data, features = rownames(data)) %>%
        RunPCA(features=VariableFeatures(object=data))

    # Clustering
    message("/// Clustering.")
    data = FindNeighbors(data, dims = 1:params$dimensions)
    data = FindClusters(data, resolution=params$resolution)
    # UMAP
    data = RunUMAP(data, dims=1:params$dimensions)
    # TSNE
    data = RunTSNE(data)

    message("/// Identifying marker genes.")
    result$clusters = unique(Idents(data))
    result$nclusters = length(result$clusters)
    result$markers = rs_findMarkers(data, minpct=params$minpct, logfc=params$logfc)
    if (save.markers) {
        result$markersfile = paste("markers", paste(params$sample, "markers", "csv", sep="."), sep="/")
        write.table(result$markers, file=result$markersfile, sep='\t', quote=FALSE, col.names=NA, row.names=TRUE)
    }

    result$seurat = data
    
#####
    ###ERRATA01
    saveRDS(result, file="results.rds")
    saveRDS(params, file='params.rds')
    ####END OF ERRATA01

    # Label clusters with SingleR, if requested
    if (params$celldex != FALSE) {
        result = run_singler(params, result, by.cell=classify.cells)
    }

    if (save.rds) {
        result$rdsfile = paste("rds", paste(params$sample, "rds", sep="."), sep="/")
        saveRDS(result, file=result$rdsfile)
    }
    
    message("/// Seurat analysis done.")
    return(result)
}

run_singler <- function(params, result, by.cell=FALSE, use.idents=TRUE, idents=FALSE, wanted.assays=c("RNA")) {
    # Run singleR on a `result' object using parameters from `params'.
    # Assignment of cells to cluster is specified by the idents vector. If use.idents is TRUE,
    # use the Idents slot instead.
    #
    # Classification is performed at the cluster level, and also at the cell level if by.cell is TRUE.
    # Produces the following:
    #   result$cluster.labels.main      - classification of each cluster
    #   result$seurat$cell.clabels.main - classification of each cell copied from the cluster it belongs to
    # If by.cell == TRUE:
    #   result$seurat$cell.labels.main  - classification of each cell (individually)
    # Similarly for fine instead of main.
    
    message(paste("/// Annotating clusters using database:", params$celldex))
    ds = DietSeurat(result$seurat, assays=wanted.assays, dimreducs=c("umap", "tsne"))
    sce = as.SingleCellExperiment(ds)
    ref = load_celldex(params$celldex)
    #ref = celldex::DatabaseImmuneCellExpressionData()
    
    
    
    if (use.idents) {
        idents=Idents(ds)
    }

    message("// Performing annotation.")
    message("/ Main")
    annot.main = SingleR(test=sce, ref=ref, labels=ref$label.main, clusters=idents)
    result$cluster.labels.main = annot.main$pruned.labels
    result$seurat$cell.clabels.main = "NA"
    for (i in result$clusters) {
        result$seurat$cell.clabels.main[Idents(result$seurat) == i] = result$cluster.labels.main[as.numeric(i)+1]
    }
    
    message("/ Fine")
    annot.fine = SingleR(test=sce, ref=ref, labels=ref$label.fine, clusters=idents)
    result$cluster.labels.fine = annot.fine$pruned.labels
    result$seurat$cell.clabels.fine = "NA"
    for (i in result$clusters) {
        result$seurat$cell.clabels.fine[Idents(result$seurat) == i] = result$cluster.labels.fine[as.numeric(i)+1]
    }

    if (by.cell) {
        message("// Performing cell-level annotation.")
        message("/ Main")
        annot.cell = SingleR(test=sce, ref=ref, labels=ref$label.main)
        result$seurat$cell.labels.main = annot.cell$pruned.labels
        message("/ Fine")
        annot.cell = SingleR(test=sce, ref=ref, labels=ref$label.fine)
        result$seurat$cell.labels.fine = annot.cell$pruned.labels
    }
    
    result$annotated = ds
    return(result)
}

load_celldex <- function(name) {
    if (name == "primary") {
        return(HumanPrimaryCellAtlasData())
    } else if (name == "blueprint") {
        return(BlueprintEncodeData())
    } else if (name == "mouse") {
        return(MouseRNAseqData())
    } else if (name == "immgen") {
        return(ImmGenData())
    } else if (name == "dice") {
        return(DatabaseImmuneCellExpressionData())
    } else if (name == "monaco") {
        return(MonacoImmuneData())
    } else if (name == "novershtern") {
        return(NovershternHematopoieticData())
    }
}

### Tests and random stuff

classification_score <- function(data, clusters, idents, fine=FALSE) {
    columns = c('Classification', 'Cells', 'Matching', 'Score')
    res = data.frame(matrix(nrow=0, ncol=length(columns)))
    colnames(res) = columns
    clusters = as.character(clusters)
    if (fine) {
        ncells = length(data$seurat$cell.labels.fine)
        for (i in 1:length(clusters)) {
            cl = clusters[i]
            res[cl,] = c(Classification=data$cluster.labels.fine[i], Cells=0, Matching=0, Score=0)
        }
    } else {
        ncells = length(data$seurat$cell.labels.main)
        for (i in 1:length(clusters)) {
            cl = clusters[i]
            res[cl,] = c(Classification=data$cluster.labels.main[i], Cells=0, Matching=0, Score=0)
        }
    }

    ids = as.character(idents)
    for (i in 1:ncells) {
        cl = ids[i]
        if (fine) {
            a = data$seurat$cell.clabels.fine[i]
            b = data$seurat$cell.labels.fine[i]
        } else {    
            a = data$seurat$cell.clabels.main[i]
            b = data$seurat$cell.labels.main[i]
        }
        res[cl,]$Cells = as.integer(res[cl,]$Cells) + 1
        if ( !is.na(a) & !is.na(b) & a == b ) {
            res[cl,]$Matching = as.integer(res[cl,]$Matching) + 1
        }
    }

    for (cl in clusters) {
        res[cl,]$Score = as.integer(res[cl,]$Matching) / as.integer(res[cl,]$Cells)
    }
    return(res)
}

test1 <- function() {
    params=list()
    params$sample = "PAD1"
    params$mingenes = 500
    params$maxgenes = 5000
    params$mtpatt = "^MT-"
    params$maxmt = 10
    params$method = "LogNormalize"
    params$nfeatures = 2000
    params$ntop = 20
    params$dimensions = 10
    params$resolution = 0.5
    params$minpct = 0.25
    params$logfc = 0.25
    # params$celldex = FALSE
    params$celldex = "primary"
    
    return(seurat_main(params))
}
