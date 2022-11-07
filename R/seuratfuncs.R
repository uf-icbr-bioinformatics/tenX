# ariva, 11/4/2019

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

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
    data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
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

# Integration and differential analysis

make_seurat_object <- function(name, mtpatt="^mt:", maxmt=5, mingenes=200, maxgenes=5000) {
    message(paste("Loading dataset", name, sep=" "))
    so = CreateSeuratObject(counts = Read10X_h5(paste(name, "/outs/filtered_feature_bc_matrix.h5", sep="")), project=name)
    so[["percent.mt"]] = PercentageFeatureSet(so, pattern=mtpatt)
    genes.no.mt = rownames(so)[grep(mtpatt, rownames(so), invert=TRUE)]
    so = subset(so, subset = percent.mt < maxmt & nFeature_RNA > mingenes & nFeature_RNA < maxgenes, features = genes.no.mt) %>%
        NormalizeData() %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
        AddMetaData(col.name = "Condition", metadata = name)
    return(so)
}

integrated_analysis <- function(smp1, smp2, plots=TRUE, resolution=0.25, logfc=0.25, mtpatt="^mt:", maxmt=5, outdir="") {
    label = paste(outdir, smp1, "_", smp2, sep="")
    so1 = make_seurat_object(smp1, mtpatt=mtpatt, maxmt=maxmt)
    so2 = make_seurat_object(smp2, mtpatt=mtpatt, maxmt=maxmt)

    solist = c(so1, so2)

    features <- SelectIntegrationFeatures(object.list = solist)
    anchors <- FindIntegrationAnchors(object.list = solist, anchor.features = features)
    combined <- IntegrateData(anchorset = anchors)

    DefaultAssay(combined) <- "integrated"

    # Run the standard workflow for visualization and clustering
    combined <- ScaleData(combined, verbose = FALSE) %>%
        RunPCA(npcs = 30, verbose = FALSE) %>%
        RunUMAP(reduction = "pca", dims = 1:30) %>%
        FindNeighbors(reduction = "pca", dims = 1:30) %>%
        FindClusters(resolution = resolution)


#    imgname = paste(label, "ident", "png", sep=".")
#    imgname = paste(label, "clusters", "png", sep=".")
                                        #dev.off()
#    }
    
    # Save marker genes
    combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = logfc)
    combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
    markersfile = paste(label, "markers", "csv", sep=".")
    write.table(combined.markers, file=markersfile, sep='\t', quote=FALSE, col.names=NA, row.names=TRUE)

    # Save all genes
    # all.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
    # all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
    # markersfile = paste(label, "full", "csv", sep=".")
    # write.table(all.markers, file=markersfile, sep='\t', quote=FALSE, col.names=NA, row.names=TRUE)


    cluster_counts = t(table(combined$orig.ident, Idents(combined)))
    # Convert to percentages
    for (row in 1:nrow(cluster_counts)) {
        cluster_counts[row,] = 100.0 * cluster_counts[row,]/sum(cluster_counts[row,])
    }
    write.table(cluster_counts, file=paste(label, "counts", "csv", sep="."), quote=FALSE, col.names=TRUE, row.names=TRUE)
    return(combined)
}

diff_analysis <- function(integ, smp1, smp2, mincells=50, outdir="") {
    filenames = c()
    markers = c()
    label = paste(outdir, smp1, "_", smp2, sep="")
    DefaultAssay(integ) = "RNA"
    clusters = sort(unique(Idents(integ)))

    # Add label consisting of condition and cluster number
    integ$clust.sample = paste("clust", Idents(integ), "_", integ$Condition, sep="")
    integ$celltype = Idents(integ)
    Idents(integ) = "clust.sample"

    for (cl in clusters) {
        id1 = paste("clust", cl, "_", smp1, sep="")
        id2 = paste("clust", cl, "_", smp2, sep="")

        # Write all fold changes
        foldchanges = FoldChange(integ, ident.1 = id1, ident.2 = id2)
        outfile = paste(label, paste("clust", cl, sep=""), "fc", "csv", sep=".")
        write.table(foldchanges, file=outfile, sep='\t', quote=FALSE, col.names=NA, row.names=TRUE)

        if ( sum(integ$clust.sample == id1) > mincells & sum(integ$clust.sample == id2) > mincells) {
            message("*** Processing cluster ", cl)
            response = FindMarkers(integ, ident.1 = id1, ident.2 = id2)
            outfile = paste(label, paste("clust", cl, sep=""), "csv", sep=".")
            filenames = c(filenames, outfile)
            markers = c(markers, response)
            write.table(response, file=outfile, sep='\t', quote=FALSE, col.names=NA, row.names=TRUE)

        }
    }
    return(list(filenames, markers))
}
