library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(harmony)
library(foreach)
library(ggplot2)

create_objects = function(names) {
    objects = foreach(i = 1:length(names), .combine='c', .init=c()) %do% {
        name = names[i]
        message(paste("Loading", name, "..."))
        path = paste(name, "/outs/filtered_feature_bc_matrix/", sep="")
        data = Read10X(data.dir = path)
        CreateSeuratObject(counts = data, project = name)
    }
    return(objects)
}

merge_objects = function(names, objects) {
    o1 = objects[[1]]
    others = objects[2:length(objects)]
    merged = merge(o1, y = others, add.cell.ids = names, project = "Merged")
    merged <- NormalizeData(merged) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
    merged <- RunUMAP(merged, dims = 1:30)
    p = DimPlot(merged, group.by = "orig.ident")
    ggsave("harmony/pre-correction.png", plot=p)
    return(merged)
}

run_harmony = function(names, merged) {       
    merged <- RunHarmony(merged, group.by.vars = "orig.ident")
    merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30)
    merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30) %>% FindClusters()
                                        # "orig.ident" = original identity
    p = DimPlot(merged, group.by = "orig.ident")
    ggsave("harmony/corrected-ident.png", plot=p)
                                        # "ident" = identity, which are clusters
    p = DimPlot(merged, group.by = "ident", split.by = 'orig.ident')
    ggsave("harmony/corrected-clusters.png", plot=p, units='in', width=30, height=6)
    return(merged)
}

export_projections = function(names, merged) {
    # split the object
    corrected.data <- SplitObject(merged, split.by = "orig.ident")
    corrected.clusters = c()
    
    umaps = foreach (i = 1:length(names)) %do% {
        name = names[i]
        barcode = rownames(Embeddings(object = corrected.data[[name]], reduction = "umap"))
        barcode = gsub(paste(name, "_", sep=""), "", barcode)
        barcode = gsub('.{2}$', '', barcode)
        barcode = paste(barcode, i, sep="-")

        clusters = Idents(corrected.data[[name]])
        cluster.data = cbind("Barcode" = barcode, data.frame("clusters" = clusters))
        corrected.clusters = rbind(corrected.clusters, cluster.data)
        
        proj = Embeddings(object = corrected.data[[name]], reduction = "umap")
        cbind("Barcode" = barcode, proj)
    }
    corrected.umap = do.call(rbind, umaps)
    write.table(corrected.umap, file="harmony/corrected_umap.csv", sep = ",", quote = F, row.names = F, col.names = T)
    write.table(corrected.clusters, file="harmony/corrected_clusters.csv", sep = ",", quote = F, row.names = F, col.names = T)
}

names = commandArgs(trailingOnly=TRUE)

message(paste("names:", names))
objects = create_objects(names)
message("Merging objects ...")
merged = merge_objects(names, objects)
message("Running harmony ...")
merged = run_harmony(names, merged)
message("Exporting projections and clusters ...")
export_projections(names, merged)
message("Done.")
