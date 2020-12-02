library(neat)
library(igraph)

# Load input data and put it into a unique data.frame
  a <- read.table("data/genes_per_patient.csv",sep="\t")
  b <- read.table("data/patient_clusters.csv",sep= ",",header=T)
  colnames(a)[1] <- "names"
  colnames(b)[1] <- "names"
  genes_per_patient_and_clusters <- merge(b,a,by="names",all = T)

# Load Reactome pathway network
  reactome <- as.data.frame(read.csv(file = "https://raw.githubusercontent.com/cirillodavide/gene_multilayer_network/master/networks/Reactome_pathways.19-10-2019.gr",sep = " ",header = F, stringsAsFactors=F))
  class(reactome[,1]) <- "character"
  class(reactome[,2]) <- "character"
  class(reactome[,3]) <- "character"
  reactome_graph <- graph_from_data_frame(reactome,directed = FALSE)
  E(reactome_graph)$reactome_ID <- reactome[,3]
 #############################################  
# 3. Reactome Data
    
    # Obtain all vertices of reactome net layer
    reactome_available <- names(V(reactome_graph))
    
    # genes per patient appearing in the pathway network
    genes_lists <- list()
    for(i in 1:nrow(genes_per_patient_and_clusters)){ 
      genes <- strsplit(as.character(genes_per_patient_and_clusters[i,4]),"_")[[1]]
      genes_lists[[i]] <- unique(genes[which(genes %in% reactome_available)])
    }
    names(genes_lists) <- genes_per_patient_and_clusters[,1] # genes_lists is an R list object with the genes associated to each patient appearing in the reactome network: our 'set_a' to test.
    
    # For each pathway, obtain genes associated because of it
    genes_per_pathway <- list()
    pathways <- unique(reactome[,3])
    for(i in 1:length(pathways)){
      pathway <- pathways[i]
      genes_per_pathway[[i]] <- unique(unlist(reactome[which(reactome[,3]==pathway),][,1:2]))
    }
    names(genes_per_pathway) <- pathways # genes_per_pathway is a list of pathways, and the genes associated to each in the reactome network: our 'set_b', the target set.
    
    message('Starting pathway enrichments. Please note this process may take a while')
    pathway_enrichments <- neat(genes_lists,blist = genes_per_pathway,network = reactome_graph,nettype = 'undirected',nodes=reactome_available,mtc.type = 'fdr',alpha = 0.05	)
    #Expected run time= 17 min. Up to 16 GB of ram occupied.
    # Save enrichment results
    write.table(pathway_enrichments,file="Output/Reactome_enrichment.tsv",sep="\t",col.names=T,row.names=F)
