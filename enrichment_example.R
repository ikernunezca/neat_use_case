# Example enrichment: Supposing we have a group of 38 patients, all of them with a particular signature of candidate genes,
# and we want to test if those genes are enriched in any of the pathways existing in reactome. 

library(neat)
library(igraph)

# Load input data and put it into a unique data.frame
  a <- read.table("data/genes_per_patient.csv",sep="\t") # An input table indicating the genes associated to each patient of the study. 
  b <- read.table("data/patient_clusters.csv",sep= ",",header=T) # An input table indicating the patient cluster of each patient. Not relevant for this example, I just keep it in order not to affect the script in any way.
  colnames(a)[1] <- "names"
  colnames(b)[1] <- "names"
  genes_per_patient_and_clusters <- merge(b,a,by="names",all = T) # Merging of a and b, a table of 4 columns: 'Patient','Original_cluster','Cluster' and 'Genes'
  # For this use case, we are only be looking to the first and fourth column of this table.

# Download Reactome pathway network and load into an igraph object
  reactome <- as.data.frame(read.csv(file = "https://raw.githubusercontent.com/cirillodavide/gene_multilayer_network/master/networks/Reactome_pathways.19-10-2019.gr",sep = " ",header = F, stringsAsFactors=F))
  class(reactome[,1]) <- "character" # 'Classic R'
  class(reactome[,2]) <- "character"
  class(reactome[,3]) <- "character"
  reactome_graph <- graph_from_data_frame(reactome,directed = FALSE) # In this network, the nodes (genes) are connected if they are part of a same pathway at the reactome database.
  E(reactome_graph)$reactome_ID <- reactome[,3] # To each of the network edges, add a feature that indicates the pathway shared by both nodes.
 #############################################  
# 3. Reactome Data
    
    # Obtain all nodes of reactome net layer
    reactome_available <- names(V(reactome_graph))
    
    # Obtain the first input for neat: a list with the genes per patient appearing in the pathway network
    genes_lists <- list()
    for(i in 1:nrow(genes_per_patient_and_clusters)){ 
      genes <- strsplit(as.character(genes_per_patient_and_clusters[i,4]),"_")[[1]]
      genes_lists[[i]] <- unique(genes[which(genes %in% reactome_available)])
    }
    names(genes_lists) <- genes_per_patient_and_clusters[,1] # genes_lists is an R list object of length 38, with the genes associated to each patient appearing in the reactome network: our 'set_a' to test.
    
    # Obtain the second input for neat: For each pathway, obtain the genes composing it
    genes_per_pathway <- list()
    pathways <- unique(reactome[,3])
    for(i in 1:length(pathways)){
      pathway <- pathways[i]
      genes_per_pathway[[i]] <- unique(unlist(reactome[which(reactome[,3]==pathway),][,1:2]))
    }
    names(genes_per_pathway) <- pathways # genes_per_pathway is an R list objects with all the Reactome pathways, and the genes associated to each in the reactome network: our 'set_b', the target set.
    
    message('Starting pathway enrichments. Please note this process may take a while')
    # Expected run time= 17 min. Up to 16 GB of ram occupied. (This in a machine similar to those at BSC)
    pathway_enrichments <- neat(genes_lists,blist = genes_per_pathway,network = reactome_graph,nettype = 'undirected',nodes=reactome_available,mtc.type = 'fdr',alpha = 0.05	)

    # Save enrichment results
    write.table(pathway_enrichments,file="Output/Reactome_enrichment.tsv",sep="\t",col.names=T,row.names=F)
