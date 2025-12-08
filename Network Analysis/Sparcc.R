conda activate fastspar

cd /mnt/80G7/zhang/Netana/

./fsparcc.sh Abundance_merge.txt  NJFL


#Correlation matrix of observed values

cor_sparcc <- read.delim('NJFL_cor.txt', row.names = 1, sep = '\t', check.names = FALSE)

#Pseudo p-value matrix

pvals <- read.delim('pvalues.tsv', row.names = 1, sep = '\t', check.names = FALSE)

#Retain values with |correlation| ≥ 0.6 and p < 0.01

cor_sparcc[abs(cor_sparcc) < 0.6] <- 0

pvals[pvals >= 0.01] <- -1
pvals[pvals < 0.01 & pvals >= 0] <- 1
pvals[pvals == -1] <- 0

#Filtered adjacency matrix

adj <- as.matrix(cor_sparcc) * as.matrix(pvals)
diag(adj) <- 0  # Convert diagonal values (indicating autocorrelation) to 0
write.table(data.frame(adj, check.names = FALSE), 'network.adj.txt', col.names = NA, sep = '\t', quote = FALSE)

#Network format conversion

library(igraph)

#Input data: adjacency matrix

neetwork_adj <- read.delim('network.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(neetwork_adj)[1:6]  # Network file in adjacency matrix format

#Adjacency matrix -> igraph adjacency list, resulting in a weighted undirected network 
 
g <- graph_from_adjacency_matrix(as.matrix(neetwork_adj), mode = 'undirected', weighted = TRUE, diag = FALSE)
g  # igraph adjacency list

#In this conversion mode, the default edge weight represents the sparcc-calculated correlation (may contain negative values)
#Since edge weights are generally positive, it's best to take the absolute value, with the correlation stored in a new column as a record 

E(g)$sparcc <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)


#Convert back from igraph adjacency list to adjacency matrix 

adj_matrix <- as.matrix(get.adjacency(g, attr = 'sparcc'))

write.table(data.frame(adj_matrix, check.names = FALSE), 'network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#GML format, which can be opened and visually edited using Gephi software
 
write.graph(g, 'network.gml', format = 'gml')

