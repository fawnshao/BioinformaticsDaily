library(igraph)
args <- commandArgs(TRUE)
# args <- c("interaction", "nodes")
# args <- c("interaction.SRR5831489.HKG.links", "interaction.SRR5831489.HKG.links.nodes")

# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
    }
    vertex.size <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
    }
    
    symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
            stars=cbind(vertex.size, vertex.size, vertex.size),
            add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)
#################################################################
# generic star vertex shape, with a parameter for number of rays
mystar <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
    }
    vertex.size  <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
    }
    norays <- params("vertex", "norays")
    if (length(norays) != 1 && !is.null(v)) {
        norays <- norays[v]
    }
    
    mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
           FUN=function(x, y, bg, size, nor) {
               symbols(x=x, y=y, bg=bg,
                       stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
                       add=TRUE, inches=FALSE)
           })
}
# no clipping, edges will be below the vertices anyway
add_shape("star", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=5))
#################################################################

links <- read.table(args[1], header = F, sep = "\t")
vertices <- read.table(args[2], header = F, sep = "\t")
network <- graph_from_data_frame(d = links, directed = F, vertices = vertices)
# vertices$group <- edge.betweenness.community(network)$membership # betweeness centrality for each node for grouping
vertices$shape <- rep("circle", nrow(vertices))
vertices$shape[vertices$V2 == "TSG"] <- "sphere"
vertices$shape[vertices$V2 == "HKG"] <- "star"
vertices$shape[vertices$V2 == "Enh"] <- "triangle"
vertices$color <- rep("grey", nrow(vertices))
vertices$color[vertices$V2 == "TSG"] <- "blue"
vertices$color[vertices$V2 == "HKG"] <- "red"
vertices$color[vertices$V2 == "Enh"] <- "green"
pdf(file = paste(args[1], "igraph.pdf", sep = "."), width = 50, height = 50)
plot(network,
     # mark.groups=vertices$group, # group vertices by betweeness indicator (redish blob background)
     layout = layout.auto, #layout_nicely, 
     # vertex.color = vertices$group, # color vertices by edge betweeness
     vertex.color = vertices$color,
     vertex.label = NA, # no vertex label (name)
     vertex.size = 5,
     # edge.arrow.size = 0.8,
     edge.width = min(links$V3, 5),
     vertex.shape = vertices$shape)
dev.off()

##### other analysis
# # g <- sample_gnp(100, 0.3)
# clique_num(network)
# cliques(network, min = 6)
# largest_cliques(network)
# 
# # To have a bit less maximal cliques, about 100-200 usually
# max_cliques(network)
# 
# a <- largest.cliques(network)
# # let's just take the first of the largest cliques
# # (in this case there's just one clique)
# clique1 <- a[[1]]
# 
# # subset the original graph by passing the clique vertices
# g2 <- induced.subgraph(graph = network, vids = clique1)
# 
# # plot the clique
# plot(g2)


################
# g <- simplify(
#     graph.compose(
#         graph.ring(10), 
#         graph.star(5, mode = "undirected")
#     )
# ) + edge("7", "8")

g <- sample_gnp(20, 1/20)
plot(g)
gc <- clusters(g)
lapply(seq_along(gc$csize)[gc$csize > 1], function(x) 
    V(g)$name[gc$membership %in% x])

set.seed(1)
g <- erdos.renyi.game(20, 1/20)
V(g)$name <- letters[1:20]
par(mar=rep(0,4))
plot(g)

dg <- decompose.graph(g) # returns a list of three graphs
plot(dg[[1]]) # plot e.g. the 1st one
plot(dg[[2]])
plot(dg[[3]])
allnodes <- length(V(graph = g))
node2community <- data.frame(Node = rep(NA, allnodes), Community = rep(NA, allnodes), CommunitySize = rep(NA, allnodes))
k <- 1
for(i in 1:length(dg)){
    a <- names(V(graph = dg[[i]]))
    size <- length(a)
    for(j in 1:size){
        node2community[k, ] <- c(a[j], paste("comp", i, sep = "."), size)
        k <- k + 1
    }
}
node2community
cl <- clusters(g)
lapply(seq_along(cl$csize)[cl$csize > 1], function(x) 
    V(g)$name[cl$membership %in% x])



components <- decompose(g, min.vertices=2)
sapply(components, diameter)
plot(components)

clu <- components(g)
groups(clu)
plot(clu)

subcomponent(g, 1, "all")
plot(subgraph.edges(g, 1:5))
