##########################################################
## Codes for building the Bayesian netowrk for HIF-1
## pathway using gene expression data from AD and 
## control subjects
##########################################################

library("BiGGR")
library("org.Hs.eg.db")
library("BayesNetBP")
library("graphite")
library("pcalg")
library("igraph")
data("Recon2")

file.name <- system.file("extdata", "brainmodel_reactions.txt", package="BiGGR")
reaction.ids <- scan(file.name, what=" ")
reaction.ids[reaction.ids=="R_SSALxm"] <- "R_r0179"
gene.info <- extractGeneAssociations(Recon2)
sbml.model <- buildSBMLFromReactionIDs(reaction.ids, Recon2)
genes.in.rx <- gene.info[reaction.ids]

###################################################
## Set the group: control or AD
###################################################

group <- "AD" #  options: "control", "AD"

##########################################################
## map genes to the reactions
##########################################################

extract.genes <- function(s){
  if (is.null(s)) return(NULL)
  s1 <- strsplit(s, split=" or | and ")[[1]]
  s2 <- sub("^ ?\\(|^ ", "", s1)
  return(unique(sub("\\.[0-9]\\)?$", "", s2))) }
genes <- lapply(genes.in.rx, extract.genes)
sym <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(sym)
mapping <- as.list(sym[mapped_genes])
symbols.rx <- lapply(genes, mp <- function(g){unlist(mapping[g])}) 

########################################################
##  Load in the gene expression data
########################################################
load("./data/brain_dats.RData")
colnames(brain_dats)[6:18] <- paste("Controls", 1:13, sep = ":")
colnames(brain_dats)[19:28] <- paste("AD", 1:10, sep = ":")

control_dats <- brain_dats[,c(1:18)]
AD_dats <- brain_dats[,c(1:5, 19:28)]

# check the fold change
means_control_dats <- apply(control_dats[,6:18], 1, "mean")
mean_AD_dats <- apply(AD_dats[,6:15], 1, "mean")
fold <-mean_AD_dats/means_control_dats

###########################################################################
##  obtain the gene ids that correspond to the metabolic model
###########################################################################

reactions <- c()
symbols <- c()
entrezs <- c()
for(i in 1:length(symbols.rx)){
    temp <- length(symbols.rx[[i]])
    reactions <- c(reactions, rep(names(symbols.rx[i]), temp))
    symbols <- c(symbols, unlist(symbols.rx[[i]]))
    entrezs <- c(entrezs, names(symbols.rx[[i]]))
}

names(symbols) <- NULL
mappings <- data.frame(reactions, symbols, entrezs)

control_sub <- c()
AD_sub <-c()
flag <- c()
for (i in 1:length(mappings[,1])){
    temp <- which(as.matrix(control_dats$Entrez) == as.matrix(mappings$entrezs)[i])
    temp <- temp[1]
    if (is.na(temp)){
        flag <- c(flag, i)
    } else {
        AD_sub <- rbind(AD_sub, AD_dats[temp,])
        control_sub <- rbind(control_sub, control_dats[temp, ])
    }
}

AD_sub <- cbind(mappings[-flag,], AD_sub)
control_sub <- cbind(mappings[-flag, ], control_sub)

########################################################################
## build the Bayesian network model
########################################################################

kegg  <- pathways("hsapiens", "kegg")
graph <- convertIdentifiers(kegg[["HIF-1 signaling pathway"]], "symbol")
graph <- pathwayGraph(graph)
graph::nodes(graph) <- substr(graph::nodes(graph), start = 8, stop = 1000)

in.graph <- nodes(graph)
in.react <- unique(unlist(symbols.rx))
intersect(in.graph, in.react)

symbols.in.data <- list()
symbols.in.data <- symbols.rx
for (i in 1:length(symbols.rx)) {
  symbols.in.data[[i]] <- symbols.rx[[i]][symbols.rx[[i]] %in% control_dats$SYMBOL]
}

symbols.in.both <- intersect(in.graph, control_dats$SYMBOL)
hif.igraph <- igraph.from.graphNEL(graph)
hif.igraph.2 <- induced_subgraph(hif.igraph, symbols.in.both)
hif.graph.2 <- igraph.to.graphNEL(hif.igraph.2)
graph.3 <- pcalg::pdag2dag(hif.graph.2, keepVstruct = FALSE)[[1]]

##########################################################

genes <- as.character(brain_dats$SYMBOL)
## get the expression data for control or AD group
if (group == "AD") {
  dat <- t(brain_dats[, 19:28])  ## 6:18 control; 19:28 AD
} else {
  dat <- t(brain_dats[, 6:18])
}

s <- c()
for (i in 1:length(symbols.in.both)) {
  ind <- which(genes == symbols.in.both[i])
  if (length(ind)==1) {
    s[i] <- ind
    next
  }
  temp <- dat[, ind]
  vars <- apply(temp, 2, var)
  pos <- which.max(vars)
  s[i] <- ind[pos]
}

# as.character(control_dats$SYMBOL[s])
# symbols.in.both
df <- data.frame(dat[, s])
colnames(df) <- symbols.in.both
rownames(df) <- NULL
node.class <- rep(FALSE, length(symbols.in.both))
names(node.class) <- symbols.in.both

if (group == "AD") {
  # initialize the AD model
  tree.init.ad <- Initializer(graph.3, df, node.class) 
} else {
  # initialize the Control model, need to update the df
  tree.init.ctrl <- Initializer(graph.3, df, node.class) 
}

SummaryMarginals(Marginals(tree.init.ad, "HIF1A"))

##########################################################
## Filtering the reactions
##########################################################

rxs <- names(symbols.rx)
constraint.rx <- list()
rx <- c()
k <- 1
for (i in 1:length(symbols.rx)) {
  this.enzymes <- symbols.rx[[i]][symbols.rx[[i]] %in% symbols.in.both]
  if (length(this.enzymes)>0) {
    constraint.rx[[k]] <- this.enzymes
    rx[k] <- rxs[i]
    k <- k + 1
  }
}
names(constraint.rx) <- rx

##########################################################
## Before saving the files, run above codes for control 
## and AD group separately,
## so that both tree.init.ctrl and tree.init.ad are generated
##########################################################

save(tree.init.ctrl, tree.init.ad, constraint.rx, file="temp.rda")
