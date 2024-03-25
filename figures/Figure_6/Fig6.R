library(ape)
library(geiger)
library(nlme)
library(phytools)
library(nlme)
library(phylosignal)
library(adephylo)
library(ape)
library(phylobase)
### phyloSignal
df <- read.delim("Traits.txt", sep = "\t", row.names = 1, header = T)
df <- as.data.frame(df)
tree <- read.tree("tree.nwk")
p4d <- phylo4d(tree, df)
barplot.phylo4d(p4d, tree.type = "phylo", tree.ladderize = TRUE)
result_df <- phyloSignal(p4d = p4d, method = "all")
H.swap <- phyloCorrelogram(p4d, trait = "Host.swap.rate")
plot(H.swap)
Co.phy <- phyloCorrelogram(p4d, trait = "Co.phylogeny.rate")
plot(Co.phy)
### PGLS
anoleData <- read.csv("data.csv", row.names = 1)
anoleTree <- read.tree("tree.nwk")
#Let’s see what this tree looks like.
plot(anoleTree)
#Geiger has a function to check that the names match between the tree and the data frame.
name.check(anoleTree, anoleData)
#Is there a correlation between awesomeness and hostility?
#plot(anoleData[, c("awesomeness", "hostility")])  
plot(anoleData[, c("Host_swap_rate", "pN_pS_ratio")])
# Extract columns
pN_pS_ratio <- anoleData[,"pN_pS_ratio"]
Host_swap_rate <- anoleData[,"Host_swap_rate"]
# Give them names
names(pN_pS_ratio) <- names(Host_swap_rate) <- rownames(anoleData)

########################################################
#方法 PGLS（Phylogenetic Generalized Least Squares）
########################################################

pglsModel <- gls(pN_pS_ratio ~ Host_swap_rate, correlation = corBrownian(phy = anoleTree, form = ~Species),
                 data = anoleData, method = "ML")

summary(pglsModel)

plot(pN_pS_ratio ~ Host_swap_rate)
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])
