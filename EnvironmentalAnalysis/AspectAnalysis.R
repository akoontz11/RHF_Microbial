library(ape)
library(PhyloMeasures)
library(pez)

# 16S----
setwd("/home/austin/Documents/PearseLab/RHFSoilSequencing/Analysis/Output/")
bact.table <- readRDS("16S_table.Rdata")
rownames(bact.table) <- gsub(pattern = "-16S", replacement = "", rownames(bact.table))
# Make a presence/absence version of the table
bact.table.PA <- bact.table
bact.table.PA[bact.table.PA != 0] <- 1
# Add together species for each site (row)
bact.sums <- apply(bact.table.PA, 1, sum)
length(bact.sums)
# Read in environmental data from Elizabeth
setwd("/home/austin/Documents/PearseLab/RHFSoilSequencing/Analysis/EnvironmentalData/")
env <- read.csv("clean_rhf_terrain.csv", as.is=TRUE)
rownames(env) <- env$plot_id
length(env$aspect)
# Aspects include one more site than bacterial OTU table does
which(!env$plot_id %in% names(bact.sums))
# Site 51. Remove from environmental data
env <- env[-51,]

# Analyze alpha diversity
cor(env$aspect, bact.sums)
plot(bact.sums ~ env$aspect, main = "Bacterial OTU richness (alpha diversity)", pch=16)
abline(lm(bact.sums ~ env$aspect), col="red", lwd=3)

# Calculate phylogenetic diversity
setwd("/home/austin/Documents/PearseLab/RHFSoilSequencing/Analysis/Output/")
bact.tree <- read.tree("RHF_16S_RootedTree.nwk")
length(bact.tree$tip.label)
ncol(bact.table)
# Phylogeny and OTU table differ by one taxa...
which(!bact.tree$tip.label %in% colnames(bact.table))
# Tip 26814. Remove that tip from the phylogeny
bact.tree <- drop.tip(bact.tree, 26814)
# Extract abundances of each OTU from bact.table
bact.abundances <- apply(bact.table, 2, sum)
# Calculate MPDs from tree, presence/absence table, and abundances
bact.MNTDs <- mntd.query(bact.tree, bact.table.PA, abundance.weights = bact.abundances)

# Analyze phylogenetic diversity
cor(env$aspect, bact.MNTDs)
plot(bact.MNTDs ~ env$aspect, main = "Bacterial OTU MNTDs (phylogenetic diversity)", pch=16)
abline(lm(bact.MNTDs ~ env$aspect), col="red", lwd=3)

# ITS----
setwd("/home/austin/Documents/PearseLab/RHFSoilSequencing/Analysis/Output/")
fung.table <- readRDS("ITS_table.Rdata")
rownames(fung.table) <- gsub(pattern = "-ITS", replacement = "", rownames(fung.table))
# Make a presence/absence version of the table
fung.table.PA <- fung.table
fung.table.PA[fung.table.PA != 0] <- 1
# Add together species for each site (row)
fung.sums <- apply(fung.table.PA, 1, sum)
length(fung.sums)
# Use same env$aspect vector (preliminary analyses showed that site 51 was missing from fungal OTUs, as well)

# Analyze alpha diversity
cor(env$aspect, fung.sums)
plot(fung.sums ~ env$aspect, main = "Fungal OTU richness (alpha diversity)", pch=16)
abline(lm(fung.sums ~ env$aspect), col="red", lwd=3)

# Calculate phylogenetic diversity
setwd("/home/austin/Documents/PearseLab/RHFSoilSequencing/Analysis/Output/")
fung.tree <- read.tree("RHF_ITS_RootedTree.nwk")
length(fung.tree$tip.label)
ncol(fung.table)
# Phylogeny and OTU table differ by one taxa...
which(!fung.tree$tip.label %in% colnames(fung.table))
# Tip 26814. Remove that tip from the phylogeny
fung.tree <- drop.tip(fung.tree, 9282)
# Extract abundances of each OTU from fung.table
fung.abundances <- apply(fung.table, 2, sum)
# Calculate MPDs from tree, presence/absence table, and abundances
fung.MNTDs <- mntd.query(fung.tree, fung.table.PA, abundance.weights = fung.abundances)

# Analyze phylogenetic diversity
cor(env$aspect, fung.MNTDs)
plot(fung.MNTDs ~ env$aspect, main = "Fungal OTU MNTDs (phylogenetic diversity)", pch=16)
abline(lm(fung.MNTDs ~ env$aspect), col="red", lwd=3)
# Plot all 4
par(mfcol=c(2,2), mar = c(4,2,0,4), oma = c(0, 4, .5, .5), mgp = c(2, 0.6, 0))

plot(bact.sums ~ env$aspect, pch=16)
abline(lm(bact.sums ~ env$aspect), col="red", lwd=3)

plot(bact.MNTDs ~ env$aspect, pch=16)
abline(lm(bact.MNTDs ~ env$aspect), col="red", lwd=3)

plot(fung.sums ~ env$aspect, pch=16)
abline(lm(fung.sums ~ env$aspect), col="red", lwd=3)

plot(fung.MNTDs ~ env$aspect, pch=16)
abline(lm(fung.MNTDs ~ env$aspect), col="red", lwd=3)

# Plant data----
setwd("/home/austin/Documents/PearseLab/RHFSoilSequencing/Analysis/EnvironmentalData/")
# Load in cover data, reformat, remove NAs, remove species essentially there as notes
data <- read.csv("rhf_2018_cover.csv", as.is=TRUE)
plant.comm <- with(data, tapply(Cover, list(Plot_id,Species), function(x) mean(x,na.rm=TRUE)))
plant.comm[is.na(plant.comm)] <- 0
plant.comm <- plant.comm[,!grepl("^[a-z]+", colnames(plant.comm))]
plant.comm <- plant.comm[,!grepl("\\(|/", colnames(plant.comm))]
plant.comm <- plant.comm[,!grepl("sp\\.", colnames(plant.comm))]
colnames(plant.comm) <- tolower(colnames(plant.comm))
# Plant sites differ from soil microbial sites
which(!rownames(plant.comm) %in% names(bact.sums))
# Remove site 51 from plant matrix
plant.comm <- plant.comm[-51,]
# Load in plant tree, reformat sp names
plant.tree <- read.tree("Vascular_Plants_rooted.dated.tre")
plant.tree$tip.label <- tolower(gsub("_", " ", plant.tree$tip.label))
#plant.tree <- congeneric.merge(plant.tree, colnames(plant.comm), split = " ")
colnames(plant.comm) %in% plant.tree$tip.label
# Not all RHF plants are in Zanne vascular plant tree
# So, remove the plants that we can't calculate phylogenetic diversity for, and prune down the large tree
plant.comm <- plant.comm[,which(colnames(plant.comm) %in% plant.tree$tip.label)]
plant.tree <- keep.tip(plant.tree, colnames(plant.comm))
# Create presence/absence matrix
plant.comm.PA <- plant.comm
plant.comm.PA[plant.comm.PA != 0] <- 1
# Add together species for each site (row)
plant.sums <- apply(plant.comm.PA, 1, sum)

# Analyze alpha diversity
# Bacteria
cor(bact.sums, plant.sums)
plot(plant.sums ~ bact.sums, main = "Plants vs. Bacteria", pch=16)
abline(lm(plant.sums ~ bact.sums), col="red", lwd=3)
# Fungi
cor(fung.sums, plant.sums)
plot(plant.sums ~ fung.sums, main = "Plants vs. Fungi", pch=16)
abline(lm(plant.sums ~ fung.sums), col="red", lwd=3)

# Calculate phylogenetic diversity
length(plant.tree$tip.label)
ncol(plant.comm)
# Extract abundances of each species from plant.table
plant.abundances <- apply(plant.comm, 2, sum)
# Calculate MPDs from tree, presence/absence table, and abundances
plant.MNTDs <- mntd.query(plant.tree, plant.comm.PA, abundance.weights = plant.abundances)

# Analyze phylogenetic diversity
# Bacteria
cor(bact.MNTDs, plant.MNTDs)
#plot(plant.MNTDs ~ bact.MNTDs, main = "Plant vs. Bacteria MNTDs", pch=16, col=4)
plot(plant.MNTDs ~ scale(bact.MNTDs), main = "Plant vs. Microbial MNTD per site", pch=16, col=4,  
     ylab = "Plant MNTD values", xlab = "Microbial MNTD values")
#abline(lm(plant.MNTDs ~ bact.MNTDs), col="red", lwd=3)
# Fungi
cor(fung.MNTDs, plant.MNTDs)
#plot(plant.MNTDs ~ fung.MNTDs, main = "Plant vs. Fungi MNTDs", pch=16)
points(plant.MNTDs ~ scale(fung.MNTDs), main = "Plant vs. Fungi MNTDs", pch=16, col=6)
legend("bottomright", c("Bacteria", "Fungi"), col = c(4,6), pch = 16)
#abline(lm(plant.MNTDs ~ fung.MNTDs), col="red", lwd=3)

# Aspect
cor(env$aspect, scale(plant.MNTDs))
plot(plant.MNTDs ~ env$aspect, main = "Plant MNTDs vs. aspect", pch=16)
abline(lm(plant.MNTDs ~ env$aspect), col="red", lwd=3)


# SES_MNTD calculations--NOT RUN----
c.bact.data <- comparative.comm(bact.tree, bact.table, env=env)
bact.ses.mntd <- .ses.mntd(c.bact.data)

c.fung.data <- comparative.comm(fung.tree, fung.table, env=env)
fung.ses.mntd <- .ses.mntd(c.fung.data)
