require(vegan)
require(dplyr)
require(reshape2)
require(grid)
require(gridExtra)
require(phyloseq)
require(RColorBrewer)
require(scales)

lis_file <- list(rplB_n="/mnt/research/tg/g/glbrc/xander/mercy_run_from_qiong/workdir_otu_nucl/rplB/otutable.tsv_0.05",
                 rplB_p="/mnt/research/tg/g/glbrc/xander/mercy_run_from_qiong/workdir_otu/rplB/otutable.tsv_0.05",
                 rpsC_n="/mnt/research/tg/g/glbrc/xander/mercy_run_from_qiong/workdir_otu_nucl/rpsC/otutable.tsv",
                 rpsC_p="/mnt/research/tg/g/glbrc/xander/mercy_run_from_qiong/workdir_otu/rpsC/otutable.tsv", 
                 SSU_sg="/mnt/research/tg/g/glbrc/ssu_analysis/clust.ssu.sub.clu/workdir_otu/otutable.tsv_0.05",
                 SSU_am="/mnt/research/tg/g/glbrc/rplB_vs_ssu/wkdir_amplicon/core_diversity_analysis.out/table_even57418.biom")


vec_gene <- names(lis_file)
print(vec_gene)

vegan_table <- function(phylo){
    otutable <- otu_table(phylo)
    if (taxa_are_rows(otutable)){
        otutable <- t(otutable)
    }
    otutable <- as.matrix(otutable)
    return(otutable)
}

file_ext <- function(path){
    parts <- strsplit(path, "\\.")[[1]]
    ext <- parts[length(parts)]
    return(ext)
}

# sample data
vec_plant <- rep(c("C","M","S"), c(7,7,7))

func <- function (x){return(paste0(x, 1:7))}
lis_names <- lapply(c("C","M","S"), func)
vec_names <- unlist(lis_names)

df_env <- data.frame(Plant=as.factor(vec_plant), row.names=vec_names)
df_env$dummy <- -1

par(oma=c(4,0,0,0),
    mar=c(2.5,2.5,1.5,0.8), 
    mgp=c(1.5,0.5,0), 
    xpd=F)

colnum <- 2
vec_bottom <- c(length(vec_gene))
vec_temp <- seq(length(vec_gene)-1, 1, -1)
for (i in vec_temp){
    if (i %% colnum == 0){
        break
    }else{
        vec_bottom <- c(vec_bottom, i)
    }
}

#layout(matrix(1:12, ncol=3, byrow=TRUE))
m <- matrix(seq(length(vec_gene)), ncol=colnum, byrow=TRUE)
layout(mat = m)
lis_df <- NULL
# score scaling in vegan: 1 <-> "sites", 2 <-> "species", 3 <-> "symmetric"
scl <- 1 ## scaling = 1
#colvec <- c("red2", "green4", "mediumblue")
colvec <- brewer.pal(n=3, name="Set1")

screen_num <- 1
for (gene in vec_gene){
    path <- lis_file[[gene]]
    letter <- LETTERS[screen_num]
    if (file_ext(path) == "biom"){
        # convert .biom for SSU_am to otu_table
        phylo <- import_biom(path)
        # remove singletons
        #phylo <- prune_taxa(taxa_sums(phylo) >= 2, phylo)
        df <- vegan_table(phylo)
        ### hardcoded; remove 2013 summer samples: T2C1, T2C2..
        df <- df[substring(rownames(df), 1, 2) == "T1",]
        rownames(df) <- substring(rownames(df), 3, 4)
        # remove OTU with no counts across all samples
        df <- df[, ! colSums(df) == 0]
        # missing C1, no good for adonis
        print(sort(rownames(df)))
        df_env_adj <- df_env[rownames(df),] # adjust metadata if some samples are missing
        
    } else {
        df <- read.table(path, sep='\t', header=T)
        rownames(df) <- df$Sample
        df$Sample <- NULL
        df_env_adj <- df_env[rownames(df),]
    }

    lis_df[[gene]] <- df
    
    # convert to relative abundance
    vec_rowsum <- rowSums(df)
    print(vec_rowsum)
    df_subsample <- df/vec_rowsum
    df_subsample_h <- decostand(df_subsample, method="hellinger")
    
    betad <- vegdist(df_subsample, method="bray")
    
    # Plant col in df_env has to be same lenght as rownames in df_subsample (# of sample names in otutable)
    ado <- adonis(betad ~ Plant, data=df_env_adj)
    print(ado)
    p <- ado$aov.tab$Pr[1]
    R2 <- ado$aov.tab$R2[1]
    
    if (p < 0.001){
        label <- "***"
    } else if (p < 0.01){
        label <- "**"
    } else if (p < 0.05){
        label <- "*"
    } else if (p < 0.1){
        label <- "."
    } else {
        label <- ""
    }
    
    # xlab only at bottom, ylab only at left
    ylab <- ""
    xlab <- ""
    if (screen_num %% colnum == 1){
        ylab <- "PC2"
    }else{
        ylab <- ""
    }
    
    if (screen_num %in% vec_bottom){
        xlab <- "PC1"
    }else{
        xlab <- ""
    }

    #title_main <- sprintf("%s(R^2=%.1f%%,%s)", gene, R2*100, label)
    R2_s <- sprintf("%.1f%%", R2*100)
    mod <- rda(df_subsample, scale=F)
    plot(mod, type="n", main="", xlab=xlab, ylab=ylab, scaling=scl)
    title(main=bquote(paste(.(letter), ". ", .(gene), .(label), "(", R^2, "=", .(R2_s), ")")), 
          line=0.8, cex.main=1, fond.main=2)
    
    # plot sites
    with(df_env, points(mod, display = "sites", 
                        col = alpha(colvec[Plant], 0.8), scaling = scl, pch = 21, cex=1.5))
    
    screen_num <- screen_num + 1
    cat("\n\n", gene, " is finished..\n")
}

## add legend if there is one pane left blank
#plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#with(df_env, legend("center", inset=0, legend = levels(Plant), bty = "n", 
#                   col = colvec, pch = 21, horiz=F, cex=2))

# add legend to bottom margin if there is no pane left blank
# inspired by: http://dr-k-lo.blogspot.com/2014/03/the-simplest-way-to-plot-legend-outside.html
par(fig=c(0,1,0,1), oma=c(0,0,0,0), mar=c(0,0,0,0), new=T)
plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n")
legend("bottom", legend = levels(df_env$Plant), xpd=T, horiz=T, inset=c(0,0), bty="n", pch=21, col=colvec, cex=2)
