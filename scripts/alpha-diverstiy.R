require(ggplot2)
require(vegan)
require(dplyr)
require(reshape2)
require(grid)
require(gridExtra)
require(ggpubr)
require(phyloseq)

lis_file <- list(
                 rplB_p="/mnt/research/tg/g/glbrc/xander/mercy_run_from_qiong/workdir_otu/rplB/otutable.tsv_0.05", 
                 rplB_n="/mnt/research/tg/g/glbrc/xander/mercy_run_from_qiong/workdir_otu_nucl/rplB/otutable.tsv_0.05",
                 rpsC_n="/mnt/research/tg/g/glbrc/xander/mercy_run_from_qiong/workdir_otu_nucl/rpsC/otutable.tsv",
                 rpsC_p="/mnt/research/tg/g/glbrc/xander/mercy_run_from_qiong/workdir_otu/rpsC/otutable.tsv", 
                 SSU_sg="/mnt/research/tg/g/glbrc/ssu_analysis/clust.ssu.sub.clu/workdir_otu/otutable.tsv_0.05",
                 SSU_am="/mnt/research/tg/g/glbrc/rplB_vs_ssu/wkdir_amplicon/core_diversity_analysis.out/table_even57418.biom"
)


even_sample_lis <- list(rplB_p=1200, rplB_n=1200, rpsC_n=1200, rpsC_p=1200, SSU_sg=1200, SSU_am=1200)

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


vec_gene <- names(lis_file)
print(vec_gene)

# sample data
vec_plant <- rep(c("C","M","S"), c(7,7,7))

func <- function (x){return(paste0(x, 1:7))}
lis_names <- lapply(c("C","M","S"), func)
vec_names <- unlist(lis_names)

df_env <- data.frame(Plant=as.factor(vec_plant), row.names=vec_names)
df_env$dummy <- -1

lis_df <- NULL
# score scaling in vegan: 1 <-> "sites", 2 <-> "species", 3 <-> "symmetric"
scl <- 1 ## scaling = 1
colvec <- c("red2", "green4", "mediumblue")
cols <- with(df_env, colvec[Plant])
print(df_env$Plant)
print(cols)

lis_df_div <- NULL
for (gene in vec_gene){
    path <- lis_file[[gene]]
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
        df_env <- df_env[rownames(df),]
        }
    else{
        df <- read.table(path, sep='\t', header=T)
        rownames(df) <- df$Sample
        df$Sample <- NULL
    }

    lis_df[[gene]] <- df

    # convert to int, remove rows and cols with all 0
    df_int <- round(df)
    
    print(rowSums(df_int))
    
    pick_row_remove0 <- !!rowSums(abs(df_int))
    df_int <- df_int[pick_row_remove0,]
    df_env_temp <- df_env[pick_row_remove0,]
    # remove col with all 0
    df_int <- df_int[,!!colSums(abs(df_int))]
    
    
    vec_rowsum <- rowSums(df_int)
    minsize <- min(vec_rowsum)
    raremax <- even_sample_lis[[gene]]
    pick_row_subsample <- vec_rowsum >= raremax
    df_int <- df_int[pick_row_subsample,]
    df_env_temp <- df_env_temp[pick_row_subsample,]
    vec_rowsum2 <- rowSums(df_int)
    
    cat(gene, ":\n")
    cat(nrow(df), " rows ==> ", nrow(df_int), " rows after convert float to int\n")
    cat(ncol(df), " cols ==> ", ncol(df_int), " cols after convert float to int\n")
    
    sample_vec_temp <- rownames(df_env_temp)
    sample_vec_temp_sort <- sample_vec_temp[order(vec_rowsum2)]
    
    cat("Sample size per site (sorted): ", paste(sample_vec_temp_sort, sort(vec_rowsum2), sep=":"), "\n")
    cat("Even sampled to", raremax, "\n")
    cat("========\n")
    
    df_int_rare <- rrarefy(df_int, raremax)
    H <- diversity(df_int_rare, "shannon")
    simp <- diversity(df_int_rare, "simpson")
    invsimp <- diversity(df_int_rare, "inv")
    otu_num <- specnumber(df_int_rare)
    chao1 <- estimateR(df_int_rare)["S.chao1",]
    
    df_div <- data.frame(otu_num=otu_num, chao1=chao1, H=H, Plant=df_env_temp$Plant)
    #df_div <- data.frame(H=H, Plant=df_env_temp$Plant)
    df_div$Gene <- gene
    lis_df_div[[gene]] <- df_div
    print(df_div)
    
    # ===================== otu_num
    oneway_test_out <- oneway.test(otu_num~Plant, data=df_div)
    p <- oneway_test_out$p.value
    print("-----------> oneway.test p value:")
    print(p)
    print(oneway_test_out)
    
    pair_t_test <- pairwise.t.test(df_div$otu_num, df_div$Plant, p.adj="none")
    print(pair_t_test)
    
    # ===================== chao1
    oneway_test_out <- oneway.test(chao1~Plant, data=df_div)
    p <- oneway_test_out$p.value
    print("-----------> oneway.test p value:")
    print(p)
    print(oneway_test_out)
    
    pair_t_test <- pairwise.t.test(df_div$chao1, df_div$Plant, p.adj="none")
    print(pair_t_test)
    
    # ======================= H
    oneway_test_out <- oneway.test(H~Plant, data=df_div)
    p <- oneway_test_out$p.value
    print("-----------> oneway.test p value:")
    print(p)
    print(oneway_test_out)
    
    pair_t_test <- pairwise.t.test(df_div$H, df_div$Plant, p.adj="none")
    print(pair_t_test)
}

df_all <- do.call("rbind", lis_df_div)
write.table(df_all, file="alpha_div_df_all.tsv", sep="\t")
melt_df <- melt(df_all, id.vars=c("Plant", "Gene"), measure.vars=c("otu_num", "chao1", "H"), 
                variable.name="Type", value.name="Value")

write.table(melt_df, file="alpha_div_melt_df.tsv", sep="\t")

print("aaaa")
theme_set(theme_bw())
gg <- (
    ggplot(data=melt_df, aes(x=Plant, y=Value, color=Gene)) + geom_boxplot() 
    + facet_wrap(~Type, ncol=1, scales="free_y", labeller=labeller(Type=c(otu_num="OTU #", chao1="Chao1", H="Shannon")))
    + labs(x="", y="Alpha diversity") + color_palette(pal="jco")
    ) #
print(gg)
gg <- (
    ggboxplot(melt_df, x="Plant", y="Value", color="Gene", palette="jco") 
    + facet_wrap(~Type, ncol=1, scales="free_y", labeller=labeller(Type=c(otu_num="OTU #", chao1="Chao1", H="Shannon")))
    + labs(x="", y="Alpha diversity")
    + stat_compare_means(aes(group="Gene"), label="p.format")
)
print(gg)

theme_set(theme_classic())
gg2 <- ggplot(data=df_all, aes(x=Plant, y=otu_num, color=Gene)) + geom_boxplot() + labs(x="", y="OTU #") + scale_color_brewer(type="qual", palette="Set1")
print(gg2)

