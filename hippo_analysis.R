library(MetaGx)

nbcore <- 8
availcore <- parallel::detectCores()
if (is.null(nbcore) || nbcore > availcore) { nbcore <- availcore }
options("mc.cores"=nbcore)


saveres <- "saveres"
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE, recursive=TRUE) }


## breast cancer datasets
datasets <- read.csv(system.file(file.path("extdata", "datasets.csv"), package="MetaGx"), stringsAsFactors=FALSE)
rownames(datasets) <- as.character(datasets[ , "Dataset.ID"])

########################
## compendium of datasets with dmfs+rfs
########################

saveres2 <- file.path(saveres, "COMPENDIUM_ALL_HIPPO_test")
if(!file.exists(saveres2)) { dir.create(saveres2, showWarnings=FALSE, recursive=TRUE) }

## download data compendium
myfn <- file.path(saveres2, "datall.RData")
if (!file.exists(myfn)) {
  ddix <- rep(TRUE, nrow(datasets))
  datall <- MetaGx::getBrCaData(datasets=datasets[ddix, , drop=FALSE], resdir="cache", remove.duplicates=FALSE, topvar.genes=1000, duplicates.cor=0.98, sbt.model="scmod1", merging.method="union", merging.std="robust.scaling", nthread=8, verbose=TRUE)
  save(list="datall", compress=TRUE, file=myfn)
} else { load(myfn) }


## load sigHippo, weights are 1!!!
EntrezGene <- c(81788,9891,10413,25937)
weights <- rep(1, length(EntrezGene))
sig <- cbind(EntrezGene, weights)

ensembl.db <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

ss <- "entrezgene"
gid <- as.character(sig[ ,"EntrezGene"])
gene.an <- biomaRt::getBM(attributes=c(ss, "ensembl_gene_id", "hgnc_symbol", "unigene", "description", "chromosome_name", "start_position", "end_position", "strand", "band"), filters=ss, values=sort(unique(gid)), mart=ensembl.db)
gene.an[gene.an == "" | gene.an == " "] <- NA
gene.an <- gene.an[!is.na(gene.an[ , ss]) & !duplicated(gene.an[ , ss]) & is.element(gene.an[ , ss], gid), , drop=FALSE]
annot <- data.frame(matrix(NA, nrow=nrow(sig), ncol=ncol(gene.an), dimnames=list(gid, colnames(gene.an))))
annot[match(gene.an[ , ss], gid), colnames(gene.an)] <- gene.an
annot <- data.frame("EntrezGene"=gid, annot, "weight"=as.numeric(sig[ ,"weights"]))
sigHippo <- annot
save(list="sigHippo", compress=TRUE, file="sigHippo.rda")

## use sigHippo
sigHippo <- sigHippo[complete.cases(sigHippo),]
sigs <- data.frame(cbind(sigHippo[,"entrezgene"], sigHippo[,"weight"]))
sigs <- data.frame(paste("geneid", sigs[,1], sep="."), sigs)
colnames(sigs) <- c("probe", "EntrezGene.ID", "coefficient")


######### removing duplicate entrezgenes, only keeping the highest coefficient (otherwise duplicate row.names error occurs due to duplicates from multiple functions)
sigs <- do.call(rbind, 
                by(sigs, INDICES=list(sigs$probe), 
                   FUN=function(x) head(x[order(x$coefficient), ], 1)))
rownames(sigs) <- sigs[,"probe"]
sigs <- c(list("DONE_SIG"=sigs), sigs <- lapply(apply(sigs, 1, list), function (x) { return (data.frame("probe"=x[[1]]["probe"], "EntrezGene.ID"=as.character(as.numeric(x[[1]]["EntrezGene.ID"])), "coefficient"=as.numeric(x[[1]]["coefficient"]))) }))

ss <- NULL
for (i in 1:length(sigs)) {
  ss <- c(ss, list(data.frame("feature"=paste("geneid", as.character(sigs[[i]][complete.cases(sigs[[i]]), "EntrezGene.ID"]), sep="."), "coefficient"=sigs[[i]][complete.cases(sigs[[i]]), "coefficient"])))
}

names(ss) <- c(names(sigs))


## association wrt subtypes
rr <- MetaGx::subtypeAssociation(eset=datall$merged, sig=ss, plot=TRUE, weighted=FALSE, condensed=TRUE, resdir=saveres2, nthread=nbcore, method="weighted.average", scaling="robust")

## subtype-scpecific co-expression
rr <- MetaGx::subtypeCorrelation(eset=datall$merged, sig=ss, weighted=FALSE, condensed=TRUE, resdir=saveres2, nthread=nbcore, method="spearman", sig.method="weighted.average", sig.scaling="robust")

## prognostic value: dmfs
rr <- MetaGx::subtypeSurvival(eset=datall$merged, sig=ss, weighted=FALSE, surv.type="dmfs", time.cens=10 * 365, condensed=TRUE, resdir=saveres2, nthread=nbcore, sig.method="weighted.average", sig.scaling="robust")

