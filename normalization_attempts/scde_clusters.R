o.ifm <- scde.error.models(counts = cd, groups = groups, n.cores = 35, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
valid.cells <- o.ifm$corr.a > 0
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
> batch<-factor(batch, levels=c("JUNE5", "JUNE29"))
> ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  35, verbose  =  1, batch=batch)
#> write.table(ediff$batch.adjusted[order(abs(ediff$batch.adjusted$Z), decreasing = TRUE), ], file = "batch_adjusted_ediff_results_batchcluster12.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

> write.table(ediff$batch.adjusted[order(2*pnorm(-abs(ediff$batch.adjusted$Z)), -ediff$batch.adjusted$ce, decreasing=FALSE), ], file="batch_adjusted_ediff_results_batchcluster12.txt", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)


genes_2fold_p_lt_0.01<-ordered_ediff[abs(ordered_ediff$ce) >= 2 & 2*pnorm(-abs(ordered_ediff$Z)) < 0.01,]

> write.table(genes_2fold_p_lt_0.01, file="batch_adjusted_ediff_results_batchcluster_12_sig_genes_p_lt_0.01_2fold.txt", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
> batch_effect_only<-ediff$batch.effect[order(2*pnorm(-abs(ediff$batch.effect$Z)), -ediff$batch.effect$ce, decreasing=FALSE),]
> write.table(batch_effect_only, file="batch_effect_only_ediff_results_batchcluster_12.txt", row.names=TRUE, col.names=TRUE, quote=FALSE)

