Sys.setlocale('LC_CTYPE','English_United States.utf8');suppressPackageStartupMessages(library(limma));set.seed(1632026)
O<-normalizePath(getwd(),winslash='/');R<-dirname(O);D<-file.path(O,'RESULTS');out<-function(x,n)write.csv(x,file.path(D,paste0(n,'.csv')),row.names=FALSE,na='NA')
z<-readRDS(file.path(D,'PREPROCESSED.rds'));stopifnot(file.exists(file.path(D,'BLIND_QC_DECISIONS_LOCKED.csv')))
m<-read.csv(file.path(R,'162_GSE85358_AGILENT_HEADER_AND_TECHNICAL_BATCH_AUDIT_2026-09-17','SAMPLE_LEVEL_TECHNICAL_METADATA_Samples.csv'),check.names=FALSE)
m<-m[match(z$sample_key$GSM,m$GSM),];stopifnot(!anyNA(m$GSM),!anyDuplicated(m$GSM),length(unique(m$donor_id))==48)
m$group<-factor(m$frozen_age_group,levels=c('young','old'));m$scan<-factor(m$scan_date_batch);m$array<-factor(m$array_or_slide_id);stopifnot(all(table(m$group)==24))
mp<-read.csv(file.path(R,'131_ANNOTATION_ADJUDICATION_FREEZE_2026-09-09','csv/GSE85358_FINAL_MAPPING.csv'),check.names=FALSE)
sel<-mp[tolower(as.character(mp$selected_representative))=='true'&nzchar(mp$final_HGNC_ID),];stopifnot(!anyDuplicated(sel$final_HGNC_ID));out(sel,'FROZEN_REPRESENTATIVES_USED')
gene<-list();cnt<-list();probeDE<-list();geneDE<-list();summ<-list();designaudit<-list()
fit<-function(E,X,coef,label,weights=NULL){
 stopifnot(qr(X)$rank==ncol(X),all(is.finite(E)));h<-rowSums(qr.Q(qr(X))^2)
 designaudit[[label]]<<-data.frame(analysis=label,n=nrow(X),p=ncol(X),rank=qr(X)$rank,residual_df=nrow(X)-ncol(X),max_leverage=max(h),condition_number=kappa(X),contrast=colnames(X)[coef])
 ff<-eBayes(lmFit(E,X,weights=weights));tt<-topTable(ff,coef=coef,number=Inf,sort.by='none',adjust.method='BH');tt$id<-rownames(tt);tt$SE<-ff$stdev.unscaled[,coef]*sqrt(ff$s2.post)
 summ[[label]]<<-data.frame(analysis=label,n_samples=ncol(E),tested=nrow(tt),FDR005=sum(tt$adj.P.Val<.05),up=sum(tt$adj.P.Val<.05&tt$logFC>0),down=sum(tt$adj.P.Val<.05&tt$logFC<0),paper_rule=sum(tt$adj.P.Val<.01&abs(tt$logFC)>.25),paper_up=sum(tt$adj.P.Val<.01&tt$logFC>.25),paper_down=sum(tt$adj.P.Val<.01&tt$logFC< -.25))
 if(grepl('_age_',label))summ[[label]][,c('paper_rule','paper_up','paper_down')]<<-NA
 tt
}
X<-model.matrix(~group,m)
for(b in names(z$probe)){
 E<-z$probe[[b]];probeDE[[b]]<-fit(E,X,2,paste0(b,'_probe_M0'));out(probeDE[[b]],paste0('PROBE_',b,'_M0'))
 ss<-sel[sel$original_id%in%rownames(E),];G<-E[match(ss$original_id,rownames(E)),,drop=FALSE];rownames(G)<-ss$final_HGNC_ID;gene[[b]]<-G
 tt<-fit(G,X,2,paste0(b,'_gene_M0'));tt$gene_symbol<-ss$final_approved_symbol;tt$representative_probe<-ss$original_id;geneDE[[b]]<-tt;out(tt,paste0('GENE_',b,'_M0'))
 cnt[[b]]<-data.frame(branch=b,retained_unique_probes=nrow(E),frozen_representatives=nrow(sel),mapped_genes=nrow(G),frozen_representatives_missing=nrow(sel)-nrow(G),nonrepresentative_or_outside_scope=nrow(E)-nrow(G))
}
G<-gene$P2;sens<-list(M0=geneDE$P2)
sens$M1<-fit(G,model.matrix(~group+scan,m),2,'P2_gene_M1')
if(all(table(m$array)>=2)&&all(table(m$array,m$group)>0)&&qr(model.matrix(~group+array,m))$rank==ncol(model.matrix(~group+array,m)))sens$M2<-fit(G,model.matrix(~group+array,m),2,'P2_gene_M2')
aw<-arrayWeights(G,design=X);out(data.frame(GSM=m$GSM,age_group=m$group,weight=aw),'MODEL_ARRAY_WEIGHTS');sens$M3<-fit(G,X,2,'P2_gene_M3',aw)
keep<-!z$qc$review_flag
if(any(!keep)&&sum(keep)>=20&&all(table(m$group[keep])>=8))sens$M4<-fit(G[,keep],X[keep,,drop=FALSE],2,'P2_gene_M4_EXPLORATORY_EXCLUSION')
for(n in names(sens))out(sens[[n]],paste0('GENE_P2_',n))
cont<-list(GEO=as.numeric(m$exact_age_GEO),SUPPLEMENT=as.numeric(m$exact_age_supplement));ages<-list()
for(n in names(cont)){xx<-model.matrix(~cont[[n]]);colnames(xx)<-c('Intercept','age_per_year');ages[[n]]<-fit(G,xx,2,paste0('P2_age_',n));out(ages[[n]],paste0('AGE_',n))}
compare<-function(a,b,label){ii<-intersect(a$id,b$id);a<-a[match(ii,a$id),];b<-b[match(ii,b$id),];sa<-a$id[a$adj.P.Val<.05];sb<-b$id[b$adj.P.Val<.05];un<-union(sa,sb)
 row<-data.frame(comparison=label,n_common=length(ii),logFC_pearson=cor(a$logFC,b$logFC),t_pearson=cor(a$t,b$t),P_rank_spearman=cor(a$P.Value,b$P.Value,method='spearman'),direction_agreement=mean(sign(a$logFC)==sign(b$logFC)),significant_overlap=length(intersect(sa,sb)),significant_union=length(un),Jaccard=if(length(un))length(intersect(sa,sb))/length(un) else NA)
 for(k in c(100,500,1000))row[[paste0('top',k,'_overlap')]]<-length(intersect(head(a$id[order(a$P.Value,a$id)],k),head(b$id[order(b$P.Value,b$id)],k)))
 row
}
comparisons<-do.call(rbind,c(lapply(setdiff(names(geneDE),'P2'),function(n)compare(geneDE$P2,geneDE[[n]],paste0('P2_vs_',n))),lapply(setdiff(names(sens),'M0'),function(n)compare(sens$M0,sens[[n]],paste0('M0_vs_',n)))))
out(comparisons,'SENSITIVITY_COMPARISON');out(compare(ages$GEO,ages$SUPPLEMENT,'GEO_vs_SUPPLEMENT_CONTINUOUS_AGE'),'AGE_SOURCE_COMPARISON')
prior<-read.csv(file.path(R,'22_phase2_5_results/GSE85358_GENE_LEVEL_ALL_RESULTS.csv'));old<-data.frame(id=prior$gene_symbol,logFC=prior$logFC_old_minus_young,t=prior$t_statistic,P.Value=prior$p_value,adj.P.Val=prior$FDR)
priorcomp<-do.call(rbind,lapply(names(geneDE),function(b){tt<-geneDE[[b]];tt$id<-tt$gene_symbol;compare(old,tt,paste0('LEGACY_vs_',b))}));out(priorcomp,'PRIOR_COMPARISON');out(data.frame(tested=nrow(old),significant=sum(old$adj.P.Val<.05)),'PRIOR_COUNTS_VERIFIED')
out(do.call(rbind,summ),'MODEL_COUNTS');out(do.call(rbind,designaudit),'DESIGN_AUDIT');out(do.call(rbind,cnt),'GENE_MAPPING_COUNTS');out(m,'SAMPLE_KEY_LABELLED')
# Descriptive PC associations; never used to select covariates.
assoc<-list();vars<-list(age_group=m$group,scan_date=m$scan,exact_scan_time=factor(m$scan_datetime),ArrayName=m$array,QC_DetectionLimit=as.numeric(m$QC__Metric_DetectionLimit),QC_NegCtrl=as.numeric(m$QC__gNegCtrlAveNetSig))
for(k in 1:5)for(v in names(vars)){xx<-vars[[v]];yy<-z$pc$x[,k];ff<-lm(yy~xx);ss<-summary(ff);fs<-ss$fstatistic;assoc[[paste(k,v)]]<-data.frame(PC=k,variable=v,R_squared=ss$r.squared,adjusted_R_squared=ss$adj.r.squared,P=pf(fs[1],fs[2],fs[3],lower.tail=FALSE),df1=fs[2],df2=fs[3])}
aa<-do.call(rbind,assoc);aa$BH<-p.adjust(aa$P,'BH');out(aa,'PCA_ASSOCIATIONS')
qc<-cbind(GSM=m$GSM,group=m$group,z$qc);qc$group_centroid_distance<-NA_real_
for(g in levels(m$group)){idx<-which(m$group==g);cent<-rowMeans(G[,idx]);qc$group_centroid_distance[idx]<-sqrt(colMeans((G[,idx]-cent)^2))};out(qc,'SAMPLE_QC_DECISIONS');out(z$raw_qc,'RAW_NORMALIZED_QC')
for(v in names(vars)[1:4]){png(file.path(O,'FIGURES',paste0('PCA_',v,'.png')),1500,1000,res=150);cc<-as.integer(factor(vars[[v]]));plot(z$pc$x[,1:2],col=cc,pch=19,xlab='PC1',ylab='PC2',main=paste('P2 PCA:',v));legend('topright',legend=levels(factor(vars[[v]])),col=seq_along(unique(cc)),pch=19,cex=.6);dev.off()}
tt<-geneDE$P2
png(file.path(O,'FIGURES','P2_M0_VOLCANO_MA.png'),1800,900,res=150);par(mfrow=c(1,2));plot(tt$logFC,-log10(tt$P.Value),pch=16,cex=.35,col=ifelse(tt$adj.P.Val<.05,'firebrick','grey'),xlab='old - young log2FC',ylab='-log10 raw P');plot(tt$AveExpr,tt$logFC,pch=16,cex=.35,col=ifelse(tt$adj.P.Val<.05,'firebrick','grey'),xlab='Average log2 expression',ylab='old - young log2FC');dev.off()
top<-head(order(tt$P.Value,tt$id),50);hm<-G[top,];rownames(hm)<-tt$gene_symbol[top];png(file.path(O,'FIGURES','P2_TOP50.png'),1600,1600,res=150);heatmap(hm,scale='row',ColSideColors=ifelse(m$group=='old','firebrick','steelblue'),cexRow=.6,cexCol=.5);dev.off();out(data.frame(gene=rownames(hm),hm),'TOP50_HEATMAP_SOURCE')
# Frozen member maps only, no new enrichment tests.
pm<-read.csv(file.path(R,'131_ANNOTATION_ADJUDICATION_FREEZE_2026-09-09/csv/PATHWAY_MEMBER_MAPPING.csv'),check.names=FALSE);pm<-pm[pm$member_mapping_status=='RESOLVED_UNIQUE',];stab<-list()
for(b in names(geneDE)){t<-geneDE[[b]];ii<-match(pm$HGNC_candidates,t$id);ok<-!is.na(ii);stab[[b]]<-data.frame(branch=b,object=pm$object[ok],pathway=pm$pathway[ok],HGNC=pm$HGNC_candidates[ok],gene=t$gene_symbol[ii[ok]],logFC=t$logFC[ii[ok]],t=t$t[ii[ok]],FDR=t$adj.P.Val[ii[ok]])}
out(do.call(rbind,stab),'PATHWAY_GENE_DIRECTIONS')
saveRDS(list(gene=gene,geneDE=geneDE,probeDE=probeDE,sensitivity=sens,ages=ages,metadata=m),file.path(D,'ANALYSIS_OBJECTS.rds'))
writeLines(capture.output(sessionInfo()),file.path(O,'SESSION_INFO.txt'));cat('MODELS_COMPLETE\n');print(do.call(rbind,summ))
