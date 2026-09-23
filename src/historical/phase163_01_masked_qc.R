Sys.setlocale('LC_CTYPE','English_United States.utf8')
suppressPackageStartupMessages(library(limma))
O<-normalizePath(getwd(),winslash='/');R<-dirname(O);set.seed(1632026)
out<-function(x,n)write.csv(x,file.path(O,'RESULTS',paste0(n,'.csv')),row.names=FALSE,na='NA')
rawdir<-file.path(R,'162_GSE85358_AGILENT_HEADER_AND_TECHNICAL_BATCH_AUDIT_2026-09-17','raw')
files<-sort(list.files(rawdir,pattern='txt.gz$',full.names=TRUE));stopifnot(length(files)==48)
want<-c('FeatureNum','ProbeName','ControlType','gMedianSignal','gBGMedianSignal','gProcessedSignal','gIsSaturated','gIsWellAboveBG','gIsFeatNonUnifOL','gIsFeatPopnOL','IsManualFlag')
con<-gzfile(files[1],open='rt');hd<-readLines(con,n=10);close(con);heads<-strsplit(hd[10],'\t')[[1]];stopifnot(heads[1]=='FEATURES',all(want%in%heads))
classes<-rep('NULL',length(heads));classes[heads%in%want]<-'numeric';classes[heads=='ProbeName']<-'character'
mats<-list();identity<-NULL;counts<-list()
for(i in seq_along(files)){
 con<-gzfile(files[i],open='rt');d<-read.delim(con,skip=9,header=TRUE,colClasses=classes,check.names=FALSE,quote='',comment.char='');close(con)
 if(i==1){identity<-d[,c('FeatureNum','ProbeName','ControlType')];for(k in setdiff(want,names(identity)))mats[[k]]<-matrix(NA_real_,nrow(d),48)}
 stopifnot(identical(identity,d[,names(identity)]));for(k in names(mats))mats[[k]][,i]<-d[[k]]
 counts[[i]]<-data.frame(sample=paste0('S',sprintf('%02d',i)),GSM=sub('_.*','',basename(files[i])),rows=nrow(d),missing=sum(is.na(d)),foreground_nonpositive=sum(d$gMedianSignal<=0),processed_nonpositive=sum(d$gProcessedSignal<=0))
 cat('read',i,basename(files[i]),'\n');flush.console()
}
sid<-paste0('S',sprintf('%02d',1:48));for(k in names(mats))colnames(mats[[k]])<-sid
stopifnot(all(vapply(mats,function(x)all(is.finite(x)),logical(1))))
bio<-identity$ControlType==0
badpaper<-mats$gIsSaturated!=0|mats$gIsWellAboveBG==0
bad<-mats$gIsFeatNonUnifOL!=0|mats$gIsFeatPopnOL!=0|mats$IsManualFlag!=0
keep0<-bio & rowSums(badpaper)<24
keep3<-bio & rowSums(mats$gIsWellAboveBG!=0)>=24 & rowSums(mats$gIsSaturated!=0)<24 & rowSums(bad)<24
E1<-log2(pmax(mats$gMedianSignal[bio,],1));Eproc<-log2(pmax(mats$gProcessedSignal[keep0,],1))
bg<-new('EListRaw',list(E=mats$gMedianSignal[bio,],Eb=mats$gBGMedianSignal[bio,]));bg<-backgroundCorrect(bg,method='normexp',offset=50)
branches<-list(P0=normalizeBetweenArrays(Eproc,method='quantile'),P1=normalizeBetweenArrays(E1,method='quantile'),P2=normalizeBetweenArrays(log2(bg$E),method='quantile'),P3=normalizeBetweenArrays(log2(pmax(mats$gMedianSignal[keep3,],1)),method='quantile'))
mask<-list(P0=keep0,P1=bio,P2=bio,P3=keep3);probe<-list();audit<-list()
for(b in names(branches)){
 ids<-identity$ProbeName[mask[[b]]];x<-branches[[b]];sx<-rowsum(x,ids,reorder=FALSE);freq<-table(ids);probe[[b]]<-sx/as.numeric(freq[rownames(sx)])
 audit[[b]]<-data.frame(branch=b,input_features=nrow(identity),removed_controls=sum(!bio),removed_quality=sum(bio)-sum(mask[[b]]),retained_spots=sum(mask[[b]]),duplicate_spots_collapsed=sum(mask[[b]])-nrow(sx),unique_probes=nrow(sx),n_samples=48,nonfinite=sum(!is.finite(sx)))
}
# Blind technical QC uses P2 normalized probe matrix and intercept-only weights.
E<-probe$P2;pc<-prcomp(t(E),center=TRUE,scale.=FALSE);C<-cor(E);medcor<-apply(C,2,function(z)median(z[z<1]));distpc<-sqrt(rowSums(pc$x[,1:5,drop=FALSE]^2));aw<-arrayWeights(E,design=matrix(1,48,1))
flag<-medcor<median(medcor)-3*mad(medcor)|distpc>median(distpc)+3*mad(distpc)
qc<-data.frame(sample=sid,median_correlation=medcor,PC_distance=distpc,blind_array_weight=aw,review_flag=flag,main_analysis_action='RETAIN')
out(qc,'BLIND_QC_DECISIONS_LOCKED');out(do.call(rbind,counts),'RAW_READ_AUDIT');out(do.call(rbind,audit),'PREPROCESSING_COUNTS');out(data.frame(sample=sid,pc$x[,1:10]),'BLIND_PCA');out(data.frame(sample=sid,C),'SAMPLE_CORRELATION')
ctrl<-do.call(rbind,lapply(sort(unique(identity$ControlType)),function(k)data.frame(sample=sid,ControlType=k,n_features=sum(identity$ControlType==k),median_raw=apply(mats$gMedianSignal[identity$ControlType==k,,drop=FALSE],2,median))))
out(ctrl,'CONTROL_PROBE_QC')
png(file.path(O,'FIGURES','BLIND_DISTRIBUTIONS.png'),1800,1400,res=160);par(mfrow=c(2,2),mar=c(6,4,3,1));boxplot(mats$gMedianSignal[bio,],outline=FALSE,las=2,main='Raw foreground (masked labels)');boxplot(E1,outline=FALSE,las=2,main='log2 foreground before normalization');boxplot(E,outline=FALSE,las=2,main='P2 normalized probes');plot(density(E[,1]),main='P2 density by masked sample');for(i in 2:48)lines(density(E[,i]),col=i);dev.off()
png(file.path(O,'FIGURES','BLIND_STRUCTURE.png'),1800,900,res=150);par(mfrow=c(1,2));plot(pc$x[,1:2],pch=19,xlab='PC1',ylab='PC2',main='P2 masked PCA');text(pc$x[,1:2],labels=sid,pos=3,cex=.55);plot(hclust(as.dist(1-C)),main='P2 correlation clustering',xlab='',sub='');dev.off()
saveRDS(list(probe=probe,identity=identity,mask=mask,sample_key=do.call(rbind,counts),qc=qc,pc=pc,cor=C,control=ctrl,raw_qc=data.frame(sample=sid,raw_median=apply(mats$gMedianSignal[bio,],2,median),background_median=apply(mats$gBGMedianSignal[bio,],2,median),log_median=apply(E1,2,median),normalized_median=apply(E,2,median))),file.path(O,'RESULTS','PREPROCESSED.rds'))
writeLines(capture.output(sessionInfo()),file.path(O,'SESSION_INFO.txt'))
cat('MASKED_QC_COMPLETE: review flags',sum(flag),'no samples excluded\n')
