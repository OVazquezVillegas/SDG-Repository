file ='NL'  
output = 'NL_SDG'
Thesauro = 'keywords_ramirez_v3-2'
cores =5

#######################################
#######################################
#######################################
#######################################
#######################################
Search2=function(j){
  if(words[j,'Mode']=='stemmed'){
    if(words[j,'type.1']==1){
      if(length(grep(paste("\\b", tolower(words[j,'stemmed']), "\\b", sep=""),datas_stemmed))!=0){
        paste(grep(paste("\\b", tolower(words[j,'stemmed']), "\\b", sep=""),datas_stemmed),j,
              unlist(lapply(datas_stemmed[grep(paste("\\b", tolower(words[j,'stemmed']), "\\b", sep=""),datas_stemmed)],
                            function(k){
                              length(unlist(gregexpr(paste("\\b", tolower(words[j,'stemmed']), "\\b", sep=""),k)))})),sep='#')}
    }
    else{
      x=table(unlist(lapply(unlist(strsplit(words[j,'stemmed'],'#')),
                            function(h){grep(paste("\\b",tolower(h),"\\b",sep=""),datas_stemmed)})))
      if(length(x[x==length(unlist(strsplit(words[j,'stemmed'],'#')))])!=0){
        paste(names(x[x==length(unlist(strsplit(words[j,'stemmed'],'#')))]),j,
              apply(matrix(unlist(lapply(unlist(strsplit(words[j,'stemmed'],'#')),
                                         function(h){unlist(lapply(datas_stemmed[as.numeric(names(x[x==length(unlist(strsplit(words[j,'stemmed'],'#')))]))],
                                                                   function(k){length(unlist(gregexpr(paste("\\b",h,"\\b",sep=""),k)))}))})),
                           byrow=T,ncol=length(unlist(strsplit(words[j,'stemmed'],'#')))),1,min),sep='#')}}}
  else{
    if(words[j,'type.1']==1){
      if(length(grep(paste('\\b',tolower(words[j,'edit_JV']),sep=""),datas))!=0){
        paste(grep(paste('\\b',tolower(words[j,'edit_JV']),sep=""),datas),j,
              unlist(lapply(datas[grep(paste('\\b',tolower(words[j,'edit_JV']),sep=""),datas)],
                            function(k){
                              length(unlist(gregexpr(paste('\\b',tolower(words[j,'edit_JV']),sep=""),k)))})),sep='#')}}
    else{
      x=table(unlist(lapply(unlist(strsplit(words[j,'edit_JV'],'#')),
                            function(h){grep(paste('\\b',tolower(h),sep=""),datas)})))
      if(length(x[x==length(unlist(strsplit(words[j,'edit_JV'],'#')))])!=0){
        paste(names(x[x==length(unlist(strsplit(words[j,'edit_JV'],'#')))]),j,
              apply(matrix(unlist(lapply(unlist(strsplit(words[j,'edit_JV'],'#')),
                                         function(h){unlist(lapply(datas[as.numeric(names(x[x==length(unlist(strsplit(words[j,'edit_JV'],'#')))]))],
                                                                   function(k){length(unlist(gregexpr(paste('\\b',h,sep=""),k)))}))})),
                           byrow=T,ncol=length(unlist(strsplit(words[j,'edit_JV'],'#')))),1,min),sep='#')}}}
}

library(textstem);library(tm);library(parallel)
################
################
################
data=as.matrix(read.csv(paste(file,".csv",sep="")))
words=as.matrix(read.csv(paste(Thesauro,'.csv',sep='')))
words[,'edit_JV']=unlist(lapply(words[,'edit_JV'], function(x) gsub('”', '\\\\b', x)))
words[,'edit_JV']=unlist(lapply(words[,'edit_JV'], function(x) gsub('\\*', '', x)))
data=data[!data[,'PY']%in%2020:2021,]
datas=datas1=data[,c('TI','AB','DE')]
data[,'UT']=substring(data[,'UT'],1,30)
datas=gsub("  "," ",gsub('[[:digit:]]+', '', gsub('  ',' ',gsub('[[:punct:] ]+',' ',tolower(unlist(mclapply(1:dim(datas)[1],function(i){paste(datas[i,1],datas[i,2],datas[i,3])},mc.cores=cores)))))))
datas=gsub("  "," ",unlist(mclapply(datas,function(i){paste(unlist(strsplit(i,' '))[!unlist(strsplit(i,' '))%in%stopwords()],collapse=" ")},mc.cores=cores)))
datas_stemmed=gsub("  "," ",unlist(mclapply(datas,function(i){paste(stem_words(unlist(strsplit(i,' '))[!unlist(strsplit(i,' '))%in%stopwords()]),collapse=" ")},mc.cores=cores)))

search=matrix(unlist(strsplit(unlist(mclapply(1:dim(words)[1],Search2,mc.cores=4)),"#")),byrow=T,ncol=3)
m=matrix(unlist(strsplit(unlist(lapply(1:dim(search)[1],function(i){paste(data[as.numeric(search[i,1]),'UT'], words[as.numeric(search[i,2]),2],words[as.numeric(search[i,2]),1],words[as.numeric(search[i,2]),3],search[i,3],sep='##')})),'##')),byrow=T,ncol=5)
colnames(m)=c('UT','KeyWord','SDG','Type','Freq')
write.csv(m,paste(output,"Words_Searched.csv",sep='_'),row.names=F)

datas=xtabs(as.numeric(m[,'Freq'])~m[,'UT']+m[,'SDG'])
datas=cbind(datas,rowSums(datas))
colnames(datas)=c(paste('SDG', 1:17, sep = '.'),'FreqT')
class(datas)='numeric'
data=data[data[,'UT']%in%rownames(datas),]
datas=datas[unlist(mclapply(1:dim(data)[1],function(i){which(rownames(datas)==data[i,'UT'])},mc.cores=cores)),]
data=data.frame(data,datas)
SDG=unlist(mclapply(1:dim(data)[1],function(i){ifelse(sum(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T)[1])>.75,                                                   paste(names(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T))[1],collapse='#'),                                                      ifelse(sum(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T)[1:2])>.60,                                                             paste(names(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T))[1:2],collapse='#'),                                                           ifelse(sum(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T)[1:3])>.50,                                                            paste(names(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T))[1:3],collapse='#'),                                                                   paste(names(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T))[1:3],collapse='#'))))},mc.cores=cores))

data=cbind(data,SDG);
write.csv(data,paste('NL_s','.csv',sep=''),row.names=F)
