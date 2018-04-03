ample directory with per-sample quantifications
expression_file_dir<-"/broad/hptmp/ht_for_beril/results/"

expression_list<-list.files(expression_file_dir,full.names=T)

datalist = lapply(expression_list, function(x){
  tmp<-read.delim(file=x,header=F,stringsAsFactors=F)
  sample_name = gsub("^/broad/hptmp/ht_for_beril/results//","",x)
  sample_name  = gsub(".expression.txt$","",sample_name)
  names(tmp)<-c("transcript",sample_name)
  return(tmp) } )


y<-datalist[[1]]
expression_matrix<-Reduce(function(x,y) {merge(x,y,by="transcript")}, datalist)

mappings = read.delim("mappingKnownIsoform.txt",stringsAsFactors = F,header=F,col.names=c("id",'enst',"ensg"))
for(i in 1:nrow(expression_matrix)){
  transcript = expression_matrix[i,"transcript"]
  one_transcript = strsplit(transcript,split="\\|")[[1]][1]
  gene = mappings[mappings$enst==one_transcript,"ensg"] 
  gene = strsplit(gene,split="\\.")[[1]][1]
  expression_matrix[i,'Gene']<-gene}


rownames(expression_matrix)<-expression_matrix$Gene
expression_matrix$Gene<-NULL
expression_matrix$transcript<-NULL

write.table(expression_matrix,file="final.gencodev10.expression.matrix.txt",col.names=T,row.names=T,quote=F,sep="\t")
