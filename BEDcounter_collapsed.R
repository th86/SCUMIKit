#Read Copy Counter and Aggregator for the Collapsed UMI-Tagged Reads in BED files
#Taihsien Ouyang
#2015 Dec. 10

fileNames=dir()

## Create the annotation vector ##
annot=NULL
for(i in 1:length(fileNames)){
		cat("Scanning file", i,"...\n" )
		bed<-read.table(fileNames[i])
		annot<-unique(c( annot, as.character(bed$V1) ))
}

## Read the UMI-collapsed BED files ##

geList<-list()
for(j in 1:length(fileNames)){ #
	cat("Processing",fileNames[j],"\n")
	bed<-read.table(fileNames[j])
  geList[[j]]<-table( as.character(bed$V1) )

}

## Aggregate the profiles ##

ge<-matrix(0,length(annot),length(fileNames))
rownames(ge)<-annot
colnames(ge)<-fileNames

for(i in 1:length(geList)){
  cat(i,"\n")
	geneList<-intersect(rownames(ge), names(geList[[i]]))
	ge[geneList,i]<-geList[[i]][geneList] 
}

cat("Writing outputs\n")
save(ge, file="count_matrix.rda")
write.table(ge, file="count_matrix.csv",sep=",")
