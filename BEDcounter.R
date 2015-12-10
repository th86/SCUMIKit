#Multi-Cored Counter and Aggregator of the Mapped Reads in BED files
#Taihsien Ouyang
#2015 Nov 9

library("doMC")

args <- commandArgs(TRUE)

cat("== Multi-Cored BED Read Counter ==\n")

## Passing arguments to R ##

if(length(args)>0){
	NUM_CORE=as.numeric(args[1])
	cat("Set number of cores =", NUM_CORE ,"\n")
}else{
	NUM_CORE=12
}

registerDoMC(NUM_CORE)

fileNames=dir()

## Create the annotation vector ##
annot=NULL
for(i in 1:length(fileNames)){
		cat("Scanning file", i,"...\n" )
		bed<-read.table(fileNames[i])
		annot<-unique(c( annot, as.character(bed$V1) ))
}

## Read the BED files in parallel ##

geList<-foreach(j = 1:length(fileNames)) %dopar% {

	cat("Processing file", j,"...\n" )
	bed<-read.table(fileNames[j])
	cat("creating UMI table", j,"...\n")
	umiList<-rep("",nrow(bed))
	for(i in 1:nrow(bed)){
		umiList[i]=strsplit(as.character(bed$V4[i]), "_")[[1]][2]
	}

	geneList<-as.character(bed$V1)
	geneNames<-unique(geneList)
	g<-rep(0,length(annot))
	names(g)<-annot
	cat("Summarizing", j,"...\n")
	for(i in 1:length(annot)){
		g[annot[i]] =length(unique( umiList[ which(geneList==geneNames[i]) ] ))
	}
	cat("Done", j,"\n")
	return(g)
}

## Aggregate the profiles ##

ge<-matrix(0,length(annot),length(fileNames))
rownames(ge)<-annot
colnames(ge)<-fileNames

cat("Aggregating expression profiles...")
for(i in 1:length(geList)){
	ge[,i]<-geList[[i]]
}

cat("Writing outputs\m")
save(ge, file="count_matrix.rda")
write.table(ge, file="count_matrix.csv",sep=",")
cat("Done\n")
