#kNN-Based Noise Correction Model for Single Cell RNA-Seq data
#Taihsien Ouyang
#2015 Dec. 9

library("impute")

args <- commandArgs(TRUE)

cat("== kNN-Based Noise Correction of Single Cell RNA-Seq data ==\n")

## Passing arguments to R ##

if(length(args)>0){
	load(args[1])
}else{
	stop("Please provide the count matrix file.\n")
}

## Estimate the parameters of Poisson Model ##

ge[is.na(ge)]<-0

mu_all<-apply(ge,1,mean)
sd_all<-apply(ge,1,sd)
mu_all_sorted<-sort(mu_all)

## Noise Filtering ##
## The Poisson model is based on doi:10.1038/nmeth.2772

noisygenes<-names( which(sd_all>(3.7*(mu_all)^0.5+0.3)) )

for( i in 1:length(noisygenes)){
	gValues<-ge[noisygenes[i],]
	while(length(gValues)>2){
		gValues=gValues[-which.max(gValues)]

		if(sd(gValues)<=(3.7*mean(gValues)^0.5+0.3) ){
			cat(noisygenes[i], ": remove" , as.numeric(ncol(ge) - length(gValues)) , "value(s) out of", ncol(ge), " values.\n") #Usually only remove 1 value
			ge[noisygenes[i], setdiff(colnames(ge),names(gValues))]=NA
			break
		}
		
	}

}

## Impute the matrix ##
cat("Imputing the matrix\n")
ge.impute<-impute.knn(ge)$data

cat("Writing results\n")
save(ge.impute, file="count_matrix_imputed.rda")

cat("== Done ==\n")