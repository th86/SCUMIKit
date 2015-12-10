# SCUMIKit
A pipeline for processing UMI-tagged single cell RNA-Seq data

## Prerequisites
Read alignment and Exact UMI collapsing: [SRAToolkit](https://github.com/ncbi/sra-tools), [SAMtools](http://samtools.sourceforge.net/), [HomerTools](http://homer.salk.edu/homer/ngs/homerTools.html) ,[BEDtools](http://bedtools.readthedocs.org/en/latest/), [UMItools](https://github.com/brwnj/umitools/), [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) and [BWA](http://bio-bwa.sourceforge.net/)

UMI sequencing error correction: [CTK](http://zhanglab.c2b2.columbia.edu/index.php/CTK_Documentation), [CZPlib](http://sourceforge.net/projects/czplib/)

Matrix Imputation: [impute](www.bioconductor.org/packages/release/bioc/html/impute.html)

Multicore Computation: [doMC](https://cran.r-project.org/web/packages/doMC/index.html) and [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/00Index.html)

## Usage

Install the prerequisite programs and set up the paths in SCUMI.sh  
Then simply run 

```r
sh ./SCUMI.sh
```

A directory with resulting BED files will be created in the designated directory. Run the following command to create the count matrix from the file using multiple processors. The values after the R script is the number of cores.

```r
Rscript /path/to/BEDcounter.R 12
```

## Reference Transcriptomes

The reference transcriptomes can be downloaded from [Ensembl](http://www.ensembl.org/info/data/ftp/index.html).

## References

Zhang, C., Darnell, R.B. 2011. Mapping in vivo protein-RNA interactions at single-nucleotide resolution from HITS-CLIP data. Nat. Biotech. 29:607-614. 

## Acknowledgements

We are grateful for the advice and help of Professor Chaolin Zhang, Professor Peter Sims, and Professor Yufeng Shen during the development of this toolkit as our final project of BINFG4017 DEEP SEQUENCING at Columbia University.