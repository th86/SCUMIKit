# SCUMIKit
A streamlined pipeline for processing UMI-tagged single cell RNA-Seq data. The pipeline does the following workflow in three simple commands.

![alt tag](https://raw.githubusercontent.com/th86/SCUMIKit/master/doc/SCUMIpipeline.png)

## Prerequisites
Read alignment and exact UMI collapsing: [SRAToolkit](https://github.com/ncbi/sra-tools), [SAMtools](http://samtools.sourceforge.net/), [HomerTools](http://homer.salk.edu/homer/ngs/homerTools.html) ,[BEDtools](http://bedtools.readthedocs.org/en/latest/), [UMItools](https://github.com/brwnj/umitools/), [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) and [BWA](http://bio-bwa.sourceforge.net/)

UMI sequencing error correction: [CTK](http://zhanglab.c2b2.columbia.edu/index.php/CTK_Documentation), [CZPlib](http://sourceforge.net/projects/czplib/)

Matrix imputation: [impute](www.bioconductor.org/packages/release/bioc/html/impute.html)

Multicore computing: [doMC](https://cran.r-project.org/web/packages/doMC/index.html) and [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/00Index.html)


Before running the pipeline, we need to prepare the index for the aligners.

For BWA, run the following command:
```r
bwa index -a bwtsw Mus_musculus.GRCm38.cdna.all.fa
```

For Bowtie:
```r
bowtie-build Mus_musculus.GRCm38.cdna.all.fa GRCm38
```

## Usage

Install the prerequisite programs and set up the paths in SCUMI.sh  
Then simply run 

```r
sh ./SCUMI.sh
```

A directory with resulting BED files will be created in the designated directory. Please run the following command to create the count matrix from the UMI-collapsed BED file in the directory. The output is an R object file (count_matrix.rda) and a CSV (count_matrix.csv) of the count matrix.

```r
Rscript /path/to/BEDcounter_collapsed.R
```

If you want to collapse the file using multiple processors without CTK or UMItools. Please use the following command. The values after the R script is the number of cores.

```r
Rscript /path/to/BEDcounter.R 12
```

To correct the count matrix using the k-NN-based single cell expression noise model, please run the follows. The output is an R object file (count_matrix_imputed.rda).

```r
Rscript /path/to/NoiseCorrection.R count_matrix.rda
```

## Reference Transcriptomes

The reference transcriptomes can be downloaded from [Ensembl](http://www.ensembl.org/info/data/ftp/index.html).

## References

Bose, S., Wan Z., Carr A., Rizvi, A.H. et al. Scalable microfluidics for single-cell RNA printing and sequencing. Genome Biol. 16:120.

Islam, S., Zeisel, A., Joost, S., La Manno G. et. al. Quantitative single-cell RNA-seq with unique molecular identifiers. Nat. Methods 11:163-166.

Zhang, C., Darnell, R.B. 2011. Mapping in vivo protein-RNA interactions at single-nucleotide resolution from HITS-CLIP data. Nat. Biotech. 29:607-614. 

## Acknowledgements

We are grateful for the advice and help of Professor Chaolin Zhang, Professor Peter Sims, and Professor Yufeng Shen during the development of this toolkit as our final project of BINFG4017 DEEP SEQUENCING at Columbia University.