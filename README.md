# SCUMIKit
A pipeline for processing UMI-tagged single cell RNA-Seq data

## Prerequisites
Exact UMI collapsing: SRAToolkit, SAMTools, BEDTools, UMItools, Bowtie or BWA
UMI sequencing error correction: CTK

## Usage

Install the prerequisite programs and set up the paths in SCUMI_Exact.sh. 
Then simply run 

```r
sh ./SCUMI_Exact.sh
```

## Reference Transcriptomes

The reference transcriptomes can be downloaded from [Ensembl](http://www.ensembl.org/info/data/ftp/index.html).
