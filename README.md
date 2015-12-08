# SCUMIKit
A pipeline for processing UMI-tagged single cell RNA-Seq data

## Prerequisites
Exact UMI collapsing: [SRAToolkit](https://github.com/ncbi/sra-tools), [SAMtools](http://samtools.sourceforge.net/), [BEDtools](http://bedtools.readthedocs.org/en/latest/), [UMItools](https://github.com/brwnj/umitools/), [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) and [BWA](http://bio-bwa.sourceforge.net/)

UMI sequencing error correction: [CTK](http://zhanglab.c2b2.columbia.edu/index.php/CTK_Documentation), [CZPlib](http://sourceforge.net/projects/czplib/)

## Usage

Install the prerequisite programs and set up the paths in SCUMI.sh  
Then simply run 

```r
sh ./SCUMI.sh
```

## Reference Transcriptomes

The reference transcriptomes can be downloaded from [Ensembl](http://www.ensembl.org/info/data/ftp/index.html).
