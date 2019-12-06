
<!-- README.md is generated from README.Rmd. Please edit that file -->



[![Build Status](https://travis-ci.org/mskilab/fragCounter.svg?branch=master)](https://travis-ci.org/mskilab/fragCounter) [![codecov.io](https://img.shields.io/codecov/c/github/mskilab/fragCounter.svg)](https://codecov.io/github/mskilab/fragCounter?branch=master)

fragCounter
===========

The goal of fragCounter is to correct Whole genome or targeted sequencing data for GC and mappability bias.                                                                                                                                                                                                                 
The GC bias curve is determined by loess regression of read count by GC and mappability scores.                                                                                                                                                                                                                             
Segmentation is done by circular binary segmentation (CBS) algorithm after getting tumor/normal ratios of corrected read counts. 

Installation (R package)                                                                                                                                                                                                                                                                                                    
------------                                                                                                                                                                                                                                                                                                                
1. Install devtools from CRAN (if you don't have it already)                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                            
```{r}                                                                                                                                                                                                                                                                                                                      
install.packages('devtools')                                                                                                                                                                                                                                                                                                
```                                                                                                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                                                                                                            
2. Install fragCounter                                                                                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                                                                                                            
```{r}                                                                                                                                                                                                                                                                                                                      
devtools::install_github('mskilab/fragCounter')                                                                                                                                                                                                                                                                              
```   

Example                                                                                                                                                                                                                                                                                                    
------------  
(after installing R package) Add fragCounter directory to PATH and test the executable 

```{bash}
$ export PATH=${PATH}:$(Rscript -e 'cat(paste0(installed.packages()["fragCounter", "LibPath"], "/fragCounter/extdata/"))')
$ frag -h ## to see the help message
```

```{bash}
$ ./frag -b inst/extdata/chr21.bam -d inst/extdata/gcMAP21/ -w 200  

Rprofile Loading                                                                                                                                                                                                                                                                                                            
Rprofile Finished Loading                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                                                            
███████╗██████╗  █████╗  ██████╗  ██████╗ ██████╗ ██╗   ██╗███╗   ██╗████████╗███████╗██████╗                                                                                                                                                                                                                               
██╔════╝██╔══██╗██╔══██╗██╔════╝ ██╔════╝██╔═══██╗██║   ██║████╗  ██║╚══██╔══╝██╔════╝██╔══██╗                                                                                                                                                                                                                              
█████╗  ██████╔╝███████║██║  ███╗██║     ██║   ██║██║   ██║██╔██╗ ██║   ██║   █████╗  ██████╔╝                                                                                                                                                                                                                              
██╔══╝  ██╔══██╗██╔══██║██║   ██║██║     ██║   ██║██║   ██║██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗                                                                                                                                                                                                                              
██║     ██║  ██║██║  ██║╚██████╔╝╚██████╗╚██████╔╝╚██████╔╝██║ ╚████║   ██║   ███████╗██║  ██║                                                                                                                                                                                                                              
╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝  ╚═════╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                                                            
Calling samtools view -f 0x02 -F 0x10 inst/extdata/chr21.bam -q 1 | cut -f "3,4,9" 
Starting fragment count on inst/extdata/chr21.bam and min mapQ 1 and   insert size limit 10000  
Finished computing coverage, and making GRanges  
Finished acquiring coverage 
Loaded GC and mappability  
length cov is 314827, length gc is 3663, length map is 3663   
Synced coverage, GC, and mappability   
Modified gc / mappability correction 
Converting to data.table 
Converting to GRanges                                                                                                                                                                                                                                                                                                       
Made GRanges

```



Usage (frag executable)                                                                                                                                                                                                                                                                                                     
------------

```{bash}
$ ./frag -h

Rprofile Loading                                                                                                                                                                                                                                                                                                            
Rprofile Finished Loading                                                                                                                                                                                                                                                                                                   
Usage: ./frag [options]                                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                                                            
Options:                                                                                                                                                                                                                                                                                                                    
        -b BAM, --bam=BAM                                                                                                                                                                                                                                                                                                   
                Path to .bam file                                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                                                            
        -c COV, --cov=COV                                                                                                                                                                                                                                                                                                   
                Path to existing coverage rds or bedgraph                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                            
        -m MIDPOINT, --midpoint=MIDPOINT                                                                                                                                                                                                                                                                                    
                If TRUE only count midpoint if FALSE then count bin footprint of every fragment interval                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                            
        -w WINDOW, --window=WINDOW                                                                                                                                                                                                                                                                                          
                Window / bin size                                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                                                            
        -d GCMAPDIR, --gcmapdir=GCMAPDIR                                                                                                                                                                                                                                                                                    
                Mappability / GC content dir                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                            
        -q MINMAPQ, --minmapq=MINMAPQ                                                                                                                                                                                                                                                                                       
                Minimal map quality                                                                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                                                                                                            
        -p PAIRED, --paired=PAIRED                                                                                                                                                                                                                                                                                          
                Is paired                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                            
        -e EXOME, --exome=EXOME                                                                                                                                                                                                                                                                                             
                Use exons as bins instead of fixed window                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                            
        -u USE.SKEL, --use.skel=USE.SKEL                                                                                                                                                                                                                                                                                    
                Use user defined regions instead of default exome skeleton                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                            
        -s SKELETON, --skeleton=SKELETON                                                                                                                                                                                                                                                                                    
                Path to skeleton file                                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                            
        -o OUTDIR, --outdir=OUTDIR                                                                                                                                                                                                                                                                                          
                Directory to dump output into                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                            
        -l LIBDIR, --libdir=LIBDIR                                                                                                                                                                                                                                                                                          
                Directory containing this R file                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                                                            
        -h, --help                                                                                                                                                                                                                                                                                                          
                Show this help message and exit
                
```



GC and Mappability files for 1 kb WGS (hg19)                                                                                                                                                                                                                                                                                                     
------------
These can be used as inputs to fragCounter to correct GC and Mappability biases. Make sure to save them in same directory and supplying the path to the directory. Retain the name of the file.

GC file: http://mskilab.com/dryclean/hg19/gc1000.rds


Mappability file: http://mskilab.com/dryclean/hg19/map1000.rds
