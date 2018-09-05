# npInv
npInv: an accurate tool for detecting and genotyping inversion using multiple alignment long reads

# Installation
Download the runnable jar such as npInv1.2.jar.

# Requirement
Java version 1.8+.

# Usage
```
Program function: Read a SE bam file and get the inversion
--output[String] file to write
--input[String] file to read
optional:
--region[String] Specify the region for running.
                 Such as chr9:1-1000 OR chr9 OR all Default[all]
				 --minAln[int] minimum size for Alignment & Inv. Default[500]
				 --IRdatabase[String] An inverted repeat file for the reference in bed format. Default[none]
				 --min[int] minimum size of inversion. Default[500]
				 --max[int] maximum size of inversion. Default[10000]
				 For example: java -jar npInv.jar --input sample.bam --output sample.VCF
```

# Contact
Email: uqhshao@uq.edu.au or haojingshao@gmail.com
