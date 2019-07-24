# LGC - a pipeline for identificating long non-coding RNAs based on feature relationship

Introduction
----------

LGC characterizes and identifies lncRNAs based on the relationship between ORF (open reading frame) Length and GC content. LGC is able to accurately distinguish lncRNAs from protein-coding RNAs in a cross-species manner without species-specific adjustments, and is robustly effective in discriminating lncRNAs from protein-coding RNAs across species that range from plants to mammals. More details on LGC methods and benchmarking for lncRNA identification are described in:

Wang, G. *et al.* (2019) Characterization and identification of long non-coding RNAs based on feature relationship.
*Bioinformatics*, [doi:10.1093/bioinformatics/btz008](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz008/5288512)

LGC web server is publicly available at [http://bigd.big.ac.cn/lgc/calculator](http://bigd.big.ac.cn/lgc/calculator).

Required packages for executing LGC
----------

Install Python 2.7 and Biopython:
1. Get Python 2.7 at https://www.python.org/download/ or install with your operating system’s package manager.
2. Get Biopython at http://biopython.org/wiki/Biopython .
3. You can follow the guide ( http://biopython.org/DIST/docs/install/Installation.html#sec4) to install Python and Biopython step by step.

User should include the PATH of above mentioned libraries / packages inside their SYSTEM PATH variable. 

Installation and Execution
----------
```
  $ tar zxf LGC-1.0.tar.gz # Depress LGC-1.0.tar.gz
  $ cd LGC-1.0 # Open the folder
  $ python LGC-1.0.py input.fasta output.txt # Run LGC
```
Successful run of LGC will print as following:
```
  $ Input: input.fasta # Input file
  $ Output: output.txt # Output file
  $ Scan ORF ... # Scan ORF and calculate coding potential score
  $ Done # LGC runs to completion
  $ Computation time XXX # Computation time of LGC
```

Input
----------

Fasta format:
Users can upload fasta-formatted file (<100 Mb) from local disk or paste fasta-formatted sequence(s) (small data set) into text area.
BED/GTF format:
Users can upload bed/gtf-formatted file (<3 Mb) or paste data into the text area.
When input file is BED/GTF format, the reference genome is required and the assembly version is important.
This web server now supports Human (GRCh 38, hg19), Mouse (mm10, mm9), Fly (dm3) and Zebrafish (Zv9).

Output
----------
After finishing calculation, results will be shown in a new page. Users can sort the results by any column by clicking on the column header. Also, LGC will assign an unique Task ID for each request. Users can also retrieve results by inputting the Task ID in the homepage. There are nine columns in the output file.

* **Sequence name**: name of transcript sequence
* **ORF Length**: length of the longest ORF
* **GC Content**: GC content of the longest ORF
* **Coding Potential Score**: coding potential score for a transcript, which is protein-coding RNA if greater than 0 or ncRNA if smaller than 0. '0' indicates that mRNA probability equals lncRNA probability. Also, if the ORF length is shorter than 100nt, '0' is output.
* **Coding Label**: "Coding" represents mRNA and “Non-coding” represents lncRNA.
* **pc**: Probability of ORF for coding sequence
* **pnc**: Probability of ORF for non-coding sequence
* **fc**: Stop-codon probability for coding sequence
* **fnc**:Stop-codon probability for non-coding sequence




