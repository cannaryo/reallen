=============
About ReALLEN
=============

ReALLEN is a software package for finding large genomic rearrangement using NGS data.
The software is developed by Genome Team in Fukushima Medical University.

============
Requirements
============

Before running ReALLEN, please confirm that the following packages are properly installed on your system.

1. Ruby 2.2.3 or lator.  
https://github.com/ruby/ruby.git
See offcial site for information.  
https://www.ruby-lang.org/

example for the installation of ruby-2.2.5 :

    wget https://cache.ruby-lang.org/pub/ruby/2.2/ruby-2.2.5.tar.bz2
    tar xvf ruby-2.2.5.tar.bz2
    cd ruby-2.2.5
    ./configure
    make
    make install

2. BWA software package.  
https://github.com/lh3/bwa.git

3. SAMtools library to process SAM/BAM files.  
https://github.com/samtools/samtools.git

I reccomend you to use git for installing ReALLEN and the depending packages.

============
Installation
============

    git clone https://github.com/cannaryo/reallen.git
    cd reallen-master
    ./setup.sh

You need network connection to download the reference and annotations during the instalation.

I recommend you to make symbolic link of 'reallen-start.sh' and 'annotate-sv.sh' into a location with a PATH (for example, /usr/local/bin)

========
Test run
========

Try following commands to check the instlation.

    cd testdata
    ../reallen-start.sh sample.bam
    ../annotate-sv.sh sample.bp

If ReALLEN properly worked, 'sample.bp' and 'sample.csv' were created.  
'sample.bp' is row output which includes the postions of breakpoints.  
'sample.csv' is annotated table.

==============================
Run ReALLEN with default value
==============================

    reallen-start.sh <your-bam-file-name.bam>

It is easy and well-oplitized run for human genome (hg19).
To make table format file, use following command.

    annotate-sv.sh <your-bam-file-name.bp>

Note that 'annotate-sv.sh' should be run in the same directory where you have run 'reallen-start.sh'.

========================
Run ReALLEN step by step
========================

Each script file is located at 'reallen/src/'  
Each script '*.rb' is simply excutable if your system has Ruby interpreter.
See run-step-by-step.txt to learn how to apply each script on your data.

If you want to use another resources (for example, reference of another species), you have to make them by yourselves, 
however, please contact me without hesitation if you have any request on your analysis about ReALLEN.


========
Citation
========

* ReALLEN: structural variation discovery in cancer genome by sensitive analysis of single-end reads.  
 Ryo Kanno, Daisuke Tanaka, Hideaki Nanamiya, Takao Isogai
