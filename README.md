# Branch for SEM

seqcode-core
============

SeqCode: Java code for the analysis of high-throughput sequencing data

Purpose
-------
This repository contains a core set of Java classes for interacting with genomic data.
Several Java projects that are released by Shaun Mahony's lab at Penn State University depend on seqcode-core for basic interactions with data. For example, anything relating to getting sequencing experiment data from databases or navigating genomic locations. This is the only repo that other SeqCode projects should depend on.


Organization
------------
* org.seqcode.data: Packages for interacting with databases & files. Command-line tools for importing/exporting are typically stored in tools subpackages.
* org.seqcode.deepseq: Classes for loading sequencing-based experimental data from ReadDB and files (e.g. SAM/BAM).
* org.seqcode.genome: Classes that represent genomes in our databases, as well as locations and sequences within those genomes.
* org.seqcode.gsebricks: Encodes the core "GSEBricks" methods as described in the GSE paper. The philosophy is to base various components for interacting with data (genome annotations, experimental data, etc) on Java's Iterator interface. The "verb" modules are mainly extensions from Java's Mapper, Filter, or Expander interfaces. Mapper produces Objects on a one-to-one basis with the input. Expander returns Iterators from each execute call. The Iterators of Iterators idea is supposed to echo Lisp-like functional composition. Once you get your head around the general idea, the code enables quite powerful and extensible systems where you can chain one iterator into another. However, this package is currently somewhat bloated since it probably contains a lot of legacy code.
* org.seqcode.gseutils: General classes from the GSE project that are used throughout SeqCode.
* org.seqcode.ml: Aims to be a store of general machine learning algorithms for the lab. Currently has GSE implementations of some clustering methods, but these methods may not be accurate or useful, so check them before use and please add more methods here!
* org.seqcode.motifs: Aims to contain a generally useful set of motif analysis classes. Not deeply/usefully populated yet.
* org.seqcode.projects.seed: SEED provides a set of classes for performing statistical-enrichment-based event detection for various types of sequencing experiments. The SEED package also aims to abstract the common functionality that most event detection approaches would require (e.g. finding enriched regions on the genome, etc). In doing so, SEED can be used to rapidly "seed" the development of event detection tools for particular purposes.
* org.seqcode.projects.seqview: The SeqView genome browser.
* org.seqcode.tools: Command line tools for various operations on genomic locations or sequences.
* org.seqcode.tutorials: Example code using various aspects of seqcode-core
* org.seqcode.viz: Some visualization code.


Relationship with GSE
---------------------
Much of this package (particularly org.seqcode.genome and org.seqcode.gsebricks) is based on GSE, the "Genomic Spatial Events" codebase from the Gifford lab at CSAIL, MIT (package: edu.mit.csail.cgs.gse). 
GSE is outlined in this publication: http://www.ncbi.nlm.nih.gov/pubmed/18229714
The core classes of GSE were built by Tim Danford and Alex Rolfe. Later additions were also contributed by Bob Altshuler, Shaun Mahony, Yuchun Guo, Chris Reeder, and Giorgos Papachristoudis. 

Major Modification History
--------------------------
Jun 6th 2016:  Major reorganization of packages. Separated this project (seqcode-code) from other seqcode projects in preparation for move to public github repo. Refactored to current namespace (org.seqcode.*). Breaking up old gse package into clearer structure. Collected all database and file interactions into org.seqcode.data. 

Jun 25th 2015: Overhauled database communications; Tomcat JDBC now used to handle connection pooling, connections now mainly opened and closed in same methods. 

May-Aug 2014: Substantial reorganization was carried out, resulting in the deletion of many outdated subpackages and classes, and the separating out of genome/region loading functionality and deep-sequencing experiment functionality into org.seqcode.genome and org.seqcode.deepseq, respectively.

Oct 24th 2012: This version of GSE was branched from the main MIT code. Initial refactoring deleted a number of outdated subpackages, updated the naming spaces for compatibility with the Mahony lab servers at PSU, and aimed to remove any interactions with Oracle databases in favor of mysql.
  
