seqcode-core
============

SeqCode: Java code for the analysis of high-throughput sequencing data

This repository contains a core set of Java classes for interacting with genomic data.
Several Java projects that are released by Shaun Mahony's lab at Penn State University depend on seqcode-core for basic interactions with data.  

Much of this package (particularly org.seqcode.genome and org.seqcode.gsebricks) is based on GSE, the "Genomic Spatial Events" codebase from the Gifford lab at CSAIL, MIT (package: edu.mit.csail.cgs.gse). 
GSE is outlined in this publication: http://www.ncbi.nlm.nih.gov/pubmed/18229714
The core classes of GSE were built by Tim Danford and Alex Rolfe. Later additions were also contributed by Bob Altshuler, Shaun Mahony, Yuchun Guo, Chris Reeder, and Giorgos Papachristoudis.

Development of seqcode-core is currently maintained by Shaun Mahony, Akshay Kakumanu, and Naomi Yamada.
	

Major History:
--------------
Jun 6th 2016:  Major reorganization of packages. Separated this project (seqcode-code) from other seqcode projects in preparation for move to public github repo. Refactored to current namespace (org.seqcode.*). Breaking up old gse package into clearer structure. Collected all database and file interactions into org.seqcode.data. 

Jun 25th 2015: Overhauled database communications; Tomcat JDBC now used to handle connection pooling, connections now mainly opened and closed in same methods. 

May-Aug 2014: Substantial reorganization was carried out, resulting in the deletion of many outdated subpackages and classes, and the separating out of genome/region loading functionality and deep-sequencing experiment functionality into org.seqcode.genome and org.seqcode.deepseq, respectively.

Oct 24th 2012: This version of GSE was branched from the main MIT code. Initial refactoring deleted a number of outdated subpackages, updated the naming spaces for compatibility with the Mahony lab servers at PSU, and aimed to remove any interactions with Oracle databases in favor of mysql.
  
