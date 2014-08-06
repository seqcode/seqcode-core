package edu.psu.compbio.seqcode.projects.akshay.Utils;

import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collection;
import java.util.LinkedList;

import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class PairedBEDExporter {

	
	private Genome gen;
	private SeqDataLoader loader;
	private boolean closeLoader;
	private SeqLocator locator;
	private Collection<SeqAlignment> alignments;
	String outName="out";
	
	public PairedBEDExporter(SeqLocator expt, Genome g) throws SQLException, IOException, NotFoundException {
		
		this.locator = expt;
		this.loader = new SeqDataLoader();
		this.closeLoader = true;
		this.gen = g;
		Collection<SeqAlignment> alignments = locator.loadAlignments(loader, gen);
	}
}	
	
//	public void execute(){
//		try{
//			FileWriter fw = new FileWriter(outName);
//			ChromRegionIterator chroms = new ChromRegionIterator(gen);
//			
//			while(chroms.hasNext()){
//			
//		}
//	}
//}
