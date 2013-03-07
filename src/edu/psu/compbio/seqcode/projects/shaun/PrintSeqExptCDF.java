package edu.psu.compbio.seqcode.projects.shaun;

import java.util.ArrayList;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class PrintSeqExptCDF {

	private static Organism org;
	private static Genome gen;
	private static int cdfThres=100;

	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("genome")||!ap.hasKey("expt") ||!ap.hasKey("out")) { 
            System.err.println("Usage:\n" +
                               "SeqExptProbLandscapeKmers " +
                               "--species <organism name> " +
                               "--genome <genome version> "+
                               "--expt <solexa expt> " +
                               "--out <out filename> ");
            return;
        }
        String species = ap.getKeyValue("species");
        String genome = ap.getKeyValue("genome");
        String exptName = ap.getKeyValue("expt");
        String outName = ap.getKeyValue("out");
        
    	try {
			org = Organism.getOrganism(species);
			gen = org.getGenome(genome);
			
			SeqDataHandler exHandle = new SeqDataHandler(org, gen, exptName);
			exHandle.compileCDF(cdfThres);
			exHandle.printCDF(outName);
			
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    

	}

}
