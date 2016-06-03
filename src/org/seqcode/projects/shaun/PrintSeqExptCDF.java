package org.seqcode.projects.shaun;

import java.util.ArrayList;

import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.NotFoundException;


public class PrintSeqExptCDF {

	private static Species org;
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
			org = Species.getSpecies(species);
			gen = new Genome(org, genome);
			
			SeqDataHandler exHandle = new SeqDataHandler(org, gen, exptName);
			exHandle.compileCDF(cdfThres);
			exHandle.printCDF(outName);
			
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    

	}

}
