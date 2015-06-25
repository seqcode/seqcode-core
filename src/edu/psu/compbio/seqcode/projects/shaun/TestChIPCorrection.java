package edu.psu.compbio.seqcode.projects.shaun;

import java.util.Iterator;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

/**
 * 
 * @author Shaun Mahony
 * 
 * TestChIPCorrection: a simple peak caller based on a ChIP-seq error model that 
 * uses the relative likelihood of observed hit counts.
 *
 */
public class TestChIPCorrection {

	private static double DEFT_T = -10.5;
	/**
	 * @param args
	 */
	public static void main(String[] args) {
//		String exptName = "YoungLab_Solexa_Oct4";
//		String exptName = "PPG_Solexa_RAR_ES+2d";
//		String exptName = "PPG_Solexa_RAR_8hr";
//		String exptName = "PPG_Solexa_RAR_2+1";
//		String WCEName = "PPG_Solexa_WCE_ES+2d";
//		String exptName = "PPG_Solexa_WCE_2+1";
		
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("genome")||!ap.hasKey("expt") || !ap.hasKey("back")) { 
            System.err.println("Usage:\n" +
                               "TestChIPCorrection" +
                               "--species <organism name> " +
                               "--genome <genome version> "+
                               "--expt <solexa expt> " +
                               "--back <background expt>"+
                               "--t <threshold on peak ll-ratio> " +
                               "--ecdf <expt CDF> " +
                               "--bcdf <back CDF> ");
            return;
        }
        String species = ap.getKeyValue("species");
        String genome = ap.getKeyValue("genome");
        String exptName = ap.getKeyValue("expt");
        String backName = ap.getKeyValue("back");
        double userT = ap.hasKey("t") ? Double.valueOf(ap.getKeyValue("t")) : DEFT_T;
        String ecdfFile = ap.hasKey("ecdf") ? ap.getKeyValue("ecdf") : "NOFILE";
        String bcdfFile = ap.hasKey("bcdf") ? ap.getKeyValue("bcdf") : "NOFILE";
        
		try {
			System.out.println(exptName+"\t"+backName);
			Species org = Species.getSpecies(species);
			Genome gen = new Genome(org, genome);
				
			SeqDataHandler IPhandle = new SeqDataHandler(org, gen, exptName);
			SeqDataHandler backhandle = new SeqDataHandler(org, gen, backName);
			
			if(!ecdfFile.equals("NOFILE")){IPhandle.loadCDFFile(ecdfFile);}
			if(!bcdfFile.equals("NOFILE")){backhandle.loadCDFFile(bcdfFile);}
			
			SeqDataLikelihoodRatio selr = new SeqDataLikelihoodRatio(IPhandle, backhandle);
			selr.calcPeaks(userT);
			selr.printPeaks(annots);
			
			
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	public static String [] annots = new String[]{
		"refGene" //, "knownGene", "mgcGenes", "ensGene"
	};

}
