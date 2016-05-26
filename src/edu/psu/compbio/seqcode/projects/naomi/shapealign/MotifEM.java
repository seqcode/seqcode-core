package edu.psu.compbio.seqcode.projects.naomi.shapealign;

import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;

/**
 * MEME algorithm to find a motif given a prior, meant to be integrated with shapeEM, which is to be developed later.
 * 
 * @author naomi yamada
 *
 */

public class MotifEM {
	protected GenomeConfig gconfig;	
	protected List<StrandedPoint> strandedPoints;
	protected int window;
	
	int q = 10 ; // number of iterations
	int N; // number of sequences
	int W = 6; // length of the motif
	int L; // length of each sequences
	double[][][] z_ijq; // estimate after q iterations of EM of the probabilities that the site begins at position j in sequence i given the model and the data.  
	double[][][] p_lkq; // estimate after q iterations of EM of the probabilities of letter l appearing in position k of the motif
	double[][] po_q; // estimate after q iterations of EM of the base frequencies of outside of motif
	
	public MotifEM(GenomeConfig gcon, List<StrandedPoint> spoints, int win){
		gconfig = gcon;	
		strandedPoints = spoints;
		window = win;		
		N = strandedPoints.size();
		L = window+1;
		
		z_ijq = new double[N][L-W+1][q];
		p_lkq = new double[4][W][q];
		po_q = new double [4][q];	
	}
	
	public void runEM(){
		// converting stranded points to stranded regions
		List<StrandedRegion> regionList = new ArrayList<StrandedRegion>();
		for(Point p: strandedPoints){		
			int start = Math.max(1, p.getLocation() - window/2 );
			int end = Math.min(p.getLocation() + window/2, p.getGenome().getChromLength(p.getChrom()));				
			StrandedRegion strandedReg = new StrandedRegion(p.getGenome(), p.getChrom(), start, end, p.getStrand());					
			regionList.add(strandedReg);
		}
		
		// get sequences from regions
		List<String> sequences = new ArrayList<String>();
		SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>(gconfig.getGenome());		
		for (StrandedRegion reg : regionList){
			String seq = seqgen.execute(reg);
			sequences.add(seq.toUpperCase());
		}	
		
		for (int itr = 0 ; itr < q; itr++){
			//Expectation
			
			//Maximization
			
			
			
			
		}
		
		
	}
	
	public void Expectation(){
		
	}
	
	public void Maximization(){
		
	}
	

	// method test
	public static void main(String[] args){
		ArgParser ap = new ArgParser(args);		
        if((!ap.hasKey("peaks") && !ap.hasKey("regions")) ) { 
            System.err.println("Usage:\n " +
                               "MotifEM " +
                               "--species <organism;genome> OR\n" +
                               "--geninfo <genome info> AND --seq <path to seqs>\n" +
                               "--peaks <file containing coordinates of peaks> \n" +
                               "--win <window of sequence to take around peaks> OR\n" +
                               "--regions \n" +
                               "");
            System.exit(0);
        }
        
		GenomeConfig gconf = new GenomeConfig(args);
		
		// parsing command line arguments	
		int win = Args.parseInteger(args, "win", 100);
		List<StrandedPoint> strandedPoints = RegionFileUtilities.loadStrandedPointsFromMotifFile(gconf.getGenome(), ap.getKeyValue("peaks"), win);
		
		MotifEM motifEM = new MotifEM(gconf, strandedPoints, win);
	}

}
