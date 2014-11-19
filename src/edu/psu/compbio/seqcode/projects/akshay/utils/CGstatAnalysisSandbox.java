package edu.psu.compbio.seqcode.projects.akshay.utils;


import java.io.File;
import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Organism;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.viz.metaprofile.BinningParameters;

public class CGstatAnalysisSandbox {
	private Genome gen;
	private int winSize;
	private int nbins;
	private SequenceGenerator seqgen;
	private List<Point> peaks;
	private List<Region> regions;
	
	
	public CGstatAnalysisSandbox(Genome g, String peakfile, int win, int nbin, String GenomeSeqFile) {
		gen=g;
		winSize = win;
		nbins = nbin;
		peaks = RegionFileUtilities.loadPeaksFromPeakFile(gen, peakfile, winSize);
		regions = RegionFileUtilities.loadRegionsFromPeakFile(gen, peakfile, winSize);
		seqgen = new SequenceGenerator(gen);
		if(GenomeSeqFile != null){
			seqgen.useCache(true);
			seqgen.useLocalFiles(true);
			seqgen.setGenomePath(GenomeSeqFile);
		}
	}
	
	
	public void printProfile(){
		
		BinningParameters params = new BinningParameters(winSize, nbins);
		double[] cumalativeProfile = new double[params.getNumBins()];
		CGScorer scorer = new CGScorer(params,seqgen);
		for(Region r : regions){
			CGScoreProfile profile = scorer.execute(r);
			double[] CGs = profile.getCGpercProfile();
			for(int i=0; i<CGs.length; i++){
				cumalativeProfile[i] = cumalativeProfile[i]+CGs[i];
			}
					
		}
		
		for(int i=0; i<cumalativeProfile.length; i++){
			cumalativeProfile[i] = cumalativeProfile[i]/regions.size();
			System.out.println(Integer.toString(i)+"\t"+Double.toString(cumalativeProfile[i]));
		}	
	}
	
	
	public void printTotals(){
		BinningParameters params = new BinningParameters(winSize, nbins);
		CGScorer scorer = new CGScorer(params,seqgen);
		int MaxCG = (int)winSize/2;
		for(Region r : regions){
			CGScoreProfile profile = scorer.execute(r);
			int total = profile.getCGcount();
			double perc = total/MaxCG;
			System.out.println(r.getMidpoint().getLocationString()+"\t"+Double.toString(perc));
		}
		
		
	}
	
	
	
	public static void main(String[] args) throws NotFoundException{
		ArgParser ap = new ArgParser(args);
		int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():200;
		int bins = ap.hasKey("bins") ? new Integer(ap.getKeyValue("bins")).intValue():10;
		Genome currgen=null;
		Organism currorg=null;
		//Load genome
		if(ap.hasKey("species")){
			Pair<Organism, Genome> pair = Args.parseGenome(args);
			if(pair != null){
				currorg = pair.car();
				currgen = pair.cdr();
			}
		}else{
			if(ap.hasKey("geninfo") || ap.hasKey("g")){
				//Make fake genome... chr lengths provided
				String fName = ap.hasKey("geninfo") ? ap.getKeyValue("geninfo") : ap.getKeyValue("g");
				currgen = new Genome("Genome", new File(fName), true);
			}else{
				currgen = null;
			}
		}
		
		
		String posFile = ap.getKeyValue("peaks");
		String genomeSequencePath = ap.hasKey("seq") ? ap.getKeyValue("seq") : null;
		
		
		boolean printprofile = ap.hasKey("printProfile");
		boolean printTotal = ap.hasKey("printTotal");
		
		CGstatAnalysisSandbox analyzer = new CGstatAnalysisSandbox(currgen, posFile, win, bins,genomeSequencePath);
		
		if(printprofile){
			analyzer.printProfile();
		}else if(printTotal){
			analyzer.printTotals();
		}
		
		
		
		
		
		
	}
	
	
	

}
