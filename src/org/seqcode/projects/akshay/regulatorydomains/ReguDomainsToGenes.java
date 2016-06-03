package org.seqcode.projects.akshay.regulatorydomains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gse.datasets.motifs.WeightMatrix;
import org.seqcode.gse.gsebricks.verbs.location.PointParser;
import org.seqcode.gse.gsebricks.verbs.location.RegionParser;
import org.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.io.RegionFileUtilities;
import org.seqcode.projects.shaun.MotifAnalysisSandbox;


/**
 * 
 * @author akshaykakumanu
 *
 */
public class ReguDomainsToGenes {
	
	private Genome gen;
	private List<RegulatoryRegion> rRegs;
	private String OutputFormat = "SUMMARY";
	private List<String> motifnames = new ArrayList<String>();
	private List<Double> motifMarkovThreshs = new ArrayList<Double>();
	private List<WeightMatrix> motiflist = new ArrayList<WeightMatrix>();
	private int numClus = 0;
	private int win;
	private double topPerc;
	// is the sorting criteria fold-change or -1*log(p/q) value
	private String sortType = "foldchange";
	
	public ReguDomainsToGenes(Genome g, int w) {
		gen=g;
		win=w;
	}
	
	//settors
	public void setRegRs(List<RegulatoryRegion> rRs){rRegs = rRs;}
	public void setOutFormat(String outf){OutputFormat = outf;}
	public void setMotifnames(List<String> mnames){motifnames = mnames;}
	public void setMotifThresholds(List<Double> mthreshs){motifMarkovThreshs = mthreshs;}
	public void setNumClus(int nC){numClus = nC;}
	public void setTopPerc(double tpc){topPerc = tpc;}
	public void serMotifList(List<WeightMatrix> mList){motiflist =mList; }
	public void setSortType(String st){sortType = st;}
	
	
	// Creates different RegulatoryClassProfiles and print them
	public void execute(){
		// Fist deal with any inf
		// No filers just usign all the regRegions
		RegulatoryClassProfile currProfile = new RegulatoryClassProfile(rRegs, 100, sortType,  null, null, null, 1, numClus, OutputFormat,"Base");
		// Now use only top-scoreing regions
		currProfile = new RegulatoryClassProfile(rRegs, topPerc, sortType, null, null, null, 2, numClus, OutputFormat,"Base+topPerc");
		
		// That have the primary motif (ASSUMING THE PRIMARY MOTIF HAS INDEX 0 ALWAYS!!!)
		if(motiflist.size()>0){ // at leas the primary motif was provided
			List<Integer> primaryMotifInd = new ArrayList<Integer>();
			primaryMotifInd.add(0);
			List<Integer> minMotifHitCount = new ArrayList<Integer>();
			minMotifHitCount.add(1);
			currProfile = new RegulatoryClassProfile(rRegs, topPerc, sortType, primaryMotifInd, minMotifHitCount, motifMarkovThreshs, 2, numClus, OutputFormat,"Base+topPerc+PrimaryMotif(1)");
			
			// That have atleast two instance of the primary motif
			minMotifHitCount = new ArrayList<Integer>();
			minMotifHitCount.add(2);
			currProfile = new RegulatoryClassProfile(rRegs, topPerc, sortType, primaryMotifInd, minMotifHitCount, motifMarkovThreshs, 2, numClus, OutputFormat,"Base+topPerc+PrimaryMotif(2)");
			
		}
		
		// Now add one motif at a time
		for(int m=1; m<motiflist.size(); m++){
			List<Integer> seclectedMotifInds = new ArrayList<Integer>();
			seclectedMotifInds.add(0);
			seclectedMotifInds.add(m);
			List<Integer> seclectedMotifMinHitCounts = new ArrayList<Integer>();
			seclectedMotifMinHitCounts.add(1);
			seclectedMotifMinHitCounts.add(1);
			currProfile = new RegulatoryClassProfile(rRegs, topPerc, sortType,  seclectedMotifInds, seclectedMotifMinHitCounts, motifMarkovThreshs, 2, numClus, OutputFormat,"Base+topPerc+PrimaryMotif(1)+"+motifnames.get(m));
		}
		// Now add all motifs
		List<Integer> allmotifInds = new ArrayList<Integer>();
		List<Integer> allmotifMinHitCounts = new ArrayList<Integer>();
		for(int m=0; m <motiflist.size(); m++){
			allmotifInds.add(m);
			allmotifMinHitCounts.add(1);
		}
		currProfile = new RegulatoryClassProfile(rRegs, topPerc, sortType,allmotifInds, allmotifMinHitCounts, motifMarkovThreshs, 2, numClus, OutputFormat,"Base+topPerc+PrimaryMotif(1)+AllOtherMotifs");
		
		// Now all motifs plus atleast 2 instance of the primany motif
		allmotifMinHitCounts = new ArrayList<Integer>();
		allmotifMinHitCounts.add(2);
		for(int m=1; m <motiflist.size(); m++){
			allmotifMinHitCounts.add(1);
		}
		currProfile = new RegulatoryClassProfile(rRegs, topPerc, sortType,allmotifInds, allmotifMinHitCounts, motifMarkovThreshs, 2, numClus, OutputFormat,"Base+topPerc+PrimaryMotif(2)+AllOtherMotifs");
		
		// Now add homotypic also
		//currProfile = new RegulatoryClassProfile(rRegs, topPerc, sortType, allmotifInds, motifMarkovThreshs, 2, numClus, OutputFormat,"Base+topPerc+Allmotifs+Homotypics");
		
	}
	
	
	
	
	
	
	
	
	public static void main(String[] args) throws IOException, ParseException{
		GenomeConfig gcon = new GenomeConfig(args);
		int win = Args.parseInteger(args, "win", 150);
		ReguDomainsToGenes runner = new ReguDomainsToGenes(gcon.getGenome(),win);
		
		// Load misc arguments 
		ArgParser ap = new ArgParser(args);
		double topPerc = Args.parseDouble(args, "topPerc", 20);
		runner.setTopPerc(topPerc);
		
		// Load peak related features
		List<Point> peaks = new ArrayList<Point>();
		List<Double> bindingstrength = new ArrayList<Double>();
		List<Double> bindingdynamics = new ArrayList<Double>();
		List<Double> homotopyicIndexes = new ArrayList<Double>();
		List<String> Seqs = new ArrayList<String>();
		// The bindingevents or the peaks file can have 3 columns
		// 1st column:- peak, or region,
		//2nd column :- binding intensity (fold-over input)
		// 3rd colum :- differential binding (fol-change in binding over two conditions)
		String peaksfilename = ap.getKeyValue("bindingevents");
		File peaksfile = new File(peaksfilename);
		if(!peaksfile.isFile()){System.err.println("Invalid file name: "+peaksfilename);System.exit(1);}
        BufferedReader reader = new BufferedReader(new FileReader(peaksfile));
        String line;
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            String[] words = line.split("\t");
            if(words.length>0 && !words[0].contains("#") && !words[0].equals("Region") && !words[0].equals("Position")){
            	if(words[0].contains(":") && words[0].contains("-")){
            		RegionParser rparser = new RegionParser(gcon.getGenome());
            		Region q = rparser.execute(words[0]);
            		peaks.add(q.getMidpoint());
            	}else if(words[0].contains(":") && !words[0].contains("-")){
            		PointParser pparser = new PointParser(gcon.getGenome());
            		Point p = pparser.execute(words[0]);
            		peaks.add(p);
            	}else{
            		System.err.println("Invalid binding events file format");System.exit(1);
            	}
            	// A very inefficient way to deal with infinities; but can't think of a better way to deal with them at then moment
            	// -Infinities are hard-coded to change to -100
            	// Infinitiy are hard-coded to chan ge to 20
            	if(words[1].equals("-Infinity") || words[1].equals("-Inf") || words[1].equals("-infinity") || words[1].equals("-inf")){
            		words[1] = "-100";
            	}
            	if(words[1].equals("Infinity") || words[1].equals("-Inf") || words[1].equals("-infinity") || words[1].equals("-inf")){
            		words[1] = "20";
            	}
            	bindingstrength.add(Double.parseDouble(words[1]));
            	if(words.length>2){
            		if(words[2].equals("-Infinity") || words[2].equals("-Inf") || words[2].equals("-infinity") || words[2].equals("-inf")){
                		words[2] = "-100";
                	}
                	if(words[2].equals("Infinity") || words[2].equals("-Inf") || words[2].equals("-infinity") || words[2].equals("-inf")){
                		words[2] = "20";
                	}
            		bindingdynamics.add(Double.parseDouble(words[2]));
            	}
            	if(words.length>3){
            		homotopyicIndexes.add(Double.parseDouble(words[3]));
            	}
            }
        }reader.close();
        
        
        String sortBY = Args.parseString(args, "sortCriteria", "BS");
        
        if(!sortBY.equals("BS") && !sortBY.equals("BD")){
        	System.err.println("Invalid sort option: Has to be either 'BS' or 'BD'");System.exit(1);
        }
        
        String sortType = Args.parseString(args, "sortType", "foldcahnge");
        runner.setSortType(sortType);
        
        
        
        // Now load motif PWMs
        String motiffile = ap.getKeyValue("motiffile");
		String backfile = ap.getKeyValue("back");
		List<WeightMatrix> matrixList = new ArrayList<WeightMatrix>();
		if(motiffile != null && backfile !=null){
			matrixList = MotifAnalysisSandbox.loadMotifFromFile(motiffile, backfile, gcon.getGenome());
		}else{
			System.err.println("No motifs provided");
		}
		runner.serMotifList(matrixList);
		
		// Now load motif names and markov thresholds..
		String thresFile;
		double thresLevel;
		if(ap.hasKey("motiffile") && !ap.hasKey("motifThresh") ){System.err.println("Provide thresholdfile for the motifs"); System.exit(1);}
		List<Double> motifthreshs = new ArrayList<Double>();
		if(ap.hasKey("motiffile") && ap.hasKey("motifThresh")){
			thresFile = ap.getKeyValue("motifThresh");
			thresLevel = Args.parseDouble(args, "threshlvl", 0.1);
			HashMap<String,Double> thres = MotifAnalysisSandbox.loadThresholdsFromFile(thresFile, thresLevel);
			List<String> motifnames = new ArrayList<String>();
			
			for(WeightMatrix wm : matrixList){
				motifnames.add(wm.getName());
				motifthreshs.add(thres.get(wm.getName()));
			}
			runner.setMotifnames(motifnames);
			runner.setMotifThresholds(motifthreshs);
				
		}
	
		// Now load seqquences
		@SuppressWarnings("rawtypes")
		SequenceGenerator seqgen = new SequenceGenerator();
		boolean useCache = ap.hasKey("seq");
		seqgen.useCache(useCache);
		if(useCache){
			seqgen.setGenomePath(ap.getKeyValue("seq"));
		}
		for(int p=0; p<peaks.size(); p++){
			Seqs.add(seqgen.execute(peaks.get(p).expand(win/2)));
		}
		SequenceGenerator.clearCache();
		
		// Load all target gene related files now
		// Assumes the file has 4 columns separated by a tab
		// 1st column is the peak location string
		//2nd column is gene name
		//3rd column is fold-change
		//4th column is cluster membership
		
		HashMap<String,String> peakGenePairs = new HashMap<String,String>();
		HashMap<String,Double> geneFCs = new HashMap<String,Double>();
		HashMap<String,Integer> geneClusIndex = new HashMap<String,Integer>();
		
		
		String genesfilename = ap.getKeyValue("targetgenes");
		int nClus = Args.parseInteger(args, "numClusters", 0);
		runner.setNumClus(nClus);
		File genefile = new File(genesfilename);
		if(!genefile.isFile()){System.err.println("Invalid file name: "+genesfilename);System.exit(1);}
        BufferedReader greader = new BufferedReader(new FileReader(genefile));
     
        while ((line = greader.readLine()) != null) {
            line = line.trim();
            String[] words = line.split("\t");
            if(words.length>0 && !words[0].contains("#") && !words[0].equals("Region") && !words[0].equals("Position")){
            	if(words[0].contains(":") && words[0].contains("-")){
            		RegionParser rparser = new RegionParser(gcon.getGenome());
            		Region q = rparser.execute(words[0]);
            		peakGenePairs.put(q.getMidpoint().getLocationString(), words[1]);
            		
            	}else if(words[0].contains(":") && !words[0].contains("-")){
            		PointParser pparser = new PointParser(gcon.getGenome());
            		Point p = pparser.execute(words[0]);
            		peakGenePairs.put(p.getLocationString(), words[1]);
            	}else{
            		System.err.println("Invalid target genes file format");System.exit(1);
            	}
            	geneFCs.put(words[1], Double.parseDouble(words[2]));
        		if(words.length>3){
        			geneClusIndex.put(words[1], Integer.parseInt(words[3]));
        		}
            }
        }greader.close();
      
        
		// Now create Regulatory regions
		List<RegulatoryRegion> regRs = new ArrayList<RegulatoryRegion>();
		for(int p=0; p<peaks.size(); p++){
			if(sortBY.equals("BS")){
				if(bindingdynamics != null){
					regRs.add(new BSRegulatoryRegion(peaks.get(p), bindingstrength.get(p) ,bindingdynamics.get(p) ,win, matrixList,motifthreshs, Seqs.get(p)));
				}else{
					regRs.add(new BSRegulatoryRegion(peaks.get(p), bindingstrength.get(p) ,win, matrixList,motifthreshs, Seqs.get(p)));
				}
			}else{
				regRs.add(new BDRegulatoryRegion(peaks.get(p), bindingstrength.get(p) ,bindingdynamics.get(p) ,win, matrixList,motifthreshs, Seqs.get(p)));
			}
		}
		
		// Now add homotypic peaks and target genes
		// First hash the peaks by chromosome for faster look up
		HashMap<String, List<Point>> byChr = new HashMap<String, List<Point>>();
		for(Point p : peaks){
			if(!byChr.containsKey(p.getChrom()))
				byChr.put(p.getChrom(), new ArrayList<Point>());
			byChr.get(p.getChrom()).add(p);
		}
		
		for(int rR=0; rR<regRs.size(); rR++){
			String currChr = regRs.get(rR).getChrom();
			// First add the target gene
			String currLocation = regRs.get(rR).getPeakLocation();
			String currTargetGene = peakGenePairs.get(currLocation);
			GeneDomain gd = new GeneDomain(currTargetGene,geneFCs.get(currTargetGene));
			if(geneClusIndex.containsKey(currTargetGene)){
				gd.setClusterInd(geneClusIndex.get(currTargetGene));
			}
			regRs.get(rR).setTargetGene(gd);
			
			if(byChr.containsKey(currChr)){
				for(Point p: byChr.get(currChr)){
					if(regRs.get(rR).coversPeak(p)){
						regRs.get(rR).addHomotypicPeak(p);
					}
				}
			}
		}
		
		
		// If the homotypic indexes were provided in the binding events file override the current indexes with those values
		if(homotopyicIndexes.size() > 0){
			for(int rR=0; rR<regRs.size(); rR++){
				regRs.get(rR).setHomotypicIndex(homotopyicIndexes.get(rR));
			}
		}
		
		
		// Now set the regulatory regions
		runner.setRegRs(regRs);
		
		// Finally, Get the type of the oputput format 
		String Output = ap.getKeyValue("OutFormat");
		if(!Output.equals("SUMMARY") && !Output.equals("PEAKLISTS") && !Output.equals("UPPERC") && !Output.equals("DOWNPERC") && !Output.equals("CLUSPERC") && !Output.equals("AVGFC")){
			System.err.println("Invalid output format type: Provide either SUMMARY; UPPERC; DOWNPERC; CLUSPERC; AVGFC");System.exit(1);
		}
		runner.setOutFormat(Output);
		
		runner.execute();
		
	}
	
	
	
	

}
