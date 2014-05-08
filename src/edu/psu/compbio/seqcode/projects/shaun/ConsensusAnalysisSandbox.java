package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedPoint;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.viz.metaprofile.BinningParameters;

public class ConsensusAnalysisSandbox {

	private Genome gen;
	private ConsensusSequence consensus;
	private ConsensusSequenceScorer scorer;
	private int misMatchThreshold = 0;
	private char searchStrand='.';
	private SequenceGenerator seqgen;
	private List<Region> regions=null;
	private List<Point> peaks=null;
	private List<String> inFileLines=null;
	
	public static void main(String[] args) throws IOException, ParseException {
		String genomeSequencePath=null;
		ConsensusAnalysisSandbox tools;
		
		ArgParser ap = new ArgParser(args);
        if((!ap.hasKey("species") && !ap.hasKey("geninfo")) || (!ap.hasKey("consensus"))) { 
            System.err.println("Usage:\n " +
                               "ConsensusAnalysisSandbox " +
                               "--species <organism;genome> OR\n" +
                               "--geninfo <genome info> AND --seq <path to seqs>\n" +
                               "--consensus <IUPAC consensus>\n" +
                               "--mismatch <mismatch limit>\n" +
                               "--peaks <file containing coordinates of peaks> \n" +
                               "--win <window of sequence to take around peaks> \n" +
                               "--strand <+/-/.>"+
                               "\nOPTIONS:\n" +
                               "--printpeakswithconsensus \n" +
                               "--printpeaksnoconsensus \n" +
                               "--printprofiles [motif-profiler style vectors] --bins <num bin for profile>\n"+
                               "");
        }
        String consensusStr = ap.getKeyValue("consensus");
        int maxMismatch = ap.hasKey("mismatch") ? new Integer(ap.getKeyValue("mismatch")).intValue():0;
        ConsensusSequence consensus = new ConsensusSequence(consensusStr);
        genomeSequencePath = ap.hasKey("seq") ? ap.getKeyValue("seq") : null;
        boolean havePeaks = ap.hasKey("peaks");
        String posFile = ap.getKeyValue("peaks");
        int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():-1;
        boolean usingWin= win>0;
        String searchStrand =  ap.hasKey("strand") ? ap.getKeyValue("strand") : ".";
        boolean printPeaksNoConsensus = ap.hasKey("printpeaksnoconsensus");
        boolean printPeaksWithConsensus = ap.hasKey("printpeakswithconsensus");
        boolean printprofiles = ap.hasKey("printprofiles");
        int bins = ap.hasKey("bins") ? new Integer(ap.getKeyValue("bins")).intValue():1; //Used by profile printing
        List<Region> posRegs = new ArrayList<Region>();
        List<Point> posPeaks = new ArrayList<Point>();
        List<String> posLines= new ArrayList<String>();
        
        try {
			Genome currgen=null;
			//Load genome
			if(ap.hasKey("species")){
				Pair<Organism, Genome> pair = Args.parseGenome(args);
				if(pair != null){
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
			
			//Load the positive hits
	        if(havePeaks){
		        File pFile = new File(posFile);
				if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	            BufferedReader reader = new BufferedReader(new FileReader(pFile));
	            String line;
	            while ((line = reader.readLine()) != null) {
	                line = line.trim();
                	Point pt=null; Region reg=null;
                	String [] curr = line.split("\\s+");
        			String coord = curr[0];
        			if(!curr[0].contains("#")){
        				if(curr.length>=3 && curr[2].contains(":")){coord = curr[2];}
	        			
	        			if(coord.contains(":")) {
	        				String [] currB = coord.split(":");
	        				String chrom = currB[0]; chrom=chrom.replaceFirst("chr", "");
	        				char strand = '?';
	        				if(currB.length==3)
	        					strand = currB[2].charAt(0);
	        				int location=-1, rstart=-1, rend=-1;
	        				if(currB[1].contains("-")){
	        					String [] currC = currB[1].split("-");
	        					location = (new Integer(currC[0])+new Integer(currC[1]))/2;
	        					rstart = new Integer(currC[0]);
	        					rend = new Integer(currC[1]);
	        				}else{
	        					location = new Integer(currB[1]);
	        					rstart = Math.max(0, location-(win/2));
	        					rend = Math.min(0, location+(win/2));
	        				}
	        				if(strand!='?')
	        					pt = new StrandedPoint(currgen, chrom, location, strand);
	        				else
	        					pt = new Point(currgen, chrom, location);
	        				
	        				if(usingWin)
	        					reg = pt.expand(win/2);
	        				else
	        					reg = new Region(currgen, chrom, rstart, rend);
	        			}
                	
	                    if(pt!=null && reg!=null){
		                	posPeaks.add(pt);
		                	posRegs.add(reg);
		                	posLines.add(line);
		                }	                    
	                }
	            }
	        }
	        
	        tools = new ConsensusAnalysisSandbox(currgen, consensus, maxMismatch, posLines, posPeaks, posRegs, genomeSequencePath, searchStrand);
	        if(printPeaksWithConsensus)
	        	tools.printPeaksWithConsensus();
	        if(printPeaksNoConsensus)
	        	tools.printPeaksWithoutConsensus();
	        if(printprofiles)
	        	tools.printConsensusProfiles(win, bins);
	        
        } catch (NotFoundException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	//Constructor
	public ConsensusAnalysisSandbox(Genome g, ConsensusSequence cons, int maxMismatch, List<String> inL, List<Point> pospeak, List<Region> posreg, String genomeSequencePath, String strand){
		gen=g;
		consensus= cons;
		inFileLines = inL;
		peaks=pospeak;
		regions=posreg;
		misMatchThreshold = maxMismatch;
		searchStrand = strand.charAt(0);
		seqgen = new SequenceGenerator(g);
		if(genomeSequencePath != null){
			seqgen.useCache(true);
			seqgen.useLocalFiles(true);
			seqgen.setGenomePath(genomeSequencePath);
		}
		scorer = new ConsensusSequenceScorer(consensus, seqgen);
	}
	
	//Print peaks that contain the consensus
	public void printPeaksWithConsensus(){
		for(int i=0; i<regions.size(); i++){
			Region r = regions.get(i);
			ConsensusSequenceScoreProfile profiler = scorer.execute(r, searchStrand);
			if(profiler.getLowestMismatch()<=misMatchThreshold)
				System.out.println(inFileLines.get(i));
		}
	}
	
	//Print peaks that don't contain the consensus
	public void printPeaksWithoutConsensus(){
		for(int i=0; i<regions.size(); i++){
			Region r = regions.get(i);
			ConsensusSequenceScoreProfile profiler = scorer.execute(r, searchStrand);
			if(profiler.getLowestMismatch()>misMatchThreshold)
				System.out.println(inFileLines.get(i));
		}
	}
	
	//Print vectors showing presence/absense of consensus matches at binwidth resolution
	public void printConsensusProfiles(int winLen, int bins){
		BinningParameters params = new BinningParameters(winLen, bins);
		for(int x=0; x<peaks.size(); x++){
			Point a = peaks.get(x);
			double[] array = new double[params.getNumBins()];
			for(int i = 0; i < array.length; i++) { array[i] = 0; }
			
			int window = params.getWindowSize();
			int left = window/2;
			int right = window-left-1;
			
			int start = Math.max(1, a.getLocation()-left);
			int end = Math.min(a.getLocation()+right, a.getGenome().getChromLength(a.getChrom()));
			Region query = new Region(gen, a.getChrom(), start, end);
			boolean strand = (a instanceof StrandedPoint) ? 
					((StrandedPoint)a).getStrand() == '+' : true;
			
			String seq = seqgen.execute(query);
			ConsensusSequenceScoreProfile profiler = scorer.execute(seq);
			for(int i=query.getStart(); i<query.getEnd(); i+=params.getBinSize()){
				for(int j=i; j<i+params.getBinSize() && j<query.getEnd(); j++){
					int offset = j-query.getStart();
					
					if(profiler.getLowestMismatch(offset)<=misMatchThreshold){
						int bin = params.findBin(offset);
						array[bin] = Math.max(array[bin],1);
					}
				}
			}
			//Print the vector
			System.out.print(a);
			for(int z=0; z<array.length; z++){
				System.out.print("\t"+array[z]);
			}System.out.print("\n");
		}
	}
}
