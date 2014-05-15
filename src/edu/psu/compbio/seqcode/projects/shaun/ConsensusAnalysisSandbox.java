package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.ScoredStrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedPoint;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.FASTALoader;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
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
                               "--strand <W/C/.>"+
                               "\nOPTIONS:\n" +
                               "--printhits \n" +
                               "--printpeakswithconsensus \n" +
                               "--printpeaksnoconsensus \n" +
                               "--printpeakswithconsensusbounds \n" +
                               "--printpeaksnoconsensusbounds \n" +
                               "--lbound <left> \n"+
                               "--rbound <right> \n" +
                               "--printseqshits \n"+
                               "--printseqswithconsensus \n"+
                               "--printseqsnoconsensus \n"+
                               "--printseqswithconsensusbounds \n"+
                               "--printseqsnoconsensusbounds \n"+
                               "--fasta <seqFile for 4 above options>\n" +
                               "--printprofiles [motif-profiler style vectors] --bins <num bin for profile>\n" +
                               "--oneperseq [flag to add one hit per sequence (profiles only)]\n" +
                               "");
            System.exit(0);
        }
        String consensusStr = ap.getKeyValue("consensus");
        int maxMismatch = ap.hasKey("mismatch") ? new Integer(ap.getKeyValue("mismatch")).intValue():0;
        ConsensusSequence consensus = new ConsensusSequence(consensusStr);
        genomeSequencePath = ap.hasKey("seq") ? ap.getKeyValue("seq") : null;
        boolean havePeaks = ap.hasKey("peaks");
        String posFile = ap.getKeyValue("peaks");
        int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():-1;
        boolean usingWin= win>0;
        boolean oneperseq = ap.hasKey("oneperseq");
        String searchStrand =  ap.hasKey("strand") ? ap.getKeyValue("strand") : ".";
        boolean printHits = ap.hasKey("printhits");
        boolean printPeaksNoConsensus = ap.hasKey("printpeaksnoconsensus");
        boolean printPeaksWithConsensus = ap.hasKey("printpeakswithconsensus");
        boolean printPeaksNoConsensusBounds = ap.hasKey("printpeaksnoconsensusbounds");
        boolean printPeaksWithConsensusBounds = ap.hasKey("printpeakswithconsensusbounds");
        int lbound = ap.hasKey("lbound") ? new Integer(ap.getKeyValue("lbound")).intValue():-1;
        int rbound = ap.hasKey("rbound") ? new Integer(ap.getKeyValue("rbound")).intValue():0;
        boolean printprofiles = ap.hasKey("printprofiles");
        int bins = ap.hasKey("bins") ? new Integer(ap.getKeyValue("bins")).intValue():1; //Used by profile printing
        boolean printSeqsHits = ap.hasKey("printseqshits");
        boolean printSeqsWithConsensus = ap.hasKey("printseqswithconsensus");
        boolean printSeqsNoConsensus = ap.hasKey("printseqsnoconsensus");
        boolean printSeqsWithConsensusBounds = ap.hasKey("printseqswithconsensusbounds");
        boolean printSeqsNoConsensusBounds = ap.hasKey("printseqsnoconsensusbounds");
        String fastaFile = ap.getKeyValue("fasta");
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
	        				if(strand!='?'){
	        					pt = new StrandedPoint(currgen, chrom, location, strand);
	        					if(usingWin)
	        						reg = new StrandedRegion(pt.expand(win/2), strand);
	        					else
	        						reg = new StrandedRegion(currgen, chrom, rstart, rend, strand);
	        				}else{
	        					pt = new Point(currgen, chrom, location);
	        					if(usingWin)
	        						reg = pt.expand(win/2);
	        					else
	        						reg = new Region(currgen, chrom, rstart, rend);
	        				}
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
	        if(printHits)
	        	tools.printConsensusHits();
	        if(printPeaksWithConsensus)
	        	tools.printPeaksWithConsensus();
	        if(printPeaksNoConsensus)
	        	tools.printPeaksWithoutConsensus();
	        if(printPeaksWithConsensusBounds)
	        	tools.printPeaksWithConsensusBoundedRegion(lbound, rbound);
	        if(printPeaksNoConsensusBounds)
	        	tools.printPeaksWithoutConsensusBoundedRegion(lbound, rbound);
	        if(printprofiles)
	        	tools.printConsensusProfiles(win, bins, oneperseq);
	        if(printSeqsHits)
	        	tools.printFastaConsensusHits(fastaFile);
	        if(printSeqsWithConsensus)
	        	tools.printFastaWithConsensus(fastaFile);
	        if(printSeqsNoConsensus)
	        	tools.printFastaWithoutConsensus(fastaFile);
	        if(printSeqsWithConsensusBounds)
	        	tools.printFastaWithConsensusBoundedRegion(fastaFile, lbound, rbound);
	        if(printSeqsNoConsensusBounds)
	        	tools.printFastaWithoutConsensusBoundedRegion(fastaFile, lbound, rbound);
	        
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
		if(searchStrand!='.' && searchStrand!='W' && searchStrand!='C'){
			System.err.println("Ignoring search strand options; should be W for Watson or C for Crick");
			searchStrand='.';
		}
		seqgen = new SequenceGenerator(g);
		if(genomeSequencePath != null){
			seqgen.useCache(true);
			seqgen.useLocalFiles(true);
			seqgen.setGenomePath(genomeSequencePath);
		}
		scorer = new ConsensusSequenceScorer(consensus, seqgen);
	}
	
	//printing the hit sequences
	public void printConsensusHits(){
		int numHits=0, peaksWithHits=0, totalPeaks=0;
		List<ScoredStrandedRegion> hits = new ArrayList<ScoredStrandedRegion>();
		List<String> hitseqs = new ArrayList<String>();
		
		for(int i=0; i<regions.size(); i++){
			totalPeaks++;
			Region r = regions.get(i);
			String seq = seqgen.execute(r);
			
			ConsensusSequenceScoreProfile profiler = scorer.execute(r, searchStrand);
			boolean goodMotif=false;
			for(int z=0; z<r.getWidth()-consensus.getLength()+1; z++){
				int currScore= profiler.getLowestMismatch(z);
				if(currScore<=misMatchThreshold){
					numHits++;
					goodMotif=true;
					String subseq = seq.substring(z, z+consensus.getLength());
					Region hitreg =new Region(gen, r.getChrom(), r.getStart()+z, r.getStart()+z+consensus.getLength()-1);
					hits.add(new ScoredStrandedRegion(gen, r.getChrom(), hitreg.getStart(), hitreg.getEnd(), currScore, profiler.getLowestMismatchStrand(z)));
					if(profiler.getLowestMismatchStrand(z)=='+'){hitseqs.add(subseq);
					}else{hitseqs.add(SequenceUtils.reverseComplement(subseq));}
				}
			}
			if(goodMotif){
				peaksWithHits++;
			}
        }
		double perc = (double)peaksWithHits/(double)totalPeaks;
		System.out.println(consensus.getSequence()+" hits: "+numHits+" hits in "+peaksWithHits+" regions from "+totalPeaks+" total peaks ("+perc+").");
		for(int h=0; h<hits.size(); h++){
			System.out.println(hits.get(h).toTabString()+"\t"+hitseqs.get(h));
		}
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
	
	//Print peaks that contain the consensus within a bounded region with respect to the region center
	public void printPeaksWithConsensusBoundedRegion(int left, int right){
		for(int x=0; x<peaks.size(); x++){
			Point a = peaks.get(x);
			int halfwin = Math.max(Math.abs(left), Math.abs(right))+consensus.getLength()+1;
			int start = Math.max(1, a.getLocation()-halfwin);
			int end = Math.min(a.getLocation()+halfwin, a.getGenome().getChromLength(a.getChrom()));
			Region query = new Region(gen, a.getChrom(), start, end);
			
			char strand = (a instanceof StrandedPoint) ? 
					((StrandedPoint)a).getStrand() : '.';
			String seq = seqgen.execute(query);
			if(strand=='-')
				seq = SequenceUtils.reverseComplement(seq);
			
			ConsensusSequenceScoreProfile profiler = scorer.execute(seq, searchStrand=='.' ? '.':(searchStrand=='W' ? '+' : '-'));
			boolean matchFound=false;
			for(int i=0; i<seq.length(); i++){
				int offset = strand=='-' ?
						(a.getLocation() - (query.getEnd()-i)) :
						(i+query.getStart() - a.getLocation());
				if(offset>=left && offset<=right)
					if(profiler.getLowestMismatch(i)<=misMatchThreshold)
						matchFound=true;
			}
			if(matchFound)
				System.out.println(inFileLines.get(x));
		}
	}
	
	//Print peaks that don't contain the consensus within a bounded region with respect to the region center
	public void printPeaksWithoutConsensusBoundedRegion(int left, int right){
		for(int x=0; x<peaks.size(); x++){
			Point a = peaks.get(x);
			int halfwin = Math.max(Math.abs(left), Math.abs(right))+consensus.getLength()+1;
			int start = Math.max(1, a.getLocation()-halfwin);
			int end = Math.min(a.getLocation()+halfwin, a.getGenome().getChromLength(a.getChrom()));
			Region query = new Region(gen, a.getChrom(), start, end);
			
			char strand = (a instanceof StrandedPoint) ? 
					((StrandedPoint)a).getStrand() : '.';
			String seq = seqgen.execute(query);
			if(strand=='-')
				seq = SequenceUtils.reverseComplement(seq);
			
			ConsensusSequenceScoreProfile profiler = scorer.execute(seq, searchStrand=='.' ? '.':(searchStrand=='W' ? '+' : '-'));
			boolean matchFound=false;
			for(int i=0; i<seq.length(); i++){
				int offset = strand=='-' ?
						(a.getLocation() - (query.getEnd()-i)) :
						(i+query.getStart() - a.getLocation());
				if(offset>=left && offset<=right)
					if(profiler.getLowestMismatch(i)<=misMatchThreshold)
						matchFound=true;
			}
			if(!matchFound)
				System.out.println(inFileLines.get(x));
		}
	}
	
	/**
	 * Print FastA entries hits that contain matches to the consensus
	 * @param inFile
	 */
	public void printFastaConsensusHits(String inFile){
		int numHits=0, peaksWithHits=0, totalPeaks=0;
		List<ScoredStrandedRegion> hits = new ArrayList<ScoredStrandedRegion>();
		List<String> hitseqs = new ArrayList<String>();
		
		FASTALoader loader = new FASTALoader();
		File f = new File(inFile);
		Iterator<Pair<String, String>> it = loader.execute(f);
		while(it.hasNext()){
			Pair<String, String> p = it.next();
			String name = p.car();
			String seq = p.cdr();
			ConsensusSequenceScoreProfile profiler = scorer.execute(seq, searchStrand=='.' ? '.':(searchStrand=='W' ? '+' : '-'));
			boolean goodMotif=false;
			for(int z=0; z<seq.length()-consensus.getLength()+1; z++){
				int currScore= profiler.getLowestMismatch(z);
				if(currScore<=misMatchThreshold){
					numHits++;
					goodMotif=true;
					char strand = profiler.getLowestMismatchStrand(z);
					String subseq = seq.substring(z, z+consensus.getLength());
					System.out.println(name+"\t"+currScore+"\t"+z+"\t"+strand+"\t"+subseq);
				}
			}
			if(goodMotif){
				peaksWithHits++;
			}
        }
		double perc = (double)peaksWithHits/(double)totalPeaks;
		System.out.println(consensus.getSequence()+" hits: "+numHits+" hits in "+peaksWithHits+" regions from "+totalPeaks+" total peaks ("+perc+").");
		
	}
	
	/**
	 * Print FastA entries that contain matches to the consensus
	 * @param inFile
	 */
	public void printFastaWithConsensus(String inFile){
		FASTALoader loader = new FASTALoader();
		File f = new File(inFile);
		Iterator<Pair<String, String>> it = loader.execute(f);
		while(it.hasNext()){
			Pair<String, String> p = it.next();
			String name = p.car();
			String sequence = p.cdr();
			ConsensusSequenceScoreProfile profiler = scorer.execute(sequence, searchStrand=='.' ? '.':(searchStrand=='W' ? '+' : '-'));
			
			if(profiler.getLowestMismatch()<=misMatchThreshold){
				System.out.println(name+"\n"+sequence);
			}
		}
	}
	/**
	 * Print FastA entries that doesn't contain matches to the consensus
	 * @param inFile
	 */
	public void printFastaWithoutConsensus(String inFile){
		FASTALoader loader = new FASTALoader();
		File f = new File(inFile);
		Iterator<Pair<String, String>> it = loader.execute(f);
		while(it.hasNext()){
			Pair<String, String> p = it.next();
			String name = p.car();
			String sequence = p.cdr();
			ConsensusSequenceScoreProfile profiler = scorer.execute(sequence, searchStrand=='.' ? '.':(searchStrand=='W' ? '+' : '-'));
			
			if(profiler.getLowestMismatch()>misMatchThreshold){
				System.out.println(name+"\n"+sequence);
			}
		}
	}
	
	/**
	 * Print FastA entries that contain matches to the consensus within bounds
	 * @param inFile
	 */
	public void printFastaWithConsensusBoundedRegion(String inFile, int left, int right){
		FASTALoader loader = new FASTALoader();
		File f = new File(inFile);
		Iterator<Pair<String, String>> it = loader.execute(f);
		while(it.hasNext()){
			Pair<String, String> p = it.next();
			String name = p.car();
			String sequence = p.cdr();
			ConsensusSequenceScoreProfile profiler = scorer.execute(sequence, searchStrand=='.' ? '.':(searchStrand=='W' ? '+' : '-'));
			boolean matchFound=false;
			for(int i=0; i<sequence.length(); i++){
				int offset = i-(sequence.length()/2);
				if(offset>=left && offset<=right)
					if(profiler.getLowestMismatch(i)<=misMatchThreshold)
						matchFound=true;
			}
			if(matchFound){
				System.out.println(name+"\n"+sequence);
			}
		}
	}
	/**
	 * Print FastA entries that doesn't contain matches to the consensus within bounds
	 * @param inFile
	 */
	public void printFastaWithoutConsensusBoundedRegion(String inFile, int left, int right){
		FASTALoader loader = new FASTALoader();
		File f = new File(inFile);
		Iterator<Pair<String, String>> it = loader.execute(f);
		while(it.hasNext()){
			Pair<String, String> p = it.next();
			String name = p.car();
			String sequence = p.cdr();
			ConsensusSequenceScoreProfile profiler = scorer.execute(sequence, searchStrand=='.' ? '.':(searchStrand=='W' ? '+' : '-'));
			boolean matchFound=false;
			for(int i=0; i<sequence.length(); i++){
				int offset = i-(sequence.length()/2);
				if(offset>=left && offset<=right)
					if(profiler.getLowestMismatch(i)<=misMatchThreshold)
						matchFound=true;
			}
			if(!matchFound){
				System.out.println(name+"\n"+sequence);
			}
		}
	}
	
	//return an indicator array for covered/uncovered by good motif
	private boolean[] motifCoverage(String seq){
		boolean [] covered = new boolean [seq.length()];
		for(int i=0; i<seq.length(); i++){covered[i]=false;}
		ConsensusSequenceScoreProfile profiler = scorer.execute(seq, searchStrand=='.' ? '.':(searchStrand=='W' ? '+' : '-'));
		for(int z=0; z<seq.length(); z++){
			double currScore= profiler.getLowestMismatch(z);
			if(currScore<=misMatchThreshold){
				for(int j=0; j<consensus.getLength() && z+j<seq.length(); j++){
					covered[z+j]=true;
				}
			}
		}
		return(covered);
	}
	/**
	 * Print vectors showing presence/absence of consensus matches at binwidth resolution
	 * @param winLen
	 * @param bins
	 */
	public void printConsensusProfiles(int winLen, int bins, boolean oneperseq){
		BinningParameters params = new BinningParameters(winLen, bins);
		for(int x=0; x<peaks.size(); x++){
			Point a = peaks.get(x);
			double[] array = new double[params.getNumBins()];
			for(int i = 0; i < array.length; i++) { array[i] = 0; }
			
			Region query = a.expand(winLen/2);
			
			char strand = (a instanceof StrandedPoint) ? 
					((StrandedPoint)a).getStrand() : '.';
			String seq = seqgen.execute(query);
			if(strand=='-')
				seq = SequenceUtils.reverseComplement(seq);
			
			ConsensusSequenceScoreProfile profiler = scorer.execute(seq, searchStrand=='.' ? '.':(searchStrand=='W' ? '+' : '-'));
			if(oneperseq){
				List<Integer> hits = new ArrayList<Integer>();
				for(int i=0; i<seq.length() && i<winLen; i++){
					if(profiler.getLowestMismatch(i)<=misMatchThreshold){
						int bin = params.findBin(i);
						hits.add(bin);
					}
				}
				if(hits.size()>0){
					if(hits.size()==1)
						array[hits.get(0)]++;
					else{
						int rand = (int)(Math.floor((Math.random() * (double)hits.size())));
						int bin = params.findBin(hits.get(rand));
						array[bin]++;
					}
				}
			}else{
				for(int i=0; i<seq.length() && i<winLen; i++){
					if(profiler.getLowestMismatch(i)<=misMatchThreshold){
						int bin = params.findBin(i);
						array[bin]++;
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
