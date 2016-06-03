package org.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.datasets.motifs.CountsBackgroundModel;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.RealValuedHistogram;
import org.seqcode.gse.utils.io.BackgroundModelIO;
import org.seqcode.gse.utils.io.RegionFileUtilities;
import org.seqcode.gse.utils.sequence.SequenceUtils;


public class KmerAnalysisSandbox {

	private GenomeConfig gConfig;
	private Genome gen;
	private int k=6;
	private int win=200;
	private List<Region> posRegions=null;
	private List<Point> posPeaks=null;
	private List<String> posSeq=null;
	private List<String> posLines=null;
	private boolean negLoaded=false;
	private List<Region> negRegions=null;
	private List<Point> negPeaks=null;
	private List<String> negSeq=null;
	private double negPseudoCount=1;
	private int numRand=100000; //Only used if simulating or randomly picking regions
	private int histoBinSize=5;
	private CountsBackgroundModel backCounts = null;
        	
	public static void main(String[] args) throws IOException, ParseException {
		ArgParser ap = new ArgParser(args);
		GenomeConfig gConfig = new GenomeConfig(args);
        if(!ap.hasKey("species") || (!ap.hasKey("k"))|| (!ap.hasKey("win"))|| (!ap.hasKey("peaks"))) { 
            System.err.println("Usage:\n " +
                               "KmerAnalysisSandbox \n" +
                               "Required:\n" +
                               "\t--species <organism;genome>\n" +
                               "\t--k <kmer length>\n" +
                               "\t--win <window around peak to examine>\n" +
                               "\t--peaks <peaks file>\n" +
                               "Options:\n" +
                               "\t--seq <genome fasta seq directory>\n" +
                               "\t--model <ranked kmer file>\n" +
                               "\t--neg <filename/random/back>\n" +
                               "\t--backcounts <counts background model that includes appropriate kmer counts>\n" +
                               "\t--numrand <number of random positions if random negative selected)\n" +
                               "\t--out <output filename>\n" +
                               "\t--rankthres <threshold for kmer matching>\n" +
                               "\t--enumerate [make a k-mer model]\n" +
                               "\t--roc [TP vs FP at the peaks level]\n" +
                               "\t--peakswithkmers [print peaks where the regions contain one of the top kmers]\n" +
                               "\t--kmerhits [print hits to the top kmers in the peaks]\n" +
                               "\t--kmerdisthist \n" +
                               "\t--recenterpeaks [recenter peaks to the nearest top kmer (where region contains kmer]\n" +
                               "\n");
            
        }else{
	        int k = ap.hasKey("k") ? new Integer(ap.getKeyValue("k")).intValue():6;
	        int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():200;
	        String peaksFile = ap.hasKey("peaks") ? ap.getKeyValue("peaks") : null;
	        String neg = ap.hasKey("neg") ? ap.getKeyValue("neg") : null;
	        String backfile = null; 
	        CountsBackgroundModel back = null;
	        String out = ap.hasKey("out") ? ap.getKeyValue("out") : "out";
	        int rankThres = ap.hasKey("rankthres") ? new Integer(ap.getKeyValue("rankthres")).intValue():-1;
	        int numrand = ap.hasKey("numrand") ? new Integer(ap.getKeyValue("numrand")).intValue():-1;
	        
	        if(ap.hasKey("backcounts")){
	        	backfile = ap.getKeyValue("backcounts");
	        	back = BackgroundModelIO.parseCountsBackgroundModel(backfile, gConfig.getGenome());
	        }
	        
	        KmerAnalysisSandbox analyzer = new KmerAnalysisSandbox(gConfig, k, win, peaksFile, back);
	        analyzer.setNumRand(numrand);
	        
	        KmerModel model = null;
	        if(ap.hasKey("enumerate")){
	        	model = analyzer.estimateKmerModel(neg);
	        	analyzer.printKmerModel(model, out);
	        }else if(ap.hasKey("model")){
	        	model = analyzer.loadKmerModel(ap.getKeyValue("model"));
	        }
	        
	        if(model != null){
	        	if(ap.hasKey("roc")){
	        		analyzer.ROC(model, neg);
	        	}
	        	if(ap.hasKey("peakswithkmers") && rankThres>0){
	        		analyzer.printPeaksWithKmers(model, rankThres);
	        	}
	        	if(ap.hasKey("kmerhits") && rankThres>0){
	        		analyzer.printKmersAtPeaks(model, rankThres);
	        	}
	        	if(ap.hasKey("recenterpeaks") && rankThres>0){
	        		analyzer.printRecenteredPeaksWithKmers(model, rankThres, out, neg);
	        	}
	        	if(ap.hasKey("kmerdisthist") && rankThres>0){
	        		analyzer.peak2motifHisto(model, rankThres);
	        	}
	        }
	   
        }                               
	}

	public KmerAnalysisSandbox(GenomeConfig g, int k, int win, String pFile, CountsBackgroundModel back){
		gConfig=g;
		gen = gConfig.getGenome();
		this.k=k;
		this.win=win;
		posRegions = RegionFileUtilities.loadRegionsFromPeakFile(gen, pFile, win);
		posPeaks = RegionFileUtilities.loadPeaksFromPeakFile(gen, pFile, win);
		posLines = RegionFileUtilities.loadLinesFromFile(pFile);
		posSeq = RegionFileUtilities.getSequencesForRegions(posRegions, gConfig.getSequenceGenerator());
		backCounts = back;
	}
	
	//Accessors
	public void setNumRand(int r){if(r>0){numRand=r;}}
	
	/**
	 * Make a new Kmer model, which is just a ranked list of Kmers.
	 * @param negative
	 * @return
	 */
	public KmerModel estimateKmerModel(String negative){
		List<Kmer> kmerList = new ArrayList<Kmer>();
		System.err.println("Estimating Kmer models");
		
		//Enumerate all k-mers in positive and negative sequences. 
		//Reuse BackgroundModel code for this purpose
		CountsBackgroundModel negCounts;
		if(negative.equals("back") && backCounts!=null){
			System.err.println("Loading negative kmers from background file");
			negCounts = backCounts;
		}else{
			loadNegativeRegions(negative);
			System.err.println("Counting kmers in negative sequences");
			negCounts = CountsBackgroundModel.modelFromSeqList(gen, negSeq, k);
			negCounts.degenerateStrands();
		}
		
		System.err.println("Counting kmers in positive sequences");
		CountsBackgroundModel posCounts = CountsBackgroundModel.modelFromSeqList(gen, posSeq, k);
		posCounts.degenerateStrands();
		
		
		//Go through all k-mers, and calculate pos/neg enrichment
		System.err.println("Positive vs negative enrichment");
		List<String> kmers = RegionFileUtilities.getAllKmers(k);
		double posTotal=0, negTotal=0;
		for(String kmer : kmers){
			String revkmer = SequenceUtils.reverseComplement(kmer);
			if(RegionFileUtilities.seq2int(kmer) <= RegionFileUtilities.seq2int(revkmer)){//Count each k-mer once only
				posTotal += posCounts.getKmerCount(kmer);
				negTotal += negCounts.getKmerCount(kmer)>0 ? negCounts.getKmerCount(kmer) : negPseudoCount;
			}
		}
		for(String kmer : kmers){
			String revkmer = SequenceUtils.reverseComplement(kmer);
			if(RegionFileUtilities.seq2int(kmer) <= RegionFileUtilities.seq2int(revkmer)){//Count each k-mer once only
			    double posFreq = (double)(posCounts.getKmerCount(kmer)/posTotal); 
			    double negFreq = (double)((negCounts.getKmerCount(kmer)>0 ? negCounts.getKmerCount(kmer) : negPseudoCount)/negTotal); 
			    double eScore = posFreq/negFreq;
			    
			    Kmer currK = new Kmer(kmer, eScore, posFreq, negFreq);
			    kmerList.add(currK);
			}
		}
		KmerModel model = new KmerModel(kmerList);
		return model;
	}
	
	/**
	 * TP vs FP at level of peaks
	 * @param model
	 */
	public void ROC(KmerModel model, String negative){
		loadNegativeRegions(negative);
		double ROCAUC=0;
		double numPosPassed=0, numNegPassed=0;
		double totalPos = (double) posRegions.size();
		double totalNeg = (double) negRegions.size();
		boolean[] posPassed = new boolean[posRegions.size()];
		for(int i=0; i<posRegions.size(); i++){posPassed[i]=false;}
		boolean[] negPassed = new boolean[negRegions.size()];
		for(int i=0; i<negRegions.size(); i++){negPassed[i]=false;}
		double TP=0, FP=0, lastTP=0, lastFP=0;
		
		System.out.println("Kmer\tRevKmer\tEnrich\tPosFreq\tNegFreq\tTP\tFP");
		//Iterate through ranked list
		for(Kmer mer : model.getKmerList()){
			//Look for matches to this mer and its reverse complement in the positive seqs
			for(int s=0; s<posSeq.size(); s++){
				String seq = posSeq.get(s);
				if(!posPassed[s]){
					if(seq.indexOf(mer.kmer)!=-1 || seq.indexOf(mer.revkmer)!=-1){
						posPassed[s]=true;
						numPosPassed++;
					}
				}
			}
			//Look for matches to this mer and its reverse complement in the negative seqs
			for(int s=0; s<negSeq.size(); s++){
				String seq = negSeq.get(s);
				if(!negPassed[s]){
					if(seq.indexOf(mer.kmer)!=-1 || seq.indexOf(mer.revkmer)!=-1){
						negPassed[s]=true;
						numNegPassed++;
					}
				}
			}
			TP = (numPosPassed/totalPos);
			FP = (numNegPassed/totalNeg);
			System.out.println(mer+"\t"+TP+"\t"+FP);

			ROCAUC += ((FP-lastFP)*lastTP) + ((TP-lastTP)*(FP-lastFP)/2);
			
			lastFP = FP; 
			lastTP = TP;
			if(TP==1 && FP==1)
				break;
		}
		System.out.println("\n\n#ROCAUC\t"+ROCAUC);
	}
	
	/**
	 * Print the peaks that contain any of the top rankThres k-mers
	 * @param model
	 * @param rankThres
	 */
	public void printPeaksWithKmers(KmerModel model, int rankThres){
		int kmerCount=0;
		boolean[] posPassed = new boolean[posRegions.size()];
		for(int i=0; i<posRegions.size(); i++){posPassed[i]=false;}
		for(Kmer mer : model.getKmerList()){
			//Look for matches to this mer and its reverse complement in the positive seqs
			for(int s=0; s<posSeq.size(); s++){
				String seq = posSeq.get(s);
				if(!posPassed[s]){
					if(seq.indexOf(mer.kmer)!=-1 || seq.indexOf(mer.revkmer)!=-1){
						posPassed[s]=true;
						System.out.println(posLines.get(s));
					}
				}
			}
			kmerCount++;
			if(kmerCount>rankThres)
				break;
		}
	}
	/**
	 * Print the peaks that contain any of the top rankThres k-mers
	 * @param model
	 * @param rankThres
	 */
	public void printKmersAtPeaks(KmerModel model, int rankThres){
		int kmerCount=0;
		boolean[] posPassed = new boolean[posRegions.size()];
		for(int i=0; i<posRegions.size(); i++){posPassed[i]=false;}
		for(Kmer mer : model.getKmerList()){
			//Look for matches to this mer and its reverse complement in the positive seqs
			for(int s=0; s<posSeq.size(); s++){
				String seq = posSeq.get(s);
				Region r = posRegions.get(s);
				Point p = posPeaks.get(s);
				
				for(int z=0; z<seq.length()-k; z++){
					String subseq = seq.substring(z, z+k);
					int start = r.getStart()+z;
					int stop = r.getStart()+z+k-1;
					if(subseq.equals(mer.kmer)){
						System.out.println(r.getChrom()+":"+start+"-"+stop+":+\t"+mer.enrichment+"\t"+mer.kmer);
					}else if(subseq.equals(mer.revkmer)){ 
						System.out.println(r.getChrom()+":"+start+"-"+stop+":-\t"+mer.enrichment+"\t"+mer.kmer);
					}
				}
			}
			kmerCount++;
			if(kmerCount>rankThres)
				break;
		}
	}
	
	//Histogram of all matches to the motif vs distance to peak 
	public void peak2motifHisto(KmerModel model, int rankThres){
		int binSize = 10, kmerCount=0;
		int numBins = (k+(win/2))/binSize;
		double [] motifDistHisto = new double [numBins+1];
		for(int h=0; h<=numBins; h++){motifDistHisto[h]=0;}
		
		for(Kmer mer : model.getKmerList()){
			//Look for matches to this mer and its reverse complement in the positive seqs
			for(int s=0; s<posSeq.size(); s++){
				String seq = posSeq.get(s);
				Region r = posRegions.get(s);
				Point p = posPeaks.get(s);
				
				for(int z=0; z<seq.length()-k; z++){
					String subseq = seq.substring(z, z+k);
					if(subseq.equals(mer.kmer) || subseq.equals(mer.revkmer)){ 
						int dist = (p.getLocation()-r.getStart())-z;
						motifDistHisto[(Math.abs(dist))/binSize]++;
					}
				}
			}
			kmerCount++;
			if(kmerCount>rankThres)
				break;
		}
		
		System.out.println("\nPeak-to-Motif Histogram\nDistance\tMotifCounts\tMotifCountsNorm");
		for(int h=0; h<=numBins; h++){
			int bin =h*binSize; 
			System.out.println(bin+"\t"+motifDistHisto[h]+"\t"+motifDistHisto[h]/(double)posRegions.size());
		}
	}
	
	/**
	 * Print the peaks that contain any of the top rankThres k-mers, recentering the peak on the closest match.
	 * If negative is null, you just don't get a histogram of expected shifts. 
	 * @param model
	 * @param rankThres
	 */
	public void printRecenteredPeaksWithKmers(KmerModel model, int rankThres, String outFile, String negative){
		HashMap<Point,Integer> posClosestKmerShift = new HashMap<Point, Integer>();
		HashMap<Point,Kmer> posClosestKmer = new HashMap<Point, Kmer>();
		HashMap<Point,Integer> negClosestKmerShift = new HashMap<Point, Integer>();
		RealValuedHistogram posShiftHisto = new RealValuedHistogram(-win/2, win/2, (win/histoBinSize));
		RealValuedHistogram negShiftHisto = new RealValuedHistogram(-win/2, win/2, (win/histoBinSize));
		
		int kmerCount=0;
		for(Kmer mer : model.getKmerList()){
			Pattern pfor = Pattern.compile(mer.kmer);
			Pattern prev = Pattern.compile(mer.revkmer);
			 
			//Look for matches to this mer and its reverse complement in the positive seqs
			for(int s=0; s<posRegions.size(); s++){
				Region r = posRegions.get(s);
				Point p = posPeaks.get(s);
				int peakOffset = p.getLocation()-r.getStart();
				int currShift = posClosestKmerShift.containsKey(p) ? posClosestKmerShift.get(p) : 10000000; 

				//Forward matches
				Matcher m = pfor.matcher(posSeq.get(s));
				while(m.find()){
					int koffset = ((m.start()+m.end())/2)-peakOffset;
					if(Math.abs(koffset)<Math.abs(currShift)){
						posClosestKmerShift.put(p, koffset);
						posClosestKmer.put(p, mer);
					}
				}
				//Reverse matches
				m = prev.matcher(posSeq.get(s));
				while(m.find()){
					int koffset = ((m.start()+m.end())/2)-peakOffset;
					if(Math.abs(koffset)<Math.abs(currShift)){
						posClosestKmerShift.put(p, koffset);
						posClosestKmer.put(p, mer);
					}
				}
			}
			kmerCount++;
			if(kmerCount>rankThres)
				break;
		}
		
		//Print the shifted positive peaks
		try{
			FileWriter fw = new FileWriter(outFile);
			for(Point p : posPeaks){
			    if(posClosestKmerShift.containsKey(p)){
				Point shifted = new Point(p.getGenome(), p.getChrom(), p.getLocation()+posClosestKmerShift.get(p));
				int shift = posClosestKmerShift.get(p);
				fw.write(shifted+"\t"+shift+"\t"+posClosestKmer.get(p).kmer+"\n");
				//Add to histogram
				posShiftHisto.addValue(shift);
			    }
			}
			System.out.println("#Positive peak kmer shift histogram:");
			posShiftHisto.printContents();
			fw.close();			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//Process negative sequences if provided
		if(negative!=null){
			loadNegativeRegions(negative);
			kmerCount=0;
			for(Kmer mer : model.getKmerList()){
				Pattern pfor = Pattern.compile(mer.kmer);
				Pattern prev = Pattern.compile(mer.revkmer);
				 
				//Look for matches to this mer and its reverse complement in the positive seqs
				for(int s=0; s<negRegions.size(); s++){
					Region r = negRegions.get(s);
					Point p = negPeaks.get(s);
					int peakOffset = p.getLocation()-r.getStart();
					int currShift = negClosestKmerShift.containsKey(p) ? negClosestKmerShift.get(p) : 10000000; 
							
					//Forward matches
					Matcher m = pfor.matcher(negSeq.get(s));
					while(m.find()){
						int koffset = ((m.start()+m.end())/2)-peakOffset;
						if(Math.abs(koffset)<Math.abs(currShift))
							negClosestKmerShift.put(p, koffset);
					}
					//Reverse matches
					m = prev.matcher(negSeq.get(s));
					while(m.find()){
						int koffset = ((m.start()+m.end())/2)-peakOffset;
						if(Math.abs(koffset)<Math.abs(currShift))
							negClosestKmerShift.put(p, koffset);
					}
				}
				kmerCount++;
				if(kmerCount>rankThres)
					break;
			}
			for(Point p : negClosestKmerShift.keySet()){
				int shift = negClosestKmerShift.get(p);
				////Add to histogram
				negShiftHisto.addValue(shift);
			}
			System.out.println("#Negative peak kmer shift histogram:");
			negShiftHisto.printContents();
		}
	}
	/**
	 * Print the kmerModel to an output file
	 * @param out
	 */
	public void printKmerModel(KmerModel model, String out){
		try{
			FileWriter fw = new FileWriter(out);
			for(Kmer curr : model.getKmerList()){
				fw.write(curr.toString()+"\n");
			}
			fw.close();			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * Load a kmer model from a file
	 * @param filename
	 * @return
	 */
	public KmerModel loadKmerModel(String filename){
		List<Kmer> kList = new ArrayList<Kmer>();
		try{
			File kFile = new File(filename);
			if(!kFile.isFile()){System.err.println("Invalid kmer model file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(kFile));
	        String line = null;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            String k = words[0];
	            Double enrich = new Double(words[2]);
	            kList.add(new Kmer(k, enrich));
	        }reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return new KmerModel(kList);
	}
	
	/**
	 * Load negative regions & sequences
	 * @param negative
	 */
	private void loadNegativeRegions(String negative){
		if(negative==null){System.err.println("Negative sequences required for this functionality."); System.exit(1);}
		System.err.println("Fetching negative sequences");
		if(!negLoaded){
			//Get the negative sequences first
			if(negative==null || negative.equals("random")){
				negRegions = RegionFileUtilities.randomRegionPick(gen, posRegions, numRand, win);
				negPeaks = RegionFileUtilities.regions2midpoints(negRegions);
				negSeq = RegionFileUtilities.getSequencesForRegions(negRegions, gConfig.getSequenceGenerator());
			}else{
				negRegions = RegionFileUtilities.loadRegionsFromPeakFile(gen, negative, win);
				negPeaks = RegionFileUtilities.loadPeaksFromPeakFile(gen, negative, win);
				negSeq = RegionFileUtilities.getSequencesForRegions(negRegions, gConfig.getSequenceGenerator());
			}
			negLoaded=true;
		}
	}
	/**
	 * Kmer: represents a scored string. 
	 * The score itself is implementation-specific, so don't rely on any given interpretation.
	 * Think of this instead as a rankable String.
	 * @author mahony
	 *
	 */
	public class Kmer implements Comparable<Kmer>{
		public String kmer;
		public String revkmer;
		public Double enrichment;
		public Double posFreq;
		public Double negFreq;

		public Kmer(String k, double e){ this(k, e, -1, -1);}
		public Kmer(String k, double e, double p, double n){
			kmer = k;
			revkmer = SequenceUtils.reverseComplement(k);
			enrichment = e;	
			posFreq=p;
			negFreq=n;
		}
		public int compareTo(Kmer k){
			if(enrichment > k.enrichment)
				return -1;
			else if (enrichment < k.enrichment)
				return 1;
			return 0;
		}
		public String toString(){
			return(new String(kmer+"\t"+revkmer+"\t"+enrichment+"\t"+posFreq+"\t"+negFreq));
		}
	}
	
	public class KmerModel{
		private List<Kmer> kmers=null;
		public KmerModel(List<Kmer> m){
			kmers = m;
			Collections.sort(kmers);
		}
		
		public List<Kmer> getKmerList(){return kmers;}
	}
}
