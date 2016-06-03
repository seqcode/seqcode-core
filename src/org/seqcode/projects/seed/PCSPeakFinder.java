package org.seqcode.projects.seed;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.datasets.motifs.WeightMatrix;
import org.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import org.seqcode.gse.utils.sequence.SequenceUtils;
import org.seqcode.projects.multigps.utilities.Utils;
import org.seqcode.projects.seed.features.EnrichedFeature;
import org.seqcode.projects.seed.features.EnrichedPCSPeakFeature;
import org.seqcode.projects.seed.features.Feature;


/**
 * Permanganate-ChIP-seq peak-finder
 * @author mahony
 *
 */
public class PCSPeakFinder extends PeakFinder{
	public static String version = "0.1";
	
	protected final float bubbleWindow=30;
	protected int tagSeqWin = 20;
	protected float[][][] tagSeqComposition; //Sequence composition around tag 5' ends; indexed by Sample, then by relative position around tag 5' end, then by base (ACGT)

	public PCSPeakFinder(GenomeConfig gcon, ExptConfig econ, SEEDConfig scon, ExperimentManager man) {
		super(gcon, econ, scon, man);
		
		tagSeqComposition = new float[manager.getSamples().size()][tagSeqWin+1][4];
		for(int s=0; s<manager.getSamples().size(); s++)
			for(int i=0; i<=tagSeqWin; i++)
				for(int j=0; j<4; j++)
					tagSeqComposition[s][i][j]=0;
		
	}
	
	/**
	 * Return the class name
	 */
	public String getProgramName(){
		return "org.seqcode.projects.seed.PCSPeakFinder";
	}
	
	/**
	 * Return a thread for this implementation
	 */
	public FeatureDetectionThread getMyThread(List<Region> regs){
		return new PCSPeakFinderThread(regs);
	}
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.setProperty("java.awt.headless", "true");
		System.err.println("Permanganate-ChIP-seq Peak Finder version "+PeakFinder.version+"\n\n");
		
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		SEEDConfig scon = new SEEDConfig(gcon, args);
		
		if(scon.helpWanted()){
			//System.out.println(PeakFinder.getPCSPeakFinderArgs());
			System.err.println(gcon.getArgsList()+
					econ.getArgsList()+
					scon.getArgsList());
		}else{
			ExperimentManager man = new ExperimentManager(econ);
			PCSPeakFinder finder = new PCSPeakFinder(gcon, econ, scon, man);
			System.err.println("\nBeginning peak finding...");
			finder.execute();
			man.close();
		}
	}
	
	/**
	 * Return a help string
	 * @return
	 */
	public static String getPCSPeakFinderArgs() {
		//TODO
		return("TODO");
	}
	
	
	
	/**
	 * Multiple-hypothesis correction. Print the output files. 
	 * Assumes that all peaks have been found.
	 * @return : final features
	 */
	public Map<ExperimentCondition, List<Feature>> postProcess() {
		System.err.println("\nPeak finding complete.");
		
		//Normalize the tag sequence compositions
		for(int s=0; s<manager.getSamples().size(); s++)
			for(int i=0; i<=tagSeqWin; i++){
				float currTot=0;
				for(int j=0; j<4; j++)
					currTot+=tagSeqComposition[s][i][j];
				for(int j=0; j<4; j++)
					tagSeqComposition[s][i][j]/=currTot;
			}
		//Print tag sequence composition motifs
		for(Sample s : manager.getSamples()){
		    String fName = sconfig.getOutputParentDir()+File.separator+sconfig.getOutBase()+s.getName()+"-tag-sequence-motif.motif";
		    String imName = sconfig.getOutputParentDir()+File.separator+sconfig.getOutBase()+s.getName()+"-tag-sequence-motif.png";
		    String motifLabel = s.getName()+" tag-sequence-motif";
		    WeightMatrix wm = makeTagSeqCompositionWeightMatrix(s);
		    try {
		    	FileWriter fout = new FileWriter(fName);
		    	fout.write(WeightMatrix.printTransfacMatrix(wm, motifLabel));
		    	fout.close();
		    } catch (IOException e) {
		    	e.printStackTrace();
		    }
		    Utils.printMotifLogo(wm, new File(imName), 75, motifLabel); //Note that the painter converts the motif to log oddds
		}
		
		
		//Correct domains for multiple testing
		for(ExperimentCondition cond : manager.getConditions()){
			stats.benjaminiHochbergCorrection(features.get(cond));
		}
		
		Map<ExperimentCondition, List<Feature>> signifFeatures = this.filter(features, sconfig.perBinBinomialPThres, true);
       	
		//All domains
		this.printEventsFile(features, ".all.pcspeaks.txt");
        
		//Filtered by q-value
		this.printEventsFile(signifFeatures, ".p"+sconfig.perBinBinomialPThres+".pcspeaks.txt");
		
		//Summarize
		for(ExperimentCondition cond : manager.getConditions())
			System.err.println(cond.getName()+"\t"+features.get(cond).size()+" peaks\t"+signifFeatures.get(cond).size()+" peaks below threshold.");
		
		return features;
	}
	
	protected WeightMatrix makeTagSeqCompositionWeightMatrix(Sample s){
		WeightMatrix matrix = new WeightMatrix(tagSeqWin+1);
	    matrix.setNameVerType(s.getName()+"-tag-sequence-motif", "freq", "CUSTOM");
	    for (int i = 0; i <= tagSeqWin; i++) {
	    	//Assume normalized already
	    	matrix.matrix[i]['A'] = tagSeqComposition[s.getIndex()][i][0];
	    	matrix.matrix[i]['C'] = tagSeqComposition[s.getIndex()][i][1];
			matrix.matrix[i]['G'] = tagSeqComposition[s.getIndex()][i][2];
			matrix.matrix[i]['T'] = tagSeqComposition[s.getIndex()][i][3];
	    }matrix.setLogOdds();
		return matrix;
	}
	
	/**
	 * PeakFinderThread: thread that searches for domains
	 * @author mahony
	 *
	 */
	public class PCSPeakFinderThread extends PeakFinderThread {

		protected char[] currRegionSeq;
		protected char[] currRegionSeqRC;
		protected SequenceGenerator<Region> seqgen;
		
		public PCSPeakFinderThread(List<Region> regs) {
			super(regs);
			seqgen = gconfig.getSequenceGenerator();
		}
		
		/**
		 * findFeatures: 
		 * 	1) Count frequencies of bases around all tag 5' positions
		 * 	2) Find potentially enriched domains on the genome by comparing all tag counts to a background model,
		 *	3) Find the most likely transcriptional bubble positions within the domain
		 * 
		 * @param subRegion : region to run analysis on
		 * @return Map of Lists of Features in the subRegion, Indexed by ExperimentCondition 
		 */
		public Map<ExperimentCondition, List<Feature>> findFeatures(Region subRegion) {
			
			//Build the sequence model in the window around the tag 5' locations
			currRegionSeq = seqgen.execute(subRegion).toCharArray();
			currRegionSeqRC = currRegionSeq.clone();
			SequenceUtils.reverseComplement(currRegionSeqRC);
			
			//Tag sequence composition
			calculateTagSequenceComposition(subRegion);
			
			//Get the enriched domains using the functionality from the superclass (PeakFinderThread) 
			//and the over-ridden processDomain methods below 
			Map<ExperimentCondition, List<Feature>> peaks = super.findFeatures(subRegion);
			
			return peaks;
		}
		
		/**
		 * ProcessDomains: 
		 *  - Trims feature coordinates back to agree with overlapping hits. 
		 *  - Counts hits in each feature, per sample 
		 *  - Finds the peak & bubble position in the domain according to a choice of methods.
		 *  
		 * @param currFeatures
		 * @param current region
		 * @return : Lists of EnrichedFeatures, indexed by condition
		 */
		protected Map<ExperimentCondition, List<EnrichedFeature>> processDomains(Map<ExperimentCondition,List<EnrichedFeature>> currFeatures, Region currSubRegion){
			Map<ExperimentCondition, List<EnrichedFeature>> peakFeatures = new HashMap<ExperimentCondition, List<EnrichedFeature>>();
			for(ExperimentCondition cond : manager.getConditions())
				peakFeatures.put(cond, new ArrayList<EnrichedFeature>());
			
			for(ExperimentCondition currCondition : manager.getConditions()){
				for(EnrichedFeature currDomain : currFeatures.get(currCondition)){
					Map<Sample, List<StrandedBaseCount>> fHitsPos = overlappingHits(hitsPos, currDomain);
					Map<Sample, List<StrandedBaseCount>> fHitsNeg = overlappingHits(hitsNeg, currDomain);
					
					//Trim the coordinates
					trimFeature(currDomain, fHitsPos, fHitsNeg, currCondition);
					
					//Quantify the feature in each Sample and in the condition in which it was found
					quantifyFeature(currDomain, fHitsPos, fHitsNeg, currCondition);
					
					//Find the peaks
					EnrichedPCSPeakFeature pcspeak= findPCSBubble(fHitsPos, fHitsNeg, (EnrichedFeature)currDomain, currCondition, currSubRegion);
	                
					if(pcspeak!=null)
						peakFeatures.get(currCondition).add(pcspeak);
					
				}
			}
			return(peakFeatures);
		}
		
		/**
		 * Calculate the sequence composition around every tag 5' position in the current region, 
		 * and add the counts to the global models.
		 * 
		 *   Assumes that currRegionSeq and currRegionSeqRC have been properly initialized
		 * @param currReg
		 */
		protected void calculateTagSequenceComposition(Region currReg){
			
			float[][][] localTagSeqComposition = new float[manager.getSamples().size()][tagSeqWin+1][4];
			for(int s=0; s<manager.getSamples().size(); s++)
				for(int i=0; i<=tagSeqWin; i++)
					for(int j=0; j<4; j++)
						localTagSeqComposition[s][i][j]=0;
			int halfSeqWin = tagSeqWin/2;
			for(Sample s : manager.getSamples()){
				for(StrandedBaseCount sbc : hitsPos.get(s)){
					int w=0;
					for(int x=sbc.getCoordinate()-halfSeqWin-currReg.getStart(); x<=sbc.getCoordinate()+halfSeqWin-currReg.getStart(); x++){
					    if(x>=0 && x<currRegionSeq.length){
					    	int y = SequenceUtils.char2int(currRegionSeq[x]);
					    	if(y>=0)
						       localTagSeqComposition[s.getIndex()][w][y]+=sbc.getCount();
					    }
						w++;
					}
				}
				for(StrandedBaseCount sbc : hitsNeg.get(s)){
					int w=0;
					for(int x=currReg.getEnd()-sbc.getCoordinate()-halfSeqWin; x<=currReg.getEnd()-sbc.getCoordinate()+halfSeqWin; x++){
					    if(x>=0 && x<currRegionSeqRC.length){
					    	int y =	SequenceUtils.char2int(currRegionSeqRC[x]);
				    		if(y>=0)	   
				    			localTagSeqComposition[s.getIndex()][w][y]+=sbc.getCount();
					    }
						w++;
					}
				}
			}
			synchronized(tagSeqComposition){
				for(int s=0; s<manager.getSamples().size(); s++)
					for(int i=0; i<=tagSeqWin; i++)
						for(int j=0; j<4; j++)
							tagSeqComposition[s][i][j]+=localTagSeqComposition[s][i][j];
			}
		}
		
		/**
		 * Get the base preceding the 5' end of the tag
		 * @param a
		 * @param queryReg
		 * @return
		 */
		protected char getPrecedingBase(StrandedBaseCount a, Region queryReg){
			char b = '.';
			int wantedPos = a.getStrand()=='+' ? 
					a.getCoordinate()-1 : 
					a.getCoordinate()+1;
			if(wantedPos>=queryReg.getStart() && wantedPos<queryReg.getEnd()){
				b = a.getStrand()=='+' ? currRegionSeq[wantedPos-queryReg.getStart()] : currRegionSeqRC[queryReg.getEnd()-wantedPos];
			}
			return b;
		}
		
		/**
		 * Find the peak locations based on maximum overlapping read counts (T-preceding tags only). 
		 * Also quantify the number of tags and bases in the central bubble region, for computing the bubble index. 
	     * 
		 * @return : Feature (EnrichedPCSPeakFeature)
		 */
		protected EnrichedPCSPeakFeature findPCSBubble(Map<Sample, List<StrandedBaseCount>> fHitsPos, Map<Sample, List<StrandedBaseCount>> fHitsNeg, 
									EnrichedFeature domain, ExperimentCondition currCondition, Region currSubRegion){
			float [][] sum = new float[domain.getCoords().getWidth()+1][4];
			for(int s=0; s<=domain.getCoords().getWidth(); s++)
				for(int b=0; b<4; b++)
					sum[s][b]=0;

			for(Sample s : currCondition.getSignalSamples()){
				if(!strandedEventDetection || domain.getCoords().getStrand()=='+')
					for(StrandedBaseCount h : fHitsPos.get(s)){
						int pbase = SequenceUtils.char2int(getPrecedingBase(h, currSubRegion));
						if(pbase>=0){
							int start = getLeft(h)-domain.getCoords().getStart(); 
				            int stop= getRight(h)-domain.getCoords().getStart();
				            for(int i=start; i<stop; i++)
				                if(i>=0 && i<=sum.length)
				                    sum[i][pbase]+=h.getCount();
						}
					}
				if(!strandedEventDetection || domain.getCoords().getStrand()=='-')
					for(StrandedBaseCount h : fHitsNeg.get(s)){
						int pbase = SequenceUtils.char2int(getPrecedingBase(h, currSubRegion));
						if(pbase>=0){
							int start = getLeft(h)-domain.getCoords().getStart(); 
							int stop= getRight(h)-domain.getCoords().getStart();
				            for(int i=start; i<=stop; i++)
				                if(i>=0 && i<sum.length)
				                    sum[i][pbase]+=h.getCount();
						}
					}
			}
			//Find the peak according to T tags
			float max = 0; int maxPos = -1;
			for(int s=0; s<sum.length; s++){
				if(sum[s][3]>max){ //peaks should be determined only using T-preceding tags in permanganate-ChIP-seq
					max= sum[s][3];
					maxPos=s;
				}
			}
			Point p = new Point(gen, domain.getCoords().getChrom(), maxPos+domain.getCoords().getStart());
			
			//Calc the counts needed for the bubble index
			float[] bubbleTags = new float[4];
			float[] bubbleBases = new float[4];
			for(int b=0; b<4; b++){bubbleTags[b]=0; bubbleBases[b]=0;}
			int bubbleStart=(int)(maxPos-(bubbleWindow/2));
			int bubbleEnd=(int)(maxPos+(bubbleWindow/2));
			for(int z=bubbleStart; z<bubbleEnd; z++){
				if(z>=0 && z<=domain.getCoords().getWidth()){
					for(int b=0; b<4; b++){
						bubbleTags[b]+=sum[z][b];
					}
					bubbleBases[SequenceUtils.char2int(currRegionSeq[z])]++;
					bubbleBases[SequenceUtils.char2int(SequenceUtils.complementChar(currRegionSeq[z]))]++;
				}
			}
			
			//Define the feature
			EnrichedPCSPeakFeature epf=null;
			try {
				epf = new EnrichedPCSPeakFeature(domain.getCoords(), p, domain.getSampleCountsPos(), domain.getSampleCountsNeg(), 
																	domain.getSignalCount(), domain.getControlCount(), domain.getScore(), 
																	max, bubbleTags, bubbleBases);
			} catch (Exception e) {}
			return(epf);
		}
	}
}
