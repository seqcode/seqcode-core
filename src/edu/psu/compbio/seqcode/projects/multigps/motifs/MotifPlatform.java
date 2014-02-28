package edu.psu.compbio.seqcode.projects.multigps.motifs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScorer;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.mixturemodel.BindingComponent;

public class MotifPlatform {
	protected Config config;
	protected ExperimentManager manager;
	protected SequenceGenerator<Region> seqgen;
	protected List<Region> randomRegions = new ArrayList<Region>(); //Randomly chosen regions for motif significance tests
	protected String[] randomSequences; //Randomly chosen sequences for motif significance tests
	protected MEMERunner meme;
		
	/**
	 * Constructor for motif platform
	 * @param c
	 * @param man
	 * @param regionsOfInterest: this list contains all possible regions that motif-finding/scanning may be run on. 
	 */
	public MotifPlatform(Config c, ExperimentManager man, List<Region> regionsOfInterest){
		config = c;
		manager=man;
		seqgen = new SequenceGenerator<Region>();
		seqgen.useCache(true);
		seqgen.useLocalFiles(true);
		seqgen.setGenomePath(config.getGenomeSequencePath());
		if (!seqgen.isRegionCached()){
			System.err.println("Caching sequences");
			randomRegions = randomRegionPick(null, config.MOTIF_FINDING_NEGSEQ, config.MOTIF_FINDING_SEQWINDOW);
			randomSequences = seqgen.setupRegionCache(regionsOfInterest, randomRegions);
			System.err.println("Caching completed");
		}
		meme = new MEMERunner(config, man);
	}

	/**
	 * Extract sequences around top BindingComponents and call the MEME runner 
	 * @param cond
	 * @param activeComponents
	 * @param trainingRound
	 */
	public void findMotifs(ExperimentCondition cond, HashMap<Region, List<List<BindingComponent>>> activeComponents, int trainingRound){
		List<BindingComponent> peaks = new ArrayList<BindingComponent>();
		//Choose which components to include
		for(Region r : activeComponents.keySet()){
			for(BindingComponent bc : activeComponents.get(r).get(cond.getIndex())){
				//Component must not be at the edge of the region 
				if((bc.getPosition()-r.getStart()>config.MOTIF_FINDING_SEQWINDOW/2) && (r.getEnd()-bc.getPosition()>config.MOTIF_FINDING_SEQWINDOW/2)){
					peaks.add(bc);
		}}}
		//Sort by responsibilities
		Collections.sort(peaks, new Comparator<BindingComponent>(){
            public int compare(BindingComponent o1, BindingComponent o2) {return o1.compareByResp(o2);}
        });
		Collections.reverse(peaks);
		
		//Get the top sequences that don't have too many lowercase letters (repeats)
		List<String> seqs = new ArrayList<String>();
		int addedSeqs=0;
		for(int b=0; b<peaks.size() && addedSeqs<config.MOTIF_FINDING_TOPSEQS; b++){
			BindingComponent bc = peaks.get(b);
			Region peakReg = new Region(bc.getCoord().getGenome(), bc.getCoord().getChrom(), bc.getPosition()-config.MOTIF_FINDING_SEQWINDOW/2, bc.getPosition()+config.MOTIF_FINDING_SEQWINDOW/2); 
			String currSeq = seqgen.execute(peakReg);
			if(lowercaseFraction(currSeq)<=config.MOTIF_FINDING_ALLOWED_REPETITIVE){
				seqs.add(currSeq);
				addedSeqs++;
			}
		}
		//Execute MEME
		Pair<List<WeightMatrix>,List<WeightMatrix>> matrices = meme.execute(seqs, new String("motif_"+cond.getName()+"_t"+trainingRound), false);
		List<WeightMatrix> wm = matrices.car();
		List<WeightMatrix> fm = matrices.cdr();
		
		if(wm.size()>0){
			//Evaluate the significance of the discovered motifs
			int bestMotif=0;
			double rocScores[] = motifROCScores(wm, seqs, randomSequences);
			double maxRoc=0;
			for(int i=0; i<rocScores.length; i++)
				if(rocScores[i]>maxRoc){
					maxRoc = rocScores[i]; 
					bestMotif=i;
				}
			//Results summary
			if(config.isVerbose()){
				System.err.println("MEME results for: "+cond.getName());
				for(int w=0; w<fm.size(); w++){
					if(fm.get(w)!=null){
						System.err.println("\t"+fm.get(w).getName()+"\t"+ WeightMatrix.getConsensus(fm.get(w))+"\tROC:"+String.format("%.2f",rocScores[w]));
					}
				}
			}
			
			//Set the condition's motif if the ROC is above threshold
			if(maxRoc >= config.MOTIF_MIN_ROC){
				System.err.println("\t"+fm.get(bestMotif).getName() + " chosen as best motif.");
				cond.setMotif(wm.get(bestMotif));
				cond.setFreqMatrix(fm.get(bestMotif));
			}else{
				System.err.println("\tNo motif passes minimum ROC score threshold.");
				cond.setMotif(null);
				cond.setFreqMatrix(null);
			}
		}else{
			cond.setMotif(null);
			cond.setFreqMatrix(null);
		}
	}
	
	/**
	 * Align motifs to get relative offsets
	 */
	public void alignMotifs(){
		SimpleMotifAligner aligner = new SimpleMotifAligner(6);
		
		//Condition 0's motif is usually the reference
		WeightMatrix refMotif=null;
		int refCondIndex=0;
		for(int e=0; e<manager.getNumConditions(); e++){
			if(manager.getExperimentSet().getIndexedCondition(e).getMotif()!=null){
				refMotif = manager.getExperimentSet().getIndexedCondition(e).getFreqMatrix();
				refCondIndex = e;
				break;
			}
		}
		if(refMotif!=null){
			//Reference offset is the center of the motif
			int refOffset = (int)(refMotif.length()/2);
			manager.getExperimentSet().getIndexedCondition(refCondIndex).setMotifOffset(refOffset);
			for(int e=0; e<manager.getNumConditions(); e++){
				if(e!=refCondIndex && manager.getExperimentSet().getIndexedCondition(e).getFreqMatrix()!=null){
					Pair<Integer,Double> forAlignment = aligner.align(refMotif, manager.getExperimentSet().getIndexedCondition(e).getFreqMatrix());
					Pair<Integer,Double> revAlignment = aligner.align(refMotif, WeightMatrix.reverseComplement(manager.getExperimentSet().getIndexedCondition(e).getFreqMatrix()));
					
					if(revAlignment.cdr()>forAlignment.cdr()){
						manager.getExperimentSet().getIndexedCondition(e).setFreqMatrix(WeightMatrix.reverseComplement(manager.getExperimentSet().getIndexedCondition(e).getFreqMatrix()));
						manager.getExperimentSet().getIndexedCondition(e).setMotif(WeightMatrix.reverseComplement(manager.getExperimentSet().getIndexedCondition(e).getMotif()));
						manager.getExperimentSet().getIndexedCondition(e).setMotifOffset(refOffset+revAlignment.car());
					}else{
						manager.getExperimentSet().getIndexedCondition(e).setMotifOffset(refOffset+forAlignment.car());
					}
				}else if (manager.getExperimentSet().getIndexedCondition(e).getFreqMatrix()==null){
					manager.getExperimentSet().getIndexedCondition(e).setMotifOffset(0);
				}
			}
		}
	}
	
	/**
	 * Scan the region with the motifs from each condition. Suitable for making a motif prior.  
	 * @param reg
	 * @return array of array of motif scores. Indexed by condition. 
	 */
	public double[][] scanRegionWithMotifs(Region reg, String regSeq){
		double[][] scanScores = new double[manager.getNumConditions()][reg.getWidth()];
		
		for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
			int e = cond.getIndex();
			for(int z=0; z<reg.getWidth(); z++)
				scanScores[e][z]=0;
			
			WeightMatrix motif = cond.getMotif();
			if(motif!=null){
				int motifOffset = cond.getMotifOffset();
				WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
				WeightMatrixScoreProfile profiler = scorer.execute(regSeq);
				for(int z=0; z<reg.getWidth()-motif.length()+1; z++){
					double currScore= profiler.getMaxScore(z);
					int zOff = z+motifOffset;
					if(currScore>0){
						if(zOff>0 && zOff<reg.getWidth())
							scanScores[e][zOff] = currScore;
					}
				}
			}
		}
		return scanScores;
	}
	
	/**
	 * Scan the region with the motifs from each condition.   
	 * @param reg
	 * @return Pair of: array of array of motif max seqs and array of array of motif scores. Indexed by condition. 
	 */
	public Pair<Double[][],String[][]> scanRegionWithMotifsGetSeqs(Region reg, String regSeq){
		Double[][] scanScores = new Double[manager.getNumConditions()][reg.getWidth()];
		String[][] scanSeqs = new String[manager.getNumConditions()][reg.getWidth()];
		
		for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
			int e = cond.getIndex();
			for(int z=0; z<reg.getWidth(); z++){
				scanScores[e][z]=0.0;
				scanSeqs[e][z]= "";
			}
			
			WeightMatrix motif = cond.getMotif();
			if(motif!=null){
				int motifOffset = cond.getMotifOffset();
				WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
				WeightMatrixScoreProfile profiler = scorer.execute(regSeq);
				for(int z=0; z<reg.getWidth()-motif.length()+1; z++){
					double currScore= profiler.getMaxScore(z);
					String currSeq = regSeq.substring(z, z+motif.length());
					if(profiler.getMaxStrand(z)=='-')
						currSeq = SequenceUtils.reverseComplement(currSeq);
					int zOff = z+motifOffset;
					if(zOff>0 && zOff<reg.getWidth()){
						if(currScore>0)
							scanScores[e][zOff] = currScore;
						scanSeqs[e][zOff] = currSeq;
					}
				}
			}
		}
		Pair<Double[][],String[][]> scoresAndSeqs = new Pair<Double[][],String[][]>(scanScores, scanSeqs);
		return scoresAndSeqs;
	}
	
	/**
	 * Get a sequence for a given region
	 */
	public String getSeq(Region reg){
		String regSeq = seqgen.execute(reg);
		return regSeq;
	}
	
	/**
	 * Scan the region with the motif from one condition. Suitable for making a motif prior.  
	 * @param cond
	 * @param reg
	 * @return array of array of motif scores. Indexed by condition. 
	 */
	public double[] scanRegionWithMotif(Region reg, ExperimentCondition cond){
		double[] scanScores = new double[reg.getWidth()];
		String regSeq = seqgen.execute(reg);
		for(int z=0; z<reg.getWidth(); z++)
			scanScores[z]=0;
		
		WeightMatrix motif = cond.getMotif();
		if(motif!=null){
			int motifOffset = cond.getMotifOffset();
			WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
			WeightMatrixScoreProfile profiler = scorer.execute(regSeq);
			for(int z=0; z<reg.getWidth()-motif.length()+1; z++){
				double currScore= profiler.getMaxScore(z);
				int zOff = z+motifOffset;
				if(currScore>0){
					if(zOff>0 && zOff<reg.getWidth())
						scanScores[zOff] = currScore;
				}
			}
		}
		return scanScores;
	}
	
	/**
	 * Compute the fraction of letters in the sequence that are lowercase or N
	 * @param seq
	 * @return
	 */
	protected double lowercaseFraction(String seq){
		double count = 0;
		for (char c:seq.toCharArray())
			if (Character.isLowerCase(c) || c=='N')
				count++;
		return count/(double)seq.length();
	}
	/**
	 * Randomly pick a set of Regions
	 * @param gen
	 * @param blackList
	 * @param numSamples
	 * @param sampleSize
	 * @return
	 */
	protected List<Region> randomRegionPick(List<Region> blackList, int numSamples, int sampleSize){
		List<Region> regs = new ArrayList<Region>();
		Random rand = new Random();
		int validSamples=0;
		
		//First see how big the genome is:
		int numChroms=0;
		long genomeSize=0;
		long [] chromoSize = new long[config.getGenome().getChromList().size()];
		String [] chromoNames = new String[config.getGenome().getChromList().size()];
		Iterator<NamedRegion> chroms = new ChromRegionIterator(config.getGenome());
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			genomeSize += (double)currentChrom.getWidth();
			chromoSize[numChroms]=currentChrom.getWidth();
			chromoNames[numChroms]=currentChrom.getChrom();
			numChroms++;				
		}

		//Now, iteratively generate random positions and check if they are valid and not overlapping repeats. 
		while(validSamples<numSamples){
			Region potential;				
			long randPos = (long)(1+(rand.nextDouble()*genomeSize));
			//find the chr
			boolean found=false;
			long total=0;
			for(int c=0; c<numChroms && !found; c++){
				if(randPos<total+chromoSize[c]){
					found=true;
					if(randPos+sampleSize<total+chromoSize[c]){
						potential = new Region(config.getGenome(), chromoNames[c], (int)(randPos-total), (int)(randPos+sampleSize-total));
						
						//is this region in the blacklist? 
						boolean valid=true;
						if(blackList!=null){
							for(Region r : blackList){
								if(potential.overlaps(r)){valid=false;}
							}
						}
						if(valid){
							validSamples++;
							regs.add(potential);
						}
					}
				}total+=chromoSize[c];
			}
		}
		return(regs);
	}
	
	/**
	 * Calculate ROC scores for all motifs. 
	 * @param matrices
	 * @param posSeqs
	 * @param negSeqs
	 * @return
	 */
	protected double[] motifROCScores(List<WeightMatrix> matrices, List<String> posSeqs, String[] negSeqs){
		double[] rocScores = new double[matrices.size()];
		int m=0;
		for(WeightMatrix motif : matrices){
			List<Double> posScores = new ArrayList<Double>();
			List<Double> negScores = new ArrayList<Double>();
			if(motif!=null){
				WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
				for(String posSeq : posSeqs){
					WeightMatrixScoreProfile profiler = scorer.execute(posSeq);
					posScores.add(profiler.getMaxScore());
				}
				for(int s=0; s<negSeqs.length; s++){
					WeightMatrixScoreProfile profiler = scorer.execute(negSeqs[s]);
					negScores.add(profiler.getMaxScore());
				}
			}
			rocScores[m] = calcROCAUC(posScores, negScores);
			m++;
		}
		return rocScores;
	}
	
	/**
	 * Calculate the area under a motif-scoring ROC
	 * @param posMaxScores
	 * @param negMaxScores
	 * @param printROC
	 * @param motif
	 * @return
	 */
	private double calcROCAUC(List<Double> posMaxScores, List<Double> negMaxScores) {
		double auc = 0;
		if(posMaxScores.size()==0)
			return 0;
		if(negMaxScores.size()==0)
			return 1;
		ArrayList<LabeledDouble> data = new ArrayList<LabeledDouble>();
		for(Double d : posMaxScores)
			data.add(new LabeledDouble(d, 1));
		for(Double d : negMaxScores)
			data.add(new LabeledDouble(d, 0));
		
		Collections.sort(data);
		double pCount = (double)posMaxScores.size();
		double nCount = (double)negMaxScores.size();
		int x=0;
		double possum=0;
		double lastsn=0;
		double lastfpr=0;
		double lastdval = 10000000;
		
		for(LabeledDouble d : data){
			possum+=d.label;
			if(d.dat!=lastdval){
				double sn = possum/pCount;
				double fp = (x+1)-possum;
				double sp = (nCount-fp)/nCount;
				double fpr=1-sp;
				if(x>0){
						    //Rectangle             //Triangle
					auc += ((fpr-lastfpr)*lastsn) + ((sn-lastsn)*(fpr-lastfpr)/2);
				}
				lastfpr=fpr;
				lastsn = sn;
			}
			lastdval = d.dat;
			x++;
		}
		return auc;
	}
	/**
	 * Simple class for ROC analysis
	 * @author Shaun Mahony
	 * @version	%I%, %G%
	 */
	protected class LabeledDouble implements Comparable<LabeledDouble>{
		public Double dat;
		public Integer label;
		public LabeledDouble(Double d, Integer i){dat=d; label=i;}
		public int compareTo(LabeledDouble ld) {
			if(dat > ld.dat){return(-1);}
			else if(dat < ld.dat){return(1);}
			else{return 0;}
		}
	}
}
