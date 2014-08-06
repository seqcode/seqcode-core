package edu.psu.compbio.seqcode.projects.akshay.bayesments.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.motifs.WeightMatrixScoreProfile;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.motifs.WeightMatrixScorer;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;


public class MotifPlatform {
	protected Config conf;
	protected List<Region> motif_search_locations = new ArrayList<Region>();
	protected List<Region> randomRegions;
	protected String[] randomSequences;
	protected MEMERunner meme;
	protected SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
	protected List<WeightMatrix> enriched_motifs;
	protected String outname;
	
	public MotifPlatform(Config c, List<Point> search_regions, String name){
		this.conf = c;

		for(Point p: search_regions){
			Region reg = p.expand(this.conf.getSeqWinSize());
			this.motif_search_locations.add(reg);
		}
		
		seqgen.useCache(true);
		seqgen.useLocalFiles(true);
		seqgen.setGenomePath(conf.getGenomeSeqPath());
		SequenceGenerator.setOffRegionCache();
		if(!seqgen.isRegionCached()){
			System.err.println("Caching sequences");
			randomRegions = randomRegionPick(null, conf.MOTIF_FINDING_NEGSEQ, conf.getSeqWinSize());
			randomSequences = seqgen.setupRegionCache(this.motif_search_locations, randomRegions);
			System.err.println("Caching completed");
		}
		
		this.meme = new MEMERunner(this.conf);
		this.outname = name;
	}
	
	//Accesors
	
	public List<WeightMatrix> getMotifs(){return this.enriched_motifs;}
	

	public void findMotifs(){
		List<String> seqs = new ArrayList<String>();
		for(int i=0; i<this.motif_search_locations.size(); i++){
			String currSeq = seqgen.execute(this.motif_search_locations.get(i));
			if(lowercaseFraction(currSeq)<=conf.MOTIF_FINDING_ALLOWED_REPETITIVE){
				seqs.add(currSeq);
			}
		}
		
		//Execute MEME
		Pair<List<WeightMatrix>,List<WeightMatrix>> matrices = meme.execute(seqs, this.outname, false);
		List<WeightMatrix> wm = matrices.car();
		List<WeightMatrix> fm = matrices.cdr();
		
		if(wm.size()>0){
			//Evaluate the significance of the discovered motifs
			double rocScores[] = motifROCScores(wm,seqs,this.randomSequences);
			System.err.println("MEME results for:" );
			for(int w=0; w<fm.size(); w++){
				if(fm.get(w)!=null){
					System.err.println("\t"+fm.get(w).getName()+"\t"+ WeightMatrix.getConsensus(fm.get(w))+"\tROC:"+String.format("%.2f",rocScores[w]));
				}
			}
			
			double maxRoc = 0.0;
			for(int i=0; i<rocScores.length; i++){
				if(rocScores[i]>maxRoc){
					maxRoc = rocScores[i]; 
				}
			}
			
			if(maxRoc >= conf.MOTIF_MIN_ROC){
				this.enriched_motifs = new ArrayList<WeightMatrix>();
				for(int m=0; m< rocScores.length; m++){
					if(rocScores[m]> conf.MOTIF_MIN_ROC){
						this.enriched_motifs.add(wm.get(m));
					}
				}
			}else{
				this.enriched_motifs = null;
			}
		}
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
	public class LabeledDouble implements Comparable<LabeledDouble>{
		public Double dat;
		public Integer label;
		public LabeledDouble(Double d, Integer i){dat=d; label=i;}
		public int compareTo(LabeledDouble ld) {
			if(dat > ld.dat){return(-1);}
			else if(dat < ld.dat){return(1);}
			else{return 0;}
		}
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
		long [] chromoSize = new long[conf.getGenome().getChromList().size()];
		String [] chromoNames = new String[conf.getGenome().getChromList().size()];
		Iterator<NamedRegion> chroms = new ChromRegionIterator(conf.getGenome());
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
						potential = new Region(conf.getGenome(), chromoNames[c], (int)(randPos-total), (int)(randPos+sampleSize-total));
						
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
}
