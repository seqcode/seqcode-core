package org.seqcode.motifs.scores2motifs;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.RepeatMaskedRegion;
import org.seqcode.genome.location.ScoredPoint;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gseutils.Pair;
import org.seqcode.motifs.MemeER;


public class Scores2Motifs {
	
	protected S2MConfig s2mConfig;
	protected List<Hill> hills;
	protected String[] randomSequences = new String[MemeER.MOTIF_FINDING_NEGSEQ];
	protected File basedir_profiles;
	protected File basedir_meme;	
	protected List<Integer> clusterAssignment = new ArrayList<Integer>();
	protected int bestNumClusters = 2;
	
	/** List of all discriminative motifs for all the labels*/ 
	protected List<WeightMatrix> discrimMotifs = new ArrayList<WeightMatrix>();
	
	/** The ROC value of the identified motifs that indicates their significance */
	protected HashMap<String, Double> discrimMotifsRocs = new HashMap<String,Double>();
	
	

	/**
	 * Constructor
	 * @param h2mcon
	 */
	public Scores2Motifs(S2MConfig h2mcon) {
		s2mConfig = h2mcon;
		basedir_profiles = new File(s2mConfig.getOutDir().getAbsoluteFile()+File.separator+"kmer_profiles");
		basedir_meme = new File(s2mConfig.getOutDir().getAbsoluteFile()+File.separator+File.separator+"meme");
		
	}
	
	/**
	 * execute: Run the method
	 * @throws IOException
	 */
	public void execute() throws IOException{
		makeMemeDirs();
		
		//Load a set of random regions
		setRandomRegs();
		
		//Define the hills according to one of the input options
		setHills();
		
		//Cluster hill profiles
		if(hills.size()>100){
			ClusterProfiles clusterManager = new ClusterProfiles(hills, S2MConfig.ITRS_CLUS,s2mConfig.getKmin(),s2mConfig.getKmax(),basedir_profiles);
			clusterAssignment = clusterManager.execute();
			bestNumClusters = clusterManager.getNumClusters();
		}
		
		//Run motif-finding
		runMEME();
		
	}

	
	/**
	 * setHills:
	 * 	There is currently only one way to define hills; by reading in a list of scored binding 
	 * 	locations from a file.
	 * 	However, in the future we will want to extend this to determining hills from a scored 
	 * 	landscape, or even from a model + sequence
	 * @return
	 */
	protected List<Hill> setHills(){
		hills = new ArrayList<Hill>();
		
		if(s2mConfig.getPeaksFile() != null){
			//Load the scored points from a file
			List<ScoredPoint> peaks = processScoredPoints();

			for(ScoredPoint sp : peaks){
				//Expand to region
				Region searchWin = sp.expand(s2mConfig.getMemeSearchWin());
				//Load sequences for each point
				String seq = s2mConfig.getSeqGen().execute(searchWin);
				//Define the Hill objects
				hills.add(new Hill(searchWin, sp.getScore(), seq));
			}
		}
		
		//Add k-mer profiles to hills
		List<int[]> profiles = getProfilesAtHills(hills);
		for(int i=0; i<profiles.size(); i++)
			hills.get(i).setKmerProfile(profiles.get(i));
		
		//sort by score
		Collections.sort(hills);		
		
		return hills;
	}
	
	
	
	/**
	 * Read the scored points from a file and optionally screen repeats. 
	 * @return
	 */
	protected List<ScoredPoint> processScoredPoints(){
		if(s2mConfig.getPeaksFile() != null){
			// Loading peaks
			List<ScoredPoint> sp = RegionFileUtilities.loadScoredPointsFromFile(s2mConfig.getGenome(), s2mConfig.getPeaksFile());
			
			//Filter repeats if appropriate
			if(s2mConfig.getScreenReps()){
				if(sp.size() > 0){

					Iterator<ScoredPoint> peakItr = sp.iterator();
					while(peakItr.hasNext()){
						ScoredPoint currPeak = peakItr.next();
						Region currReg = currPeak.expand(s2mConfig.getRepMaskWin());
						double repLen = 0;
						Iterator<RepeatMaskedRegion> repItr = s2mConfig.getRepMask().execute(currReg);
						while (repItr.hasNext()) {
							RepeatMaskedRegion currRep = repItr.next();
							if (currRep.overlaps(currReg)) {
								repLen += (double) currRep.getOverlapSize(currReg);
							}
						}
						if (repLen / (double) currReg.getWidth() > s2mConfig.getRepPropLimit()){
							peakItr.remove();
						}
					}
				}
			}
			return sp;
		}
		return null;
	}
	
	/**
	 *  Generate random sequences needed for MEME analysis to assess motifs
	 */
	protected void setRandomRegs(){
		List<Region> randomRegions = MemeER.randomRegionPick(s2mConfig.getGenome(), null, MemeER.MOTIF_FINDING_NEGSEQ,s2mConfig.getMemeSearchWin());
		for(int i=0; i<randomRegions.size(); i++){
			randomSequences[i] = s2mConfig.getSeqGen().execute(randomRegions.get(i));
		}
	}
	

	/**
	 * Make k-mer profiles from hill regions
	 * @param rs
	 * @return
	 */
	protected List<int[]> getProfilesAtHills(List<Hill> hills){
		List<int[]> ret = new ArrayList<int[]>();

		for(Hill h: hills){
			int[] pfl = new int[s2mConfig.getNumK()];
			String seq = h.getSeq().toUpperCase();
			for(int k=s2mConfig.getKmin(); k<=s2mConfig.getKmax(); k++){
				for(int i=0; i<(seq.length()-k+1); i++){
					String currk = seq.substring(i, i+k);
					String revcurrk = SequenceUtils.reverseComplement(currk);
					int currKInt = RegionFileUtilities.seq2int(currk);
					int revcurrKInt = RegionFileUtilities.seq2int(revcurrk);
					int kmer = currKInt<revcurrKInt ? currKInt : revcurrKInt;
					int baseind = s2mConfig.getKmerBaseInd(currk);	
					pfl[baseind+kmer]++;
				}
			}
			ret.add(pfl);
		}
		return ret;
	}
	
	
	/**
	 * Make MEME directories
	 */
	protected void makeMemeDirs(){
		File memeDir = new File(s2mConfig.getOutDir().getAbsoluteFile()+File.separator+"meme");
				memeDir.mkdirs();
	}
	
	/**
	 * Run motif-finding with MEME
	 */
	protected void runMEME(){
		System.err.println("Finished clustering K-mer profiles; now running MEME");

		// Now do meme search on each clusters separately
		String memeargs = s2mConfig.getMemeArgs();
		MemeER meme = new MemeER(s2mConfig.getMemePath(), memeargs);
		meme.setMotifMinROC(s2mConfig.getMotifMinROC());
		for(int c=0; c<bestNumClusters; c++){ // Over each cluster
			System.err.println("Loading sequences for meme analysis : Cluster"+c);
			int numHillsLoaded = 0;
			List<String> seqs = new ArrayList<String>();
			for(int p=0; p<hills.size(); p++){
				if(clusterAssignment.get(p) == c && numHillsLoaded < S2MConfig.NUM_HILLS){
					numHillsLoaded++;
					String s = hills.get(p).getSeq();
					if(MemeER.lowercaseFraction(s)<= S2MConfig.MOTIF_FINDING_ALLOWED_REPETITIVE){
						seqs.add(s);
					}
				}
				if(numHillsLoaded >= S2MConfig.NUM_HILLS)
					break;
			}
			//List<WeightMatrix> selectedMotifs = new ArrayList<WeightMatrix>();
			File meme_outFile = new File(basedir_meme+File.separator+"Cluster-"+Integer.toString(c)+"_meme");

			System.err.println("Running meme now");
			Pair<List<WeightMatrix>,List<WeightMatrix>> matrices = meme.execute(seqs, meme_outFile, false);
			System.err.println("Finished running meme, Now evaluating motif significance");
			List<WeightMatrix> wm = matrices.car();
			List<WeightMatrix> fm = matrices.cdr();
			// Index for all the selected motifs
			int motInd = 0;
			
			if(wm.size()>0){
				//Evaluate the significance of the discovered motifs
				double rocScores[] = meme.motifROCScores(wm,seqs,randomSequences);
				System.err.println("MEME results for:" );
				for(int w=0; w<fm.size(); w++){
					if(fm.get(w)!=null){
						System.err.println("\t"+fm.get(w).getName()+"\t"+ WeightMatrix.getConsensus(fm.get(w))+"\tROC:"+String.format("%.2f",rocScores[w]));
					}
					if(rocScores[w] > meme.getMotifMinROC()){
						//selectedMotifs.add(fm.get(w));
						fm.get(w).setName("c"+Integer.toString(c)+"_"+Integer.toString(motInd));
						discrimMotifs.add(fm.get(w));
						discrimMotifsRocs.put("c"+Integer.toString(c)+"_"+Integer.toString(motInd), rocScores[w]);
					}
					motInd++;
				}
			}
		}
	}

}
