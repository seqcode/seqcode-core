package edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder.motifs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder.clusterkmerprofile.ClusterProfiles;
import edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder.framework.SeqUnwinderConfig;
import edu.psu.compbio.seqcode.projects.akshay.utils.MemeER;
import edu.psu.compbio.seqcode.projects.multigps.utilities.Utils;

public class Discrim {
	
	protected SeqUnwinderConfig seqConfig;
	/** List of all discriminative motifs for all the lables*/ 
	protected List<WeightMatrix> discrimMotifs = new ArrayList<WeightMatrix>();
	
	/** The ROC value of the identified motifs that indicates their significance */
	protected HashMap<String, Double> discrimMotifsRocs = new HashMap<String,Double>();
	protected String[] randomSequences = new String[5000];
	
	public void setRandomRegs(){
		// Generate random sequences later need for meme analysis to assess motifs
		
		List<Region> randomRegions = MemeER.randomRegionPick(seqConfig.getGenome(), null, MemeER.MOTIF_FINDING_NEGSEQ,seqConfig.getWin());
		for(int i=0; i<randomRegions.size(); i++){
			randomSequences[i] = seqConfig.getSeqGen().execute(randomRegions.get(i));
		}
	}

	public Discrim(SeqUnwinderConfig seqcon) {
		seqConfig = seqcon;
	}
	
	public void execute() throws IOException{
		setRandomRegs();
		for(String s : seqConfig.getKmerWeights().keySet()){
			if(!s.equals("RootRandom")){
				KmerModelScanner scanner = new KmerModelScanner(s);
				scanner.execute();
			}
		}
		// Update Config
		seqConfig.setDiscrimMotifs(discrimMotifs);
		seqConfig.setDiscrimMotifsRocs(discrimMotifsRocs);
	}

	
	//Slave methods
	public void makeMemeDirs(){
		for(String s : seqConfig.getModelNames()){
			if(!s.equals("RootRandom")){
				File memeDir = new File(seqConfig.getOutDir().getAbsoluteFile()+File.separator+s+File.separator+"meme");
				memeDir.mkdirs();
				File kmerProfDir = new File(seqConfig.getOutDir().getAbsoluteFile()+File.separator+s+File.separator+"kmer_profiles");
				kmerProfDir.mkdirs();
			}
		}

	}

	
	public class KmerModelScanner {
		//Model specific features
		protected String kmerModelName;
		protected List<Region> modelRegions = new ArrayList<Region>();
		
		protected List<Region> posHills = new ArrayList<Region>();
		protected ArrayList<int[]> profiles = new ArrayList<int[]>();
		protected List<Integer> clusterAssignment = new ArrayList<Integer>();
		protected HashMap<Integer,String> posHillsToIndex =new HashMap<Integer,String>();
		protected HashMap<Integer,Double> posHillScores =new HashMap<Integer,Double>();
		protected File basedir_profiles;
		protected File basedir_meme;
		
		/** 
		 * Constructor
		 */
		public KmerModelScanner(String modName) {
			kmerModelName = modName;
			basedir_profiles = new File(seqConfig.getOutDir().getAbsoluteFile()+File.separator+kmerModelName+File.separator+"kmer_profiles");
			basedir_meme = new File(seqConfig.getOutDir().getAbsoluteFile()+File.separator+kmerModelName+File.separator+"meme");
			
			// Now add all the training examples that are linked to the current K-mer model
			// First; split to get subgroup names
			String[] subGsMod = kmerModelName.split("#");
			
			for(int p=0; p< seqConfig.getRegions().size(); p++){
				String[] pieces = seqConfig.getPeakAnnotations().get(p).split(";");
				boolean addReg = true;
				
				for(String subGsInMod : subGsMod){
					boolean SGinPeak = false;
					for(String piece : pieces){
						if(subGsInMod.equals(piece)){
							SGinPeak = true;
						}
					}
					if(!SGinPeak)
						addReg = false;
				}
				
				if(addReg)
					modelRegions.add(seqConfig.getRegions().get(p));
			}
			
		}
		
		
		public void execute() throws IOException{
			// First, find, cluster and print the hills
			clusterKmerProfilesAtMountains();
			if(posHills.size() > 100){
				System.err.println("Finished clustering K-mer profiles for model: "+kmerModelName);

				// Now do meme search on each clusters separately
				String memeargs = seqConfig.getMemeArgs();
				MemeER meme = new MemeER(seqConfig.getMemePath(), memeargs);
				for(int c=0; c<seqConfig.getNumDiscrimClusters(); c++){ // Over each cluster
					System.err.println("Loading sequences for meme analysis : "+kmerModelName+ "Cluster"+c);
					List<String> seqs = new ArrayList<String>();
					for(int p=0; p<posHills.size(); p++){
						if(clusterAssignment.get(p) == c){
							Region tmpR = posHills.get(p).getMidpoint().expand(seqConfig.getMemeSearchWin()/2);
							String s = seqConfig.getSeqGen().execute(tmpR);
							if(MemeER.lowercaseFraction(s)<= SeqUnwinderConfig.MOTIF_FINDING_ALLOWED_REPETITIVE){
								seqs.add(s);
							}
						}
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
							if(rocScores[w] > MemeER.MOTIF_MIN_ROC){
								//selectedMotifs.add(fm.get(w));
								fm.get(w).setName(kmerModelName+"_c"+Integer.toString(c)+"_"+Integer.toString(motInd));
								discrimMotifs.add(fm.get(w));
								discrimMotifsRocs.put(kmerModelName+"_c"+Integer.toString(c)+"_"+Integer.toString(motInd), rocScores[w]);
							}
							motInd++;
						}
					}
				}
			}
		
		}
		
		
		//Settors
		
		
		//Fillers
		private void fillHills() throws IOException{
			List<Pair<Region,Double>> hills= new ArrayList<Pair<Region,Double>>();
			hills= findHills();
			List<Pair<Region,Double>> top_hills = sorthills(hills);
			int index=0;
			HashMap<String,Integer> addedHills = new HashMap<String,Integer>();
			for(Pair<Region,Double> pr : top_hills){
				if(!addedHills.containsKey(pr.car().getLocationString())){
					posHills.add(pr.car());
					posHillsToIndex.put(index,pr.car().getLocationString() );
					posHillScores.put(index++, pr.cdr());
					addedHills.put(pr.car().getLocationString(), 1);
				}
			}
			profiles = getProfilesAtPeaks(posHills);
		}
		
		private List<Pair<Region,Double>> sorthills(List<Pair<Region,Double>> hlls){
			List<Pair<Region,Double>> ret = new ArrayList<Pair<Region,Double>>();
			HillsIndexComparator comp = new HillsIndexComparator(hlls);
			Integer[] indicies = comp.createIndexArray();
			Arrays.sort(indicies, comp);
			
			for(int i=0; i<Math.min(hlls.size(),SeqUnwinderConfig.NUM_HILLS); i++){
				ret.add(hlls.get(indicies[i]));
			}
			
			return ret;
		}
		
		/**
		 * Prints the mountain composition
		 * @param useCache
		 * @param genpath
		 * @param oddsThresh
		 */
		public void clusterKmerProfilesAtMountains() throws IOException{
			fillHills();
			System.err.println("No of hills for K-mer model "+kmerModelName+" are :"+posHills.size());
			if(posHills.size() > 100){
				ClusterProfiles clusterManager = new ClusterProfiles(SeqUnwinderConfig.ITRS_CLUS,seqConfig.getNumDiscrimClusters(),profiles,posHillsToIndex,seqConfig.getMinM(),seqConfig.getMaxM(),posHillScores,basedir_profiles);
				clusterAssignment = clusterManager.execute();
			}
		}
		
		private ArrayList<int[]> getProfilesAtPeaks(List<Region> rs){
			ArrayList<int[]> ret = new ArrayList<int[]>();
			
			for(Region r: rs){
				int numK = 0;
				for(int k=seqConfig.getKmin(); k<seqConfig.getKmax(); k++){
					numK += (int)Math.pow(4, k);
				}
				int[] pfl = new int[numK];
				String seq = seqConfig.getSeqGen().execute(r).toUpperCase();
				for(int k=seqConfig.getKmin(); k<=seqConfig.getKmax(); k++){
					for(int i=0; i<(seq.length()-k+1); i++){
						String currk = seq.substring(i, i+k);
						String revcurrk = SequenceUtils.reverseComplement(currk);
						int currKInt = RegionFileUtilities.seq2int(currk);
						int revcurrKInt = RegionFileUtilities.seq2int(revcurrk);
						int kmer = currKInt<revcurrKInt ? currKInt : revcurrKInt;
						int baseind = seqConfig.getKmerBaseInd(currk);	
						pfl[baseind+kmer]++;
					}
				}
				ret.add(pfl);
			}
			
			return ret;
		}
		
		
		
		/**
		 * Finds mountains for a given list of regions and given scoring threshold
		 * Odds threshold can be -ve when scanning the model at negPeaks. +ve when scanning the model at posPeaks
		 * Also populates the 
		 * @param useCache
		 * @param genpath
		 * @param oddsThresh
		 */
		private List<Pair<Region,Double>> findHills(){
			
			// Detects if scanning is done for the positive or neg set from the sign of oddsThresh
			//int classDetector = oddsThresh>0 ? 1:-1;
			
			List<Pair<Region,Double>> ret = new ArrayList<Pair<Region,Double>>();
			
			for(Region r : modelRegions){
				String seq = seqConfig.getSeqGen().execute(r).toUpperCase();
				if(seq.contains("N"))
					continue;
				List<Pair<Region,Double>> mountains = new ArrayList<Pair<Region,Double>>();
				for(int l=seqConfig.getMinM(); l<=seqConfig.getMaxM(); l++){
					for(int i=0; i<(seq.length()-l+1); i++){
						String motif = seq.substring(i, i+l);
						double score=0.0;
						for(int k=seqConfig.getKmin(); k<=seqConfig.getKmax(); k++){
							for(int j=0; j<motif.length()-k+1; j++){
								String currk = motif.substring(j, j+k);
								String revcurrk = SequenceUtils.reverseComplement(currk);
								int  currKInt = RegionFileUtilities.seq2int(currk);
								int  revCurrKInt = RegionFileUtilities.seq2int(revcurrk);
								int kmer = currKInt<revCurrKInt ? currKInt : revCurrKInt;
								int baseInd = seqConfig.getKmerBaseInd(currk);
								score = score+seqConfig.getKmerWeights().get(kmerModelName)[baseInd+kmer];
							}
						}
						
						Region hill = new Region(seqConfig.getGenome(),r.getChrom(),r.getStart()+i,r.getStart()+i+l-1);
						
						Iterator<Pair<Region,Double>> it = mountains.iterator();
						boolean add=true;
						while(it.hasNext() && add){
							Pair<Region,Double> pr = it.next();
							Region currHill = pr.car();
							Double currScore = pr.cdr();
							if(currHill.overlaps(hill) && currScore<=score){
								it.remove();
								add=true;
							}else if(currHill.overlaps(hill) && currScore> score){
								add=false;
							}
						}
						if(add && score > seqConfig.getHillsThresh()){
							mountains.add(new Pair<Region,Double>(hill,score));
						}
						
					}
				}
				ret.addAll(mountains);	
			}
			return ret;
		}
		
		public class HillsIndexComparator implements Comparator<Integer>{
			
			List<Pair<Region,Double>> hill;
			
			public HillsIndexComparator(List<Pair<Region,Double>> h) {
				hill = h;
			}
			
			public Integer[] createIndexArray(){
				Integer[] indexes = new Integer[hill.size()];
				for(int i=0; i<indexes.length; i++){
					indexes[i] = i;
				}
				return indexes;
			}

			@Override
			public int compare(Integer o1, Integer o2) {
				return -1*(hill.get(o1).cdr().compareTo(hill.get(o2).cdr()));
			}
			
		}
	}

}
