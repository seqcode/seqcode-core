package edu.psu.compbio.seqcode.projects.naomi.multiscalesignal;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromosomeGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.projects.seed.SEEDConfig;

/**
 * MultiScaleSignalRepresentation
 *
 * Methods refer to two papers
 * Probabilistic Multiscale Image Segmentation, Vincken et al. IEEE (1997)
 * 
 * @author naomi yamada
 *
 **/

public class MultiScaleSignalRepresentation {

	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected SEEDConfig sconfig;
	protected int numScale;

	//external parameters
	protected int threePrimReadExt = 200;
	protected int binWidth = 1;
	
	protected Map<Region, HashMap<Integer,Set<Integer>>> segmentationTree = new HashMap<Region, HashMap<Integer, Set<Integer>>>();

	public MultiScaleSignalRepresentation(GenomeConfig gcon, ExptConfig econ, SEEDConfig scon, int scale){	
		gconfig = gcon;
		econfig = econ;
		sconfig = scon;
		numScale = scale;
	}
	 
	public void runMSR(){
		
		ExperimentManager manager = new ExperimentManager(econfig);
		Genome genome = gconfig.getGenome();
		
		//fix here to get parameters only if they are specified
		binWidth = sconfig.getBinWidth();
		threePrimReadExt = sconfig.getTag3PrimeExtension();
		
		//test to print binWidth and threePrimReadExt
		System.out.println("binWidth is: "+binWidth);
		System.out.println("threePrimReadExt is: "+threePrimReadExt);
		
		//get scaling ratio
		double scaling = 1;		
		for (ControlledExperiment rep: manager.getReplicates()){
			scaling = rep.getControlScaling();
			System.out.println("Condition: "+rep.getCondName()+"\tScalingFactor: "+scaling);
		}
		
		Iterator<Region> chroms = new ChromosomeGenerator<Genome>().execute(genome);
		//iterating each chromosome (each chromosome is a region).
		while (chroms.hasNext()) {
			
			Region currChrom = chroms.next();	
			int currchromSize = currChrom.getWidth();
			int currchromBinSize = (int) Math.ceil(currchromSize/binWidth);
			
			Map<Sample,float[]> condiGaussianBlur = new HashMap<Sample,float[]>();
			float[] sampleCounts = new float[currchromBinSize];
			//primitive array to store signal and the subsequent convolved signals
			//its index correspond to the coordinates
			float[][] gaussianBlur = new float[currchromBinSize][2];
			for (int i = 0; i<currchromBinSize; i++){
				for (int j = 0; j<2; j++)
					gaussianBlur[i][j] = 0;
			}
			
			//get StrandedBaseCount list for signal and control
			Map<Sample, List<StrandedBaseCount>> condiCountsMap = new HashMap<Sample, List<StrandedBaseCount>>();
			
			for (ControlledExperiment rep: manager.getReplicates()){
				condiCountsMap.put(rep.getSignal(), rep.getSignal().getBases(currChrom));
				if (rep.hasControl())
					condiCountsMap.put(rep.getControl(), rep.getControl().getBases(currChrom));			
			}
			
			//StrandedBasedCount object contains positive and negative strand separately
			//store all base counts indexed by positions at column[1]
			//extend reads to 3' end and bin according to bin size			
			for (Sample sample : manager.getSamples()){	
				for (int i = 0;i<currchromBinSize;i++)
						sampleCounts[i] = 0;
				
				List<StrandedBaseCount> currentCounts = condiCountsMap.get(sample);
				for (StrandedBaseCount hits: currentCounts){
					for (int i = 0; i<threePrimReadExt+1; i++){
						if (hits.getStrand()=='+' && hits.getCoordinate()+i<currchromSize){
							sampleCounts[(int) Math.ceil((hits.getCoordinate()+i)/binWidth)]+=hits.getCount();
						}else if (hits.getStrand()=='+' && hits.getCoordinate()-i >=0){
							sampleCounts[(int) Math.ceil((hits.getCoordinate()-i)/binWidth)]+=hits.getCount();
						}
					}
				}
				condiGaussianBlur.put(sample, sampleCounts);
				currentCounts = null;
			}
			
			//if there is a single sample, build SegmentaionTree using single sample
			for (ControlledExperiment rep : manager.getReplicates()){
				if (!rep.hasControl()){
					for (Sample sample : manager.getSamples()){
						float[] counts = condiGaussianBlur.get(sample);
						for (int i = 0; i<currchromBinSize; i++)
							gaussianBlur[i][1] = counts[i];					
					}
				}else{ //for control and signal; construct a gaussianBlur by taking signal*scaling-control
					float[] signalCounts = condiGaussianBlur.get(rep.getSignal());
					float[] controlCounts = condiGaussianBlur.get(rep.getControl());
					for (int i = 0; i<currchromBinSize; i++)
						gaussianBlur[i][1] = (float) (signalCounts[i]*scaling)-controlCounts[i];
				}
			}
			
			
			/*********************
			 * Starting nodes
			 */					
			//linkageMap contains index of kids and parents
			Map <Integer, Integer> linkageMap = new HashMap<Integer, Integer>();
			//adding starting nodes; to qualify for the starting nodes the signal intensity needs to be different from the subsequent signal intensity
			//adding the starting and end positions in the kids at start and end positions  
			//setting max & min signal intensity  
			List <Integer> nonzeroList = new ArrayList<Integer>();
			linkageMap.put(0,0);
			float DImax = 0;
			float DImin = (float) Integer.MAX_VALUE;
			for (int i = 0 ; i< gaussianBlur.length-1; i++){ 
				if (gaussianBlur[i][1] != gaussianBlur[i+1][1])
					linkageMap.put(i,i);
				if (gaussianBlur[i][1] > DImax)
					DImax = gaussianBlur[i][1];
				if (gaussianBlur[i][1] < DImin)
					DImin = gaussianBlur[i][1];		
				if (gaussianBlur[i][1]!=0)
					nonzeroList.add(i);
			}
			linkageMap.put(gaussianBlur.length-1,gaussianBlur.length-1);
			
			//determine the first nonzero and last nonzero from signal	
			int trailingZero = 0;
			int zeroEnd = 0;
			if (!nonzeroList.isEmpty()){
				trailingZero = Collections.min(nonzeroList)-1;
				zeroEnd = Collections.max(nonzeroList)+1;
			}
			if (trailingZero == -1)
				trailingZero = 0;
			
			System.out.println("DImax is: "+DImax+"\t"+"DImin is: "+DImin+
					"\t"+"trailingZero: "+trailingZero+"\t"+"zeroEnd"+"\t"+zeroEnd);	
			
			//build segmentationTree
			SegmentationTree segtree = new SegmentationTree(gconfig, econfig, sconfig, numScale);	

			Map<Region, HashMap<Integer,Set<Integer>>> segmentationTree = segtree.buildTree(currChrom, currchromBinSize, gaussianBlur, linkageMap, DImax, DImin, trailingZero, zeroEnd);
			
			System.out.println("from returned values from segmentationTree");
			
			for (Region chrom : segmentationTree.keySet()){
				System.out.println("current chrom is: "+chrom);
				HashMap<Integer,Set<Integer>> chromTree = segmentationTree.get(chrom);
				for (Integer scale : chromTree.keySet()){
					System.out.println("current scale is:"+scale);
					Set<Integer> segmentation = chromTree.get(scale);
					System.out.println("current size is : "+segmentation.size());
					for (Integer coord : segmentation){
						System.out.println(coord);
					}
				}
			}			
		}// end of chromosome iteration		
		manager.close();
	}
		
	public static void main(String[] args) {
		
		/***
		 * Need to specify --tag3ext & --binwidth --scale
		 ***/
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
		SEEDConfig sconf = new SEEDConfig(gconf, args);
		ArgParser ap = new ArgParser(args);
		int numScale = 20;
		if (ap.hasKey("scale")){
			numScale = Args.parseInteger(args,"scale",3);
		}		
		MultiScaleSignalRepresentation msr = new MultiScaleSignalRepresentation (gconf, econf, sconf,numScale);	
		msr.runMSR();		
	}	
}
