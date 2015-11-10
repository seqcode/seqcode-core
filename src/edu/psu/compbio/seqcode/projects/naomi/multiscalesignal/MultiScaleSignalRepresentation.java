package edu.psu.compbio.seqcode.projects.naomi.multiscalesignal;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.StringEscapeUtils;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
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
import edu.psu.compbio.seqcode.projects.seed.stats.FeatureStatistics;

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
	protected ExperimentManager manager;
	protected int numScale;

	//external parameters
	protected int threePrimReadExt = 200;
	protected int binWidth = 1;
	
	protected double scaling = 1;	
	
	
	
	protected Map<Integer,List<Region>> segRegionTree = new HashMap<Integer, List<Region>>();
	
	public MultiScaleSignalRepresentation(GenomeConfig gcon, ExptConfig econ, SEEDConfig scon, ExperimentManager man, int scale){	
		gconfig = gcon;
		econfig = econ;
		sconfig = scon;
		manager = man;
		numScale = scale;
	}
	 
	public void runMSR(){

		Genome genome = gconfig.getGenome();
		
		//fix here to get parameters only if they are specified
		binWidth = sconfig.getBinWidth();
		threePrimReadExt = sconfig.getTag3PrimeExtension();
		
		//test to print binWidth and threePrimReadExt
		System.out.println("binWidth is: "+binWidth);
		System.out.println("threePrimReadExt is: "+threePrimReadExt);
		
		//get scaling ratio	
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
			
			System.out.println("current chrom is "+currChrom.getChrom());
			
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
						if (hits.getStrand()=='+' && (int) Math.ceil((hits.getCoordinate()+i)/binWidth)<currchromBinSize){
							sampleCounts[(int) Math.ceil((hits.getCoordinate()+i)/binWidth)]+=hits.getCount();
						}else if (hits.getStrand()=='-' && (int) Math.ceil((hits.getCoordinate()-i)/binWidth) >=0 
								&& (int) Math.ceil((hits.getCoordinate())/binWidth)<currchromBinSize){								
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
						gaussianBlur[i][1] = (float) (signalCounts[i]-scaling*controlCounts[i]);
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
			float DImax = (float) Integer.MIN_VALUE;
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
			
//			System.out.println("DImax is: "+DImax+"\t"+"DImin is: "+DImin+
//					"\t"+"trailingZero: "+trailingZero+"\t"+"zeroEnd"+"\t"+zeroEnd);	
			float maxInt = DImax-DImin;
//			System.out.println("max intensity is: "+maxInt);
			
			final long startTime = System.currentTimeMillis();

			//build segmentationTree
			SegmentationTree segtree = new SegmentationTree(gconfig, econfig, sconfig, numScale);	
			
//			Map<Integer,Set<Integer>> segmentationTree = segtree.buildTree(currchromBinSize, gaussianBlur, linkageMap, maxInt, trailingZero, zeroEnd);
			
			final long endTime = System.currentTimeMillis();
			
			System.out.println("building segTree for chrom "+ currChrom+ " took "+(endTime-startTime));
			System.out.println(StringEscapeUtils.escapeJava(currChrom.getChrom()));

			
			// I'm testing with small chrom21 ; beginning of test
			Map<Integer,Set<Integer>> segmentationTree = null;	
			
			if (currChrom.getChrom().contains("13") && currchromBinSize > 10000000 && currchromBinSize <20000000){			
			segmentationTree = segtree.buildTree(currchromBinSize, gaussianBlur, linkageMap, maxInt, trailingZero, zeroEnd);

			
			//printing segmenationTree for test
			System.out.println("from returned values from segmentationTree");
			for (Integer scale : segmentationTree.keySet()){
					System.out.println("current scale is:"+scale);
					Set<Integer> segmentation = segmentationTree.get(scale);
					System.out.println("current size is : "+segmentation.size());
//					System.out.println("current nodes are is : "+segmentation);	
			}	
			
			//converting coordinates to regions
			for (Integer scale : segmentationTree.keySet()){
				Integer prevCoord = 0;
				for (Integer coord : segmentationTree.get(scale)){
					if (coord>0){
						Region segRegion = new Region(genome,currChrom.getChrom(),prevCoord,coord);
						if (segRegionTree.containsKey(scale)){
							segRegionTree.get(scale).add(segRegion);	
						}else{
							List<Region> regionList = new ArrayList<Region>();
							regionList.add(segRegion);
							segRegionTree.put(scale, regionList);
						}					
						prevCoord = coord;
					}
				}
			}
			
			}//end of test with small chromosome
			
		}// end of chromosome iteration					
		manager.close();
	}
	
	// for now I am performing binomial test; later change to edgeR
	public void computeSFC(){
		
		for (Integer scale : segRegionTree.keySet()){
			System.out.println("scale is: "+scale+"size is "+segRegionTree.get(scale).size());
		}
		
		Map<Integer,List<Region>> segSFC = new HashMap<Integer, List<Region>>();
		
		Sample signal = null;
		Sample control = null;
		for (ControlledExperiment rep: manager.getReplicates()){
			signal = rep.getSignal();
			control = rep.getControl();
		}
		
		System.out.println("signal is :"+signal.getName());
		System.out.println("controi is: "+control.getName());
		
		double signalCounts = 0;
		double controlCounts = 0;
		double pval = 1;
		
		for (Integer scale : segRegionTree.keySet()){
//		List<Integer> scaleList = new ArrayList<Integer>();
//		scaleList.add(4);
//		scaleList.add(9);
//		scaleList.add(14);
//		scaleList.add(19);
		
//		for (Integer scale : scaleList){
			List<Region> rSFC = new ArrayList<Region>();
			List<Region> regList = segRegionTree.get(scale);
			System.out.println("size of region "+regList.size());
			for (Region reg : regList){

				signalCounts = signal.countHits(reg);
				controlCounts = control.countHits(reg)*scaling;
				
				FeatureStatistics stat = new FeatureStatistics();
				if (signalCounts+controlCounts>0){
					
					if (signalCounts>controlCounts){
						pval = stat.binomialPValue(controlCounts, signalCounts+controlCounts, sconfig.getMinSigCtrlFoldDifference());
					}else{
						pval = stat.binomialPValue(signalCounts, signalCounts+controlCounts, sconfig.getMinSigCtrlFoldDifference());
					}	
				}
				if (pval<0.01){rSFC.add(reg);}
				pval =1;
				
			}
			segSFC.put(scale,rSFC);
		}
		
		for (Integer scale : segRegionTree.keySet()){
			System.out.println("scale: "+scale+"size of original region "+segRegionTree.get(scale).size());
			System.out.println("size  of SFC region "+segSFC.get(scale).size());
			System.out.println("content "+segSFC.get(scale));
		}		
		manager.close();		
	}
		
	public void printCounts() throws FileNotFoundException, UnsupportedEncodingException{
		
		Sample signal = null;
		Sample control = null;
		for (ControlledExperiment rep: manager.getReplicates()){
			signal = rep.getSignal();
			control = rep.getControl();
		}		
		
//		for (Integer scale : segRegionTree.keySet()){
		
		List<Integer> scaleKeyList = new ArrayList<Integer>(segRegionTree.keySet());
		
		// selecting a set of scales for output
		List<Integer> scaleList = new ArrayList<Integer>();
		for (int i = Math.round(scaleKeyList.size()/5); i<scaleKeyList.size(); i = i+Math.round(scaleKeyList.size()/5))
			scaleList.add(i);		
		
		for (Integer scale : scaleList){
			
			String outName = scale+"_scaleCounts.txt";
			
			PrintWriter writer = new PrintWriter(outName,"UTF-8");

			List<Region> regList = segRegionTree.get(scale);
				
			writer.println("#scale is: "+scale+"size is "+segRegionTree.get(scale).size());
			
			System.out.println("total control signal is "+control.getHitCount()+" signal is "+signal.getHitCount());
			
			for (Region reg : regList){
				System.out.println(reg.getChrom()+"\t"+reg+"\t"+(control.countHits(reg)+signal.countHits(reg)));
//				if (signal.countHits(reg)+control.countHits(reg) > 0)
//					writer.println(reg+"\t"+signal.countHits(reg)+"\t"+control.countHits(reg)*scaling);
			}
			writer.close();
		}
		manager.close();				
	}
	
		
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
		
		/***
		 * Need to specify --tag3ext & --binwidth --scale
		 ***/
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
		SEEDConfig sconf = new SEEDConfig(gconf, args);
		ArgParser ap = new ArgParser(args);
		ExperimentManager man = new ExperimentManager(econf); 
		int numScale = 20;
		if (ap.hasKey("scale")){
			numScale = Args.parseInteger(args,"scale",3);
		}		
		MultiScaleSignalRepresentation msr = new MultiScaleSignalRepresentation (gconf, econf, sconf, man, numScale);	
		msr.runMSR();	
		msr.computeSFC();
//		msr.printCounts();
		
		
	}	
}
