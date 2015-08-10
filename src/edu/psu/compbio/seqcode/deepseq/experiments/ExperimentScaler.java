package edu.psu.compbio.seqcode.deepseq.experiments;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import Jama.Matrix;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.utils.models.Model;
import edu.psu.compbio.seqcode.gse.utils.models.data.DataFrame;
import edu.psu.compbio.seqcode.gse.utils.models.data.DataRegression;
import edu.psu.compbio.seqcode.gse.viz.scatter.ScatterPlot;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingManager;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingModel;
import edu.psu.compbio.seqcode.projects.multigps.framework.MultiGPSConfig;
import edu.psu.compbio.seqcode.projects.multigps.framework.PotentialRegionFilter;

/**
 * ExperimentScaler: calculate a scaling transformation between all Sample pairs in an ExperimentCondition
 * This is performed on the condition level so that a scaling can also be defined between pooled hits from all signals & controls
 *  
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExperimentScaler {
	
	public ExperimentScaler(){}
	
	/**
	 * Calculate a scaling ratio by fitting a line through the hit count pairs.
	 * Using a 10Kbp window, this is the same as PeakSeq with Pf=0
	 * @return double 
	 */
	public double scalingRatioByRegression(List<Float> setA, List<Float> setB){
		double scalingRatio=1;
		if(setA.size()!=setB.size()){
			System.err.println("ExperimentScaler is trying to scale lists of two different lengths");
			System.exit(1);
		}
			
		List<PairedCounts> scalingData = new ArrayList<PairedCounts>();
		for(int x=0; x<setA.size(); x++)
			scalingData.add(new PairedCounts(setA.get(x), setB.get(x)));                

		//Scaling ratio via Tim's regression                                                                                                                               
        DataFrame df = new DataFrame(PairedCounts.class, scalingData.iterator());                                                                      
        DataRegression r = new DataRegression(df, "x~y - 1");                                                                                                               
        r.calculate();                                                                                                                                                      
        Map<String, Double> map = r.collectCoefficients();                                                                                                                  
        scalingRatio = map.get("y");                                                                                                                                        
        return(scalingRatio);
	}
	
	/**
	 * Find the median hit count ratio in bins that have non-zero counts
	 * @return
	 */
	public double scalingRatioByMedian(List<Float> setA, List<Float> setB){
		double scalingRatio=1;
		if(setA.size()!=setB.size()){
			System.err.println("ExperimentScaler is trying to scale lists of two different lengths");
			System.exit(1);
		}
			
		ArrayList<Float> ratios = new ArrayList<Float>();
	    for(int x=0; x<setA.size(); x++){
			if(setA.get(x)>0 && setB.get(x)>0)
				ratios.add((float)(setA.get(x) / setB.get(x)));
        }
        Collections.sort(ratios);
		scalingRatio = ratios.get(ratios.size() / 2);
        return(scalingRatio);
	}
	
	/**
	 * Find the scaling ratio according to the SES method from Diaz, et al. Stat Appl Genet Mol Biol. 2012.
	 * Also sets a background proportion estimate for the signal channel.  
	 * @return
	 */
	public double scalingRatioBySES(List<Float> setA, List<Float> setB){
		double scalingRatio=1;
		if(setA.size()!=setB.size()){
			System.err.println("ExperimentScaler is trying to scale lists of two different lengths");
			System.exit(1);
		}
		
		float totalA=0, totalB=0;
		List<PairedCounts> counts = new ArrayList<PairedCounts>();
		for(int x=0; x<setA.size(); x++){
			totalA += setA.get(x);
			totalB += setB.get(x);
			counts.add(new PairedCounts(setA.get(x), setB.get(x)));                
		}
		
		Collections.sort(counts);
        
        //SES procedure
        double cumulA=0, cumulB=0, maxDiffAB=0, maxDiffAprop=0, currDiff=0;
        int maxDiffIndex=0, i=0;
        for(PairedCounts pc : counts){
        	cumulA+=pc.x;
        	cumulB+=pc.y;
        	currDiff = (cumulB/totalB)-(cumulA/totalA);
        	if(currDiff>maxDiffAB && cumulA>0 && cumulB>0){
        		maxDiffAB=currDiff;
        		maxDiffIndex=i;
        		maxDiffAprop=(cumulA/totalA);
        		scalingRatio = cumulA/cumulB;
        	}
        	i++;
        }
		return(scalingRatio);
	}
	
	/**
	 * Find the scaling ratio according to the NCIS method from Liang & Keles (BMC Bioinf 2012).
	 * Also sets a background proportion estimate for the signal channel.  
	 * Should be run using *all* genomic windows in the Lists. 
	 * Uses ratios that are based on at least 75% of genomic regions by default. 
	 * @param setA : signal list
	 * @param setB : control list
	 * @param outputFile : optional file that will contain the data 
	 * @return
	 */
	public double scalingRatioByNCIS(List<Float> setA, List<Float> setB, String outputFile){
		double scalingRatio=1;
		if(setA.size()!=setB.size()){
			System.err.println("ExperimentScaler is trying to scale lists of two different lengths");
			System.exit(1);
		}
		
		float numPairs = (float)setA.size();
		List<PairedCounts> counts = new ArrayList<PairedCounts>();
		for(int x=0; x<setA.size(); x++)
			counts.add(new PairedCounts(setA.get(x), setB.get(x))); 
		
		//NCIS uses increasing total tag counts versus enrichment ratio
		Collections.sort(counts, new Comparator<PairedCounts>(){
            public int compare(PairedCounts o1, PairedCounts o2) {return o1.compareByTotal(o2);}
        });
        
        //NCIS procedure
        double cumulA=0, cumulB=0, currRatio=0, lastRatio=-1;
        float i=0;
        List<Double> totalCounts=new ArrayList<Double>();
        List<Double> ratios=new ArrayList<Double>();
        for(PairedCounts pc : counts){
        	cumulA+=pc.x;
        	cumulB+=pc.y;
        	
        	i++;
        	if(i/numPairs > 0.75 && cumulA>0 && cumulB>0){ //NCIS estimates begin using the lower 3 quartiles of the genome (based on total tags)
	        	currRatio = (cumulA/cumulB);
	        	if(lastRatio==-1 || currRatio<lastRatio){
	        		lastRatio = currRatio;
	        	}else{
	        		break;
	        	}
        	}
        	if(pc.x>0 && pc.y>0){
        		Double ratio  = (pc.x/pc.y); 
        		totalCounts.add(pc.x+pc.y);
        		ratios.add(ratio);
        	}
        }
        scalingRatio = currRatio;
        
        /*Scatter plot generation*/
        if(outputFile!=null){
	        Matrix dataToPlot = new Matrix(totalCounts.size(),2);
	        int count=0;
			for(int d=0; d<totalCounts.size(); d++){
				dataToPlot.set(count, 0, totalCounts.get(d));
				dataToPlot.set(count, 1, ratios.get(d));
				count++;
			}
			//Generate image
			ScalingPlotter plotter = new ScalingPlotter(outputFile+" NCIS plot");
			plotter.saveXYplot(dataToPlot, "totalTagCount", "Signal/Control", outputFile+".NCIS_scaling.png", true);
        }
        
		return(scalingRatio);
	}
	
	/**
	 * Calculate the background proportion of an IP experiment by correcting the scaling ratio by the read count ratio.
	 * Be careful with this method, there are a couple of assumptions:
	 *  - The scaling ratio was calculated between IP and control experiments
	 *  - The method used to calculate the scaling ratio attempted to normalize to background and not all regions (e.g. SES method attempts background normalization) 
	 * @return
	 */
	public Double calculateBackgroundFromScalingRatio(ControlledExperiment expt){
		if(expt.getControlScaling()==-1)
			return(-1.0); //scaling not yet performed
		double ctrlCount = expt.getControl()==null ? expt.getSignal().getHitCount() : expt.getControl().getHitCount();
		return(expt.getControlScaling() / (expt.getSignal().getHitCount()/ctrlCount));
	}
	
	
	
	/**
	 * Main for testing
	 * @param args
	 */
	public static void main(String[] args){
		
		GenomeConfig gconfig = new GenomeConfig(args);
		ExptConfig econfig = new ExptConfig(gconfig.getGenome(), args);
		MultiGPSConfig mgpsconfig = new MultiGPSConfig(gconfig, args, false);
		
		
		if(gconfig.helpWanted()){
			System.err.println("ExperimentScaler:");
			System.err.println(gconfig.getArgsList()+"\n"+econfig.getArgsList());
		}else{
			ExperimentManager exptMan = new ExperimentManager(econfig);
			BindingManager bindingManager = new BindingManager(exptMan);
			//Initialize binding models & binding model record
			Map<ControlledExperiment, List<BindingModel>>  repBindingModels = new HashMap<ControlledExperiment, List<BindingModel>>();
			for(ControlledExperiment rep : exptMan.getReplicates()){
				if(mgpsconfig.getDefaultBindingModel()!=null)
					bindingManager.setBindingModel(rep, mgpsconfig.getDefaultBindingModel());
				else if(rep.getExptType()!=null && rep.getExptType().getName().toLowerCase().equals("chipexo"))
					bindingManager.setBindingModel(rep, new BindingModel(BindingModel.defaultChipExoEmpiricalDistribution));
				else
					bindingManager.setBindingModel(rep, new BindingModel(BindingModel.defaultChipSeqEmpiricalDistribution));
				repBindingModels.put(rep, new ArrayList<BindingModel>());
				repBindingModels.get(rep).add(bindingManager.getBindingModel(rep));
			}
			for(ExperimentCondition cond : exptMan.getConditions())
				bindingManager.updateMaxInfluenceRange(cond);
			
			
			//Test
			System.err.println("Conditions:\t"+exptMan.getConditions().size());
			for(ExperimentCondition c : exptMan.getConditions()){
				System.err.println("Condition "+c.getName()+":\t#Replicates:\t"+c.getReplicates().size());
			}
			for(ExperimentCondition c : exptMan.getConditions()){
				for(ControlledExperiment r : c.getReplicates()){
					System.err.println("Condition "+c.getName()+":\tRep "+r.getName());
					if(r.getControl()==null)
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount());
					else
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount()+"\tControl:\t"+r.getControl().getHitCount());
				}
			}
			
			ExperimentScaler scaler = new ExperimentScaler();
			
			//Potential regions reqd by PeakSeq method
			PotentialRegionFilter potentialFilter = new PotentialRegionFilter(mgpsconfig, econfig, exptMan, bindingManager);
			List<Region> potentials = potentialFilter.execute();
			System.out.println("\t"+potentials.size()+" potential regions.\n");
			
			//Generate the data structures for calculating scaling factors
			//Window size loaded by ExptConfig option --scalewin
			Genome genome = econfig.getGenome();
			Map<Sample, List<Float>> sampleWindowCounts = new HashMap<Sample, List<Float>>();
			Map<Sample, List<Float>> noPotSampleWindowCounts = new HashMap<Sample, List<Float>>();
			for(Sample samp : exptMan.getSamples()){
				List<Float> currSampCounts = new ArrayList<Float>();
				List<Float> noPotCurrSampCounts = new ArrayList<Float>();
				for(String chrom:genome.getChromList()) {
		            int chrlen = genome.getChromLength(chrom);
		            for (int start = 1; start  < chrlen - econfig.getScalingSlidingWindow(); start += econfig.getScalingSlidingWindow()) {
		                Region r = new Region(genome, chrom, start, start + econfig.getScalingSlidingWindow());
		                currSampCounts.add(samp.countHits(r));
		                
		                boolean overlapsPotentials=false;
		                for(Region p : potentials){
		                	if(r.overlaps(p)){
		                		overlapsPotentials=true; break;
		                	}
		                }
		                if(!overlapsPotentials)
		                	noPotCurrSampCounts.add(samp.countHits(r));
		            }
		        }
				sampleWindowCounts.put(samp, currSampCounts);
				noPotSampleWindowCounts.put(samp, noPotCurrSampCounts);
			}
			System.out.println("Sliding window size for scaling methods: "+econfig.getScalingSlidingWindow());
			System.out.println("\tNumbers of windows:\tAll="+sampleWindowCounts.get(exptMan.getSamples().get(0)).size()+"\tnoPotenials="+noPotSampleWindowCounts.get(exptMan.getSamples().get(0)).size()+"\n");
			
			
			//Hit ratios
			for(Sample sampA : exptMan.getSamples()){ 
				if(sampA.isSignal()){
					for(Sample sampB : exptMan.getSamples())
						if(sampA!=null && sampB!=null && sampA.getIndex() != sampB.getIndex()){
							double hitRatio = sampA.getHitCount()/sampB.getHitCount();
							System.out.println("HitRatio\t"+sampA.getName()+" vs "+sampB.getName()+"\t"+hitRatio);
						}
				}
			}
			
			//Median
			for(Sample sampA : exptMan.getSamples()){ 
				if(sampA.isSignal()){
					for(Sample sampB : exptMan.getSamples())
						if(sampA!=null && sampB!=null && sampA.getIndex() != sampB.getIndex())
							System.out.println("Median\t"+sampA.getName()+" vs "+sampB.getName()+"\t"+scaler.scalingRatioByMedian(sampleWindowCounts.get(sampA), sampleWindowCounts.get(sampB)));
				}
			}
			
			//Regression on full dataset (i.e. PeakSeq using Pf=0)
			for(Sample sampA : exptMan.getSamples()){ 
				if(sampA.isSignal()){
					for(Sample sampB : exptMan.getSamples())
						if(sampA!=null && sampB!=null && sampA.getIndex() != sampB.getIndex())
							System.out.println("Regression\t"+sampA.getName()+" vs "+sampB.getName()+"\t"+scaler.scalingRatioByRegression(sampleWindowCounts.get(sampA), sampleWindowCounts.get(sampB)));
				}
			}
			
			//SES
			for(Sample sampA : exptMan.getSamples()){ 
				if(sampA.isSignal()){
					for(Sample sampB : exptMan.getSamples())
						if(sampA!=null && sampB!=null && sampA.getIndex() != sampB.getIndex())
							System.out.println("SES\t"+sampA.getName()+" vs "+sampB.getName()+"\t"+scaler.scalingRatioBySES(sampleWindowCounts.get(sampA), sampleWindowCounts.get(sampB)));
				}
			}
			
			//NCIS
			for(Sample sampA : exptMan.getSamples()){ 
				if(sampA.isSignal()){
					for(Sample sampB : exptMan.getSamples())
						if(sampA!=null && sampB!=null && sampA.getIndex() != sampB.getIndex())
							System.out.println("NCIS\t"+sampA.getName()+" vs "+sampB.getName()+"\t"+scaler.scalingRatioByNCIS(sampleWindowCounts.get(sampA), sampleWindowCounts.get(sampB), null));
				}
			}
			
			//Regression after filtering out potential regions (i.e. PeakSeq using Pf=1)
			for(Sample sampA : exptMan.getSamples()){ 
				if(sampA.isSignal()){
					for(Sample sampB : exptMan.getSamples())
						if(sampA!=null && sampB!=null && sampA.getIndex() != sampB.getIndex())
							System.out.println("PeakSeq\t"+sampA.getName()+" vs "+sampB.getName()+"\t"+scaler.scalingRatioByRegression(noPotSampleWindowCounts.get(sampA), noPotSampleWindowCounts.get(sampB)));
				}
			}
			
			exptMan.close();
		}
	}
	
	
	/**
	 * Simple class for storing paired counts that are sortable in first dimension
	 * @author mahony
	 *
	 */
	public class PairedCounts extends Model implements Comparable<PairedCounts>{
		public Double x,y;
		public PairedCounts(double a, double b){
			x=a;
			y=b;
		}
		/**
		 * Sort on increasing X variables
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		public int compareTo(PairedCounts pc) {
			if(x<pc.x){return -1;}
			if(x>pc.x){return 1;}
			return 0;
		}
		
		/**
		 * Compare based on the sum of both paired counts
		 * @param pc
		 * @return
		 */
		public int compareByTotal(PairedCounts pc){
			if((x+y)<(pc.x+pc.y)){return -1;}
			if((x+y)>(pc.x+pc.y)){return 1;}
			return 0;
		}
		
	}
	
	public class ScalingPlotter extends ScatterPlot{
		public ScalingPlotter(String title) {
			super(title);
		}
		/**
		 * Make an XY scatter plot from a 2-D dataset and save image.
		 * This configuration is used by deepseq.stats.Normalization classes
		 *  
		 * @param datapoints - 2D dataset (colored grey)
		 * @param datapoints_highlight - 2D dataset (colored blue), can be null
		 * @param yLine Double - data coordinates of line to be drawn parallel to x axis (used to show scaling line) 
		 * @param outFilename - String
		 * @param rasterImage - boolean
		 */
		private void saveXYplot(Matrix datapoints, String xName, String yName, String outFilename, boolean rasterImage){
			this.setWidth(800);
			this.setHeight(800);
			this.setXLogScale(false);
			this.setYLogScale(false);
			this.addDataset("other", datapoints, new Color(75,75,75,80), 3);
			this.setXAxisLabel(xName);
			this.setYAxisLabel(yName);
			this.setXRangeFromData();
			this.setYRangeFromData();
			
			//Set the tick units according to the range
			double xUpper = daxis.getRange().getUpperBound();
			double xLower = daxis.getRange().getLowerBound();
	    	double yUpper = raxis.getRange().getUpperBound();
			double yLower = raxis.getRange().getLowerBound();
			
			try {
				this.saveImage(new File(outFilename), width, height, rasterImage);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
}
