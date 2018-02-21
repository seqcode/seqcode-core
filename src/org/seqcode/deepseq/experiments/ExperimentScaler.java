package org.seqcode.deepseq.experiments;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.NumberAxis;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.gseutils.models.Model;
import org.seqcode.ml.regression.DataFrame;
import org.seqcode.ml.regression.DataRegression;
import org.seqcode.projects.seed.DomainFinder;
import org.seqcode.projects.seed.SEEDConfig;
import org.seqcode.projects.seed.features.Feature;
import org.seqcode.viz.scatter.ScatterPlot;

import Jama.Matrix;


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
	public double scalingRatioByNCIS(List<Float> setA, List<Float> setB, String outputFile, double minFrac){
		double scalingRatio=1;
		double totalAtScaling=0;
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
        for(PairedCounts pc : counts){
        	cumulA+=pc.x;
        	cumulB+=pc.y;
        	totalAtScaling = pc.x+pc.y;
        	
        	i++;
        	if(i/numPairs > minFrac && cumulA>0 && cumulB>0){ //NCIS estimates begin using the lower 3 quartiles of the genome (based on total tags)
	        	currRatio = (cumulA/cumulB);
	        	if(lastRatio==-1 || currRatio<lastRatio){
	        		lastRatio = currRatio;
	        	}else{
	        		break;
	        	}
        	}
        }
        scalingRatio = currRatio;
        
        
        /*Scaling plot generation*/
        if(outputFile!=null){
        	//Cumulative ratio vs bin total
        	List<Double> bintotals=new ArrayList<Double>();
            List<Double> ratios=new ArrayList<Double>();
            cumulA=0; cumulB=0;
        	for(PairedCounts pc : counts){
            	cumulA+=pc.x;
            	cumulB+=pc.y;
            	if(cumulA>0 && cumulB>0){
            		Double ratio  = (cumulA / cumulB); 
            		bintotals.add(pc.x+pc.y);
            		ratios.add(ratio);
            	}
        	}
	        Matrix dataToPlot = new Matrix(bintotals.size(),2);
	        int count=0;
			for(int d=0; d<bintotals.size(); d++){
				dataToPlot.set(count, 0, bintotals.get(d));
				dataToPlot.set(count, 1, ratios.get(d));
				count++;
			}
			//Marginal ratios vs bin totals
			List<Double> bintot=new ArrayList<Double>();
			List<Double> mratios=new ArrayList<Double>();
			for(int x=0; x<counts.size(); x++){
				PairedCounts pc = counts.get(x);
				if(pc.x>0 && pc.y>0){
					double currA=pc.x, currB=pc.y;
					double currTot=pc.x+pc.y;
					while(x<counts.size()-1 && (counts.get(x+1).x + counts.get(x+1).y)==currTot){
						x++;
						pc = counts.get(x);
						currA+=pc.x; 
						currB+=pc.y;
					}
					bintot.add(currTot);
					mratios.add(currA/currB);
				}
			}
			Matrix dataToPlot2 = new Matrix(bintot.size(),2);
	        count=0;
			for(int d=0; d<bintot.size(); d++){
				dataToPlot2.set(count, 0, bintot.get(d));
				dataToPlot2.set(count, 1, mratios.get(d));
				count++;
				
			}

			//Generate images
			ScalingPlotter plotter = new ScalingPlotter(outputFile+" NCIS plot");
			plotter.saveXYplot(dataToPlot, totalAtScaling, scalingRatio, "Binned Total Tag Count", "Cumulative Count Scaling Ratio", outputFile+".NCIS_scaling-ccr.png", true);
			ScalingPlotter plotter2 = new ScalingPlotter(outputFile+" NCIS plot");
			plotter2.saveXYplot(dataToPlot2, totalAtScaling, scalingRatio, "Binned Total Tag Count", "Marginal Signal/Control Ratio", outputFile+".NCIS_scaling-marginal.png", true);
			
			//Print data points to files
			try {
				FileWriter fout = new FileWriter(outputFile+".NCIS_scaling-ccr.count");
				for(int d=0; d<bintotals.size(); d++)
					fout.write(bintotals.get(d)+"\t"+ratios.get(d)+"\n");
				fout.close();
				FileWriter fout2 = new FileWriter(outputFile+".NCIS_scaling-marginal.count");
				for(int d=0; d<bintot.size(); d++)
					fout2.write(bintot.get(d)+"\t"+mratios.get(d)+"\n");
				fout2.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        }
        
		return(scalingRatio);
	}
	
	/**
	 * Find the scaling ratio according to the total tag normalization followed by NCIS method from Liang & Keles (BMC Bioinf 2012).
	 * Also sets a background proportion estimate for the signal channel.  
	 * Should be run using *all* genomic windows in the Lists. 
	 * Uses ratios that are based on at least 75% of genomic regions by default. 
	 * @param setA : signal list
	 * @param setB : control list
	 * @param outputFile : optional file that will contain the data 
	 * @return
	 */
	public double scalingRatioByHitRatioAndNCIS(List<Float> setA, List<Float> setB, double totalA, double totalB, String outputFile, double minFrac){
		double scalingRatio=1;
		double totalAtScaling=0;
		if(setA.size()!=setB.size()){
			System.err.println("ExperimentScaler is trying to scale lists of two different lengths");
			System.exit(1);
		}
		
		//First normalize using total reads
		float tRatio = (float) (totalA/totalB);
		List<Float> setnB = new ArrayList<Float>();
		for (int x=0; x< setnB.size();x++)
			setnB.add(setnB.get(x)/tRatio);
		
		float numPairs = (float)setA.size();
		List<PairedCounts> counts = new ArrayList<PairedCounts>();
		for(int x=0; x<setA.size(); x++)
			counts.add(new PairedCounts(setA.get(x), setnB.get(x))); 
		
		//NCIS uses increasing total tag counts versus enrichment ratio
		Collections.sort(counts, new Comparator<PairedCounts>(){
            public int compare(PairedCounts o1, PairedCounts o2) {return o1.compareByTotal(o2);}
        });
        
        //NCIS procedure
        double cumulA=0, cumulB=0, currRatio=0, lastRatio=-1;
        float i=0;
        for(PairedCounts pc : counts){
        	cumulA+=pc.x;
        	cumulB+=pc.y;
        	totalAtScaling = pc.x+pc.y;
        	
        	i++;
        	if(i/numPairs > minFrac && cumulA>0 && cumulB>0){ //NCIS estimates begin using the lower 3 quartiles of the genome (based on total tags)
	        	currRatio = (cumulA/cumulB);
	        	if(lastRatio==-1 || currRatio<lastRatio){
	        		lastRatio = currRatio;
	        	}else{
	        		break;
	        	}
        	}
        }
        scalingRatio = currRatio*tRatio; //Multiply by the total tag normalization
        
        
        /*Scaling plot generation*/
        if(outputFile!=null){
        	//Cumulative ratio vs bin total
        	List<Double> bintotals=new ArrayList<Double>();
            List<Double> ratios=new ArrayList<Double>();
            cumulA=0; cumulB=0;
        	for(PairedCounts pc : counts){
            	cumulA+=pc.x;
            	cumulB+=pc.y;
            	if(cumulA>0 && cumulB>0){
            		Double ratio  = (cumulA / cumulB); 
            		bintotals.add(pc.x+pc.y);
            		ratios.add(ratio);
            	}
        	}
	        Matrix dataToPlot = new Matrix(bintotals.size(),2);
	        int count=0;
			for(int d=0; d<bintotals.size(); d++){
				dataToPlot.set(count, 0, bintotals.get(d));
				dataToPlot.set(count, 1, ratios.get(d));
				count++;
			}
			//Marginal ratios vs bin totals
			List<Double> bintot=new ArrayList<Double>();
			List<Double> mratios=new ArrayList<Double>();
			for(int x=0; x<counts.size(); x++){
				PairedCounts pc = counts.get(x);
				if(pc.x>0 && pc.y>0){
					double currA=pc.x, currB=pc.y;
					double currTot=pc.x+pc.y;
					while(x<counts.size()-1 && (counts.get(x+1).x + counts.get(x+1).y)==currTot){
						x++;
						pc = counts.get(x);
						currA+=pc.x; 
						currB+=pc.y;
					}
					bintot.add(currTot);
					mratios.add(currA/currB);
				}
			}
			Matrix dataToPlot2 = new Matrix(bintot.size(),2);
	        count=0;
			for(int d=0; d<bintot.size(); d++){
				dataToPlot2.set(count, 0, bintot.get(d));
				dataToPlot2.set(count, 1, mratios.get(d));
				count++;
				
			}

			//Generate images
			ScalingPlotter plotter = new ScalingPlotter(outputFile+" NCIS plot");
			plotter.saveXYplot(dataToPlot, totalAtScaling, scalingRatio, "Binned Total Tag Count", "Cumulative Count Scaling Ratio", outputFile+".TotalReadsAndNCIS_scaling-ccr.png", true);
			ScalingPlotter plotter2 = new ScalingPlotter(outputFile+" NCIS plot");
			plotter2.saveXYplot(dataToPlot2, totalAtScaling, scalingRatio, "Binned Total Tag Count", "Marginal Signal/Control Ratio", outputFile+".TotalReadsAndNCIS_scaling-marginal.png", true);
			
			//Print data points to files
			try {
				FileWriter fout = new FileWriter(outputFile+".TotalReadsAndNCIS_scaling-ccr.count");
				for(int d=0; d<bintotals.size(); d++)
					fout.write(bintotals.get(d)+"\t"+ratios.get(d)+"\n");
				fout.close();
				FileWriter fout2 = new FileWriter(outputFile+".TotalReadsAndNCIS_scaling-marginal.count");
				for(int d=0; d<bintot.size(); d++)
					fout2.write(bintot.get(d)+"\t"+mratios.get(d)+"\n");
				fout2.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
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
		SEEDConfig sconfig = new SEEDConfig(gconfig, args);
		
		if(gconfig.helpWanted()){
			System.err.println("ExperimentScaler:");
			System.err.println(gconfig.getArgsList()+"\n"+econfig.getArgsList());
		}else{
			ExperimentManager exptMan = new ExperimentManager(econfig);
			
			
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
			DomainFinder potentialFilter = new DomainFinder(gconfig, econfig, sconfig, exptMan);
			Map<ExperimentCondition, List<Feature>> featuresByCond = potentialFilter.execute();
			List<Feature> potentials = new ArrayList<Feature>();
			for(ExperimentCondition ec : exptMan.getConditions())
				potentials.addAll(featuresByCond.get(ec));
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
		                for(Feature f : potentials){
		                	Region p = f.getCoords();
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
							System.out.println("NCIS\t"+sampA.getName()+" vs "+sampB.getName()+"\t"+scaler.scalingRatioByNCIS(sampleWindowCounts.get(sampA), sampleWindowCounts.get(sampB), null,econfig.getNCISMinBinFrac()));
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
			
			//Total tag normalization followed by NCIS
			for(Sample sampA : exptMan.getSamples()){ 
				if(sampA.isSignal()){
					for(Sample sampB : exptMan.getSamples())
						if(sampA!=null && sampB!=null && sampA.getIndex() != sampB.getIndex())
							System.out.println("HitRatioAndNCIS\t"+sampA.getName()+" vs "+sampB.getName()+"\t"+scaler.scalingRatioByHitRatioAndNCIS(sampleWindowCounts.get(sampA), sampleWindowCounts.get(sampB),
									sampA.getHitCount(), sampB.getHitCount(), null,econfig.getNCISMinBinFrac()));
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
		private void saveXYplot(Matrix datapoints, double scalingTotal, double scalingRatio, String xName, String yName, String outFilename, boolean rasterImage){
			this.setWidth(800);
			this.setHeight(800);
			this.addDataset("other", datapoints, new Color(75,75,75,80), 3);
			this.setXAxisLabel(xName);
			this.setYAxisLabel(yName);
			this.setXLogScale(true);
			this.setYLogScale(true);
			this.setXRangeFromData();
			this.setYRangeFromData();
			if(raxis.getRange().getLowerBound() >0.1){
			    raxis.setLowerBound(0.1);
			}
			this.addDomainMarker(scalingTotal);			
			this.addRangeMarker(scalingRatio);
			
			//Set the tick units according to the range
			double xUpper = daxis.getRange().getUpperBound();
			double xLower = daxis.getRange().getLowerBound();
			if(daxis instanceof org.jfree.chart.axis.NumberAxis)
	    		((NumberAxis)daxis).setTickUnit(new NumberTickUnit(5));
	    	double yUpper = raxis.getRange().getUpperBound();
			double yLower = raxis.getRange().getLowerBound();

	    	if(raxis instanceof org.jfree.chart.axis.NumberAxis)
	    		((NumberAxis)raxis).setTickUnit(new NumberTickUnit(5));
	    	
			try {
				this.saveImage(new File(outFilename), width, height, rasterImage);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
}
