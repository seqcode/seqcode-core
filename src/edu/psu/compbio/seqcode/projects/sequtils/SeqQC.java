package edu.psu.compbio.seqcode.projects.sequtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromosomeGenerator;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.RealValuedHistogram;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.framework.StrandedBaseCount;

public class SeqQC {

	Config config;
	Genome gen;
	ExperimentManager manager;
	float[] startcounts =null;
	float densityWindow = 500;
	double testQuantile = 0.05; //Base per-read count distribution on this proportion of the lowest-ranked densities
	int histoMax = 100;
	int sesScalingWin=1000;
	boolean verbose=false;
	boolean printNonVerboseHeader=true;
	HashMap<ControlledExperiment, String> infoStrings=new HashMap<ControlledExperiment, String>();
	
	/**
	 * Constructor
	 * @param c
	 * @param man
	 */
	public SeqQC(Config c){
		config = c;
		gen = config.getGenome();
		config.setPerBaseReadFiltering(false);
		config.setSESScaling(true);
		config.setScalingSlidingWindow(sesScalingWin);
		manager = new ExperimentManager(config);
		
		for(ControlledExperiment expt : manager.getExperimentSet().getReplicates()){
			String name = expt.getSignal().getName().startsWith("EXPERIMENT") ? expt.getSignal().getSourceName() : expt.getSignal().getName();
			infoStrings.put(expt, new String(name));
		}
	}
	
	//Accessors & settors
	public void setTestQuantile(double tq){
		if(tq>0 && tq<=1.0)
			testQuantile=tq;
	}
	public void setVerbose(boolean v){verbose=v;}
	public void setPrintHeader(boolean v){printNonVerboseHeader=v;}
	
	/**
	 * Run all implemented QC checks
	 */
	public void execute(){
		if(!verbose && printNonVerboseHeader)
			printHeader();
		
		//Library size calculation
		Map<ControlledExperiment, Double> meanFragmentCoverage = estimateLibrarySize();
		
		//Estimate signal proportion
		Map<ControlledExperiment, Double> sigProps = estimateSignalProportions(meanFragmentCoverage);
		
		if(!verbose)
			for(ControlledExperiment expt : manager.getExperimentSet().getReplicates())
				System.out.println(infoStrings.get(expt));
	}
	
	/**
	 * estimateLibrarySize
	 * returns the mean coverage of each fragment
	 */
	public Map<ControlledExperiment,Double> estimateLibrarySize(){
		Map<ControlledExperiment, Double> meanFragmentCoverage = new HashMap<ControlledExperiment, Double>();
		
		//If we have multiple experiments, process one at a time
		for(ControlledExperiment expt : manager.getExperimentSet().getReplicates()){
			// 1: Find the density of covered bases surrounding each read
			List<DensityCountPair> densities = new ArrayList<DensityCountPair>();
			Iterator<Region> chroms = new ChromosomeGenerator().execute(gen);
			while (chroms.hasNext()) {
				Region currentRegion = chroms.next();
				//Split the job up into large chunks
	            for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=config.MAXSECTION){
	                int y = x+config.MAXSECTION; 
	                if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
	                Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
				
					List<StrandedBaseCount> ipHits = expt.getSignal().getUnstrandedBases(currSubRegion);
					makeHitStartArray(ipHits, currSubRegion, '+');
					float posCounts[] = startcounts.clone();
					makeHitStartArray(ipHits, currSubRegion, '-');
					float negCounts[] = startcounts.clone();
		            
					int halfDWin = (int)densityWindow/2;
					for(int index=currSubRegion.getStart()+halfDWin; index<currSubRegion.getEnd()-halfDWin; index++){  //Note that this means a tiny fraction of reads that are located on the large block boundaries will be ignored
						int i = index-currSubRegion.getStart();
						if(posCounts[i]>0){ //Treat fragments on opposite strands separately
							float dens =0; 
							for(int j=i-halfDWin; j<i+halfDWin; j++){
								if(posCounts[j]>0 && j!=i)
									dens++;
								if(negCounts[j]>0)
									dens++;
							}
							dens /= (densityWindow*2);
							densities.add(new DensityCountPair(dens, posCounts[i]));
						}
						if(negCounts[i]>0){ //Treat fragments on opposite strands separately
							float dens =0; 
							for(int j=i-halfDWin; j<i+halfDWin; j++){
								if(posCounts[j]>0)	
									dens++;
								if(negCounts[j]>0 && j!=i)
									dens++;
							}
							dens /= (densityWindow*2);
							densities.add(new DensityCountPair(dens, negCounts[i]));
						}
					}
	            }
			}
			Collections.sort(densities); //Sort the density pairs in increasing order
			
			//2: Generate a read count per base distribution for the lowest density sites. 
			double currWeight=0, fracWeight=0, quantileLimit=testQuantile*expt.getSignal().getHitCount();
			RealValuedHistogram histo = new RealValuedHistogram(0,histoMax,histoMax);
			RealValuedHistogram fullHisto = new RealValuedHistogram(0,histoMax,histoMax);
			for(DensityCountPair dcp : densities){
				fullHisto.addValue(dcp.getCount());
				if(currWeight<quantileLimit){
					histo.addValue(dcp.getCount());
					fracWeight+=dcp.getCount();
				}
				currWeight+=dcp.getCount();
			}
			
			CensusLibraryComplexity census = new CensusLibraryComplexity(histo, 1, 20);
			census.setVerbose(verbose);
			census.execute();
			
			CensusLibraryComplexity fullCensus = new CensusLibraryComplexity(fullHisto, 1, 10);
			fullCensus.setVerbose(verbose);
			fullCensus.execute();
			
			meanFragmentCoverage.put(expt,  census.getEstimatedNegBinomialMean());
			
			//Get statistics on Fragments (i.e. partial histogram)
			double pLibrarySize = census.getEstimatedPoissonLibSize() / testQuantile;
			double pMean = census.getEstimatedPoissonLambda();
			double pLibraryCoverage = census.getEstimatedPoissonObservedCoverage();
			double pLibraryCoverageDoubledSeq = census.getEstimatedPoissonObservedGivenCount(2*fracWeight);
			double pNovelFragsUnderDoubleSeq = (pLibraryCoverageDoubledSeq-pLibraryCoverage)*pLibrarySize;
			double nbLibrarySize = census.getEstimatedNegBinomialLibSize() / testQuantile;
			double nbLibraryCoverage = census.getEstimatedNegBinomialObservedCoverage();
			double nbLibraryCoverageDoubledSeq = census.getEstimatedNegBinomialObservedGivenCount(2*fracWeight);
			double nbMean = census.getEstimatedNegBinomialMean();
			double nbVar = census.getEstimatedNegBinomialVariance();
			double nbNovelFragsUnderDoubleSeq = (nbLibraryCoverageDoubledSeq-nbLibraryCoverage)*nbLibrarySize;
			double lnLibrarySize = census.getEstimatedLogNormLibSize() / testQuantile;
			double lnLibraryCoverage = census.getEstimatedLogNormObservedCoverage();
			double lnLibraryCoverageDoubledSeq = census.getEstimatedLogNormObservedGivenCount(2*fracWeight);
			double lnMean = census.getEstimatedLogNormMean();
			//double lnNovelFragsUnderDoubleSeq = (lnLibraryCoverageDoubledSeq-lnLibraryCoverage)*lnLibrarySize;
			
			//Get statistics on Positions (i.e. full histogram)
			double pPosLibrarySize = fullCensus.getEstimatedPoissonLibSize();
			double pPosMean = fullCensus.getEstimatedPoissonLambda();
			double pPosLibraryCoverage = fullCensus.getEstimatedPoissonObservedCoverage();
			double pPosLibraryCoverageDoubledSeq = fullCensus.getEstimatedPoissonObservedGivenCount(2*currWeight);
			double pPosNovelFragsUnderDoubleSeq = (pPosLibraryCoverageDoubledSeq-pPosLibraryCoverage)*pPosLibrarySize;
			double nbPosLibrarySize = fullCensus.getEstimatedNegBinomialLibSize();
			double nbPosLibraryCoverage = fullCensus.getEstimatedNegBinomialObservedCoverage();
			double nbPosLibraryCoverageDoubledSeq = fullCensus.getEstimatedNegBinomialObservedGivenCount(2*currWeight);
			double nbPosMean = fullCensus.getEstimatedNegBinomialMean();
			double nbPosVar = fullCensus.getEstimatedNegBinomialVariance();
			double nbPosNovelFragsUnderDoubleSeq = (nbPosLibraryCoverageDoubledSeq-nbPosLibraryCoverage)*nbPosLibrarySize;
			double lnPosLibrarySize = fullCensus.getEstimatedLogNormLibSize();
			double lnPosLibraryCoverage = fullCensus.getEstimatedLogNormObservedCoverage();
			double lnPosLibraryCoverageDoubledSeq = fullCensus.getEstimatedLogNormObservedGivenCount(2*currWeight);
			double lnPosMean = fullCensus.getEstimatedLogNormMean();
			
			
			if(verbose){
				String name = expt.getSignal().getName().startsWith("EXPERIMENT") ? expt.getSignal().getSourceName() : expt.getSignal().getName();
				System.out.println("Experiment: "+name+" = "+String.format("%.1f mapped tags at %.0f unique positions", expt.getSignal().getHitCount(),expt.getSignal().getHitPositionCount()));
				
				System.out.println("\nPer-Fragment Estimated Poisson statistics:");
				System.out.println("Estimated initial mappable library size = "+String.format("%.1f", pLibrarySize) +" fragments");
				System.out.println(String.format("Each fragment was sequenced on average %.5f times", pMean));
				System.out.println(String.format("Proportion of library sequenced: %.3f", pLibraryCoverage));
				System.out.println(String.format("Proportion of library sequenced if you double the number of reads: %.3f (approx %.0f previously unsequenced fragments)", pLibraryCoverageDoubledSeq, pNovelFragsUnderDoubleSeq ));
				
				System.out.println("\nPer-Fragment Estimated Negative Binomial statistics:");
				System.out.println("Estimated initial mappable library size = "+String.format("%.1f", nbLibrarySize) +" fragments");
				System.out.println(String.format("Each fragment was sequenced on average %.5f times (%.5f variance, r = %.5f, p = %.5f, k = %.5f)", nbMean, nbVar, census.getEstimatedNegBinomialRP()[0], census.getEstimatedNegBinomialRP()[1], census.getEstimatedNegBinomialGammaK()));
				System.out.println(String.format("Proportion of library sequenced: %.3f", nbLibraryCoverage));
				System.out.println(String.format("Proportion of library sequenced if you double the number of reads: %.3f (approx %.0f previously unsequenced fragments)", nbLibraryCoverageDoubledSeq, nbNovelFragsUnderDoubleSeq ));
				
				System.out.println("\nPer-Fragment Estimated LogNormal statistics:");
				System.out.println("Estimated initial mappable library size = "+String.format("%.1f", lnLibrarySize) +" fragments");
				System.out.println(String.format("Each fragment was sequenced on average %.5f times ", lnMean));
				System.out.println(String.format("Proportion of library sequenced: %.3f", lnLibraryCoverage));
				//System.out.println(String.format("Proportion of library sequenced if you double the number of reads: %.3f (approx %.0f previously unsequenced fragments)", lnLibraryCoverageDoubledSeq, lnNovelFragsUnderDoubleSeq ));
				
				System.out.println("\nPer-Position Estimated Poisson statistics:");
				System.out.println("Estimated initial mappable positions = "+String.format("%.1f", pPosLibrarySize) +" fragments");
				System.out.println(String.format("Each fragment was sequenced on average %.5f times", pPosMean));
				System.out.println(String.format("Proportion of positions sequenced: %.3f", pPosLibraryCoverage));
				System.out.println(String.format("Proportion of positions sequenced if you double the number of reads: %.3f (approx %.0f previously unsequenced fragments)", pPosLibraryCoverageDoubledSeq, pPosNovelFragsUnderDoubleSeq ));
				
				System.out.println("\nPer-Position Estimated Negative Binomial statistics:");
				System.out.println("Estimated initial mappable positions = "+String.format("%.1f", nbPosLibrarySize) +" positions");
				System.out.println(String.format("Each position was sequenced on average %.5f times (%.5f variance, r = %.5f, p = %.5f, k = %.5f)", nbPosMean, nbPosVar, fullCensus.getEstimatedNegBinomialRP()[0], fullCensus.getEstimatedNegBinomialRP()[1], fullCensus.getEstimatedNegBinomialGammaK()));
				System.out.println(String.format("Proportion of positions sequenced: %.3f", nbPosLibraryCoverage));
				System.out.println(String.format("Proportion of positions sequenced if you double the number of reads: %.3f (approx %.0f previously unsequenced fragments)", nbPosLibraryCoverageDoubledSeq, nbPosNovelFragsUnderDoubleSeq ));
				
				System.out.println("\nPer-Position Estimated LogNormal statistics:");
				System.out.println("Estimated initial mappable positions = "+String.format("%.1f", lnPosLibrarySize) +" positions");
				System.out.println(String.format("Each position was sequenced on average %.5f times", lnPosMean));
				System.out.println(String.format("Proportion of positions sequenced: %.3f", lnPosLibraryCoverage));
				//System.out.println(String.format("Proportion of positions sequenced if you double the number of reads: %.3f (approx %.0f previously unsequenced fragments)", lnPosLibraryCoverageDoubledSeq, lnPosNovelFragsUnderDoubleSeq ));
				
			}else{
				String currInfo = infoStrings.get(expt);
				currInfo = currInfo + String.format("\t%.1f\t%.0f\t%.5f\t%.5f\t%.1f\t%.1f\t%.3f\t%.3f\t%.3f\t%.3f", 
						expt.getSignal().getHitCount(),
						expt.getSignal().getHitPositionCount(),
						nbMean,
						nbVar,
						nbLibrarySize,
						nbPosLibrarySize,
						nbLibraryCoverage,
						nbPosLibraryCoverage,
						nbLibraryCoverageDoubledSeq,
						nbPosLibraryCoverageDoubledSeq
						);
				infoStrings.put(expt, currInfo);
			}
		}
		return meanFragmentCoverage;
	}
	
	/**
	 * estimateSignalProportions
	 * @param meanFragmentCoverage
	 */
	public Map<ControlledExperiment,Double> estimateSignalProportions(Map<ControlledExperiment,Double> meanFragmentCoverage){
		Map<ControlledExperiment, Double> sigProps = new HashMap<ControlledExperiment, Double>();
		
		if(verbose)
			System.out.println("\nEstimating signal proportion");
		
		//Estimate proportion via scaling
		for(ControlledExperiment expt : manager.getExperimentSet().getReplicates()){
			
			float scaling = meanFragmentCoverage.get(expt).floatValue();
			if(scaling>1)
				expt.correctSignalCounts(scaling);
			
			if(expt.hasControl()){
				double readRatio = expt.getSignal().getHitCount()/expt.getControl().getHitCount();
				sigProps.put(expt, expt.getSigProp());
				
				if(verbose){
					System.out.println(String.format("Read count ratio for: %s vs %s = %.5f", expt.getSignal().getName(), expt.getControl().getName(), readRatio));
					System.out.println(String.format("Scaling factor via SES scaling for: %s = %.5f", expt.getName(), expt.getControlScaling()));
					System.out.println(String.format("Signal proportion via SES scaling for: %s = %.5f", expt.getName(), expt.getSigProp()));
				}else{
					String currInfo = infoStrings.get(expt);
					currInfo = currInfo + String.format("\t%.5f\t%.5f",
							expt.getControlScaling(), expt.getSigProp());
					infoStrings.put(expt, currInfo);
				}
				//ExperimentScaler scaler = new ExperimentScaler(expt.getSignal(), expt.getControl());
				//scaler.scalingRatioByMedian(1000);
				//scaler.scalingRatioByRegression(1000);
			}else{
				String currInfo = infoStrings.get(expt);
				currInfo = currInfo + "\tNA\tNA";
				infoStrings.put(expt, currInfo);
			}
		}
		
		/*//Estimate proportion via distribution-fitting
		BackgroundDetector detector = new BackgroundDetector(config, manager, 200, 100);
		HashMap<Sample, Double> detectedBP = detector.execute();
		for(ControlledExperiment expt : manager.getExperimentSet().getReplicates()){
			System.out.println(String.format("Signal proportion via fitting for: %s = %.5f", expt.getName(), 1-detectedBP.get(expt.getSignal())));
		}*/
		
		return sigProps;
	}
	
	
	//Makes integer arrays corresponding to the read starts over the current region
	protected void makeHitStartArray(List<StrandedBaseCount> hits, Region currReg, char strand){
		startcounts = new float[currReg.getWidth()+1];
		for(int i=0; i<=currReg.getWidth(); i++){startcounts[i]=0;}
		for(StrandedBaseCount r : hits){
			if(strand=='.' || r.getStrand()==strand){
				int offset=inBounds(r.getCoordinate()-currReg.getStart(),0,currReg.getWidth());
				startcounts[offset]+=r.getCount();
			}
		}
	}
	protected final int inBounds(int x, int min, int max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
	
	/**
	 * Main driver method
	 */
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
		if(args.length==0 || ap.hasKey("h")){
			System.err.println("SeqQC:\n" +
					"\t--geninfo <genome info file>\n" +
					"\t--expt <experiment to test QC>\n" +
					"\t--ctrl <control for this experiment>\n" +
					"\t--format <format of experiment files (SAM/BED/IDX)>\n" +
					"\t--testquantile <proportion of lowest ranked density reads to build distribution from>\n" +
					"\t--verbose\n" +
					"\t--noheader\n");
		}else{
			Config con = new Config(args, false);
			
			SeqQC qc = new SeqQC(con);
			if(ap.hasKey("testquantile"))
				qc.setTestQuantile(new Double(ap.getKeyValue("testquantile")));
			if(ap.hasKey("verbose"))
				qc.setVerbose(true);
			if(ap.hasKey("noheader"))
				qc.setPrintHeader(false);
			qc.execute();
		}
	}
	
	private void printHeader(){
		System.out.println("Dataset\t" +
				"MappedTags\t" +
				"MappedPositions\t" +
				"EstNBMean\t" +
				"EstNBVar\t" +
				"EstFragmentLibSize\t" +
				"EstTotalPositions\t" +
				"EstFragmentCoverage\t" +
				"EstPositionCoverage\t" +
				"EstFragmentCoverageDoubleSeq\t" +
				"EstPositionCoverageDoubleSeq\t" +
				"");
	}
	
	private class DensityCountPair implements Comparable<DensityCountPair>{
		private float density;
		private float count;
		
		public DensityCountPair(float d, float c){
			density = d;
			count =c;
		}
		
		public void setCount(float count) {
			this.count = count;
		}
		public void setDensity(float density) {
			this.density = density;
		}
		public float getCount() {
			return count;
		}
		public float getDensity() {
			return density;
		}
		
		// sort according to density
		public int compareTo(DensityCountPair dc) {
			if(density < dc.density)
				return -1;
			else if (density > dc.density)
				return 1;
			else return 0;
		}

	}
}
