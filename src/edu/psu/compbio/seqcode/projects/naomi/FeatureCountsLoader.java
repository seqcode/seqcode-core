package edu.psu.compbio.seqcode.projects.naomi;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;

public class FeatureCountsLoader {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ExperimentManager manager;
		
	protected List<StrandedPoint> strandedPoints;
	protected List<StrandedRegion> strandedRegions;
	
	protected int edge = 40; // Because shift will shorten the array, edge will ensure that windows are covered with reads
	protected int fivePrimeShift = 0;
	protected int window = 1000;	
	
	// sample counts array for each stranded region for each sample
	protected Map<Sample, Map<StrandedRegion,double[][]>> strandedRegionSampleCounts = new HashMap<Sample, Map<StrandedRegion,double[][]>>();
	protected Map<Sample,double[]> sampleComposite = new HashMap<Sample,double[]>();
	
	public FeatureCountsLoader(GenomeConfig gcon, ExptConfig econ, ExperimentManager man){	
		gconfig = gcon;
		econfig = econ;
		manager = man;
	}
	
	// setters
	public void setStrandedPoints(List<StrandedPoint> p){strandedPoints = p;}
	public void setStrandedRegions(List<StrandedRegion> reg){strandedRegions = reg;} 
	public void setWindowSize(int w){window = w;}
	public void setFivePrimeShift(int s){fivePrimeShift = s;}
	
	public Map<Sample, Map<StrandedRegion,double[][]>> strandedRegionSampleCounts(){
		
		// StrandedBaseCount list for each stranded regions for each sample
		Map<Sample, Map<StrandedRegion,List<StrandedBaseCount>>> sampleCountsMap = new HashMap<Sample, Map<StrandedRegion,List<StrandedBaseCount>>>();
				
		List<StrandedRegion> regionList = new ArrayList<StrandedRegion>();
		for(Point p: strandedPoints){		
			int start = Math.max(1, p.getLocation() - (window+edge)/2 );
			int end = Math.min(p.getLocation() + (window+edge)/2, p.getGenome().getChromLength(p.getChrom()));				
			StrandedRegion strandedReg = new StrandedRegion(p.getGenome(), p.getChrom(), start, end, p.getStrand());					
			regionList.add(strandedReg);
		}
		setStrandedRegions(regionList);
		
		for (ExperimentCondition condition : manager.getConditions()){		
			for (ControlledExperiment rep: condition.getReplicates()){				
				Map<StrandedRegion,List<StrandedBaseCount>> regionCounts =  new HashMap<StrandedRegion,List<StrandedBaseCount>>();				
				for (StrandedRegion reg : strandedRegions){
					regionCounts.put(reg, rep.getSignal().getBases(reg));
				}
				sampleCountsMap.put(rep.getSignal(),regionCounts);
			}					
		}
		
		//StrandedBasedCount object contains positive and negative strand separately
		// Reverse the array depending of strand of features	
		for (Sample sample : sampleCountsMap.keySet()){
			
			Map<StrandedRegion,double[][]> regionCounts = new HashMap<StrandedRegion,double[][]>();
			
			for (StrandedRegion reg : sampleCountsMap.get(sample).keySet()){			
				double[][] sampleCounts = new double[window+edge+1][2];
				for (int i = 0;i <= window+edge;i++){
					for (int s = 0; s<2; s++)
						sampleCounts[i][s] = 0;
				}	
				if (reg.getStrand() == '+'){ // regions(features) are positive strand
					reg.getMidpoint().getLocation();
					
					for (StrandedBaseCount hits: sampleCountsMap.get(sample).get(reg)){	
						if (hits.getStrand()=='+'){
							sampleCounts[hits.getCoordinate()-reg.getMidpoint().getLocation()+(window+edge)/2][0] = hits.getCount();
						}else{
							sampleCounts[hits.getCoordinate()-reg.getMidpoint().getLocation()+(window+edge)/2][1] = hits.getCount();
						}					
					}
				}else{ // if regions (features) are reverse strand, I need to flip the strands and locations
					for (StrandedBaseCount hits: sampleCountsMap.get(sample).get(reg)){	
						if (hits.getStrand()=='+'){
							sampleCounts[reg.getMidpoint().getLocation()-hits.getCoordinate()+(window+edge)/2][1] = hits.getCount();
						}else{
							sampleCounts[reg.getMidpoint().getLocation()-hits.getCoordinate()+(window+edge)/2][0] = hits.getCount();	
						}			
					}		
				}
				regionCounts.put(reg, sampleCounts);
			}
			strandedRegionSampleCounts.put(sample, regionCounts);
		}
		return strandedRegionSampleCounts;	
	}
	
	public Map<Sample,double[]> sampleComposite(){
		
		Map<Sample, Map<StrandedRegion,double[][]>> sRegionSampleCounts = strandedRegionSampleCounts();
					
		// nonstrandedComposite is a shifted and strand merged composite to be used to measure standard deviation
		for (Sample sample : sRegionSampleCounts.keySet()){
			
			double [] nonstrandedComposite = new double[window+1];
			for (int i = 0; i <=window ; i ++)
				nonstrandedComposite[i] = 0;
			
			for (Region reg : sRegionSampleCounts.get(sample).keySet()){				
				double[][] regionCounts = sRegionSampleCounts.get(sample).get(reg);		

				// get shifted composite for forward and reverse strands
				for (int j = 0 ; j <=window ; j++){
					nonstrandedComposite[j] += regionCounts[j-fivePrimeShift+edge/2][0]; 
					nonstrandedComposite[j] += regionCounts[j+fivePrimeShift+edge/2][1];
				}
			}
			sampleComposite.put(sample, nonstrandedComposite);		
		}
		return sampleComposite;
	}	

}
