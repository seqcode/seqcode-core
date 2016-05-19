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

/**
 * FeatureCountsLoader : load counts from stranded feature
 * 
 * input : reference point
 * 
 * @author naomi yamada
 */

public class FeatureCountsLoader {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ExperimentManager manager;
		
	protected List<StrandedPoint> strandedPoints;
	protected List<StrandedRegion> strandedRegions;
	
	protected int edge = 0; // Because shift will shorten the array, edge will ensure that windows are covered with reads
	protected int fivePrimeShift = 0;
	protected int window = 1000;	
	
	// sample counts array for each stranded region for each sample
//	protected Map<Sample, Map<StrandedRegion,double[][]>> strandedRegionSampleCounts = new HashMap<Sample, Map<StrandedRegion,double[][]>>();
	protected Map<ControlledExperiment, Map<StrandedRegion,double[][]>> strandedRegionControlCounts = new HashMap<ControlledExperiment, Map<StrandedRegion,double[][]>>();
	protected Map<ControlledExperiment, double[]> controlCompositeMap = new HashMap<ControlledExperiment,double[]>();
	
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
	
	public Map<ControlledExperiment, Map<StrandedRegion,double[][]>> getControlStrandedRegionCounts(){return strandedRegionControlCounts;}
	public Map<ControlledExperiment, double[]> controlComposite(){return controlCompositeMap;}
	
	public Map<ControlledExperiment, Map<StrandedRegion,double[][]>> strandedRegionSampleCounts(){
		
		// if there tag shift is needed, set the edge to 40 to obtain counts from extra region
		if (fivePrimeShift > 0){edge=40;}
		
		// StrandedBaseCount list for each stranded regions for each sample
		Map<ControlledExperiment, Map<StrandedRegion,List<StrandedBaseCount>>> sampleCountsMap = new HashMap<ControlledExperiment, Map<StrandedRegion,List<StrandedBaseCount>>>();
		Map<ControlledExperiment, Map<StrandedRegion,List<StrandedBaseCount>>> controlCountsMap = new HashMap<ControlledExperiment, Map<StrandedRegion,List<StrandedBaseCount>>>();
		Map<ControlledExperiment, Map<StrandedRegion,double[][]>> strandedRegionSampleCounts = new HashMap<ControlledExperiment, Map<StrandedRegion,double[][]>>();	
		//<Sample, Map<StrandedRegion,double[][]>> strandedRegionControlCounts = new HashMap<Sample, Map<StrandedRegion,double[][]>>();	
		
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
				Map<StrandedRegion,List<StrandedBaseCount>> sampleRegionCounts =  new HashMap<StrandedRegion,List<StrandedBaseCount>>();
				Map<StrandedRegion,List<StrandedBaseCount>> controlRegionCounts =  new HashMap<StrandedRegion,List<StrandedBaseCount>>();	
				for (StrandedRegion reg : strandedRegions){
					sampleRegionCounts.put(reg, rep.getSignal().getBases(reg));
					if (rep.hasControl()){
						controlRegionCounts.put(reg, rep.getControl().getBases(reg));
					}
					controlRegionCounts.put(reg, rep.getControl().getBases(reg));
				}
				sampleCountsMap.put(rep,sampleRegionCounts);
				controlCountsMap.put(rep,controlRegionCounts);
			}					
		}
		
		//StrandedBasedCount object contains positive and negative strand separately
		// Reverse the array depending of strand of features	
		for (ControlledExperiment replicates : sampleCountsMap.keySet()){
			
			Map<StrandedRegion,double[][]> sampleRegionCounts = new HashMap<StrandedRegion,double[][]>();
			Map<StrandedRegion,double[][]> controlRegionCounts = new HashMap<StrandedRegion,double[][]>();
			
			for (StrandedRegion reg : sampleCountsMap.get(replicates).keySet()){			
				double[][] sampleCounts = new double[window+edge+1][2];
				double[][] controlCounts = new double[window+edge+1][2];
				for (int i = 0;i <= window+edge;i++){
					for (int s = 0; s<2; s++){
						sampleCounts[i][s] = 0;
						controlCounts[i][s] = 0;
					}
				}	
				if (reg.getStrand() == '+'){ // regions(features) are positive strand
					reg.getMidpoint().getLocation();
					
					for (StrandedBaseCount hits: sampleCountsMap.get(replicates).get(reg)){	
						if (hits.getStrand()=='+'){
							sampleCounts[hits.getCoordinate()-reg.getMidpoint().getLocation()+(window+edge)/2][0] = hits.getCount();
						}else{
							sampleCounts[hits.getCoordinate()-reg.getMidpoint().getLocation()+(window+edge)/2][1] = hits.getCount();
						}					
					}				
					for (StrandedBaseCount hits: controlCountsMap.get(replicates).get(reg)){	
						if (hits.getStrand()=='+'){
							controlCounts[hits.getCoordinate()-reg.getMidpoint().getLocation()+(window+edge)/2][0] = hits.getCount();
						}else{
							controlCounts[hits.getCoordinate()-reg.getMidpoint().getLocation()+(window+edge)/2][1] = hits.getCount();
						}					
					}				
				}else{ // if regions (features) are reverse strand, I need to flip the strands and locations
					for (StrandedBaseCount hits: sampleCountsMap.get(replicates).get(reg)){	
						if (hits.getStrand()=='+'){
							sampleCounts[reg.getMidpoint().getLocation()-hits.getCoordinate()+(window+edge)/2][1] = hits.getCount();
						}else{
							sampleCounts[reg.getMidpoint().getLocation()-hits.getCoordinate()+(window+edge)/2][0] = hits.getCount();	
						}			
					}	
					for (StrandedBaseCount hits: controlCountsMap.get(replicates).get(reg)){	
						if (hits.getStrand()=='+'){
							controlCounts[reg.getMidpoint().getLocation()-hits.getCoordinate()+(window+edge)/2][1] = hits.getCount();
						}else{
							controlCounts[reg.getMidpoint().getLocation()-hits.getCoordinate()+(window+edge)/2][0] = hits.getCount();	
						}			
					}	
				}				
				sampleRegionCounts.put(reg, sampleCounts);
				controlRegionCounts.put(reg, controlCounts);
			}
			strandedRegionSampleCounts.put(replicates, sampleRegionCounts);
			strandedRegionControlCounts.put(replicates, controlRegionCounts);
		}
		return strandedRegionSampleCounts;
	}
	
	public Map<ControlledExperiment,double[]> sampleComposite(){
		
		Map<ControlledExperiment, Map<StrandedRegion,double[][]>> sampleRegionComposite = strandedRegionSampleCounts();
		Map<ControlledExperiment, Map<StrandedRegion,double[][]>> controlRegionComposite = getControlStrandedRegionCounts();
		
		Map<ControlledExperiment,double[]> sampleCompositeMap = new HashMap<ControlledExperiment,double[]>();
					
		// nonstrandedComposite is a shifted and strand merged composite to be used to measure standard deviation
		for (ControlledExperiment rep : sampleRegionComposite.keySet()){
			
			double [] sampleComposite = new double[window+1];
			double [] controlComposite = new double[window+1];
			for (int i = 0; i <=window ; i ++){
				sampleComposite[i] = 0;
				controlComposite[i] = 0;
			}
			
			for (Region reg : sampleRegionComposite.get(rep).keySet()){				
				double[][] sampleRegionCounts = sampleRegionComposite.get(rep).get(reg);		
				double[][] controlRegionCounts = controlRegionComposite.get(rep).get(reg);	

				// get shifted composite for forward and reverse strands
				for (int j = 0 ; j <=window ; j++){
					sampleComposite[j] += sampleRegionCounts[j-fivePrimeShift+edge/2][0]; 
					sampleComposite[j] += sampleRegionCounts[j+fivePrimeShift+edge/2][1];
					controlComposite[j] += controlRegionCounts[j-fivePrimeShift+edge/2][0]; 
					controlComposite[j] += controlRegionCounts[j+fivePrimeShift+edge/2][1];
				}
			}
			sampleCompositeMap.put(rep, sampleComposite);	
			controlCompositeMap.put(rep, controlComposite);
		}
		return sampleCompositeMap;
	}	

}
