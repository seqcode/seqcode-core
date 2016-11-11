package org.seqcode.projects.galaxyexo;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;


/**
 * FeatureCountsLoader : load counts from stranded feature
 * 
 * input : reference point
 * 
 * @author naomi yamada
 */

public class FeatureCountsLoader {	
	protected GenomeConfig gconfig;		
	protected List<StrandedPoint> strandedPoints;
	protected List<StrandedRegion> strandedRegions;
	
	protected int edge = 0; // Because shift will shorten the array, edge will ensure that windows are covered with reads
	protected int fivePrimeShift = 0;
	protected int windowSize = 1000;	
	
	// sample counts array for each stranded region for each sample
	protected Map<StrandedRegion,double[][]> controlRegionCounts = new HashMap<StrandedRegion,double[][]>();
	protected double [] controlComposite;
	
	public FeatureCountsLoader(GenomeConfig gcon, List<StrandedPoint> p, int win){	
		gconfig = gcon;
		strandedPoints = p;
		windowSize = win;
	}
	public FeatureCountsLoader(GenomeConfig gcon, List<StrandedPoint> p){	
		gconfig = gcon;
		strandedPoints = p;
	}
	
	// setters
	public void setStrandedRegions(List<StrandedRegion> reg){strandedRegions = reg;} 
	public void setFivePrimeShift(int s){fivePrimeShift = s;}
	public void setWindowSize(int win) {windowSize = win;}
	// getters
	public List<StrandedRegion> getStrandedRegions(){return strandedRegions;}
	public Map<StrandedRegion,double[][]> getControlStrandedRegionCounts(){return controlRegionCounts;}
	public double[] controlComposite(){return controlComposite;}
	
	public Map<StrandedRegion,double[][]> strandedRegionSampleCounts(ControlledExperiment rep){
		
		// if there tag shift is needed, set the edge to 40 to obtain counts from extra region
		if (fivePrimeShift > 0){edge=40;}	
		
		List<StrandedRegion> regionList = new ArrayList<StrandedRegion>();
		for(Point p: strandedPoints){		
			int start = Math.max(1, p.getLocation() - (windowSize+edge)/2 );
			int end = Math.min(p.getLocation() + (windowSize+edge)/2, p.getGenome().getChromLength(p.getChrom()));				
			StrandedRegion strandedReg = new StrandedRegion(p.getGenome(), p.getChrom(), start, end, p.getStrand());					
			regionList.add(strandedReg);
		}
		setStrandedRegions(regionList);
					
		Map<StrandedRegion,List<StrandedBaseCount>> sampleCountsMap =  new HashMap<StrandedRegion,List<StrandedBaseCount>>();
		Map<StrandedRegion,List<StrandedBaseCount>> controlCountsMap =  new HashMap<StrandedRegion,List<StrandedBaseCount>>();	
		for (StrandedRegion reg : strandedRegions){
			sampleCountsMap.put(reg, rep.getSignal().getBases(reg));
			if (rep.hasControl()){
				controlCountsMap.put(reg, rep.getControl().getBases(reg));
			}
		}
		
		//StrandedBasedCount object contains positive and negative strand separately
		// Reverse the array depending of strand of features			
		Map<StrandedRegion,double[][]> sampleRegionCounts = new HashMap<StrandedRegion,double[][]>();			
		for (StrandedRegion reg : sampleCountsMap.keySet()){			
			double[][] sampleCounts = new double[windowSize+edge+1][2];
			double[][] controlCounts = new double[windowSize+edge+1][2];
			for (int i = 0;i <= windowSize+edge;i++){
				for (int s = 0; s<2; s++){
					sampleCounts[i][s] = 0;
					controlCounts[i][s] = 0;
				}
			}	
			if (reg.getStrand() == '+'){ // regions(features) are positive strand					
				for (StrandedBaseCount hits: sampleCountsMap.get(reg)){	
					if (hits.getStrand()=='+'){
						sampleCounts[hits.getCoordinate()-reg.getMidpoint().getLocation()+(windowSize+edge)/2][0] = hits.getCount();
					}else{
						sampleCounts[hits.getCoordinate()-reg.getMidpoint().getLocation()+(windowSize+edge)/2][1] = hits.getCount();
					}					
				}
			}else{ // if regions (features) are reverse strand, I need to flip the strands and locations
				for (StrandedBaseCount hits: sampleCountsMap.get(reg)){	
					if (hits.getStrand()=='+'){
						sampleCounts[reg.getMidpoint().getLocation()-hits.getCoordinate()+(windowSize+edge)/2][1] = hits.getCount();
					}else{
						sampleCounts[reg.getMidpoint().getLocation()-hits.getCoordinate()+(windowSize+edge)/2][0] = hits.getCount();	
					}			
				}	
			}
				
			// only execute if controls are loaded
			if (rep.hasControl()){
				if (reg.getStrand() == '+'){ 				
					for (StrandedBaseCount hits: controlCountsMap.get(reg)){	
						if (hits.getStrand()=='+'){
							controlCounts[hits.getCoordinate()-reg.getMidpoint().getLocation()+(windowSize+edge)/2][0] = hits.getCount();
						}else{
							controlCounts[hits.getCoordinate()-reg.getMidpoint().getLocation()+(windowSize+edge)/2][1] = hits.getCount();
						}					
					}	
				}else{
					for (StrandedBaseCount hits: controlCountsMap.get(reg)){	
						if (hits.getStrand()=='+'){
							controlCounts[reg.getMidpoint().getLocation()-hits.getCoordinate()+(windowSize+edge)/2][1] = hits.getCount();
						}else{
							controlCounts[reg.getMidpoint().getLocation()-hits.getCoordinate()+(windowSize+edge)/2][0] = hits.getCount();	
						}			
					}	
				}
			} // end of control				
			sampleRegionCounts.put(reg, sampleCounts);
			controlRegionCounts.put(reg, controlCounts);
		}
		return sampleRegionCounts;
	}
	
	public double[] sampleComposite(ControlledExperiment rep){
		
		Map<StrandedRegion,double[][]> sampleRegionComposite = strandedRegionSampleCounts(rep);
		Map<StrandedRegion,double[][]> controlRegionComposite = getControlStrandedRegionCounts();
					
		// nonstrandedComposite is a shifted and strand merged composite to be used to measure standard deviation
		double [] sampleComposite = new double[windowSize+1];
		controlComposite = new double[windowSize+1];
		for (int i = 0; i <=windowSize ; i ++){
			sampleComposite[i] = 0;
			controlComposite[i] = 0;
		}
			
		for (Region reg : sampleRegionComposite.keySet()){				
			double[][] sampleRegionCounts = sampleRegionComposite.get(reg);		
			double[][] controlRegionCounts = controlRegionComposite.get(reg);	

			// get shifted composite for forward and reverse strands
			for (int j = 0 ; j <=windowSize ; j++){
				sampleComposite[j] += sampleRegionCounts[j-fivePrimeShift+edge/2][0]; 
				sampleComposite[j] += sampleRegionCounts[j+fivePrimeShift+edge/2][1];
				controlComposite[j] += controlRegionCounts[j-fivePrimeShift+edge/2][0]; 
				controlComposite[j] += controlRegionCounts[j+fivePrimeShift+edge/2][1];
			}
		}
		return sampleComposite;
	}
}
