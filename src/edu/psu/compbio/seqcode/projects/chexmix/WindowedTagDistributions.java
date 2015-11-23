package edu.psu.compbio.seqcode.projects.chexmix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;

public class WindowedTagDistributions {
	protected ExperimentManager exptMan;	
	protected List<Region> regions = new ArrayList<Region>();
	protected int win;
	protected int numConditions;
	protected int numRegions;
	protected double[][][] perRegionPlus; //per-point, per-condition plus strand tags  {point, condition, location}
	protected double[][][] perRegionMinus;  //per-point, per-condition minus strand tags   {point, condition, location}
	protected HashMap<Region,Integer> regionIndex = new HashMap<Region,Integer>();
	protected boolean isSignal;
	
	public WindowedTagDistributions(List<Point> points, ExperimentManager eMan, int win, boolean loadSignal){
		exptMan = eMan;
		this.win = win;
		this.numConditions=exptMan.getNumConditions();
		isSignal = loadSignal;
	
		for(Point pt : points)
			regions.add(pt.expand(win/2));
		numRegions = regions.size();
		
		perRegionPlus = new double[numRegions][numConditions][win];
		perRegionMinus = new double[numRegions][numConditions][win];
		
		for(int p=0; p<numRegions; p++)
			regionIndex.put(regions.get(p), p);

		//Reset
		for(int c=0; c<numConditions; c++){
			for(int p=0; p<numRegions; p++)
				for(int w=0; w<win; w++){
					perRegionPlus[p][c][w]=0; perRegionMinus[p][c][w]=0;
				}
		}
		
		for(ExperimentCondition cond : exptMan.getConditions()){
			for(ControlledExperiment rep : cond.getReplicates()){
				
				if(loadSignal || rep.hasControl()){
					//Iterate through regions
					int p=0;
					for(Region reg : regions){
						//Load reads
						List<StrandedBaseCount> pReads = loadSignal ? 
								rep.getSignal().getStrandedBases(reg, '+') : 
									rep.getControl().getStrandedBases(reg, '+');
						List<StrandedBaseCount> mReads = loadSignal ? 
								rep.getSignal().getStrandedBases(reg, '-') :
									rep.getControl().getStrandedBases(reg, '-');
						
						
						for(StrandedBaseCount sbc : pReads){
							int sdist = sbc.getCoordinate()-reg.getStart();
							if(sdist>=0 && sdist<win){
								perRegionPlus[p][cond.getIndex()][sdist]+=sbc.getCount();
							}
						}
						for(StrandedBaseCount sbc : mReads){
							int sdist = sbc.getCoordinate()-reg.getStart();
							if(sdist>=0 && sdist<win){
								perRegionMinus[p][cond.getIndex()][sdist]+=sbc.getCount();
							}
						}
						
						p++;
					}
				}
			}
		}
	}
	
	//Accessors
	public int getWinSize(){return win;}
	public int getNumConditions(){return numConditions;}
	public double[] getRegionPlus(Region r, ExperimentCondition c){
		return perRegionPlus[regionIndex.get(r)][c.getIndex()];
	}
	public double[] getRegionMinus(Region r, ExperimentCondition c){
		return perRegionMinus[regionIndex.get(r)][c.getIndex()];
	}
	public double[][] getRegionPlus(int index){return perRegionPlus[index];}
	public double[][] getRegionMinus(int index){return perRegionMinus[index];}
	public List<Region> getRegions(){return regions;}
	public Region getRegion(int i){return regions.get(i);}

}
