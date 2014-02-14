package edu.psu.compbio.seqcode.projects.akshay.bayesments.features;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.utils.CoutPlotter;
import edu.psu.compbio.seqcode.projects.shaun.EventMetaMaker;

/**
 * Class the reads and stores the training data for the Bayesian network
 * @author akshaykakumanu
 *
 */

public class GenomicLocations {
	
	protected Genome gen;
	protected ExperimentSet experiments;
	// List of points, which are used to train the model
	protected List<Point> locations = new ArrayList<Point>();
	//2-d array of counts, rows as points and columns as chromatin condition
	protected float[][] chromatinCounts;
	//2-d array of counts, rows as points and columns as factor conditions 
	protected float[][] factorCounts;
	// 2-d array of regions, rows as points and columns as conditions
	//(note: factor regions and chromatin regions may have different win sizes)
	protected Region[][] chromatinRegions;
	protected Region[][] factorRegions;
	
	// number of chromatin and factor conditions
	protected int numChromCons;
	protected int numFacCons;
	
	/**
	 *
	 * @param manager
	 * @param config
	 */
	public GenomicLocations(ExperimentManager manager, Config config) {
		try{
			File peaksFile = config.getPeaksFile();
			locations = EventMetaMaker.loadPoints(peaksFile, this.gen);
			this.experiments = manager.getExperimentSet();
			this.gen = manager.getGenome();
			
			//Filling chromatinRegions
			List<ExperimentCondition> chromConds = manager.getChromatinConditionList();
			this.numChromCons = chromConds.size();
			chromatinRegions = new Region[locations.size()][chromConds.size()];
			chromatinCounts = new float[locations.size()][chromConds.size()];
			int pointCount=0;
			for(Point p : locations){
				int conCount=0;
				for(ExperimentCondition ec : chromConds){
					chromatinRegions[pointCount][conCount] = new Region(this.gen,p.getChrom()
							,p.getLocation()-ec.getWinSize(),p.getLocation()+ec.getWinSize());
					conCount++;
				}
				pointCount++;
			}
			
			//Filling factorRegions
			pointCount=0;
			List<ExperimentCondition> facConds = manager.getFacConditionList();
			this.numFacCons = facConds.size();
			factorRegions = new Region[locations.size()][facConds.size()];
			factorCounts = new float[locations.size()][facConds.size()];
			for(Point p : locations){
				int conCount=0;
				for(ExperimentCondition ec : facConds){
					factorRegions[pointCount][conCount] = new Region(this.gen,p.getChrom()
							, p.getLocation()-ec.getWinSize(), p.getLocation()+ec.getWinSize());
					conCount++;
				}
				pointCount++;
			}
			
			//Filling chromatinCounts
			pointCount=0;
			for(Point p : locations){
				int conCount=0;
				for(ExperimentCondition ec : chromConds){
					Region r = chromatinRegions[pointCount][conCount];
					chromatinCounts[pointCount][conCount] = ec.getTotalSignalCountInARegion(r);
					conCount++;
				}
				pointCount++;
			}
			
			//Filling factorCounts
			pointCount=0;
			for(Point p: locations){
				int conCount=0;
				for(ExperimentCondition ec : facConds){
					Region r = factorRegions[pointCount][conCount];
					factorCounts[pointCount][conCount] = ec.getTotalSignalCountInARegion(r);
					conCount++;
				}
				pointCount++;
			}
			
			// If flag is on, do asinh transformation
			if(config.doAffine()){
				updataAsineTransformation();
			}
		}
		catch(IOException e){
			e.printStackTrace();
		}
	}
	/**
	 * Does asinh transformation of the counts data and updates it to the transformed values
	 */
	public void updataAsineTransformation(){
		for(int i=0; i<locations.size(); i++){
			for(int c=0; c<numChromCons; c++){
				//System.out.println("Before transformation "+Float.toString(chromatinCounts[i][c]));
				chromatinCounts[i][c] = (float) Math.log((double)chromatinCounts[i][c]+(double) Math.sqrt(Math.pow((double)chromatinCounts[i][c], 2.0)+1.0));
				//System.out.println("after transformation "+Float.toString(chromatinCounts[i][c]));
			}
		}
		
		for(int i=0; i<locations.size(); i++){
			for(int f=0; f<numFacCons; f++){
				factorCounts[i][f] = (float) Math.log((double)factorCounts[i][f]+(double) Math.sqrt(Math.pow((double)factorCounts[i][f], 2.0)+1.0));
			}
		}
	}
	
	//Accessors
	public int getNumChromatinCons(){return this.numChromCons;}
	public int getNumFacCons(){return this.numFacCons;}
	public int getNumTrainingExamples(){return this.locations.size();}
	public float[][] getChromatinCounts(){return this.chromatinCounts;}
	public float[][] getFactorCounts(){return this.factorCounts;}
	public List<Point> getLocations(){return this.locations;}
	
	
	/**
	 * Calls the CoutPlotter class and plots the Observed experimental tracks
	 * @param conf
	 * @param manager
	 */
	public void plotData(Config conf,ExperimentManager manager){
		for(int c=0; c<numChromCons; c++){
			float[] counts = new float[locations.size()];
			for(int i=0; i<locations.size(); i++){
				counts[i] = chromatinCounts[i][c];
			}
			String name_tag = manager.getChromatinConditionList().get(c).getName();
			CoutPlotter cp = new CoutPlotter(counts, conf, name_tag);
			cp.plot();
		}
		for(int f=0; f<numFacCons; f++){
			float[] counts = new float[locations.size()];
			for(int i=0; i<locations.size(); i++){
				counts[i] = factorCounts[i][f];
			}
			String name_tag = manager.getFacConditionList().get(f).getName();
			CoutPlotter cp = new CoutPlotter(counts, conf, name_tag);
			cp.plot();
		}
		
	}
	
	

}
