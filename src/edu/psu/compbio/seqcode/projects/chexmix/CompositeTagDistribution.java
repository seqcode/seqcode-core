package edu.psu.compbio.seqcode.projects.chexmix;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.projects.multigps.utilities.Utils;

/**
 * CompositeTagDistribution: watson/crick tag distributions for a collection of aligned points and the resulting composite
 *
 * 	Coordinates are 0-based, and the center of the distribution/alignment is defined by the centerOffset variable	
 * @author mahony
 *
 */
public class CompositeTagDistribution {
	protected ExperimentManager exptMan;	
	protected List<StrandedPoint> points;
	protected int win;
	protected int centerOffset;
	protected int numConditions;
	protected int numPoints;
	protected double[][] watson; //per-condition watson tags
	protected double[][] crick;  //per-condition crick tags
	protected double[][][] perPointWatson; //per-point, per-condition watson tags
	protected double[][][] perPointCrick;  //per-point, per-condition crick tags
	protected HashMap<StrandedPoint,Integer> pointIndex = new HashMap<StrandedPoint,Integer>();
	protected boolean isSignal;
	
	public CompositeTagDistribution(List<StrandedPoint> points, ExperimentManager eMan, int win, boolean loadSignal){
		exptMan = eMan;
		this.win = win;
		centerOffset = win/2;
		this.numConditions=exptMan.getNumConditions();
		this.points = points;
		numPoints = points.size();
		isSignal = loadSignal;
	
		watson = new double[numConditions][win];
		crick = new double[numConditions][win];
		perPointWatson = new double[numPoints][numConditions][win];
		perPointCrick = new double[numPoints][numConditions][win];
		
		for(int p=0; p<numPoints; p++)
			pointIndex.put(points.get(p), p);

		//Reset
		for(int c=0; c<numConditions; c++){
			for(int w=0; w<win; w++){watson[c][w]=0; crick[c][w]=0;}
			for(int p=0; p<numPoints; p++)
				for(int w=0; w<win; w++){
					perPointWatson[p][c][w]=0; perPointCrick[p][c][w]=0;
				}
		}
		
		for(ExperimentCondition cond : exptMan.getConditions()){
			for(ControlledExperiment rep : cond.getReplicates()){
				
				if(loadSignal || rep.hasControl()){
					//Iterate through points
					int p=0;
					for(StrandedPoint pt : points){
						//Load reads
						List<StrandedBaseCount> wReads = loadSignal ? 
								rep.getSignal().getStrandedBases(pt.expand(win), pt.getStrand()) : 
									rep.getControl().getStrandedBases(pt.expand(win), pt.getStrand());
						List<StrandedBaseCount> cReads = loadSignal ? 
								rep.getSignal().getStrandedBases(pt.expand(win), pt.getStrand()=='+' ? '-' : '+') :
									rep.getControl().getStrandedBases(pt.expand(win), pt.getStrand()=='+' ? '-' : '+');
						
						
						if(pt.getStrand()=='+'){
							for(StrandedBaseCount sbc : wReads){
								int sdist = sbc.getCoordinate()-pt.getLocation()+(win/2);
								if(sdist>=0 && sdist<win){
									watson[cond.getIndex()][sdist]+=sbc.getCount();
									perPointWatson[p][cond.getIndex()][sdist]+=sbc.getCount();
								}
							}
							for(StrandedBaseCount sbc : cReads){
								int sdist = sbc.getCoordinate()-pt.getLocation()+(win/2);
								if(sdist>=0 && sdist<win){
									crick[cond.getIndex()][sdist]+=sbc.getCount();
									perPointCrick[p][cond.getIndex()][sdist]+=sbc.getCount();
								}
							}
						}else{
							for(StrandedBaseCount sbc : wReads){
								int sdist = pt.getLocation()-sbc.getCoordinate()+(win/2);
								if(sdist>=0 && sdist<win){
									watson[cond.getIndex()][sdist]+=sbc.getCount();
									perPointWatson[p][cond.getIndex()][sdist]+=sbc.getCount();
								}
							}
							for(StrandedBaseCount sbc : cReads){
								int sdist = pt.getLocation()-sbc.getCoordinate()+(win/2);
								if(sdist>=0 && sdist<win){
									crick[cond.getIndex()][sdist]+=sbc.getCount();
									perPointCrick[p][cond.getIndex()][sdist]+=sbc.getCount();
								}
							}
						}
						
						p++;
					}
				}
			}
			//Normalize
			double wsum=0, csum=0;
			for(int w=0; w<win; w++){
				wsum+=watson[cond.getIndex()][w]; csum+=crick[cond.getIndex()][w];
			}for(int w=0; w<win; w++){
				watson[cond.getIndex()][w]/=wsum; crick[cond.getIndex()][w]/=csum;
			}
		}
	}
	
	//Accessors
	public int getWinSize(){return win;}
	public int getCenterOffset(){return centerOffset;}
	public double[] getCompositeWatson(ExperimentCondition c){return watson[c.getIndex()];}
	public double[] getCompositeCrick(ExperimentCondition c){return crick[c.getIndex()];}
	public double[] getPointWatson(StrandedPoint p, ExperimentCondition c){return perPointWatson[pointIndex.get(p)][c.getIndex()];}
	public double[] getPointCrick(StrandedPoint p, ExperimentCondition c){return perPointCrick[pointIndex.get(p)][c.getIndex()];}
	
	/**
	 * Per-condition sum of tags in composites
	 * @return
	 */
	public double[] getCompositeSums(){
		double[] sums = new double[numConditions];
		for(int c=0; c<numConditions; c++){
			for(int i=0; i<win; i++)
				sums[c]=0;
			for(int i=0; i<win; i++){
				sums[c]+=watson[c][i]+crick[c][i];
			}
		}
		return sums;
	}
	
	public String toString(ExperimentCondition cond){
		String out="";
		for(int w=0; w<win; w++){
			int pos = (w-centerOffset);
			out = out + pos+"\t"+watson[cond.getIndex()][w]+"\t"+crick[cond.getIndex()][w]+"\n";
		}
		return out;
	}
	
	//Print probs to a file
	public void printToFile(ExperimentCondition cond, String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			for(int w=0; w<win; w++){
				int pos = (w-centerOffset);
				fout.write(pos+"\t"+watson[cond.getIndex()][w]+"\t"+crick[cond.getIndex()][w]+"\n");
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	
	//Main method to make new composite distributions
	public static void main(String[] args){
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		ChExMixConfig config = new ChExMixConfig(gcon, args);
		if(config.helpWanted()){
			System.err.println("CompositeTagDistribution:");
			System.err.println("\t--points <stranded point file>");
			System.err.println("\t--win <window around points>");
			System.err.println(config.getArgsList());			
		}else{
			ExperimentManager manager = new ExperimentManager(econ);
			
			int w = Args.parseInteger(args, "win", 400);
			String pFile = Args.parseString(args, "points", null);
			List<StrandedPoint> pts = Utils.loadStrandedPointsFromFile(config.getGenome(), pFile);
			
			CompositeTagDistribution maker = new CompositeTagDistribution(pts, manager, w, true);
			manager.close();
		}
	}
	
}

