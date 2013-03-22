package edu.psu.compbio.seqcode.projects.multigps.utilities;

import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.StrandedPoint;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.framework.StrandedBaseCount;

public class BindingModelMaker {
	protected ExperimentSet experiments;
	protected Config config;
	protected List<StrandedPoint> points;
	protected Integer win;
	
	
	public BindingModelMaker(List<StrandedPoint> p, ExperimentSet e, Config c, int w){
		experiments = e;
		config = c;
		win = w;
		points = p;
	}
	
	public void execute(){
		double[] watson = new double[win];
		double[] crick = new double[win];
		
		for(ExperimentCondition c : experiments.getConditions()){
			for(ControlledExperiment r : c.getReplicates()){
				//Reset
				for(int w=0; w<win; w++){
					watson[w]=0; crick[w]=0;
				}
				
				//Iterate through points
				for(StrandedPoint pt : points){
					//Load reads
					List<StrandedBaseCount> wReads = r.getSignal().getStrandedBases(pt.expand(win), pt.getStrand());
					List<StrandedBaseCount> cReads = r.getSignal().getStrandedBases(pt.expand(win), pt.getStrand()=='+' ? '-' : '+');
					
					if(pt.getStrand()=='+'){
						for(StrandedBaseCount sbc : wReads){
							int sdist = sbc.getCoordinate()-pt.getLocation()+(win/2);
							if(sdist>=0 && sdist<win)
								watson[sdist]++;
						}
						for(StrandedBaseCount sbc : cReads){
							int sdist = pt.getLocation()-sbc.getCoordinate()+(win/2);
							if(sdist>=0 && sdist<win)
								crick[sdist]++;
						}
					}else{
						for(StrandedBaseCount sbc : wReads){
							int sdist = pt.getLocation()-sbc.getCoordinate()+(win/2);
							if(sdist>=0 && sdist<win)
								watson[sdist]++;
						}
						for(StrandedBaseCount sbc : cReads){
							int sdist = sbc.getCoordinate()-pt.getLocation()+(win/2);
							if(sdist>=0 && sdist<win)
								crick[sdist]++;
						}
					}
				}
				
				//Print
				for(int w=0; w<win; w++){
					System.out.println(w-(win/2)+"\t"+watson[w]+"\t"+crick[w]);
				}
			}
		}
	}
	
	//Main
	public static void main(String[] args){
		Config config = new Config(args); 
		if(config.helpWanted()){
			System.err.println("BindingModelMaker:");
			System.err.println("\t--points <stranded point file>");
			System.err.println("\t--win <window around points>");
			System.err.println(config.getArgsList());			
		}else{
			ExperimentManager manager = new ExperimentManager(config);
			ExperimentSet eset = manager.getExperimentSet();
			
			int w = Args.parseInteger(args, "win", 400);
			String pFile = Args.parseString(args, "points", null);
			List<StrandedPoint> pts = Utils.loadStrandedPointsFromFile(config.getGenome(), pFile);
			BindingModelMaker maker = new BindingModelMaker(pts, eset, config, w);
			maker.execute();
			
			manager.close();
		}
	}
}
