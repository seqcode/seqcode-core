package edu.psu.compbio.seqcode.projects.multigps.utilities;

import java.util.List;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.projects.multigps.framework.MultiGPSConfig;

public class BindingModelMaker {
	protected ExperimentManager experiments;
	protected MultiGPSConfig config;
	protected List<StrandedPoint> points;
	protected Integer win;
	
	
	public BindingModelMaker(List<StrandedPoint> p, ExperimentManager e, MultiGPSConfig c, int w){
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
								watson[sdist]+=sbc.getCount();
						}
						for(StrandedBaseCount sbc : cReads){
							int sdist = pt.getLocation()-sbc.getCoordinate()+(win/2);
							if(sdist>=0 && sdist<win)
								crick[sdist]+=sbc.getCount();
						}
					}else{
						for(StrandedBaseCount sbc : wReads){
							int sdist = pt.getLocation()-sbc.getCoordinate()+(win/2);
							if(sdist>=0 && sdist<win)
								watson[sdist]+=sbc.getCount();
						}
						for(StrandedBaseCount sbc : cReads){
							int sdist = sbc.getCoordinate()-pt.getLocation()+(win/2);
							if(sdist>=0 && sdist<win)
								crick[sdist]+=sbc.getCount();
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
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		MultiGPSConfig config = new MultiGPSConfig(gcon, args);
		if(config.helpWanted()){
			System.err.println("BindingModelMaker:");
			System.err.println("\t--points <stranded point file>");
			System.err.println("\t--win <window around points>");
			System.err.println(config.getArgsList());			
		}else{
			ExperimentManager manager = new ExperimentManager(econ);
			
			int w = Args.parseInteger(args, "win", 400);
			String pFile = Args.parseString(args, "points", null);
			List<StrandedPoint> pts = Utils.loadStrandedPointsFromFile(config.getGenome(), pFile);
			BindingModelMaker maker = new BindingModelMaker(pts, manager, config, w);
			maker.execute();
			
			manager.close();
		}
	}
}
