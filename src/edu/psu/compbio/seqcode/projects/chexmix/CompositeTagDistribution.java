package edu.psu.compbio.seqcode.projects.chexmix;

import java.util.List;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;

public class CompositeTagDistribution {
	protected ExperimentCondition experiment;	
	protected List<StrandedPoint> points;
	protected Integer win;
	protected double[] watson;
	protected double[] crick;
	
	
	public CompositeTagDistribution(List<StrandedPoint> points, ExperimentCondition expt, int win){
		experiment = expt;
		this.win = win;
		this.points = points;
	
		watson = new double[win];
		crick = new double[win];
		
		for(ControlledExperiment r : experiment.getReplicates()){
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

