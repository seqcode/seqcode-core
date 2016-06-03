package org.seqcode.projects.shaun.dataexplorer;

import java.util.Collection;

import org.seqcode.genome.location.Region;


public class DataSourceCombiner extends DataSource{

	private DataCollection data;
	private Collection<DataSource> sources;
	private String combinationType="SUM";
	
	public DataSourceCombiner(String combinationType, DataCollection dat, Collection<DataSource> srcs,String name, double threshold, double weight){
		super(name, threshold, weight);
		data=dat;
		sources=srcs;
		this.combinationType=combinationType;
	}
	
	public double genData(Region r){
		double score=0;
		for(DataSource d : sources){
			if(combinationType.equals("SUM"))
				score+=data.getDataVal(r, d)[0]*d.getWeight();
			else if(combinationType.equals("MAX")){
				score = data.getDataVal(r, d)[0]>score ? data.getDataVal(r, d)[0] : score;
			}
		}
		return(score);
	}

	public double[] genDataVector(Region r, int binSize) {
		int numBins = (r.getWidth()/binSize)+1;
		double [] vec = new double[numBins];
		for(int i=0; i<numBins; i++){
			double score=0;
			for(DataSource d : sources){
				if(combinationType.equals("SUM"))
					score+=data.getDataVal(r, d)[i]*d.getWeight();
				else if(combinationType.equals("MAX")){
					score = data.getDataVal(r, d)[i]>score ? data.getDataVal(r, d)[0] : score;
				}
			}
		}
		return vec;
	}
	public void cleanup(){}
}
