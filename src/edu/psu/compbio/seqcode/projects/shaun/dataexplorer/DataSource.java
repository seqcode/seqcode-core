package edu.psu.compbio.seqcode.projects.shaun.dataexplorer;

import edu.psu.compbio.seqcode.genome.location.Region;

public abstract class DataSource {

		protected String name;
		protected double threshold;
		protected double weight;
		
		public DataSource(String name, double threshold, double weight){this.threshold=threshold; this.name=name;; this.weight=weight;}
				
		public double getThreshold(){return threshold;}
		public double getWeight(){return weight;}
		public String getName(){return name;}
		public abstract double genData(Region r);
		public abstract double[] genDataVector(Region r, int binSize);	
		public abstract void cleanup();
		
		protected double[] reverseVec(double[] vec){
			double[] newVec = new double[vec.length];
			for(int i=0; i<newVec.length; i++)
				newVec[i] = vec[vec.length-i-1];
			return newVec;
		}
}
