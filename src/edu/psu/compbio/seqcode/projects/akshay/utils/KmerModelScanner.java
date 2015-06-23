package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.util.List;

import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;

public class KmerModelScanner {
	
	
	protected double[][] kmerwights;
	protected List<Point> peaks;
	//protected SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
	//SequenceGenerator seqgen = new SequenceGenerator();
	
	public KmerModelScanner() {
		
	}
	
	
	public int compareKmerModels(int[][] a, int[][] b) throws IncorrectComparision{
		int score=0;
		if(a.length != b.length){
			throw new IncorrectComparision("Unequal model lenghts");
		}else{
			for(int r=0; r<a.length; r++){
				for(int c=0; c<=a.length; c++){
					if(a[r][c] != b[r][c])
					score = score++;
				}
			}
		}		
		return score;
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	public class IncorrectComparision extends Exception{
		IncorrectComparision(String s){
			super(s);
		}
	}
	
	
}
