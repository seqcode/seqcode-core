package edu.psu.compbio.seqcode.projects.seed.features;

import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;

/**
 * Feature: Any genomic feature that we wish to discover in a deep-sequencing experiment. 
 * 
 * @author mahony
 */
public abstract class Feature implements Comparable<Feature>{
	protected static int regCount=0;
	public static boolean scoreIsAPValue=true;
	
	protected int ID; //Unique ID for this feature
	protected Region coords; //Full extent of genomic feature
	protected double score; //Score for this region
	
    public Feature(Region c){
		ID = regCount; regCount++;
		coords=c;
	}
	/**
	 * Returns the score for this feature
	 * @return double
	 */
    public double getScore(){return score;}
    
    /**
     * Update the score for this feature
     * @param s
     */
    public void setScore(double s){score =s;}
    
	/**
	 * Returns coordinates of this Feature
	 */
	public Region getCoords(){return coords;}
	
	/**
	 * Update the feature coordinates
	 * @param r
	 */
	public void setCoords(Region r){coords=r;}
	
	/**
	 * Returns String describing the genomic feature (custom format)
	 */
    public abstract String toString();
    /**
     * Returns GFF-formatted String describing the feature 
     */
    public abstract String toGFF();
    /**
     * Returns a line describing each field returned by the toString method
     */
    public abstract String headString();
    /**
     * Returns the genomic sequence at/around feature
     * @param extension : extends window by half of extension either site of window (0 provides sequence from feature coords only)
     * @return
     */
    public abstract String toSequence(SequenceGenerator seqgen, int extension);
    
    /**
     * Rank according to increasing score
     */
  	public int compareTo(Feature p) {
  		if(score<p.score){return(-1);}
  		else if(score>p.score){return(1);}
  		else{return(0);}
  	}
}
