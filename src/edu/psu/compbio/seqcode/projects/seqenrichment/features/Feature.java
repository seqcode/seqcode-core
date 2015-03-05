package edu.psu.compbio.seqcode.projects.seqenrichment.features;

import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;

/**
 * Feature: Any genomic feature that we wish to discover in a deep-sequencing experiment. 
 * 
 * @author mahony
 */
public abstract class Feature{
	protected static int regCount=0;
	
	protected int ID; //Unique ID for this feature
	protected Region coords; //Full extent of genomic feature 
	
    public Feature(Region c){
		ID = regCount; regCount++;
		coords=c;
	}
	
	/**
	 * Returns coordinates of this Feature
	 */
	public Region getCoords(){return coords;}
	
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
}
