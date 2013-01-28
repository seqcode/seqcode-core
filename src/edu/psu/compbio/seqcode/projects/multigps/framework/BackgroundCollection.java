package edu.psu.compbio.seqcode.projects.multigps.framework;

import java.util.ArrayList;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
/**
 * BackgroundCollection: A container to hold diverse types of BackgroundModel and handle the various ways we may use them
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class BackgroundCollection {
	
	protected ArrayList<BackgroundModel> models = new ArrayList<BackgroundModel>();
	protected double highestHigh=0, lowestLow=0;
	
	public BackgroundCollection(){}
	
	/**
	 * Get the threshold for the genomic model
	 * @param tType
	 * @return
	 */
	public int getGenomicModelThreshold(){
		for(BackgroundModel m : models)
			if(m.isGenomeWide())
				return m.getThreshold();
		return -1;
	}
	/**
	 * Get the maximum threshold
	 * @param str strand
	 * @return
	 */
	public int getMaxThreshold(char str){
		int max=0;
		for(BackgroundModel m : models)
			if((m.strand==str ||str=='.' ) && m.getThreshold()>max)
				max = m.getThreshold();
		return max;
	}
	
	public void addBackgroundModel(BackgroundModel m){models.add(m);}
	
	/**
	 * Refresh the constituent background models using the local context
	 * @param currReg
	 * @param currOffset
	 * @param thisExptHitCounts
	 * @param otherExptHitCounts
	 */
	public void updateModels(Region currReg, int currOffset, double [] thisExptHitCounts, double [] otherExptHitCounts, float hitCountBin){
		for(BackgroundModel m : models)
			m.updateModel(currReg, currOffset, thisExptHitCounts, otherExptHitCounts, hitCountBin);		
	}
	
	/**
	 * Print thresholds
	 */
	public void printThresholds(){
		for(BackgroundModel m : models){
			System.out.println(m.modelType+"\t"+m.getThreshold()+"\t"+m.getStrand());
		}
	}
	
	/**
	 * Does the count pass all thresholds?
	 * @param count
	 * @return
	 */
	public boolean passesAllThresholds(int count){return passesAllThresholds(count, '.');}
	public boolean passesAllThresholds(int count, char strand){
		boolean pass=true;
		for(BackgroundModel m : models){
			if(m.getStrand()==strand)
				if(!m.passesThreshold(count))
					pass=false;
		}
		return pass;
	}
	public boolean passesGenomicThreshold(int count){return passesGenomicThreshold(count, '.');}
	public boolean passesGenomicThreshold(int count, char strand){
		boolean pass=true;
		for(BackgroundModel m : models){
			if(m.isGenomeWide() && m.getStrand()==strand)
				if(!m.passesThreshold(count))
					pass=false;
		}
		return pass;
	}
	
	/**
	 * Is the count under all thresholds?
	 * @param count
	 * @return
	 */
	public boolean underAllThresholds(int count){return underAllThresholds(count, '.');}
	public boolean underAllThresholds(int count, char strand){
		boolean pass=true;
		for(BackgroundModel m : models){
			if(m.getStrand()==strand)
				if(!m.underThreshold(count))
					pass=false;
		}
		return pass;
	}
	public boolean underGenomicThreshold(int count){return underGenomicThreshold(count, '.');}
	public boolean underGenomicThreshold(int count, char strand){
		boolean pass=true;
		for(BackgroundModel m : models){
			if(m.isGenomeWide() && m.getStrand()==strand)
				if(!m.underThreshold(count))
					pass=false;
		}
		return pass;
	}
	
}
