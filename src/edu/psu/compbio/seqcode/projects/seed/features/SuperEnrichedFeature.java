package edu.psu.compbio.seqcode.projects.seed.features;

import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.genome.location.Region;

public class SuperEnrichedFeature extends EnrichedFeature {
	
	protected List<EnrichedFeature> typicalEnrichedFeatures;
	
	// Overriding value set in inherited feature 
	public static boolean scoreIsAPValue=false;
	

	public SuperEnrichedFeature(EnrichedFeature ef) throws Exception {
		super(ef.coords, ef.sampleCountsPos, ef.getSampleCountsNeg(), ef.getSignalCount(), ef.getControlCount(), ef.getScore());
		// TODO Auto-generated constructor stub
		typicalEnrichedFeatures = new ArrayList<EnrichedFeature>();
		typicalEnrichedFeatures.add(ef);
		
	}
	
	
	public void addTypicalFeature(EnrichedFeature tf){
		typicalEnrichedFeatures.add(tf);
		setCoords(new Region(coords,tf.getCoords()));
	}
	
	public EnrichedFeature getFirstTEF(){return typicalEnrichedFeatures.get(0);}
	
	public EnrichedFeature getLastTEF(){return typicalEnrichedFeatures.get(typicalEnrichedFeatures.size()-1);}
	
	public List<EnrichedFeature> getTEFs(){return typicalEnrichedFeatures;}
	
	

}
