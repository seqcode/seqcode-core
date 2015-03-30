package edu.psu.compbio.seqcode.projects.seed.features;

import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.genome.location.Region;

public class SuperEnrichedFeature extends EnrichedFeature {
	
	protected List<EnrichedFeature> typicalEnrichedFeatures;
	
	
	// Is the super enriched feature a super enhancer
	public boolean isSuperEnhancer = false;
	
	// Slope of the super enhaner feature (should be calculated on a sorted list of SFE's (Sorting based on SFE read counts))
	public double slope = -1.0;
	

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
	
	
	// Get methods
	public EnrichedFeature getFirstTEF(){return typicalEnrichedFeatures.get(0);}
	public EnrichedFeature getLastTEF(){return typicalEnrichedFeatures.get(typicalEnrichedFeatures.size()-1);}
	public List<EnrichedFeature> getTEFs(){return typicalEnrichedFeatures;}
	public boolean isSEnh(){return isSuperEnhancer;}
	
	
	
	
	//Set methods
	public void setSuperEnahncer(boolean isSE){isSuperEnhancer = isSE;}
	public void setSlope(double s){slope =s;}
	
	public  String toString() {
		return new String(coords.toString()+"\t"+String.format("%.1f", signalCount)+"\t"+String.format("%.1f", controlCount)+"\t"+String.format("%.5e", score)+"\t"+String.format("%.1f", slope)+"\t"+Integer.toString(typicalEnrichedFeatures.size()));
	}
	
	public String toGFF() {
		String TEs ="";
		for(EnrichedFeature ef : typicalEnrichedFeatures){
			TEs.concat(ef.getCoords().getLocationString()+",");
		}
		return new String(coords.getChrom()+"\tSEEDS\tSuperfeature\t"+coords.getStart()+"\t"+coords.getEnd()+"\t.\t"+coords.getStrand()+"\t.\t"+
	"Note="+"Score:"+String.format("%.5e", score)+",Signal="+signalCount+",Control="+controlCount+"; StichedEnhancersCoords="+TEs+"; Slope="+String.format("%.1f", slope));
	}
	
	public String headString() {
		return new String("#Coordinates\tSignal\tControl\tscore\tSlope\tNoTEs");
	}
	
	

}
