package edu.psu.compbio.seqcode.projects.seed.features;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;

/**
 * EnrichedFeature: a type of Feature that is associated with tag enrichment values
 * @author mahony
 */
public class EnrichedFeature extends Feature{
	/**
	 * ExperimentManager is required (to define enrichment values)
	 */
	public static ExperimentManager manager = null;
	public static boolean scoreIsAPValue=true;
	
	/**
	 * It may seem wasteful to save counts for all Samples here in contexts where we only care about per-condition values. 
	 * However, this way allows some flexibility in the future for features defined across multiple experiments.  
	 */
	protected float[] sampleCountsPos; //Positive counts, indexed by Sample id 
	protected float[] sampleCountsNeg; //Negative counts, indexed by Sample id
	protected float signalCount; //Total count for signal (usually defined by condition)
	protected float controlCount; //Total count for control (usually defined by condition)
	
	public EnrichedFeature(Region c, float[] perSamplePosCounts, float[] perSampleNegCounts, float sigCount, float ctrlCount, double score) throws Exception {
		super(c);
		if(manager==null){
			throw new Exception("EnrichedFeature: must instantiate experiment manager before using this class");
		}
		
		this.sampleCountsPos =perSamplePosCounts;
		this.sampleCountsNeg =perSampleNegCounts;
		this.signalCount = sigCount;
		this.controlCount = ctrlCount;
		this.score = score;
	}

	//Accessors
	public float[] getSampleCountsPos(){return sampleCountsPos;}
	public float[] getSampleCountsNeg(){return sampleCountsNeg;}
	public float getSignalCount(){return signalCount;}
	public float getControlCount(){return controlCount;}
	
	public String toString() {
		return new String(coords.toString()+"\t"+String.format("%.1f", signalCount)+"\t"+String.format("%.1f", controlCount)+"\t"+String.format("%.5e", score));
	}

	public String toGFF() {
		return new String(coords.getChrom()+"\tSEEDS\tfeature\t"+coords.getStart()+"\t"+coords.getEnd()+"\t.\t"+coords.getStrand()+"\t.\t"+"Note="+"Score:"+String.format("%.5e", score)+",Signal="+signalCount+",Control="+controlCount);
	}

	public String headString() {
		return new String("#Coordinates\tSignal\tControl\tp-value");
	}

	public String toSequence(SequenceGenerator seqgen, int extension) {
		return seqgen.execute(coords.expand(extension/2, extension/2));
	}

	//Settors
	public void setSampleCountsPos(float[] scp){sampleCountsPos = scp;}
	public void setSampleCountsNeg(float[] scn){sampleCountsNeg = scn;}
	public void setSignalCount(float c){signalCount=c;}
	public void setControlCount(float c){controlCount=c;}
	
	/**
     * Rank according to increasing score, then by decreasing signal
     */
  	public int compareTo(Feature p) {
  		if(score<p.score){return(-1);}
  		else if(score>p.score){return(1);}
  		else{
  			if(p instanceof EnrichedFeature){
  				EnrichedFeature ef = (EnrichedFeature)p;
  				if(this.getSignalCount()>ef.getSignalCount()){return(-1);}
  				else if(this.getSignalCount()<ef.getSignalCount()){return(1);}
  				else{return(0);}
  			}else
  				return(0);
  		}
  	}
}
