package edu.psu.compbio.seqcode.projects.seqenrichment.features;

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
	protected static ExperimentManager manager = null;
	
	/**
	 * It may seem wasteful to save counts for all Samples here in contexts where we only care about per-condition values. 
	 * However, this way allows some flexibility in the future for features defined across multiple experiments.  
	 */
	protected float[] sampleCountsPos; //Positive counts, indexed by Sample id 
	protected float[] sampleCountsNeg; //Negative counts, indexed by Sample id
	protected float signalCount; //Total count for signal (usually defined by condition)
	protected float controlCount; //Total count for control (usually defined by condition)
	protected float score; //Score is defined implementation-specifically. 
	
	public EnrichedFeature(Region c, float[] perSamplePosCounts, float[] perSampleNegCounts, float sigCount, float ctrlCount, float score) throws Exception {
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

	
	public String toString() {
		return new String(coords.toString()+"\t"+String.format("%.1f", signalCount)+"\t"+String.format("%.1f", controlCount)+"\t"+String.format("%.5e", score));
	}

	public String toGFF() {
		return new String(coords.getChrom()+"\tSEEDS\tfeature\t"+coords.getStart()+"\t"+coords.getEnd()+"\t.\t"+coords.getStrand()+"\t.\t"+"Note="+"Score:"+String.format("%.5e", score)+",Signal="+signalCount+",Control="+controlCount);
	}

	public String headString() {
		return new String("Coordinates\tSignal\tControl\tScore");
	}

	public String toSequence(SequenceGenerator seqgen, int extension) {
		return seqgen.execute(coords.expand(extension/2, extension/2));
	}

}
