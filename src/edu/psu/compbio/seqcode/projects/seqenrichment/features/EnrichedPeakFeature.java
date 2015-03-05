package edu.psu.compbio.seqcode.projects.seqenrichment.features;

import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;

/**
 * PeakFeature: a punctate feature within a broader enriched region
 * @author mahony
 *
 */
public class EnrichedPeakFeature extends EnrichedFeature{

	protected Point peak;
	
	public EnrichedPeakFeature(Region c, Point p, float[] perSamplePosCounts,
			float[] perSampleNegCounts, float sigCount, float ctrlCount,
			float score) throws Exception {
		super(c, perSamplePosCounts, perSampleNegCounts, sigCount, ctrlCount, score);
		
		peak=p;
	}

	public String toString() {
		return new String(peak.toString()+"\t"+coords.toString()+"\t"+String.format("%.1f", signalCount)+"\t"+String.format("%.1f", controlCount)+"\t"+String.format("%.5e", score));
	}

	public String toGFF() {
		return new String(peak.getChrom()+"\tSEEDS\tfeature\t"+peak.getLocation()+"\t"+peak.getLocation()+1+"\t.\t"+coords.getStrand()+"\t.\t"+"Note="+"Score:"+String.format("%.5e", score)+",Signal="+signalCount+",Control="+controlCount);
	}

	public String headString() {
		return new String("Peak\tRegion\tSignal\tControl\tScore");
	}

	public String toSequence(SequenceGenerator seqgen, int extension) {
		return seqgen.execute(peak.expand(extension));
	}

}
