package org.seqcode.projects.seed.features;

import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;

/**
 * PeakFeature: a punctate feature within a broader enriched region
 * @author mahony
 *
 */
public class EnrichedPeakFeature extends EnrichedFeature{

	protected Point peak;
	protected double peakScore;
	
	public EnrichedPeakFeature(Region c, Point p, float[] perSamplePosCounts,
			float[] perSampleNegCounts, float sigCount, float ctrlCount,
			double score, float peakScore) throws Exception {
		super(c, perSamplePosCounts, perSampleNegCounts, sigCount, ctrlCount, score);
		
		this.peak=p;
		this.peakScore = peakScore;
	}

	public String toString() {
		return new String(peak.toString()+"\t"+coords.toString()+"\t"+String.format("%.1f", signalCount)+"\t"+String.format("%.1f", controlCount)+"\t"+String.format("%.5e", score)+"\t"+String.format("%.3f", peakScore));
	}

	public String toGFF() {
		return new String(peak.getChrom()+"\tSEEDS\tfeature\t"+peak.getLocation()+"\t"+peak.getLocation()+1+"\t.\t"+coords.getStrand()+"\t.\t"+"Note="+"Score:"+String.format("%.5e", score)+",PeakScore:"+String.format("%.3f", peakScore)+",Signal="+signalCount+",Control="+controlCount);
	}

	public String headString() {
		return new String("Peak\tDomainRegion\tSignal\tControl\tDomainScore\tPeakScore");
	}

	public String toSequence(SequenceGenerator seqgen, int extension) {
		return seqgen.execute(peak.expand(extension));
	}
	
	public double getPeakScore(){return peakScore;}
	public void setPeakScore(double p){peakScore = p;}

}
