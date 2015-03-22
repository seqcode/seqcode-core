package edu.psu.compbio.seqcode.projects.seed.features;

import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;

public class EnrichedPCSPeakFeature extends EnrichedPeakFeature{

	protected float[][] sampleCountsByBase; //Counts in enriched region, indexed by Sample id, then indexed by ACGT 
	
	public EnrichedPCSPeakFeature(Region c, Point p, float[] perSamplePosCounts,
			float[] perSampleNegCounts, float sigCount, float ctrlCount,
			double score, float peakScore) throws Exception {
		super(c,p, perSamplePosCounts, perSampleNegCounts, sigCount, ctrlCount, score, peakScore);
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

}
