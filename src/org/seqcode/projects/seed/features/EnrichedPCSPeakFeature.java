package org.seqcode.projects.seed.features;

import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;

public class EnrichedPCSPeakFeature extends EnrichedPeakFeature{

	protected float[] tagCountsByBase; //Per condition tag counts in enriched region, indexed by ACGT 
	protected float[] bubbleSequenceCounts; //Counts of the 4 bases in the bubble region, indexed by ACGT
	
	public EnrichedPCSPeakFeature(Region c, Point p, float[] perSamplePosCounts,
			float[] perSampleNegCounts, float sigCount, float ctrlCount,
			double score, float peakScore,
			float[] bubbleTags, float[] bubbleBases) throws Exception {
		super(c,p, perSamplePosCounts, perSampleNegCounts, sigCount, ctrlCount, score, peakScore);
		tagCountsByBase = bubbleTags;
		bubbleSequenceCounts = bubbleBases;
	}

	public String toString() {
		return new String(peak.toString()+"\t"+coords.toString()+"\t"+String.format("%.1f", signalCount)+"\t"+String.format("%.1f", controlCount)+"\t"+String.format("%.5e", score)+"\t"+String.format("%.3f", peakScore)+"\t"+String.format("%.5f", bubbleIndex()));
	}

	public String toGFF() {
		return new String(peak.getChrom()+"\tSEEDS\tfeature\t"+peak.getLocation()+"\t"+peak.getLocation()+1+"\t.\t"+coords.getStrand()+"\t.\t"+"Note="+"Score:"+String.format("%.5e", score)+",PeakScore:"+String.format("%.3f", peakScore)+",Signal="+signalCount+",Control="+controlCount+",BubbleIndex="+bubbleIndex());
	}

	public String headString() {
		return new String("Peak\tDomainRegion\tSignal\tControl\tDomainScore\tPeakScore\tBubbleIndex");
	}

	public void setTagCountsByBase(float[] counts){
		for(int b=0; b<4; b++)
			tagCountsByBase[b]=counts[b];
	}
	
	public void setBubbleSequenceCounts(float[] counts){
		for(int b=0; b<4; b++)
			bubbleSequenceCounts[b]=counts[b];
	}
	
	/**
	 * Bubble index = expected T / expected A+C+G
	 * @return
	 */
	public float bubbleIndex(){
		float bI = -1;
		float expectedT=0, expectedACG=0;
		if(bubbleSequenceCounts[3]>0){
			expectedT = tagCountsByBase[3] / bubbleSequenceCounts[3];
			float baseACG = bubbleSequenceCounts[0]+bubbleSequenceCounts[1]+bubbleSequenceCounts[2];
			if(baseACG>0)
				expectedACG = (tagCountsByBase[0]+tagCountsByBase[1]+tagCountsByBase[2])/baseACG;
			if(expectedACG>0)
				bI = expectedT/expectedACG;
		}
		return bI;
	}
}
