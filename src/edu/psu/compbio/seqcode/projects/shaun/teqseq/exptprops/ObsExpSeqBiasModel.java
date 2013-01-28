package edu.psu.compbio.seqcode.projects.shaun.teqseq.exptprops;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.AlignHit;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.GenomeLoader;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.ReadLoader;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.ReadUtils;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.RegionHitsExtractor;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.AGene;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.AnnotationFilter;

/**
 * ObsExpSeqBiasModel: Estimates sequencing bias using a custom procedure.
 * In this approach, observed and expected counts are only calculated over non-overlapping, non-UTR, known exons.
 * If the exon of length L contains R reads, the assumption of no sequencing bias would suggest that each of the 
 * L-k+1 k-mers in the exon should have R/(L-k+1) reads starting at that position. 
 * The observed k-mer counts are therefore estimated from the first x bases in the mapped 5' read position.
 * The expected k-mer counts are calculated from the exon length and read counts as above. 
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 * @see SeqBiasModel
 */
public class ObsExpSeqBiasModel extends SeqBiasModel{

	protected Map<String, List<Region>> regionsOfInterest;
	
	public ObsExpSeqBiasModel(GenomeLoader g, int k, Map<String, List<Region>> regionsOfInterest){
		super(g, k);
		this.regionsOfInterest = regionsOfInterest;
		fileExt = "."+kmerLength+"mer.obsexp.seqbias";
	}

	@Override
	public void execute(ReadLoader reads, Collection<AGene> genes) {
		Map<String, List<Region>> testExons = AnnotationFilter.getNonOverlappingMiddleExons(genes);
		RegionHitsExtractor extractor = new RegionHitsExtractor(genLoader, reads, testExons, regionsOfInterest);
		
		int numExons=0;
		int readCount=0;
		int exonSumLen=0;
		int regionsFound=0, regionsNotFound=0;
		while(extractor.hasNextRegion()){
			//Get the next Exon & AlignHit list pair
			Pair<Region, List<AlignHit>> p = extractor.getNext();
			if(p!=null){
				regionsFound++;
				Region exon = p.car();
				List<AlignHit> hits = p.cdr();
				//Get the corresponding sequence
				String seq = genLoader.getSequence(exon);
				
				//Count the read starts
				double[] posStarts = ReadUtils.makeReadStartsArray(hits, exon, '+');
				double[] negStarts = ReadUtils.makeReadStartsArray(hits, exon, '-');
				
				//Count pos & neg strand hits
				double posHitCount=0, negHitCount=0;
				for(int z=0; z<posStarts.length; z++){ posHitCount+=posStarts[z];}
				for(int z=0; z<negStarts.length; z++){ negHitCount+=negStarts[z];}
				
				//Calculate the expected count
				double posExpected = posHitCount/(double)exon.getWidth();
				double negExpected = negHitCount/(double)exon.getWidth();
				
				//Increment Obs & Exp models
				for(int x=kmerLength-1; x<=exon.getWidth()-kmerLength; x++){
					//Add this k-mer even if the hit count is zero
					String posKmer = seq.substring(x, x+kmerLength);
					int posIndex = seq2int(posKmer);
					if(posIndex!=-1){
						kmerObsCounts[posIndex]+=posStarts[x];
						contributingHits+=posStarts[x];
						kmerExpCounts[posIndex]+=posExpected;
					}
					//System.out.println(exon.getLocationString()+"\t"+x+"+\t"+posKmer+"\t"+posIndex+"\t"+posStarts[x]+"\t"+posExpected);
					
					String negKmer = SequenceUtils.reverseComplement(seq.substring(x-kmerLength+1, x+1));
					int negIndex = seq2int(negKmer);
					if(negIndex!=-1){
						kmerObsCounts[negIndex]+=negStarts[x];
						contributingHits+=negStarts[x];
						kmerExpCounts[negIndex]+=negExpected;
					}
					//System.out.println(exon.getLocationString()+"\t"+x+"-\t"+negKmer+"\t"+negIndex+"\t"+negStarts[x]+"\t"+negExpected);					
				}
				numExons++;
				exonSumLen+= exon.getWidth();
				readCount += hits.size();
			}else{
				regionsNotFound++;
			}
		}
		System.err.println("# Exons: "+numExons+", bp: "+exonSumLen+", Reads: "+readCount);
		System.err.println("# Exons found: "+regionsFound+", Exons not found: "+regionsNotFound);
		
		//Normalize & weight
		double sumStart=0, sumEnd=0;
		for(int i=0; i<numKmers; i++){
			sumStart += kmerObsCounts[i];
			sumEnd += kmerExpCounts[i];
		}
		for(int i=0; i<numKmers; i++){
			kmerObsProbs[i] = kmerObsCounts[i]/sumStart;
			kmerExpProbs[i] = kmerExpCounts[i]/sumEnd;
			kmerWeights[i] = kmerExpProbs[i]/kmerObsProbs[i];
		}
		
		//Save the model
		String outName = reads.getSourcePath()+fileExt;
		save(outName);
		System.err.println("Saved SeqBiasModel to "+outName);
	}

}

