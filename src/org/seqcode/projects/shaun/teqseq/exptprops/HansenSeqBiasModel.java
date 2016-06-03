package org.seqcode.projects.shaun.teqseq.exptprops;

import java.util.Collection;

import org.seqcode.genome.Genome;
import org.seqcode.projects.shaun.teqseq.core.GenomeLoader;
import org.seqcode.projects.shaun.teqseq.core.ReadLoader;
import org.seqcode.projects.shaun.teqseq.geneloaders.AGene;


/**
 * HansenSeqBiasModel: Estimates sequencing bias using the procedure described by Hansen, Brenner, & Dudoit, NAR, 2010.
 * In this approach, the observed k-mer counts are averaged over the first x positions in the reads (x=2 by default),
 * and the expected k-mer counts are averaged over the last y positions (y=6 by default). 
 * One problem with this approach is that it is too susceptible to bias introduced by expression. In other words, 
 * the sequencing bias estimated using this approach is not constant over different experimental conditions, given 
 * the same preparation methods.  
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 * @see SeqBiasModel
 */
public class HansenSeqBiasModel extends SeqBiasModel{

	protected int fivePrimeSpread=2;
	protected int threePrimeSpread=6;
	
	public HansenSeqBiasModel(GenomeLoader g, int k){
		super(g, k);
		fileExt = "."+kmerLength+"mer.hansen.seqbias";
	}
	
	@Override
	public void execute(ReadLoader reads, Collection<AGene> genes){
		reads.initializeReadIterator();
		int [] startIndex = new int[fivePrimeSpread];
		int [] endIndex = new int[threePrimeSpread];
		
		//Count all appropriate read k-mers
		while(reads.nextRead()){
			byte[] currSeq = reads.getCurrReadSeq();
					//System.out.println(bytes2String(currSeq));
			//Start
			int numStarts=0;
			for(int e =0; e<fivePrimeSpread; e++){
				for(int c=0; c<kmerLength; c++)
					kmer[c] = currSeq[c];
				startIndex[e] = seq2int(kmer);
				if(startIndex[e]!=-1)
					numStarts++;
			}
					//System.out.println("s\t"+bytes2String(kmer));
			//Ends 
			int numEnds = 0;
			for(int e =0; e<threePrimeSpread; e++){
				for(int c=0; c<kmerLength; c++)
					kmer[c] = currSeq[currSeq.length-e-kmerLength+c];
				endIndex[e] = seq2int(kmer);
				if(endIndex[e]!=-1)
					numEnds++;
					//System.out.println("e"+e+"\t"+bytes2String(kmer));
			}
			
			//Add the read here
			if(numStarts>0 && numEnds>0){
				contributingHits++;
				for(int e =0; e<fivePrimeSpread; e++){
					if(startIndex[e]!=-1)
						kmerObsCounts[startIndex[e]]+=reads.getCurrHit().getWeight();
				}
				for(int e =0; e<threePrimeSpread; e++){
					if(endIndex[e]!=-1)
						kmerExpCounts[endIndex[e]]+=reads.getCurrHit().getWeight();
				}
			}
		}
		
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
