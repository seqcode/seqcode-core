package edu.psu.compbio.seqcode.projects.shaun.teqseq.core;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.AGene;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.AnnotationFilter;

/**
 * GenomeSegmenter: segment a genome into manageable regions using either a set of genes or read coverage directly.
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class GenomeSegmenter {

	protected GenomeLoader genLoader;
	protected Integer minSubChrSize;
	protected Integer minSubChrSpace;
	protected Map<String, List<Region>> segments;
	protected int maxHitSpan = 50000; //A maximum just in case there are some long distance mapping artifacts
	protected int emptyCoverageLevel=0; //Maximum coverage in an 'empty' region
	
	public GenomeSegmenter(GenomeLoader gLoad, Integer minSubChrSize, Integer minSubChrSpace) {
		genLoader = gLoad;
		this.minSubChrSize = minSubChrSize;
		this.minSubChrSpace = minSubChrSpace;
		segments = new HashMap<String, List<Region>>();
	}

	
	/**
	 * Segment the genome using all reads to be analyzed. 
	 *  i.e. paint all reads and read pairs onto the genome and define the segments by searching for empty spaces.
	 * @param l
	 * @return
	 */
	public Map<String, List<Region>> segmentWithReads(Collection<ArrayList<ReadLoader>> l) {
		List<ReadLoader> loaders = new ArrayList<ReadLoader>();
		for(List<ReadLoader> rl : l)
			loaders.addAll(rl);
		return segmentWithReads(loaders);
	}
	
	/**
	 * Segment the genome using all reads to be analyzed. 
	 *  i.e. paint all reads and read pairs onto the genome and define the segments by searching for empty spaces. 
	 * @param loaders
	 * @return
	 */
	public Map<String, List<Region>> segmentWithReads(List<ReadLoader> loaders) {
		Genome gen = genLoader.getGenome();
		int numSeg=0;
		for(String c : gen.getChromList()){
			Region chrom = new Region(gen, c, 1, gen.getChromLength(c));
			System.err.println("GenomeSegmenter: Processing chromosome "+chrom.getChrom());
			segments.put(chrom.getChrom(), new ArrayList<Region>());

			//Initialize painting array
			int[] coverageArray = new int[chrom.getWidth()+1];
			for(int i=0; i<coverageArray.length; i++){ coverageArray[i]=0;}
			
			//Iterate through loaders, getting reads for the current chromosome
			for(ReadLoader rl : loaders){
				List<AlignHit> currHits = rl.getOverlappingHits(chrom);
				
				//Paint coverage onto the array. 
				//Remember that coverage here refers to all bases between paired reads, not just the actual reads themselves.
				for(AlignHit h : currHits){
					//Calculate spanning region of the read(s)
					Region hitSpan=null;
					if(h.getPairedRead()!=null && h.getChrom().equals(h.getPairedRead().getChrom())){//If this is a paired read
						//Only count each pair once
						if(h.getFivePrime()<h.getPairedRead().getLocation() && h.getPairedRead().getLocation()-h.getFivePrime()<maxHitSpan)
							hitSpan = new Region(gen, chrom.getChrom(), h.getFivePrime(), h.getPairedRead().getLocation());
					}else{//Otherwise span across all align blocks						
						hitSpan = h;
					}
					if(hitSpan!=null){
					    //Paint coverage
					    for(int j=hitSpan.getStart(); j<=hitSpan.getEnd(); j++){
						coverageArray[j]++;
					    }
					}
				}
			}
			//Segment in coverage array
			int currStart=0; int currEmpty=0;
			for(int x=0; x<coverageArray.length; x++){
				if(coverageArray[x]<=emptyCoverageLevel)
					currEmpty++;
				else
					currEmpty=0;
					
				if((x-currStart)>minSubChrSize && currEmpty>minSubChrSpace){
					segments.get(chrom.getChrom()).add(new Region(gen, chrom.getChrom(), currStart, x));
					currEmpty=0;
					currStart = x+1;
				}
			}
			//Last element
			segments.get(chrom.getChrom()).add(new Region(gen, chrom.getChrom(), currStart, chrom.getEnd()));
			numSeg+=segments.get(chrom.getChrom()).size();
		}
		System.err.println("GenomeSegmenter: Genome segmented into "+numSeg+" segments.");
		return(segments);
	}	
	
	
	/**
	 * Segment the genome using gene annotations.
	 *  i.e. find spaces between genes. 
	 * @param genes
	 * @return
	 */
	public Map<String, List<Region>> segmentWithGenes(Collection<AGene> genes){
		Map<String, List<Region>> geneRegions = AnnotationFilter.getAllGeneRegions(genes);
		
		int numSeg=0;
		for(String chr : geneRegions.keySet()){
			segments.put(chr, new ArrayList<Region>());
			//int currStart = cRegs.get(chr).get(0).getStart();
			int currStart = 0;
			for(int r=0; r<geneRegions.get(chr).size()-1; r++){
				int potEnd = geneRegions.get(chr).get(r).getEnd()+((geneRegions.get(chr).get(r+1).getStart()-geneRegions.get(chr).get(r).getEnd())/2);
				if(geneRegions.get(chr).get(r+1).getStart()-geneRegions.get(chr).get(r).getEnd() > minSubChrSpace &&
						potEnd-currStart > minSubChrSize){
					segments.get(chr).add(new Region(genLoader.getGenome(), chr, currStart, potEnd));
					currStart = potEnd + 1;
				}
			}
			//Add the last region for this chr
			int chrEnd = genLoader.getChrMap().get(chr);
			segments.get(chr).add(new Region(genLoader.getGenome(), chr, currStart, chrEnd));
			numSeg+=segments.get(chr).size();
		}
		System.err.println("GenomeSegmenter: Initialized with "+numSeg+" segments.");
		return(segments);
	}
}
