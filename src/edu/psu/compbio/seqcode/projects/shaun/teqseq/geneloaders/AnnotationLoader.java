package edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.gse.utils.RealValuedHistogram;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.GenomeLoader;

public abstract class AnnotationLoader {

	protected List<AGene> genes;
	protected GenomeLoader gLoad;
	protected Genome gen;
	
	public AnnotationLoader(GenomeLoader gl){
		gLoad = gl;
		gen = gLoad.getGenome();
		genes = new ArrayList<AGene>();
	}
	
	public abstract Collection<AGene>  loadGenes();
	public List<AGene> getGenes() {return genes;}
	
	//Simple method to analyze gene lengths in known annotations
	public void geneLengths(){
		RealValuedHistogram histo;
		ArrayList<Double> lengths=new ArrayList<Double>();
		double minLen=100000000, maxLen=0;
		
		for(AGene g : genes){
			double len = g.getCoords().getEnd() - g.getCoords().getStart();
			lengths.add(len);
			if(len>maxLen)
				maxLen = len;;
			if(len<minLen)
				minLen = len;
		}
		System.out.println("\n\nGene Lengths\nMin: "+minLen+"\tMax: "+maxLen);
		histo = new RealValuedHistogram(0, maxLen, 50);
		histo.addValues(lengths);
		histo.printContents();
	}
	//Simple method to analyze exon-exon distances in known annotations
	public void exonExonDistances(){
		RealValuedHistogram histo;
		ArrayList<Double> distances=new ArrayList<Double>();
		double minDist=100000000, maxDist=0;
		
		for(AGene g : genes){
			//Only single-isoform genes
			if(g.getIsoforms().size()==1){
				AIsoform i = g.getIsoforms().get(0);
				ArrayList<ARegion> exons = i.getExons();
				for(int x=1; x<exons.size(); x++){
					double dist;
					if(exons.get(x).getCoords().getStart() > exons.get(x-1).getCoords().getEnd())
						dist= exons.get(x).getCoords().getStart() - exons.get(x-1).getCoords().getEnd();
					else
						dist= exons.get(x-1).getCoords().getEnd()-exons.get(x).getCoords().getStart();
					distances.add(dist);
					if(dist>maxDist)
						maxDist=dist;
					if(dist<minDist)
						minDist=dist;
				}
			}
		}
		System.out.println("\n\nExon-Exon Distances\nMin: "+minDist+"\tMax: "+maxDist);
		histo = new RealValuedHistogram(0, maxDist, 50);
		histo.addValues(distances);
		histo.printContents();
	}
}
