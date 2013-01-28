package edu.psu.compbio.seqcode.projects.shaun.rnaseq;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.deepseq.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.deepseq.ReadHit;
import edu.psu.compbio.seqcode.projects.shaun.rnaseq.genemodels.GeneTUnit;
import edu.psu.compbio.seqcode.projects.shaun.rnaseq.genemodels.TUnit;

public class RNASeqRegionCounter extends RNASeqAnalyzer{

	protected HashMap<Region, Double> regionCounts = new HashMap<Region, Double>();
	
	public RNASeqRegionCounter(String[] args) {
		super(args);
	}

	public void execute(){
		Double totalHitCount = 0.0;
		for(DeepSeqExpt e : experiments){
			totalHitCount += e.getHitCount();
		}
		
		System.err.println("Counting reads in "+knownGenes.size()+" known genes");
		
		//Dumb iteration over all regions for now
		for(Region reg : getRegionsOfInterest()){
			regionCounts.put(reg, 0.0);
			
			//If we have multiple experiments, just pool all data
			for(DeepSeqExpt e : experiments){
				for(ReadHit hit : e.loadHits(reg)){
					regionCounts.put(reg, regionCounts.get(reg)+hit.getWeight());
				}
			}
		}
		
		//Print counts
		System.out.println(String.format("Coord\tHits\tRegionLength\tFPKM"));
		for(Region reg : getRegionsOfInterest()){
			double count = regionCounts.get(reg);
			double fpkm = (count/((double)reg.getWidth()/1000.0))/(totalHitCount/1000000);
			System.out.println(String.format("%s\t%.2f\t%d\t%.2f", reg.getLocationString(),count,reg.getWidth(), fpkm));
		}
		
		//Close loaders
		this.cleanup();
	}
	
	public static void main(String[] args) {
		RNASeqRegionCounter counter = new RNASeqRegionCounter(args);
		counter.execute();
	}
	
	public void printError() {
		printError("RNASeqRegionCounter");
	}	
}
