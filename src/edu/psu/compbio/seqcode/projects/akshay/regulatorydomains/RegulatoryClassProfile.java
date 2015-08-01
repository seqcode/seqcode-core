package edu.psu.compbio.seqcode.projects.akshay.regulatorydomains;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import edu.psu.compbio.seqcode.genome.location.Point;

/**
 * 
 * @author akshaykakumanu
 *
 */
public class RegulatoryClassProfile {
	
	
	// Percentage of top binding events based on specified criteria to consider
	private double topPerc;
	
	// The following two lists should contain the following
	// The motif indices that you want the regulatory regions to have
	// The no.pf motif instances that they should have (same indexing as above)
	private List<Integer> motifindexes=new ArrayList<Integer>();
	private List<Integer> minMotifHit = new ArrayList<Integer>();
	
	private List<Double> motifLogOddThrehs = new ArrayList<Double>(); // This list contains all the motifs
	private double homotypicCutoff;
	// SUMMARY; PEAKLISTS; UPPERC; DOWNPERC; CLUSPERC; AVGFC; Each of the options explained below
	private String outputFormat;
	private List<RegulatoryRegion> bindingEvents;
	private RegulatoryClass regClass;
	private String regClassName;
	private int numClusters=0; // number of clusters of genes
	private String sortType = "foldchange";
	
	/**
	 * 
	 * @param regRs
	 * @param topPC
	 * @param sType
	 * @param motifsInds
	 * @param motifThres
	 * @param homoCutoff
	 * @param nC
	 * @param Output
	 * @param ClassName
	 */
	public RegulatoryClassProfile(List<RegulatoryRegion> regRs, double topPC, String sType, List<Integer> motifsInds, List<Integer> MinMotifHits, List<Double> motifThres, double homoCutoff, int nC, String Output, String ClassName) {
		bindingEvents = regRs;
		topPerc = topPC;
		if(motifsInds != null)
			motifindexes = motifsInds;
		if(motifThres != null)
			motifLogOddThrehs = motifThres;
		if(MinMotifHits != null)
			minMotifHit = MinMotifHits;
		homotypicCutoff = homoCutoff;
		numClusters=nC;
		outputFormat=Output;
		regClassName =ClassName;
		sortType = sType;
		// The following code will generate the regulatory class
		
		// Sort the regulatory regions
		Collections.sort(bindingEvents);
		if(sortType.equals("foldchange"))
			Collections.reverse(bindingEvents);
		
		//Create a regulatory class
		int indexCutoff = (int) topPerc*bindingEvents.size()/100;
		for(int rR=0; rR<indexCutoff; rR++){
			boolean addBE = true;
			RegulatoryRegion currRr= bindingEvents.get(rR);
			//now check for motifs
			for(int m=0; m<motifindexes.size(); m++){
				if(currRr.getBestMotifScore(motifindexes.get(m)) < motifLogOddThrehs.get(motifindexes.get(m))){
					addBE=false;
				}
				if(currRr.getMotifHitCount(motifindexes.get(m)) < minMotifHit.get(m)){
					addBE=false;
				}
			}
			if(addBE){
				if(currRr.getHomotypicIndex() < homotypicCutoff){
					addBE=false;
				}
			}
			if(addBE){
				if(regClass == null){
					regClass= new RegulatoryClass(currRr,regClassName);
				}else{
					regClass.addRegR(bindingEvents.get(rR));
				}
			}
		}
		
		if(outputFormat.equals("SUMMARY")){ printSUMMARY();}
		if(outputFormat.equals("PEAKLISTS")){ printPEAKLISTS();}
		if(outputFormat.equals("UPPERC")){ printUPPERC();}
		if(outputFormat.equals("DOWNPERC")){ printDOWNPERC();}
		if(outputFormat.equals("CLUSPERC")){ printCLUSPERC(numClusters);}
		if(outputFormat.equals("AVGFC")){ printAVGFC();}

	}
	
	
	
	
	// Output printing methods
	
	// SUMMARY
	/**
	 * Prints the summary of a reg class
	 */
	private void printSUMMARY(){
		System.out.print(regClass.getSummaryString());
	}
	
	//PEAKLISTS
	/**
	 * Prints Peak-location(\t)Target-gene=foldchange(\t)Target-gene-cluster-index
	 */
	private void printPEAKLISTS(){
		System.out.println("#"+regClass.getClassName());
		List<RegulatoryRegion> rC = regClass.getRegRegs();
		for(int rR=0; rR<rC.size();rR++){
			System.out.println(rC.get(rR).getPeakLocation()+
					"\t"+rC.get(rR).getTargetGeneFoldChange()+
					"\t"+Integer.toString(rC.get(rR).getTargetGeneClusterIndex()));
		}
	}
	
	//UPPERC
	/**
	 * Just prints the percentage of target genes that are upregulated
	 */
	private void printUPPERC(){
		System.out.print("#"+regClass.getClassName()+"\t");
		System.out.println(regClass.getPercUpRegGenes());
	}

	//DOWNPERC
	/**
	 * Just prints the percentage of target genes that are downregualted
	 */
	private void printDOWNPERC(){
		System.out.print("#"+regClass.getClassName()+"\t");
		System.out.println(regClass.getPercDownRegGenes());
	}
	
	//CLUSPERC
	/*
	 * Prints the percentage of target genes that belong to each cluster class
	 * Cluster 0 are target genes which are assigned to no cluster (non-DE genes)
	 */
	private void printCLUSPERC(int nClusters){
		System.out.print("#"+regClass.getClassName()+"\t");
		StringBuilder sb = new StringBuilder();
		for(int c=0; c<=nClusters; c++){
			sb.append(regClass.getPercGenesInClusterInd(c)+"\t");
		}
		sb.deleteCharAt(sb.length()-1);
		System.out.println(sb.toString());
	}
	
	//AVGFC
	/**
	 * Prints the avg fold-change of the target genes
	 */
	private void printAVGFC(){
		System.out.print("#"+regClass.getClassName()+"\t");
		System.out.println(regClass.getAvgTrgetGeneFoldChange());
	}
	
	

}
