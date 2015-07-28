package edu.psu.compbio.seqcode.projects.akshay.regulatorydomains;

import java.util.Collections;
import java.util.List;

import edu.psu.compbio.seqcode.genome.location.Point;

public class RegulatoryClassProfile {
	
	// Has to be either BS (binding strength) or BD (binding dynamics)
	private String SORT_CRITERIA="BS";
	// Percentage of top binding events based on specified criteria to consider
	private double topPerc;
	private int[] motifindexes;
	private double homotypicCutoff;
	// SUMMARY; PEAKLISTS; DEPERC; CLUSPERC; AVGFC; Each of the options explained below
	private String outputFormat;
	private List<RegulatoryRegion> bindingEvents;
	private RegulatoryClass regClass;
	private String regClassName;
	
	public RegulatoryClassProfile(List<RegulatoryRegion> regRs, String sortBY, double topPC, int[] motifs, double homoCutoff) {
		bindingEvents = regRs;
		topPerc = topPC;
		motifindexes = motifs;
		homotypicCutoff = homoCutoff;
		try{
			if(!sortBY.equals("BS") && !sortBY.equals("BD")){
				throw new IllegalSortCriteria();
			}else{
				SORT_CRITERIA = sortBY;
			}
			
		}catch(Exception e){
			System.out.println("Invalid Sort Criteria: Privide either BS or BD");
			e.printStackTrace();
			
		}
		
		// The following code will generate the regulatory class
		
		// Sort the regulatory regions
		Collections.sort(bindingEvents);
		Collections.reverse(bindingEvents);
		
		//Creat a regulatory class
		int indexCutoff = (int) topPerc*bindingEvents.size()/100;
		regClass = new RegulatoryClass(bindingEvents.get(0),regClassName);
		for(int rR=1; rR<indexCutoff; rR++){
			regClass.addRegR(bindingEvents.get(rR));
		}
		
		
		
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	@SuppressWarnings("serial")
	public class IllegalSortCriteria extends Exception{
	}
	

}
