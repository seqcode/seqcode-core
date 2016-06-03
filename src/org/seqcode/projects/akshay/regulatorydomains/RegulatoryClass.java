package org.seqcode.projects.akshay.regulatorydomains;

import java.util.ArrayList;
import java.util.List;

/**
 * 
 * @author akshaykakumanu
 *
 */
public class RegulatoryClass {
	private List<RegulatoryRegion> regRegs;
	private String className;
	private int numRegs;
	
	public RegulatoryClass(RegulatoryRegion regR, String name) {
		regRegs = new ArrayList<RegulatoryRegion>();
		regRegs.add(regR);
		className = name;
		numRegs = 1;
	}
	
	
	public void addRegR(RegulatoryRegion regR){
		regRegs.add(regR);
		numRegs++;
	}
	
	// Gettors
	public int getNumRegRegions(){return numRegs;}
	public String getClassName(){return className;}
	public List<RegulatoryRegion> getRegRegs(){return regRegs;}
	public String getSummaryString(){
		StringBuilder sb = new StringBuilder();
		sb.append("#"+className);sb.append("\n");
		sb.append("Avg-Binding-Intensity:"+"\t");sb.append(this.getAvgBindingIntensity());sb.append("\n");
		sb.append("Average-binding-dynamics:"+"\t");sb.append(this.getAvgBindingDynamicsIndex());sb.append("\n");
		sb.append("Perc-upregulated:"+"\t");sb.append(this.getPercUpRegGenes());sb.append("\n");
		sb.append("Perc-down-regulated:"+"\t");sb.append(this.getPercDownRegGenes());sb.append("\n");
		sb.append("Class-Size"+"\t");sb.append(numRegs);sb.append("\n");
		return sb.toString();
	}

	public double getAvgBindingIntensity(){
		double avgBI=0;
		for(int rR=0; rR<regRegs.size(); rR++){
			avgBI = avgBI + regRegs.get(rR).getBindingIntensity();
		}
		avgBI = avgBI/numRegs;
		return avgBI;
	}
	
	public double getAvgBindingDynamicsIndex(){
		double avgBDI=0;
		for(int rR=0; rR<regRegs.size(); rR++){
			avgBDI = avgBDI + regRegs.get(rR).getBindingDynamicsIndex();
		}
		avgBDI = avgBDI/numRegs;
		return avgBDI;
	}
	
	public double getPercUpRegGenes(){
		double upRegPerc=0;
		for(int rR=0; rR<regRegs.size(); rR++){
			if(regRegs.get(rR).getTargetGeneFoldChange() > 0){
				upRegPerc++;
			}
		}
		upRegPerc = upRegPerc*100/numRegs;
		return upRegPerc;
	}
	
	public double getPercDownRegGenes(){
		double downRegPerc=0;
		for(int rR=0; rR<regRegs.size(); rR++){
			if(regRegs.get(rR).getTargetGeneFoldChange() < 0){
				downRegPerc++;
			}
		}
		downRegPerc = downRegPerc*100/numRegs;
		return downRegPerc;
	}
	
	public double getPercGenesInClusterInd(int c){
		double clusPerc=0;
		for(int rR=0; rR<regRegs.size(); rR++){
			if(regRegs.get(rR).getTargetGeneClusterIndex() == c){
				clusPerc++;
			}
		}
		clusPerc = clusPerc*100/numRegs;
		return clusPerc;
	}
	
	public double getAvgTrgetGeneFoldChange(){
		double avgFC=0;
		for(int rR=0; rR<regRegs.size(); rR++){
			avgFC = avgFC+regRegs.get(rR).getTargetGeneFoldChange();
		}
		avgFC =avgFC/numRegs;
		return avgFC;
	}
	
	public double getMotifhitRate(int m, double[] motifThreshLvls){
		double motHitRate = 0;
		for(int rR=0; rR<regRegs.size(); rR++){
			if(regRegs.get(rR).getBestMotifScore(m) >= motifThreshLvls[m]){
				motHitRate++;
			}
		}
		motHitRate = motHitRate*100/numRegs;
		return motHitRate;
	}
	
	public List<String> getAllTargetGeneList(){
		List<String> geneList = new ArrayList<String>();
		for(int rR=0; rR<regRegs.size(); rR++){
			geneList.add(regRegs.get(rR).getTargetGeneName());
		}
		return geneList;
	}
	
	public List<String> getUpRegTargeGeneLists(){
		List<String> geneList = new ArrayList<String>();
		for(int rR=0; rR<regRegs.size(); rR++){
			if(regRegs.get(rR).getTargetGeneFoldChange() > 0){
				geneList.add(regRegs.get(rR).getTargetGeneName());
			}
		}
		return geneList;
	}
	
	public List<String> getDownRegTargeGeneLists(){
		List<String> geneList = new ArrayList<String>();
		for(int rR=0; rR<regRegs.size(); rR++){
			if(regRegs.get(rR).getTargetGeneFoldChange() < 0){
				geneList.add(regRegs.get(rR).getTargetGeneName());
			}
		}
		return geneList;
	}
	
	public List<String> getTargetGeneListwithClusterInd(int c){
		List<String> geneList = new ArrayList<String>();
		for(int rR=0; rR<regRegs.size(); rR++){
			if(regRegs.get(rR).getTargetGeneClusterIndex() == c){
				geneList.add(regRegs.get(rR).getTargetGeneName());
			}
		}
		return geneList;
	}
	
	
}
