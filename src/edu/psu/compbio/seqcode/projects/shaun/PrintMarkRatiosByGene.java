package edu.psu.compbio.seqcode.projects.shaun;
import java.util.*;
import java.io.*;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipData;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.locators.ChipChipLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;


public class PrintMarkRatiosByGene {
	
	public static void main(String[] args) { 
		
		String exptLabel1 = "Suz12@8hr";
		String [] experiment1 =suz12Versions[1];
		String exptLabel2 = "H3K27me3@8hr";
		String [] experiment2 =h3k27Versions[2];
		int TSSUpWindow = 1000;
		int TSSDownWindow = 1000;
		boolean restrictRegions = true;
		
		Set<Region> regSet=null;
		
		try {
			//ExpressionLoader exprLoader = new ExpressionLoader();
			
			if(restrictRegions){
				File tilingList = new File("Hox_tiled_regions.txt");
				//File tilingList = new File("testRegions2.coords");
				WellTiledRegionParser tiled = new WellTiledRegionParser(tilingList, 10000);
				regSet = tiled.getRegions();
			}
			
			Genome g = Organism.findGenome("mm8");
			int numRegions=regSet.size();
			
			//Find the genes in the current regions
			List<String> coveredGenes = new ArrayList<String>();
			List<String> geneCoords = new ArrayList<String>();
			double[][] geneVals = new double[2][];
			int numCoveredGenes=0;
			int r=0;
			RefGeneGenerator<NamedRegion> geneGen = new RefGeneGenerator<NamedRegion>(g, "refGene");
			ChromRegionIterator chroms = new ChromRegionIterator(g);
			Iterator<Gene> genes = new ExpanderIterator<NamedRegion,Gene>(geneGen, chroms);
			
			while(genes.hasNext()) { 
				Gene gene = genes.next();
				if(!restrictRegions || isInRegions(gene, regSet)){
					coveredGenes.add(gene.getName());
					if(gene.getStrand() == '+'){
						geneCoords.add(String.format("%s:%d:+", gene.getChrom(), gene.getStart()));
					}else{
						geneCoords.add(String.format("%s:%d:-", gene.getChrom(), gene.getEnd()));
					}
					numCoveredGenes++;
				}
			}System.out.println("Gene Names Loaded");
			//Populate value array 1
			geneVals[0] = new double[numCoveredGenes];
			ChipChipLocator loc = 
				new ChipChipLocator(g, experiment1[0], experiment1[1]);
			ChipChipData data = loc.createObject();
			geneGen = new RefGeneGenerator<NamedRegion>(g, "refGene");
			chroms = new ChromRegionIterator(g);
			genes = new ExpanderIterator<NamedRegion,Gene>(geneGen, chroms);
			r=0;
			while(genes.hasNext()) { 
				Gene gene = genes.next();
				if(!restrictRegions || isInRegions(gene, regSet)){
					String c = gene.getChrom();
					int st, ed;
					if(gene.getStrand() == '+'){
						st= gene.getStart()-(TSSUpWindow/2);
						ed = gene.getStart()+(TSSDownWindow/2);
					}else{
						st= gene.getEnd()-(TSSDownWindow/2);
						ed = gene.getEnd()+(TSSUpWindow/2);
					}
					
					data.window(c, st, ed);//System.out.println(gene.getName()+"\t"+c+"\t"+st+"\t"+ed);
					double total=0;
					for(int i = 0; i < data.getCount(); i++) {
						total+= getMeanRatio(data, i);
					}
					double numD = (double)data.getCount();
				//	System.out.println(total+"\t"+numD);
					if(numD>0){
						geneVals[0][r] = total/numD;
					}else{
						geneVals[0][r] =0;
					}
					r++;
				}
			}System.out.println("Expt 1 Loaded");
			
			//Populate value array 2
			geneVals[1] = new double[numCoveredGenes];
			loc = 
				new ChipChipLocator(g, experiment2[0], experiment2[1]);
			data = loc.createObject();
			geneGen = new RefGeneGenerator<NamedRegion>(g, "refGene");
			chroms = new ChromRegionIterator(g);
			genes = new ExpanderIterator<NamedRegion,Gene>(geneGen, chroms);
			r=0;
			while(genes.hasNext()) { 
				Gene gene = genes.next();
				if(!restrictRegions || isInRegions(gene, regSet)){
					String c = gene.getChrom();
					int st, ed;
					if(gene.getStrand() == '+'){
						st= gene.getStart()-(TSSUpWindow/2);
						ed = gene.getStart()+(TSSDownWindow/2);
					}else{
						st= gene.getEnd()-(TSSDownWindow/2);
						ed = gene.getEnd()+(TSSUpWindow/2);
					}
					
					data.window(c, st, ed);//System.out.println(gene.getName()+"\t"+c+"\t"+st+"\t"+ed);
					double total=0;
					for(int i = 0; i < data.getCount(); i++) {
						total+= getMeanRatio(data, i);
					}
					double numD = (double)data.getCount();
				//	System.out.println(total+"\t"+numD);
					if(numD>0){
						geneVals[1][r] = total/numD;
					}else{
						geneVals[1][r] =0;
					}
					r++;
				}
			}System.out.println("Expt 2 Loaded");

			
			
			//Finally, print the rations
			geneGen = new RefGeneGenerator<NamedRegion>(g, "refGene");
			chroms = new ChromRegionIterator(g);
			genes = new ExpanderIterator<NamedRegion,Gene>(geneGen, chroms);
			//Header
			System.out.println("GeneName\tCoord\t"+exptLabel1+"\t"+exptLabel2+"\t"+exptLabel1+":"+exptLabel2+" Ratio");
			for(r=0; r<numCoveredGenes; r++){
				System.out.print(String.format("%s\t%s\t", coveredGenes.get(r), geneCoords.get(r)));
				if(!Double.isNaN(geneVals[0][r]) && !Double.isInfinite(geneVals[0][r]) && !Double.isNaN(geneVals[1][r]) && !Double.isInfinite(geneVals[1][r])){
					double rat = geneVals[0][r]/geneVals[1][r];
					System.out.print(String.format("%f\t\t%f\t%f",geneVals[0][r], geneVals[1][r], rat));
				}else{
					System.out.print("0.0000\t");
				}
				System.out.print("\n");
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (NotFoundException e) {
			e.printStackTrace();
		} catch (UnknownRoleException e) {
			e.printStackTrace();
		} 		
	}
	
	public static double getMeanRatio(ChipChipData d, int i) { 
		double sum = 0.0;
		double count = 0;
		for(int j = 0; j < d.getReplicates(i); j++) { 
			if(!Double.isInfinite(d.getRatio(i, j)) && !Double.isNaN(d.getRatio(i, j))){
				sum += d.getRatio(i, j);
				count += 1;
			}
		}
		return count > 0 ? sum / count : 0.0;
	}
	public static boolean isInRegions(StrandedRegion gene, Set<Region> regions){
		boolean covered=false;
		for(Region reg : regions){
			if(reg.contains(gene)){
				covered=true;
			}
		}		
		return(covered);
	}
	
	//public static double[] exprCutoffs = { 2.0, 2.0, 2.0 };

	public static boolean hoxArray;
	public static boolean hoxRegions;
	public static String[] names, times;
	public static int[] blocks;
	public static String[][] h3k4Versions, h3k27Versions, h3k79Versions, suz12Versions;
	public static String[] exprExptNames;
		
	
	public static File gene2Affy = new File("Mouse430_2_clean_annos_unique.txt");
	public static File affy2Gene = new File("Mouse430_2_clean_annos_probe.txt");

	static { 
	    names = new String[] { "H3K4me3", "H3K27me3", "H3K79me2" };
        times = new String[] { "ES", "ES+2d", "RA+8h", "RA+1d", "Olig2", "Hb9"};
        blocks = new int[] { 5, 2, 4, 6, 6, 3 };

        exprExptNames = new String[] { 
        		"ES Stage PMA",
        		"ES+2d Stage PMA",
        		"ES+2d Stage, post RA PMA",
        		"Olig2 Stage PMA",
        		"Hb9 Stage PMA"
        };
        

        h3k4Versions = new String[][] { 
                { "Mm H3K4me3:HBG3:ES Stage vs H3:HBG3:ES Stage", "median linefit, quantile norm" },
                { "Mm H3K4me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA", "median linefit, quantile norm" },
                { "Mm H3K4me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm" },
                { "Mm H3K4me3:HBG3:2+1 day vs H3:HBG3:2+1 day", "median linefit, quantile norm" },
                { "Mm H3K4me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm" },
                { "Mm H3K4me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm" }
        };

        h3k27Versions = new String[][] { 
                { "Mm H3K27me3:HBG3:ES Stage vs H3:HBG3:ES Stage", "median linefit, quantile norm" },
                { "Mm H3K27me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA", "median linefit, quantile norm" },
                { "Mm H3K27me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm" },
                { "Mm H3K27me3:HBG3:2+1 day vs H3:HBG3:2+1 day", "median linefit, quantile norm" },
                { "Mm H3K27me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm" },
                { "Mm H3K27me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm" }
        }; 
        
        h3k79Versions = new String[][] { 
                { "Mm H3K79me2:HBG3:mES vs H3:HBG3:mES", "median linefit, quantile norm 6tp" },
                { "Mm H3K79me2:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage", "median linefit, quantile norm 6tp" },
                { "Mm H3K79me2:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm 6tp" },
                { "Mm H3K79me2:HBG3:2+1 day vs WCE:HBG3:2+1 day", "median linefit, quantile norm 6tp" },
                { "Mm H3K79me2:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm 6tp" },
                { "Mm H3K79me2:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm 6tp" }
        };
        
        suz12Versions = new String[][] { 
        		{"Mm Suz12:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage","median linefit, quantile norm 2tp"},	
        		{"Mm Suz12:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA","median linefit, quantile norm 2tp"}	
        };
	}

		
	
	public static String getNameString(Set<String> names) { 
		StringBuilder sb = new StringBuilder();
		for(String n : names) { 
			if(sb.length() > 0) { sb.append(","); }
			sb.append(n);
		}
		return sb.toString();
	}
}
