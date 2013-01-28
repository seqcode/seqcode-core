package edu.psu.compbio.seqcode.projects.shaun;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipData;
import edu.psu.compbio.seqcode.gse.datasets.expression.ExpressionLoader;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.locators.ChipChipLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;

/* PrintSmoothMarksByRegion: 
 * Reworked by Shaun 25-May-2008
 * Real simple method... just smoothes and outputs some experimental data from the PPG
 */
public class PrintSmoothMarksByGene {
	
	public static void main(String[] args) { 
		
		String exptLabel = "H3K4me3";
		String [][] experiment =h3k4Versions;  
		int TSSWindowUpstream =250;
		int TSSWindowDownstream =250;
	    int plotType =0; //0=Heat-Map, 1=Difference from first, 2=Difference between consecutive times
		boolean restrictRegions = true;
		boolean loadGenesFromFile=true;
		String geneFile = "HoxGenes.coords";
		
		Set<Region> regSet=null;
		
		try {
			if(restrictRegions){
				File tilingList = new File("Hox_tiled_regions.txt");
				//File tilingList = new File("testRegions2.coords");
				WellTiledRegionParser tiled = new WellTiledRegionParser(tilingList, 10000);
				regSet = tiled.getRegions();
			}
			
			Genome g = Organism.findGenome("mm8");
			int numRegions=regSet.size();
			List<Gene> coveredGenes = new ArrayList<Gene>();
			List<String> geneCoords = new ArrayList<String>();
			int vLen = experiment.length;
			double[][] geneVals = new double[vLen][];
			int numCoveredGenes=0;
			int r=0;
			ChromRegionIterator chroms = new ChromRegionIterator(g);
			Iterator<Gene> genes = null;
			//Find the genes in the current regions
			if(!loadGenesFromFile){
				RefGeneGenerator<NamedRegion> geneGen = new RefGeneGenerator<NamedRegion>(g, "refGene");
				genes = new ExpanderIterator<NamedRegion,Gene>(geneGen, chroms);
			}else{
				File gFile = new File(geneFile);
				ArrayList<Gene> tmpGenes = new ArrayList<Gene>();
				if(gFile.isFile()){
					BufferedReader reader;
					reader = new BufferedReader(new FileReader(gFile));
					String line;
					while((line= reader.readLine())!=null){
						String [] tokens = line.split("[\\s*\\t\\r\\n\\f]");
						String chr = tokens[0];
						chr = chr.replaceAll("chr", "");
						int start = new Integer(tokens[1]).intValue();
						int end = new Integer(tokens[2]).intValue();
						char str = tokens[3].charAt(0);
						//String probe = tokens[4];
						//String name = tokens[5];
						//Gene currg = new Gene(g, chr, start, end, name, probe,str, "NONE");
						String name = tokens[4];
						Gene currg = new Gene(g, chr, start, end, name, name,str, "NONE");
						tmpGenes.add(currg);						
					}
					reader.close();
				}
				genes = tmpGenes.iterator();
			}
			
			while(genes.hasNext()) { 
				Gene gene = genes.next();
				if(!restrictRegions || isInRegions(gene, regSet)){
					coveredGenes.add(gene);
					if(gene.getStrand() == '+'){
						geneCoords.add(String.format("%s:%d-%d", gene.getChrom(), gene.getStart()-TSSWindowUpstream, gene.getStart()+TSSWindowDownstream));
					}else{
						geneCoords.add(String.format("%s:%d-%d", gene.getChrom(), gene.getEnd()-TSSWindowDownstream, gene.getEnd()+TSSWindowUpstream));
					}
					numCoveredGenes++;
				}
			}
			//Populate the smoothed arrays
			for(int x=0; x<vLen; x++){	
				geneVals[x] = new double[numCoveredGenes];
				ChipChipLocator loc=null;
				if(experiment[x].length==2)
					loc=new ChipChipLocator(g, experiment[x][0], experiment[x][1]);
				else if(experiment[x].length==3)
					loc=new ChipChipLocator(g, experiment[x][0], experiment[x][1], experiment[x][2]);
				ChipChipData data = loc.createObject();
				r=0;
				for(Gene gene : coveredGenes) {
					String c = gene.getChrom();
					int st, ed;
					if(gene.getStrand() == '+'){
						st= gene.getStart()-TSSWindowUpstream;
						ed = gene.getStart()+TSSWindowDownstream;
					}else{
						st= gene.getEnd()-TSSWindowDownstream;
						ed = gene.getEnd()+TSSWindowUpstream;
					}
					
					data.window(c, st, ed);
					double total=0;
					for(int i = 0; i < data.getCount(); i++) {
						total+= getMeanRatio(data, i);
					}
					double numD = (double)data.getCount();
					if(numD>0){
						geneVals[x][r] = total/numD;
					}else{
						geneVals[x][r] =0;
					}
					r++;
				}
			}
			
			//Header
			System.out.println(String.format(exptLabel));
			System.out.print("GeneName\tCoord\t");
			for(int x=0; x<vLen; x++){
				System.out.print(String.format("%s\t", times[x]));
			}System.out.print("\n");
			for(r=0; r<numCoveredGenes; r++){
				System.out.print(String.format("%s\t%s\t%s\t", coveredGenes.get(r).getID(), coveredGenes.get(r).getName(), geneCoords.get(r)));
				for(int x=0; x<vLen; x++){
					if(!Double.isNaN(geneVals[x][r]) && !Double.isInfinite(geneVals[x][r])){
						if(plotType==0)
							System.out.print(String.format("%f\t",geneVals[x][r]));
						else if(plotType==1)
							System.out.print(String.format("%f\t",geneVals[x][r]-geneVals[0][r]));
						else{
							if(x==0)
								System.out.print("0.000\t");
							else
								System.out.print(String.format("%f\t",geneVals[x][r]-geneVals[x-1][r]));
						}
					}else{
						System.out.print("0.0000\t");
					}
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
        times = new String[] { "0", "2", "2+8hRA", "3", "4", "7"};
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
                { "Mm H3K27me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm", "2 (2/1/08)"  },
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
        		{"Mm Suz12:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage",	"median linefit, quantile norm 2tp"},	
        		{"Mm Suz12:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA",	"median linefit, quantile norm 2tp"}
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
