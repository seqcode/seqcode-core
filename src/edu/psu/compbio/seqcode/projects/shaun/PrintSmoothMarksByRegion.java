package edu.psu.compbio.seqcode.projects.shaun;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingScan;
import edu.psu.compbio.seqcode.gse.datasets.binding.BindingScanLoader;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipData;
import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.expression.ExpressionLoader;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.locators.ChipChipLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.binding.BindingExpander;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;

/* PrintSmoothMarksByRegion: 
 * Started by Shaun 30-Oct-2007
 * Real simple method... just smoothes and outputs some experimental data from the PPG
 * 
 * Settings used in PPG Hox paper heat-maps: window=500bp, step=250bp
 * 
 */
public class PrintSmoothMarksByRegion {
	
	public static void main(String[] args) { 
		
		boolean chipChip=true;
		String exptLabel = "H3K27me3";
		String [][] experiment =h3k27Versions;
		String [][] background =backChipSeq;
		String [] timeLabels = times;
		int smoothWin=500;
		int smoothOff=250;
	    int plotType =0; //0=Heat-Map, 1=Difference from first, 2=Difference between consecutive times
		
		File tilingList = new File("testRegions3.coords");

		try {
			//ExpressionLoader exprLoader = new ExpressionLoader();
			
			WellTiledRegionParser tiled = new WellTiledRegionParser(tilingList, 10000);
			Set<Region> tiledRegions = tiled.getRegions();
			Genome g = Organism.findGenome("mm8");
			int numRegions=tiledRegions.size();
			
			//Find the genes in the current regions
			List<String> [] coveredGenes = new ArrayList [numRegions];
			List<String> [] geneCoords = new ArrayList [numRegions];
			int numCoveredGenes=0;
			int r=0;
			for(Region tiledRegion : tiledRegions) {
				RefGeneGenerator<NamedRegion> geneGen = new RefGeneGenerator<NamedRegion>(g, "refGene");
				ChromRegionIterator chroms = new ChromRegionIterator(g);
				Iterator<Gene> genes = new ExpanderIterator<NamedRegion,Gene>(geneGen, chroms);
				coveredGenes[r]=new ArrayList<String>();
				geneCoords[r]=new ArrayList<String>();
				while(genes.hasNext()) { 
					Gene gene = genes.next();
					if(isInRegion(gene, tiledRegion)){
						coveredGenes[r].add(String.format("%s", gene.getName()));
						if(gene.getStrand() == '+'){
							geneCoords[r].add(String.format("%s:%d:+", gene.getChrom(), gene.getStart()));
						}else{
							geneCoords[r].add(String.format("%s:%d:-", gene.getChrom(), gene.getEnd()));
						}
						numCoveredGenes++;
					}
				}
				r++;
			}
		
			 
			//Set up the 3-d arrays
			int vLen = experiment.length;
			double[][][] smoothArray = new double[vLen][][];
			String[][][] smoothLabels = new String[vLen][][];
			
			for(int x=0; x<vLen; x++){
				
				smoothArray[x] = new double [tiledRegions.size()][];
				smoothLabels[x] = new String [tiledRegions.size()][];
				r=0;
				for(Region tiledRegion : tiledRegions) { 
					String c = tiledRegion.getChrom();
					int st = tiledRegion.getStart(), ed = tiledRegion.getEnd();
					smoothArray[x][r] = new double [((ed-st)/smoothOff)+1];
					smoothLabels[x][r] = new String [((ed-st)/smoothOff)+1];
					for(int l =0; l<(ed-st)/smoothOff; l++){
						//Label is the center of the current window
						smoothLabels[x][r][l]=String.format("%s:%d", c, (st+((smoothOff)*l)+(smoothWin/2)));
					}
					r++;
				}
			}
			
			//Populate the smoothed arrays
			for(int x=0; x<vLen; x++){
				if(chipChip){
					ChipChipLocator loc=null;
					if(experiment[x].length==2)
						loc=new ChipChipLocator(g, experiment[x][0], experiment[x][1]);
					else if(experiment[x].length==3)
						loc=new ChipChipLocator(g, experiment[x][0], experiment[x][1], experiment[x][2]);
					ChipChipData data = loc.createObject();
					System.out.println(data.getCount()+" probes in array");

					r=0;
					for(Region tiledRegion : tiledRegions) { 
						
						String c = tiledRegion.getChrom();
						int st = tiledRegion.getStart(), ed = tiledRegion.getEnd();
						
						data.window(c, st, ed);
						//System.out.println(data.getCount()+" probes in "+tiledRegion.getWidth()+"bp");
						
						double[] array = new double[data.getCount()];
						int[] positions = new int[data.getCount()];
						
						for(int i = 0; i < data.getCount(); i++) {
							positions[i] = data.getPos(i);
							array[i] = getMeanRatio(data, i);
						}
						
						//Data loaded, now smooth over overlapping windows of sequence
						int lastStartIndex=0; int l=0;
						for(int i = st; i < ed; i+=smoothOff) {
							double count =0;
							boolean iFound=false;
							smoothArray[x][r][l]=0;
							for(int j=lastStartIndex; j<array.length && positions[j]<i+smoothWin; j++){
								if(!iFound && positions[j]>=i){
									lastStartIndex=j; iFound=true;
								}
								if(iFound){
									smoothArray[x][r][l]+=array[j];
									count+=1;
								}							
							}
							if(count>=1){
								smoothArray[x][r][l]=smoothArray[x][r][l]/count;
							}else{smoothArray[x][r][l]=1;}
							l++;
						}
						r++;
					}
				}else{//ChIP-Seq
					double readLength=26;
					double readExtension = 174;
					double iphittot = 0, backhittot=0;
					Organism org = Organism.getOrganism("Mus musculus");
					Genome gen = org.getGenome("mm8");
					ArrayList<ChipSeqLocator> iplocs = new ArrayList<ChipSeqLocator>();
					for(String[] l : experiment){
						iplocs.add(new ChipSeqLocator(l[0], l[1]));						
					}
					ArrayList<ChipSeqLocator> backlocs = new ArrayList<ChipSeqLocator>();
					for(String[] l : background){
						backlocs.add(new ChipSeqLocator(l[0], l[1]));
					}
					
					ArrayList<SeqExpt> IPhandles = new ArrayList<SeqExpt>();
					ArrayList<SeqExpt> backhandles = new ArrayList<SeqExpt>();
					try {
						double genomeLen = gen.getGenomeLength();
						//Load background
						for(ChipSeqLocator back : backlocs){
							System.err.print(String.format("%s\t", back.getExptName()));
							SeqExpt curr = new SeqExpt(gen, back);
			                curr.setReadLength(readLength);
			                curr.setReadExtension(readExtension);
							backhandles.add(curr);
							backhittot += curr.getHitCount();
						}
			            System.err.print(String.format("%.0f reads loaded\n", backhittot));
						//Load experiments
			            x=0;
						for(ChipSeqLocator ip : iplocs){
							IPhandles = new ArrayList<SeqExpt>();
							System.err.print(String.format("%s\t", ip.getExptName()));
							SeqExpt curr = new SeqExpt(gen, ip);
						    curr.setReadLength(readLength);
			                curr.setReadExtension(readExtension);
							IPhandles.add(curr);
							iphittot = curr.getHitCount();
						
							System.err.print(String.format("%.0f reads loaded\n", iphittot));
						
				            
				            //Iterate through regions
							r=0;
				            for(Region currRegion : tiledRegions) { 
				            	String c = currRegion.getChrom();
								int st = currRegion.getStart(), ed = currRegion.getEnd();							
								LinkedList<StrandedRegion> ipHits = new LinkedList<StrandedRegion>();
								LinkedList<StrandedRegion> backHits = new LinkedList<StrandedRegion>();
								for(SeqExpt IP: IPhandles){
									ipHits.addAll(IP.loadExtendedHits(currRegion));									
								}
								for(SeqExpt back: backhandles){
									backHits.addAll(back.loadExtendedHits(currRegion));
								}
								int [] ipHitCounts = makeHitLandscape(ipHits, currRegion, smoothWin, smoothOff);
				                int [] backHitCounts = makeHitLandscape(backHits, currRegion, smoothWin, smoothOff);
				                int l=0;
				                for(int i = st; i < ed; i+=smoothOff) {
				                	double ipWinHits=(double)ipHitCounts[l];
				                	double backWinHits=(double)backHitCounts[l];
				                	//smoothArray[x][r][l]= binomialSampleEquality(ipWinHits, backWinHits, iphittot, backhittot);
				                	smoothArray[x][r][l]=ipWinHits/iphittot;
									l++;
				                }
				                r++;
				            }
				            x++;
						}
					} catch (SQLException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				
				}
			}
			
			//Finally, print the array
			r=0;
			for(Region tiledRegion : tiledRegions) { 
				//Header
				System.out.println(String.format(exptLabel));
				System.out.print("Genes\t");
				for(int q=0; q<coveredGenes[r].size(); q++){
					System.out.print(coveredGenes[r].get(q)+"\t");	
				}System.out.print("\nGeneCoords\t");
				for(int q=0; q<geneCoords[r].size(); q++){
					System.out.print(geneCoords[r].get(q)+"\t");	
				}
				System.out.print("\nPosition\t");
				for(int x=0; x<vLen; x++){
					System.out.print(String.format("%s\t", timeLabels[x]));
				}System.out.print("\n");
				
				int st = tiledRegion.getStart(), ed = tiledRegion.getEnd();
				for(int l =0; l<(ed-st)/smoothOff; l++){
					System.out.print(String.format("%s\t", smoothLabels[0][r][l]));
					for(int x=0; x<vLen; x++){
						if(!Double.isNaN(smoothArray[x][r][l]) && !Double.isInfinite(smoothArray[x][r][l])){
							if(plotType==0)
								System.out.print(String.format("%.10f\t",smoothArray[x][r][l]));
							else if(plotType==1)
								System.out.print(String.format("%f\t",smoothArray[x][r][l]-smoothArray[0][r][l]));
							else{
								if(x==0)
									System.out.print("0.000\t");
								else
									System.out.print(String.format("%f\t",smoothArray[x][r][l]-smoothArray[x-1][r][l]));
							}
						}else{
							System.out.print("0.0000\t");
						}
					}
					System.out.print("\n");
				}
				r++;
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
	public static boolean isInRegion(StrandedRegion gene, Region region){
		boolean covered=false;
		if(region.contains(gene)){
			covered=true;
		}		
		return(covered);
	}
	
	//public static double[] exprCutoffs = { 2.0, 2.0, 2.0 };

	public static boolean hoxArray;
	public static boolean hoxRegions;
	public static String[] names, times, times3;
	public static int[] blocks;
	public static String[][] h3k4Versions, h3k27Versions, h3k27subVersions, h3k79Versions, ring1bVersions, suz12Versions, rarChipSeq, backChipSeq;
	public static String[] exprExptNames;
		
	
	public static File gene2Affy = new File("Mouse430_2_clean_annos_unique.txt");
	public static File affy2Gene = new File("Mouse430_2_clean_annos_probe.txt");

	static { 
	    names = new String[] { "H3K4me3", "H3K27me3", "H3K79me2" };
        times = new String[] { "0", "2", "2+8hRA", "3", "4", "7"};
        times3 = new String[] { "2", "2+RA+8h", "3"};
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
                { "Mm H3K27me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm", "2 (2/1/08)" },
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
        h3k27subVersions = new String[][] {
        	       { "Mm H3K27me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA", "median linefit, quantile norm" },
        	       { "Mm H3K27me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm", "2 (2/1/08)" },
                   { "Mm H3K27me3:HBG3:2+1 day vs H3:HBG3:2+1 day", "median linefit, quantile norm" },
        };
        ring1bVersions = new String[][] {
        		{"Mm Ring1b:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage", "median linefit, quantile norm"},	
        		{"Mm Ring1b:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA","median linefit, quantile norm"},	
        		{"Mm Ring1b:HBG3:2+1 day vs WCE:HBG3:2+1 day", "median linefit, quantile norm"}
        };
        suz12Versions = new String[][] {
        		{"Mm Suz12:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage","median linefit"},	
        		{"Mm Suz12:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA","median linefit"},	
        		{"Mm Suz12:HBG3:2+1 day vs WCE:HBG3:2+1 day","median linefit", "2 (7/24/08)"}	
        };
        rarChipSeq = new String[][]{
        		{"PPG_Solexa_RAR_ES+2d", "ELAND_unique"},
        		{"PPG_Solexa_RAR_8hr", "ELAND_unique"}        		
        };
        backChipSeq = new String[][]{
        		{"PPG_Solexa_WCE_2+1", "ELAND_unique"},
        		{"PPG_Solexa_WCE_ES+2d", "ELAND_unique"}
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
	private static int [] makeHitLandscape(LinkedList<StrandedRegion> hits, Region currReg, int winWidth, int winStep){
		int numBins = (int)(currReg.getWidth()/winStep);
		int [] land = new int[numBins+1];
		for(int i=0; i<=numBins; i++){land[i]=0;}
		for(StrandedRegion r : hits){
			int binstart = (int)Math.max(0, ((double)((r.getStart()-currReg.getStart())/winStep)-(Math.floor(winWidth/winStep)-1)));
			int binend = (int)(Math.min((double)(r.getEnd()-currReg.getStart()), (double)currReg.getWidth())/winStep);
			for(int i=binstart; i<=binend; i++){
				land[i]++; 
			}
		}
		return(land);
	}
	/* Binomial test for differences between two population proportions */
	private static double binomialSampleEquality(double X1, double X2, double n1, double n2){
		double P1 = X1/n1;
		double P2 = X2/n2;
		double P = (X1+X2)/(n1+n2);
		double Z0 = (P1-P2)/(Math.sqrt(P*(1-P)*((1/n1)+(1/n2))));
		return(Z0);
	}
}
