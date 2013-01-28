package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipData;
import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.locators.ChipChipLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Gene;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ExpanderIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.PointParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RefGeneGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RegionParser;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;

public class CleanAndNormArrays {

	public static Genome genome;
	public static int spikeWin=2;
	public static double spikeFilter=4; //Filter at this number of times the average of the surroundings. 
	public static int smoothWin=4000, smoothOff=500;
	public static String mode="smooth";
	public static int geneTSSWin=4000;
	
	public static void main(String[] args) {
		try{
			CleanAndNormArrays cna = new CleanAndNormArrays();
			ArrayList <ArrayList<QuantElement>> arrayValues = new ArrayList<ArrayList<QuantElement>>();
			ArrayList <ArrayList<Point>> arrayPositions = new ArrayList<ArrayList<Point>>();
			HashMap<Region, ArrayList<Point>> regionProbes = new HashMap<Region,ArrayList<Point>>();
			String regFileName = Args.parseString(args,"regions",null);
			String outRoot = Args.parseString(args,"out","out");
			FileWriter dataFW = new FileWriter(outRoot+".data");
			String exptName = Args.parseString(args,"expt","K27"); //K27, K4, K79, Suz12, Ring1b
			String [][] currExpt = exptNames_K27;
			if(exptName.equals("K4"))
				currExpt = exptNames_K4;
			if(exptName.equals("K79"))
				currExpt = exptNames_K79;
			if(exptName.equals("Suz12"))
				currExpt = exptNames_Suz12;
			if(exptName.equals("Ring1b"))
				currExpt = exptNames_Ring1b;
			mode = Args.parseString(args,"mode","smooth"); //smooth, probes, genes
			
			int spikesFiltered=0;
			int totalProbes=0;
			int arrayProbes=0;
			
			if(regFileName !=null){
				genome = Organism.findGenome("mm8");
	        	ArrayList<Region> regions = loadRegionsFromFile(regFileName, -1);
	        
				for(int b=0; b<currExpt.length; b++){
					arrayValues.add(new ArrayList<QuantElement>());
					arrayPositions.add(new ArrayList<Point>());
					
					String expt = currExpt[b][0];
			        String version = currExpt[b][1];
			        String replicate = currExpt[b][2];
			        	
			       	//Load the ChIP-chip data
					ChipChipLocator loc = null;
			        if(replicate == null || replicate.equals("all")){
			        	loc = new ChipChipLocator(genome, expt, version);	        	
			        }else{
			        	loc = new ChipChipLocator(genome, expt, version, replicate);
			        }
					ChipChipData data = loc.createObject();
	
					/////////////////////////////////////////////////////////
					/////// Loading & de-spiking
					/////////////////////////////////////////////////////////
					for(Region reg : regions){
						if(b==0)
							regionProbes.put(reg, new ArrayList<Point>());
						//if(b==0){System.out.println("#"+reg);}
						String chrName = reg.getChrom();
						int cstart = reg.getStart();
						int cend = reg.getEnd();
						data.window(chrName, cstart, cend);
						
						for(int i=0; i<data.getCount(); i++){
							double currProbe = getMeanRatio(data, i);
							Point pos = new Point(genome, chrName, data.getPos(i));
							//if(b==0){System.out.println(pos);}
							if(b==0){
								regionProbes.get(reg).add(pos);
								arrayProbes++;
							}
							
							//Filter spurious values
							double left=0, right=0, avg=0, lcount=0, rcount=0, tcount=0;
							for(int j=i-spikeWin; j<=i+spikeWin; j++){
								if(j>=0 && j<i){
									left+=getMeanRatio(data, j); lcount++; 
									avg+=getMeanRatio(data, j); tcount++;
								}else if(j<data.getCount() && j>i){
									right+=getMeanRatio(data, j); rcount++;
									avg+=getMeanRatio(data, j); tcount++;
								}
							}
							if(lcount>0){left/=lcount;}
							if(rcount>0){right/=rcount;}
							if(tcount>0){avg/=tcount;}
							if((left>0 && currProbe/left>spikeFilter) && (right>0 && currProbe/right>spikeFilter)){
								currProbe=avg; spikesFiltered++; 
								/*System.out.print(b+":"+pos);
								for(int j=i-spikeWin; j<=i+spikeWin; j++){
									if(j>=0 && j<data.getCount()){
										System.out.print("\t"+getMeanRatio(data, j));
										if(j==i){System.out.print("*");}
									}
								}System.out.println("");
								*/
							}
							
							QuantElement qe = cna.new QuantElement(pos, currProbe);
							arrayValues.get(b).add(qe);							
							arrayPositions.get(b).add(pos);
							totalProbes++; 
						}
					}
			        
				}
				System.out.println("Spikes Filtered:\t"+spikesFiltered);
				System.out.println("Total Probes:\t"+totalProbes);
				System.out.println("Probes in tiled regions: "+arrayProbes);
	
				/////////////////////////////////////////////////////////
				/////// Quantile normalization
				/////////////////////////////////////////////////////////
				//Sort the HashMaps
				for(int a=0; a<arrayValues.size(); a++){
					Collections.sort(arrayValues.get(a));
				}
				
				//Arrays should now be sorted; do a quick quality check
				int numProbes = arrayValues.get(0).size();
				for(int a=1; a<arrayValues.size(); a++){
					if(numProbes!=arrayValues.get(a).size()){
						String exptN = currExpt[a][0]+";"+currExpt[a][1]+";"+currExpt[a][2];
						System.out.println(exptN+"\t"+arrayValues.get(a).size());
					}
				}
				
				//Quantile normalize & store in a Hash
				ArrayList <HashMap<Point, Double>> mapNorm = new ArrayList<HashMap<Point, Double>>();
				for(int a=0; a<arrayValues.size(); a++){
					mapNorm.add(new HashMap<Point, Double>());
				}
				double numArrays = (double)arrayValues.size();
				for(int v=0; v<numProbes; v++){
					double sum = 0;
					for(int a=0; a<arrayValues.size(); a++){
						sum+=arrayValues.get(a).get(v).value;
					}
					for(int a=0; a<arrayValues.size(); a++){
						arrayValues.get(a).get(v).value = sum/numArrays;
						mapNorm.get(a).put(arrayValues.get(a).get(v).probe, arrayValues.get(a).get(v).value);
					}
				}
				
				
				//Print the array headings
				for(int b=0; b<currExpt.length; b++){
					String expt = currExpt[b][0];
			        String version = currExpt[b][1];
			        String replicate = currExpt[b][2];
			        dataFW.write("\t"+expt+";"+version+";"+replicate);
				}dataFW.write("\n");
				
				if(mode.equals("smooth")){
					/////////////////////////////////////////////////////////
					/////// Smooth data using a sliding window and print
					/////////////////////////////////////////////////////////
					for(Region reg : regions){
						dataFW.write("#"+reg+"\n");
						for(int winStart=reg.getStart(); winStart<(reg.getEnd()-smoothWin); winStart+=smoothOff){
							//int winEnd = winStart+smoothWin > reg.getEnd() ? reg.getEnd() : winStart+smoothWin;
							int winEnd = winStart+smoothWin+1;
							if(winEnd-winStart > smoothWin){
								dataFW.write(reg.getChrom()+":"+(winStart+winEnd)/2);
								for(int b=0; b<currExpt.length; b++){
									double val=0, count=0;
									for(Point p : regionProbes.get(reg)){
										if(p.getLocation()>=winStart && p.getLocation()<=winEnd){
											val+=mapNorm.get(b).get(p);
											count++;
										}
									}
									val/=count;
									dataFW.write("\t"+val);
								}dataFW.write("\n");
							}
						}					
					}
				}else if(mode.equals("probes")){
					/////////////////////////////////////////////////////////
					/////// Print the probes themselves
					/////////////////////////////////////////////////////////
					for(Region reg : regions){
						for(Point p : regionProbes.get(reg)){
							dataFW.write(p.getLocationString());
							for(int b=0; b<currExpt.length; b++){
								dataFW.write("\t"+mapNorm.get(b).get(p));
							}dataFW.write("\n");
						}
					}
				}else if(mode.equals("genes")){
					/////////////////////////////////////////////////////////
					/////// Average in a window around gene TSSs
					/////////////////////////////////////////////////////////
					ChromRegionIterator chroms = new ChromRegionIterator(genome);
					RefGeneGenerator<NamedRegion> geneGen = new RefGeneGenerator<NamedRegion>(genome, "refGene");
					Iterator<Gene> genes = new ExpanderIterator<NamedRegion,Gene>(geneGen, chroms);
					while(genes.hasNext()) {							
						Gene gene = genes.next();
						int st, ed;
						if(gene.getStrand() == '+'){
							st= gene.getStart()-(geneTSSWin/2);
							ed = gene.getStart()+(geneTSSWin/2);
						}else{
							st= gene.getEnd()-(geneTSSWin/2);
							ed = gene.getEnd()+(geneTSSWin/2);
						}
						Region TSS = new Region(genome, gene.getChrom(), st, ed);
						for(Region reg : regions){
							if(reg.contains(TSS)){
								dataFW.write(gene.getID()+"\t"+gene.getName());
								for(int b=0; b<currExpt.length; b++){
									double val=0, count=0;
									for(Point p : regionProbes.get(reg)){
										if(p.getLocation()>=st && p.getLocation()<=ed){
											val+=mapNorm.get(b).get(p);
											count++;
										}
									}
									val/=count;
									dataFW.write("\t"+val);
								}dataFW.write("\n");
							}
						}
					}
				}
			}
			
			dataFW.close();
	    } catch (UnknownRoleException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
    public static double getMeanRatio(ChipChipData d, int i) { 
		double sum = 0.0;
		int count = 0;
		for(int j = 0; j < d.getReplicates(i); j++) { 
			Double curr = d.getRatio(i, j);
			if(!curr.isNaN() && !curr.isInfinite()){
				sum+=curr;
				count += 1;
			}			
		}
		return count > 0 ? sum / (double)count : 0.0;
	}
	
	public static String [][] exptNames_K27 = new String[][]{
		//Hox Arrays
		{"Mm H3K27me3:HBG3:ES Stage vs H3:HBG3:ES Stage",	"median linefit",	"1"},
		{"Mm H3K27me3:HBG3:ES Stage vs H3:HBG3:ES Stage",	"median linefit",	"3"},
		{"Mm H3K27me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA",	"median linefit",	"1"},
		{"Mm H3K27me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA",	"median linefit",	"2"},
		{"Mm H3K27me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA",	"median linefit",	"2 (2/1/08)"},
		{"Mm H3K27me3:HBG3:2+1 day vs H3:HBG3:2+1 day",	"median linefit",	"1"},
		{"Mm H3K27me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage",	"median linefit",	"chip1"},
		{"Mm H3K27me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage",	"median linefit",	"chip2"},
		{"Mm H3K27me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage",	"median linefit",	"1"},
	};
	public static String [][] exptNames_K4 = new String[][]{
		//Hox Arrays
		{"Mm H3K4me3:HBG3:ES Stage vs H3:HBG3:ES Stage",	"median linefit",	"1"},
		{"Mm H3K4me3:HBG3:ES Stage vs H3:HBG3:ES Stage",	"median linefit",	"3"},
		{"Mm H3K4me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA",	"median linefit",	"1"},
		{"Mm H3K4me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA",	"median linefit",	"2"},
		{"Mm H3K4me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA",	"median linefit",	"1"},
		{"Mm H3K4me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA",	"median linefit",	"2 (2/1/08)"},
		{"Mm H3K4me3:HBG3:2+1 day vs H3:HBG3:2+1 day",	"median linefit",	"1"},
		{"Mm H3K4me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage",	"median linefit",	"chip1"},
		{"Mm H3K4me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage",	"median linefit",	"chip2"},
		{"Mm H3K4me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage",	"median linefit",	"1"},
		{"Mm H3K4me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage",	"median linefit",	"2"},
	};
	public static String [][] exptNames_K79 = new String[][]{
		//Hox Arrays
		{ "Mm H3K79me2:HBG3:mES vs H3:HBG3:mES", "median linefit", "2 (10/26/07, Hox Array)"},
	    { "Mm H3K79me2:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage", "median linefit", "1 (5/20/08)"},
	    { "Mm H3K79me2:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA", "median linefit", "1 (5/20/08)" },
	    { "Mm H3K79me2:HBG3:2+1 day vs WCE:HBG3:2+1 day", "median linefit", "1 (5/20/08)" },
	    { "Mm H3K79me2:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit", "2 (10/26/07, Hox Array)"},
	    { "Mm H3K79me2:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit", "2 (10/26/07, Hox Array)"}
	};
	public static String [][] exptNames_Ring1b = new String[][] {
    		{"Mm Ring1b:HBG3:2+1 day vs WCE:HBG3:2+1 day","median linefit","1 (7/24/08)"},
    		{"Mm Ring1b:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage","median linefit","1 (7/24/08)"},
    		{"Mm Ring1b:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA","median linefit","1 (7/24/08)"}
    };
	public static String [][] exptNames_Suz12 = new String[][] {
    		{"Mm Suz12:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage","median linefit","1 (5/20/08)"},
    		{"Mm Suz12:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA","median linefit","1 (5/20/08)"},
    		{"Mm Suz12:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA","median linefit","2 (7/24/08)"},
    		{"Mm Suz12:HBG3:2+1 day vs WCE:HBG3:2+1 day","median linefit", "2 (7/24/08)"}	
    };
	
	class QuantElement implements Comparable<QuantElement>{
		public Point probe;
		public Double value;
		public QuantElement(Point p, Double v){probe=p; value=v;}
		public int compareTo(QuantElement t1) {
			Double t1val = ((QuantElement) t1).value;  
			if(this.value > t1val){return(1);}
			else if(this.value < t1val){return(-1);}
			else{return(0);}
		}		
	}
	
	//Load a set of regions from a peak file
    public static ArrayList<Region> loadRegionsFromFile(String filename, int win){
	ArrayList<Region> regs = new ArrayList<Region>();
		try{
			File pFile = new File(filename);
		    if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
		    BufferedReader reader = new BufferedReader(new FileReader(pFile));
		    String line;// = reader.readLine(); //Ignore first line
		    while ((line = reader.readLine()) != null) {
		    	line = line.trim();
		    	String[] words = line.split("\\s+");
			    if(words.length>=3 && win!=-1){
				    PointParser pparser = new PointParser(genome);
				    Point p = pparser.execute(words[2]);
				    int rstart = p.getLocation()-(win/2)<1 ? 1:p.getLocation()-(win/2);
				    int rend = p.getLocation()+(win/2)>genome.getChromLength(p.getChrom()) ? genome.getChromLength(p.getChrom()):p.getLocation()+(win/2)-1;
				    Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
				    regs.add(r);
			    }else if(words.length>=1){	
				    RegionParser parser = new RegionParser(genome);
				    Region q = parser.execute(words[0]);
				    if(win!=-1){
					int rstart = q.getMidpoint().getLocation()-(win/2)<1 ? 1:q.getMidpoint().getLocation()-(win/2);
					int rend = q.getMidpoint().getLocation()+(win/2)>genome.getChromLength(q.getChrom()) ? genome.getChromLength(q.getChrom()):q.getMidpoint().getLocation()+(win/2)-1;
					Region r = new Region(q.getGenome(), q.getChrom(), rstart, rend);
					if(r!=null){regs.add(r);}
				    }else{
					if(q!=null){regs.add(q);}
				    }
				}
		    }reader.close();
		} catch (FileNotFoundException e) {
		    // TODO Auto-generated catch block
		    e.printStackTrace();
		} catch (IOException e) {
		    // TODO Auto-generated catch block
		    e.printStackTrace();
		}
		return(regs);
    }
}
