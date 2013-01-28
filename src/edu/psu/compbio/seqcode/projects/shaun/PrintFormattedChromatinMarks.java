package edu.psu.compbio.seqcode.projects.shaun;

import java.util.*;
import java.io.*;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipData;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.locators.ChipChipLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.PointParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RegionParser;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;

/* PrintFormattedChromatinMarks 
 * Shaun 19th April 2010 
 */
public class PrintFormattedChromatinMarks {
    private static Genome gen;
    public static void main(String[] args) { 
	HashMap<Region, ChromBlock> blocks = new HashMap<Region, ChromBlock>();
	String regFileName = Args.parseString(args,"regions",null);
	String outRoot = Args.parseString(args,"out","out");
	
	if(regFileName !=null){			
	    try {
		gen = Organism.findGenome("mm8");
		ArrayList<Region> regions = loadRegionsFromFile(regFileName, -1);
		int numRegions=regions.size();
		int numTimes = times.length;
		
		//Construct the data matrix
		//Iterate through experiments and times
		int numExperiments = experiments.length;
		int arrayIndex=0; System.out.println(numExperiments+"\t"+numTimes);
		for(int t=0; t<numTimes; t++){	
		    for(int exp=0; exp<numExperiments; exp++){System.out.println(t+"\t"+exp+"\t"+arrayIndex);
			//Load experiment
			ChipChipLocator loc=null;
			if(experiments[exp][t].length==2)
			    loc=new ChipChipLocator(gen, experiments[exp][t][0], experiments[exp][t][1]);
			else if(experiments[exp][t].length==3)
			    loc=new ChipChipLocator(gen, experiments[exp][t][0], experiments[exp][t][1], experiments[exp][t][2]);
			ChipChipData data = loc.createObject();
			
			for(Region reg : regions) { System.out.println(reg);
			    if(!blocks.containsKey(reg)){
				blocks.put(reg, new ChromBlock(reg));
			    }
			    ChromBlock block = blocks.get(reg);
			    
			    String c = reg.getChrom();
			    int st = reg.getStart(), ed = reg.getEnd();
			    data.window(c, st, ed);
			    //System.out.println(data.getCount()+" probes in "+tiledRegion.getWidth()+"bp");
			    
			    for(int i = 0; i < data.getCount(); i++) {
				Point p = new Point(reg.getGenome(), reg.getChrom(), data.getPos(i));
				double ratio = getMeanRatio(data, i);
				if(!block.pointIndex.containsKey(p)){
				    block.pointIndex.put(p, block.numPoints);
				    block.data.add(block.numPoints, new ArrayList<Double>());
				    block.points.add(block.numPoints, p);
				    block.numPoints++;
				}int bindex = block.pointIndex.get(p); //System.out.println(p.toString()+" BIndex "+bindex);
				block.data.get(bindex).add(arrayIndex, ratio);
			    }
			}
			arrayIndex++;
		    }
		}
			
		//Print the matrix
		FileWriter dataFW = new FileWriter(outRoot+".data");
		FileWriter indexFW = new FileWriter(outRoot+".index");
		for(Region reg : regions) {
		    dataFW.write("#"+reg.toString()+"\n");
		    indexFW.write("#"+reg.toString()+"\n");
		    if(blocks.containsKey(reg)){
			ChromBlock block = blocks.get(reg);
			//Print data
			for(int p=0; p<block.numPoints; p++){
			    indexFW.write(block.points.get(p)+"\n");
			    //Print filler fields
			    for(int t=0; t<numTimes; t++){
				dataFW.write("0 ");
			    }
			    for(int d=0; d<block.data.get(p).size(); d++){
				dataFW.write(block.data.get(p).get(d).toString());
				if(d<block.data.get(p).size()-1){
				    dataFW.write(" ");
				}
			    }
			    dataFW.write("\n");
			}
		    } 
		    dataFW.write("\n");
		    indexFW.write("\n");
		}
		dataFW.close();
		indexFW.close();
	    } catch (IOException e) {
		e.printStackTrace();
	    } catch (NotFoundException e) {
		e.printStackTrace();
	    } catch (UnknownRoleException e) {
		e.printStackTrace();
	    } 		
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
    
    public static String[] names, times;
    public static String[][] h3k4Versions, h3k27Versions, h3k27subVersions, h3k79Versions, ring1bVersions, suz12Versions, rarChipSeq, backChipSeq;
    public static String[][][] experiments;
    
    static { 
	names = new String[] { "H3K4me3", "H3K27me3", "H3K79me2" };
	times = new String[] { "0", "2", "2+8hRA", "3", "4", "7"};
        
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
	    { "Mm H3K79me2:HBG3:mES vs H3:HBG3:mES", "median linefit, quantile norm 6tp", "2 (10/26/07, Hox Array)"},
	    { "Mm H3K79me2:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage", "median linefit, quantile norm 6tp" },
	    { "Mm H3K79me2:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm 6tp" },
	    { "Mm H3K79me2:HBG3:2+1 day vs WCE:HBG3:2+1 day", "median linefit, quantile norm 6tp" },
	    { "Mm H3K79me2:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm 6tp", "2 (10/26/07, Hox Array)"},
	    { "Mm H3K79me2:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm 6tp", "2 (10/26/07, Hox Array)"}
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

        experiments = new String[][][]{h3k4Versions, h3k27Versions, h3k79Versions};
        
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
		    PointParser pparser = new PointParser(gen);
		    Point p = pparser.execute(words[2]);
		    int rstart = p.getLocation()-(win/2)<1 ? 1:p.getLocation()-(win/2);
		    int rend = p.getLocation()+(win/2)>gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+(win/2)-1;
		    Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
		    regs.add(r);
                }else if(words.length>=1){
		    RegionParser parser = new RegionParser(gen);
		    Region q = parser.execute(words[0]);
		    if(win!=-1){
			int rstart = q.getMidpoint().getLocation()-(win/2)<1 ? 1:q.getMidpoint().getLocation()-(win/2);
			int rend = q.getMidpoint().getLocation()+(win/2)>gen.getChromLength(q.getChrom()) ? gen.getChromLength(q.getChrom()):q.getMidpoint().getLocation()+(win/2)-1;
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

