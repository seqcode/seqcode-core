package edu.psu.compbio.seqcode.projects.shaun;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.binding.BindingScan;
import edu.psu.compbio.seqcode.gse.datasets.binding.BindingScanLoader;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.binding.BindingScanGenerator;
import edu.psu.compbio.seqcode.gse.utils.iterators.SingleIterator;

//Adapted from Yuchun's ChromatinMarksAnalysis.java
public class MakeSTDBinaryData {
	
	private static Map<String, List<Integer>> geneExpr;
	private static String [] expressionTimes; 
	private final HashMap<String, String> HMMScans;
	private final String typeHMM = "BindingEventWrapper";
	
	private static boolean wellTiledOnly=true;
	private static boolean usingExpression=false;
	private static boolean limitGeneList=false;
	
	private static int windowUpstream =500;	//The amount of upstream sequence to include in the window
	private static int windowDownstream =500; //The amount of downstream sequence to include in the window
	private static double coverageThreshold = 0.5; //Threshold for how much of window has to be covered by a HMM domain
	
	private static String expressionFilename = "timeseries_all_geneProgramDisctete2.txt";
	
	private static String outputFilename = "binary_timeseries_all.tab";
	
    private static Genome genome;
    List<Gene> allGenes;
    List<Point> allPoints;

    public static void main(String[] args) {
    	boolean useGenes=true;
    	try{
    		Genome g = Organism.findGenome("mm8");
    		MakeSTDBinaryData analysis = new MakeSTDBinaryData(g, wellTiledOnly);
    		
    		if(useGenes){
    			ArrayList<String> validGeneIDs =new ArrayList<String>();
        		if(limitGeneList){
        			validGeneIDs = analysis.loadGeneNamesFromFile("geneList.txt");
        		}
	    		ArrayList<Gene> allGenes = new ArrayList<Gene>();
	    		WellTiledRegionParser wtrp = new WellTiledRegionParser();
	    		Expander<NamedRegion,Gene> refgenes = 
		            new RefGeneGenerator<NamedRegion>(g, "refGene");
		        Iterator<NamedRegion> itr = new ChromRegionIterator(g);
		        while(itr.hasNext()) { 
		        	NamedRegion chromosome = itr.next();
		            Iterator<Gene> genes = refgenes.execute(chromosome);
		            while(genes.hasNext()) { 
		                Gene x = genes.next();
		                if(x != null){
			                Region r = expandPoint(new Point(x.getGenome(), x.getChrom(), x.getFivePrime()));
			                if(!wellTiledOnly || wtrp.isWellTiled(r)){
			                	if(!limitGeneList || validGeneIDs.contains(x.getName())){
			                		allGenes.add(x);
			                	}
			                }
		                }
		            }
		        }
		        System.out.println("Loaded gene count: \t" + allGenes.size());
	    		analysis.setAllGenes(allGenes);
	    		if(usingExpression){parseExpression(new File(expressionFilename));}
	    		analysis.getHMMMarkCoverageByGene();
    		}else{
    			ArrayList<Point> allPoints = analysis.loadPeaksFromPeakFile("all.withmotif.tiledpeaks");
    			analysis.setAllPoints(allPoints);
    			analysis.getHMMMarkCoverageByPoint(allPoints);
    		}
    	}catch (Exception e){
    		e.printStackTrace();
    	}
    }
        
    public MakeSTDBinaryData(Genome g, boolean useHoxArrayGenesOnly) {
    	genome = g;        
        
    	HMMScans = new HashMap<String, String>();
    	HMMScans.put("H3K4me3_0_ES0", "Mm H3K4me3:HBG3:ES Stage vs H3:HBG3:ES Stage,median linefit, quantile norm,{config=edu.psu.compbio.seqcode.gse.projects.ppg.quantile_norm_HMM}");
    	HMMScans.put("H3K4me3_2.0_ES2", "Mm H3K4me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA,median linefit, quantile norm,{config=edu.psu.compbio.seqcode.gse.projects.ppg.quantile_norm_HMM}");
    	HMMScans.put("H3K4me3_2.3_ES2p8h", "Mm H3K4me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA,median linefit, quantile norm,{config=edu.psu.compbio.seqcode.gse.projects.ppg.quantile_norm_HMM}");
    	HMMScans.put("H3K4me3_3_ES2p1d", "Mm H3K4me3:HBG3:2+1 day vs H3:HBG3:2+1 day,median linefit, quantile norm,{config=edu.psu.compbio.seqcode.gse.projects.ppg.quantile_norm_HMM}");
    	HMMScans.put("H3K4me3_4_Olig2", "Mm H3K4me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage,median linefit, quantile norm,{config=edu.psu.compbio.seqcode.gse.projects.ppg.quantile_norm_HMM}");
    	HMMScans.put("H3K4me3_7_Hb9", "Mm H3K4me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage,median linefit, quantile norm,{config=edu.psu.compbio.seqcode.gse.projects.ppg.quantile_norm_HMM}");

    	HMMScans.put("H3K27me3_0_ES0", "Mm H3K27me3:HBG3:ES Stage vs H3:HBG3:ES Stage,median linefit, quantile norm,{config=edu.psu.compbio.seqcode.gse.projects.ppg.quantile_norm_HMM}");
    	HMMScans.put("H3K27me3_2.0_ES2", "Mm H3K27me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA,median linefit, quantile norm,{config=edu.psu.compbio.seqcode.gse.projects.ppg.quantile_norm_HMM}");
    	HMMScans.put("H3K27me3_2.3_ES2p8h", "Mm H3K27me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA,median linefit, quantile norm,{config=edu.psu.compbio.seqcode.gse.projects.ppg.quantile_norm_HMM}");
    	HMMScans.put("H3K27me3_3_ES2p1d", "Mm H3K27me3:HBG3:2+1 day vs H3:HBG3:2+1 day,median linefit, quantile norm,{config=edu.psu.compbio.seqcode.gse.projects.ppg.quantile_norm_HMM}");
    	HMMScans.put("H3K27me3_4_Olig2", "Mm H3K27me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage,median linefit, quantile norm,{config=edu.psu.compbio.seqcode.gse.projects.ppg.quantile_norm_HMM}");
    	HMMScans.put("H3K27me3_7_Hb9", "Mm H3K27me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage,median linefit, quantile norm,{config=edu.psu.compbio.seqcode.gse.projects.ppg.quantile_norm_HMM}");

    	HMMScans.put("H3K79me2_0_ES0", "Mm H3K79me2:HBG3:mES vs H3:HBG3:mES,median linefit, quantile norm 6tp,{config=edu.psu.compbio.seqcode.gse.shaun.yyc_optimized_HMM}");
    	HMMScans.put("H3K79me2_2.0_ES2", "Mm H3K79me2:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage,median linefit, quantile norm 6tp,{config=edu.psu.compbio.seqcode.gse.shaun.yyc_optimized_HMM}");
    	HMMScans.put("H3K79me2_2.3_ES2p8h", "Mm H3K79me2:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA,median linefit, quantile norm 6tp,{config=edu.psu.compbio.seqcode.gse.shaun.yyc_optimized_HMM}");
    	HMMScans.put("H3K79me2_3_ES2p1d", "Mm H3K79me2:HBG3:2+1 day vs WCE:HBG3:2+1 day,median linefit, quantile norm 6tp,{config=edu.psu.compbio.seqcode.gse.shaun.yyc_optimized_HMM}");
    	HMMScans.put("H3K79me2_4_Olig2", "Mm H3K79me2:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage,median linefit, quantile norm 6tp,{config=edu.psu.compbio.seqcode.gse.shaun.yyc_optimized_HMM}");
    	HMMScans.put("H3K79me2_7_Hb9", "Mm H3K79me2:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage,median linefit, quantile norm 6tp,{config=edu.psu.compbio.seqcode.gse.shaun.yyc_optimized_HMM}");
    }
    public void setAllGenes(List<Gene> genes){allGenes = genes;}
    public void setAllPoints(List<Point> points){allPoints = points;}
    
    public static Region expandPoint(Point p) {
    	int start=0, end=0;
    	start = p.getLocation()-windowUpstream;
    	end=p.getLocation()+windowDownstream;
    	
    	Region reg = new Region(genome, p.getChrom(), start, end);
    	
    	return(reg);
    }
    
	public static void parseExpression(File f) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(f));
		geneExpr = new HashMap<String,List<Integer>>();
		String line = br.readLine();
		String[] split = line.split("\t");
		int numTimes = split.length-1;
		expressionTimes = new String [numTimes];		
		for (int t=1; t<split.length; t++) {expressionTimes[t-1]=split[t];}
		
		while((line = br.readLine()) != null) { 
			split = line.split("\t");
			String tmpGene = split[0];
			ArrayList<Integer> tmpList = new ArrayList<Integer>();
			for (int time=1; time<split.length; time++) {
				tmpList.add(Integer.parseInt(split[time]));
			}
			geneExpr.put(tmpGene, tmpList);
		}
	}

    
    public void getHMMMarkCoverageByGene(){
    	BindingScanGenerator generator=null;
		HashMap<Gene, HashMap<String, Integer>> geneBindingEventCounts = new HashMap<Gene, HashMap<String, Integer>>();
		for (Gene gene: allGenes){
			geneBindingEventCounts.put(gene, new HashMap<String, Integer>());
		}
					
		// for each marks, each gene, calculate the length of HMM called domain coverage
		for (String key: HMMScans.keySet()){
			try{
				BindingScanLoader loader = new BindingScanLoader();
				generator = new BindingScanGenerator(loader);
				Collection<BindingScan> scans;
				scans = loader.loadScans(genome, HMMScans.get(key), typeHMM); 
					
				Iterator<BindingScan> it = scans.iterator();
				while (it.hasNext()){
					generator.addBindingScan(it.next());				
				}
			}
			catch (SQLException e){
				e.printStackTrace();
			}
			
			// each gene, sum the length of overlapping domains within range
			for (Gene gene: allGenes){
				Region range=expandPoint(new Point(gene.getGenome(), gene.getChrom(), gene.getFivePrime()));
				
				Iterator<BindingEvent> domains = generator.execute(range);	
				int length = 0;
				while(domains.hasNext()){
					BindingEvent domain = domains.next();
					// generator returns the whole domain. only count the portion inside the promoter
					length+=domain.getOverlapSize(range);
				}
				geneBindingEventCounts.get(gene).put(key, new Integer(length));
			}
		}
		String filename = outputFilename;
		writeToFile(filename, makeGeneOutputString(geneBindingEventCounts), false);
    }
    
    public void getHMMMarkCoverageByPoint(ArrayList<Point> points){
    	BindingScanGenerator generator=null;
		HashMap<Point, HashMap<String, Integer>> pointBindingEventCounts = new HashMap<Point, HashMap<String, Integer>>();
		for (Point p : points){
			pointBindingEventCounts.put(p, new HashMap<String, Integer>());
		}
					
		// for each marks, each gene, calculate the length of HMM called domain coverage
		for (String key: HMMScans.keySet()){
			try{
				BindingScanLoader loader = new BindingScanLoader();
				generator = new BindingScanGenerator(loader);
				Collection<BindingScan> scans;
				scans = loader.loadScans(genome, HMMScans.get(key), typeHMM); 
					
				Iterator<BindingScan> it = scans.iterator();
				while (it.hasNext()){
					generator.addBindingScan(it.next());				
				}
			}
			catch (SQLException e){
				e.printStackTrace();
			}
			
			// each gene, sum the length of overlapping domains within range
			for (Point p: points){
				Region range=expandPoint(p);
				
				Iterator<BindingEvent> domains = generator.execute(range);	
				int length = 0;
				while(domains.hasNext()){
					BindingEvent domain = domains.next();
					// generator returns the whole domain. only count the portion inside the promoter
					length+=domain.getOverlapSize(range);
				}
				pointBindingEventCounts.get(p).put(key, new Integer(length));
			}
		}
		String filename = outputFilename;
		writeToFile(filename, makeOutputString(pointBindingEventCounts), false);
    }
    
 
    
    private String makeGeneOutputString(HashMap<Gene, HashMap<String, Integer>> dataTable){
		// output as data table, row-gene, column-mark count
		StringBuffer sb = new StringBuffer();
		sb.append("Gene ID\t");
		sb.append("Gene Name\t");
		sb.append("Promoter Range\t");
	
		double promLen = windowUpstream+windowDownstream+1;
		
		if(usingExpression){
			for(int i=0; i<expressionTimes.length; i++){
				sb.append(expressionTimes[i]).append("\t");
			}
		}
		
		TreeSet<String> sortedKeys = new TreeSet<String>(HMMScans.keySet());
		for (String key: sortedKeys){
			sb.append(key).append("\t");
		}
		sb.append("\n");
		for (Gene gene: allGenes){
			if(!usingExpression || geneExpr.containsKey(gene.getID())){
				sb.append(gene.getID()).append("\t");
				sb.append(gene.getName()).append("\t");
				Region range=expandPoint(new Point(gene.getGenome(), gene.getChrom(), gene.getFivePrime()));
				sb.append(range.regionString()).append("\t");
				
				if(usingExpression){
					List<Integer> gl = geneExpr.get(gene.getID());
					Iterator<Integer> gi = gl.iterator();				
					while(gi.hasNext()){
						sb.append(gi.next()).append("\t");
					}
				}
							
				for (String key: sortedKeys){
					double coverage = (dataTable.get(gene).get(key)).doubleValue()/promLen;
					//sb.append(coverage).append("\t");
					
					if(coverage>=coverageThreshold){
						sb.append("1\t");
					}else{
						sb.append("0\t");
					}
				}
				sb.append("\n");
			}
		}
		return sb.toString();
    }
    private String makeOutputString(HashMap<Point, HashMap<String, Integer>> dataTable){
		// output as data table, row-gene, column-mark count
		StringBuffer sb = new StringBuffer();
		sb.append("Point\t");
		sb.append("Point\t");
		sb.append("RegionRange\t");
	
		double promLen = windowUpstream+windowDownstream+1;
				
		TreeSet<String> sortedKeys = new TreeSet<String>(HMMScans.keySet());
		for (String key: sortedKeys){
			sb.append(key).append("\t");
		}
		sb.append("\n");
		for (Point p: allPoints){
			sb.append(p.getLocationString()).append("\t");
			sb.append(p.getLocationString()).append("\t");
			Region range=expandPoint(p);
			sb.append(range.regionString()).append("\t");
			for (String key: sortedKeys){
				double coverage = (dataTable.get(p).get(key)).doubleValue()/promLen;
				if(coverage>=coverageThreshold){
					sb.append("1\t");
				}else{
					sb.append("0\t");
				}
			}
			sb.append("\n");
		}
		
		return sb.toString();
    }
    //Load a set of regions from a peak file
	public ArrayList<Point> loadPeaksFromPeakFile(String filename){
		ArrayList<Point> peaks = new ArrayList<Point>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            
	            PointParser pparser = new PointParser(genome);
	            if(words.length>=3){
	                Point p = pparser.execute(words[2]);
	            	peaks.add(p);
                }else{
                	Point p = pparser.execute(words[0]);
	            	peaks.add(p);
                }
	        }reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(peaks);
	}
	//Load a list of gene names
	public ArrayList<String> loadGeneNamesFromFile(String filename){
		ArrayList<String> ids = new ArrayList<String>();
		try{
			File gFile = new File(filename);
			if(!gFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(gFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            ids.add(words[0]);
	        }
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}return(ids);
	}
    private void writeToFile(String filename, String s, boolean append) {
    	try {
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(filename, append)));
			out.print(s);
			out.flush();
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
    }   

}

