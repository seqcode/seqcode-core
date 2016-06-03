package org.seqcode.projects.shaun.teqseq.core;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.gsebricks.verbs.location.RegionParser;


/**
 * ReadUtils: A collection of utilites for processing reads
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ReadUtils {

	/**
	 * Make an array representing the (weighted) counts of read 5' ends at each position in a region
	 * @param hits	A List of alignment hits 
	 * @param r		A Region containing the hits
	 * @param strand	A character representing the strand wanted (+/-/.)
	 * @return A double array
	 */
	public static double[] makeReadStartsArray(List<AlignHit> hits, Region r, char strand){
		double[] starts = new double[r.getWidth()+1];
		for(int i=0; i<=r.getWidth(); i++){starts[i]=0.0;}
		
		for(AlignHit h : hits){
			int offset = h.getFivePrime() - r.getStart();
			if(offset>=0 && offset<=r.getWidth()){
				if(strand == '.' || strand == h.getStrand()){
					starts[offset]+=h.getWeight();
				}
			}
		}
		
		return(starts);
	}
	
	/**
	 * Sum the weighted starts appearing in a region
	 * @param hits	A List of alignment hits 
	 * @param r		A Region containing the hits
	 * @param strand	A character representing the strand wanted (+/-/.)
	 * @return A double array
	 */
	public static double sumHitStartWeights(List<AlignHit> hits, Region r, char strand){
		double startSum=0;
		for(AlignHit h : hits){
			int offset = h.getFivePrime() - r.getStart();
			if(offset>=0 && offset<=r.getWidth()){
				if(strand == '.' || strand == h.getStrand()){
					startSum += h.getWeight();					
				}
			}
		}
		
		return(startSum);
	}
	
	/**
	 * Simple method to load a list of regions from a file, returning them organized by chromosome
	 * @param gen
	 * @param filename
	 * @return
	 */
	public static Map<String, List<Region>> loadRegionsFromFile(Genome gen, String filename){
		Map<String, List<Region>> regs = new HashMap<String, List<Region>>();
		RegionParser rparser = new RegionParser(gen);
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line = reader.readLine(); //Ignore first line
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(words.length>=1 && words[0].contains(":")){
	            	if(words[0].contains("-")){
	                	Region q = rparser.execute(words[0]);
	                	if(!regs.containsKey(q.getChrom()))
	                		regs.put(q.getChrom(), new ArrayList<Region>());
		            	regs.get(q.getChrom()).add(q);
	            	}
		        }
	        }reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return(regs);
	}

	/** 
	 * Simple method to write a list of regions to a file
	 * @param regionsOfInterest
	 * @param string
	 */
	public static void writeRegionsToFile(Map<String, List<Region>> regions, String outFileName) {
		try {
			FileWriter outFW = new FileWriter(outFileName);
			for(String chr : regions.keySet()){
				for(Region r : regions.get(chr))
					outFW.write(r.getLocationString()+"\n");
			}
			outFW.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
