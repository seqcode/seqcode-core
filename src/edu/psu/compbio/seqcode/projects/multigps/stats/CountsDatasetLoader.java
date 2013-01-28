package edu.psu.compbio.seqcode.projects.multigps.stats;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import Jama.Matrix;

import edu.psu.compbio.seqcode.gse.utils.Pair;

public class CountsDatasetLoader {

	protected static Matrix counts;
	protected static String [] units;
	protected static int [] design;
	protected static int sampleCount=0;
	protected static HashMap<Integer, Pair<String,String>> sampleToName = new HashMap<Integer, Pair<String,String>>();
	
	public static CountsDataset loadCountsDataFile(String filename){
		int numConditions=0;
		HashMap<String, Integer> condition2Index = new HashMap<String, Integer>();
		HashMap<Integer, String> index2Condition = new HashMap<Integer, String>();
		ArrayList<Integer> designList = new ArrayList<Integer>();
		
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        
	        //First line should have experiment labels (condition:replicate)
	        String line= reader.readLine();
	        String[] labels = line.split("\\s+");
	        boolean first=true;
	        for(String s : labels){
	        	if(!first){
	        		String cond="COND", rep="REP";
	        		String[] parts = s.split(":");
	        		cond = parts[0];
	        		if(parts.length>1)
	        			rep = parts[1];
	        		sampleToName.put(sampleCount, new Pair<String, String>(cond, rep));
	        		sampleCount++;
	        		
	        		if(!condition2Index.containsKey(cond)){
	        			condition2Index.put(cond, numConditions);
	        			index2Condition.put(numConditions, cond);
	        			numConditions++;
	        		}
	        		designList.add(condition2Index.get(cond));
	        	}else{
	        		first=false;
	        	}
	        }
	        
	        //Set up design array
	        design = new int[sampleCount];
	        int x=0;
	        for(Integer i : designList){
	        	design[x]=i;
	        	x++;
	        }
	        
	        //Load rest of file
	        ArrayList<Double[]> countArray = new ArrayList<Double[]>();
	        ArrayList<String> unitArray = new ArrayList<String>();
	        int numPoints = 0;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            unitArray.add(words[0]);
	            Double[] currCounts = new Double[sampleCount];
	            for(int w=1; w<words.length; w++)
	            	currCounts[w-1] = new Double(words[w]);
	            countArray.add(currCounts);
	            numPoints++;
	        }reader.close();
	        
	        //Set up counts array
	        counts = new Matrix(numPoints,sampleCount);
	        units = new String[numPoints];
	        for(int d=0; d<numPoints; d++){
	        	Double[] currCounts = countArray.get(d);
	        	for(int s=0; s<sampleCount; s++)
	        		counts.set(d,s,currCounts[s]);
	        	units[d]=unitArray.get(d);
	        }
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		CountsDataset data = new CountsDataset(counts, units, design, sampleToName, index2Condition);
		
		return data;
	}
}
