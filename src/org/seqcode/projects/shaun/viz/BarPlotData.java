package org.seqcode.projects.shaun.viz;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class BarPlotData {
	private int cols=0, rows=0;
	public ArrayList<ArrayList<Double>> data= new ArrayList<ArrayList<Double>>();
	public ArrayList<String> colNames =new ArrayList<String>();
	public ArrayList<String> rowNames =new ArrayList<String>();
	public ArrayList<Double> colMaxVals =new ArrayList<Double>();
	public ArrayList<Double> colMinVals =new ArrayList<Double>();
	public ArrayList<Color> colMaxColors =new ArrayList<Color>();
	public ArrayList<Color> colMinColors =new ArrayList<Color>();

	public BarPlotData(){
		
	}
	
	//Accessor 
	public int getRows(){return rows;}
	public int getCols(){return cols;}
	
	//Read the data to be plotted from a tab-delimited file
	public void readData(String filename){

		boolean first=true;
		BufferedReader reader;
		try {
			File dFile = new File(filename);
			if(!dFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
        	reader = new BufferedReader(new FileReader(dFile));
		
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(first){
	            	for(int w=1; w<words.length; w++){
	            		colNames.add(words[w]);
	            		data.add(new ArrayList<Double>());
	            		if(w==1){
	            			colMaxColors.add(Color.yellow);
	            			colMinColors.add(Color.blue);
	            		}else{
	            			colMaxColors.add(new Color((int)(Math.random()*256),(int)(Math.random()*256),(int)(Math.random()*256)));
	            			colMinColors.add(Color.white);
	            		}
	            		cols++;
	            	}
	            	first=false;
	            }else{
	            	if(words.length>0){
	            		rowNames.add(words[0]);
	            		for(int w=1; w<words.length; w++){
	            			data.get(w-1).add(new Double(words[w]));
	            		}
	            		rows++;
	            	}
	            }
	        }
	        for(int c=0; c<cols; c++){
	        	colMaxVals.add(findMaxVal(data.get(c)));
	        	colMinVals.add(findMinVal(data.get(c)));System.out.println(findMaxVal(data.get(c))+"\t"+findMinVal(data.get(c)));
	        }
        } catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private double findMaxVal(ArrayList<Double> list){
		Double max = -1000000.0;
		for(Double d : list){
			if(d>max){max=d;}
		}
		return(max.doubleValue());
	}private double findMinVal(ArrayList<Double> list){
		Double min = +1000000.0;
		for(Double d : list){
			if(d<min){min=d;}
		}
		return(min.doubleValue());
	}
	
}
