package org.seqcode.projects.akshay.multigpstesting;

import java.io.*;
import java.util.*;

import org.seqcode.gse.tools.utils.Args;



public class CalculateMutualInfoContent {
	
	private List<Double> mi = new ArrayList<Double>();
	
	public CalculateMutualInfoContent(String filename) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String currentline = br.readLine();
		while(currentline != null && currentline.startsWith("D")){
			double tempMi=0;
			currentline = br.readLine();
			do{
				String[] pieces = currentline.split("\t");
				for(int i=1; i<pieces.length-1; i++){
					tempMi += (Double.parseDouble(pieces[i]) == 0.0 ? 0.0 : -1.0*Double.parseDouble(pieces[i])*Math.log(Double.parseDouble(pieces[i])));
				}
				currentline = br.readLine();
			}while(!currentline.startsWith("X"));
			currentline = br.readLine();
			mi.add(tempMi);
		}
		br.close();
	}

	public double getMax(){
		double ret = 0.0;
		for (double current: mi){
			if(current > ret){
				ret = current;
			}
		}
		return ret;
	}
	
	
	public static void main(String[] args) throws IOException{
		Set<String> commandlineArgs = Args.parseArgs(args);
		if(commandlineArgs.contains("help") || commandlineArgs.size() == 0){
			System.out.println("CalculateMutualInfoContent usage:\n" +
					"\t--motif <motif pwm file>\n");
		}
		else{
			Collection<String> motiffiles =  Args.parseStrings(args, "motif");
			String motiffilename = (String) motiffiles.toArray()[0];
			CalculateMutualInfoContent driver = new CalculateMutualInfoContent(motiffilename);
			System.out.println(driver.getMax());
		}
	}
	
	
}
