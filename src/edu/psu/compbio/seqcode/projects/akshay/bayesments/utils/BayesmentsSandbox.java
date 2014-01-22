package edu.psu.compbio.seqcode.projects.akshay.bayesments.utils;

import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentManager;

public class BayesmentsSandbox {
	
	public static void printArray(double[] vec, String header_tag){
		String out ="";
		String header = "";
		for(int i=0; i<vec.length; i++){
			header = i==0 ? header_tag+"_"+Integer.toString(i) : header + "\t" + header_tag+"_"+Integer.toString(i);
		}
		System.out.println(header);
		for(int i=0; i<vec.length; i++){
			out = i==0? Double.toString(vec[i]):out+"\t"+Double.toString(vec[i]);
		}
		System.out.println(out);
	}
	
	public static void printArray(double[][] vec, String row_tag, String column_tag, ExperimentManager manager){
		String out ="";
		String column_header = "";
		int numRows = vec.length;
		int numCols = vec[0].length;
		
		for(int i=0; i< numCols; i++){
			if(column_tag == "MUc" || column_tag == "SIGMAc"){
				column_header = column_header + "\t"+manager.getChromatinConditionList().get(i).getName();
			}
			else if(column_tag == "MUf" || column_tag == "SIGMAf"){
				column_header = column_header + "\t"+manager.getFacConditionList().get(i).getName();
			}
			else{
				column_header = column_header +"\t"+ column_tag+"_"+Integer.toString(i);
			}
		}
		System.out.println(column_header);
		for(int i=0; i< numRows; i++){
			for(int j=0; j< numCols; j++){
				out = j==0 ? row_tag + "_"+Integer.toString(i)+"\t"+Double.toString(vec[i][j]): out+"\t"+Double.toString(vec[i][j]);
			}
			System.out.println(out);
		}
	}
	

}
