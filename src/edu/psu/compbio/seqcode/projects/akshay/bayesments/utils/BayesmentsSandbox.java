package edu.psu.compbio.seqcode.projects.akshay.bayesments.utils;

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
	
	public static void printArray(double[][] vec, String row_tag, String column_tag){
		String out ="";
		String row_header = "";
		String column_header = "";
		int numRows = vec.length;
		int numCols = vec[0].length;
		
		for(int i=0; i< numCols; i++){
			column_header = column_header + "\t"+column_tag+"_"+Integer.toString(i);
		}
		for(int i=0; i< numCols; i++){
			for(int j=0; j< numRows; j++){
				out = j==0 ? row_tag + "_"+Integer.toString(i)+"\t"+Double.toString(vec[i][j]): out+"\t"+Double.toString(vec[i][j]);
			}
			System.out.println(out);
		}
	}
	

}
