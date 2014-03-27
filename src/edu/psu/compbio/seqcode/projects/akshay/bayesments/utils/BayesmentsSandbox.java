package edu.psu.compbio.seqcode.projects.akshay.bayesments.utils;

import java.io.File;

import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentManager;

/**
 * Utility class that contains some useful static methods
 * @author akshaykakumanu
 *
 */

public class BayesmentsSandbox {
	
	/**
	 * Prints a vector of doubles and also adds a header tag
	 * @param vec
	 * @param header_tag
	 */
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
	
	/**
	 * 
	 * @param vec
	 * @param row_tag
	 * @param column_tag
	 * @param column_names
	 */
	public static void printArray(double[][] vec, String row_tag, String column_tag, String[] column_names){
		String out ="";
		String column_header = "";
		int numRows = vec.length;
		int numCols = vec[0].length;
		
		for(int i=0; i< numCols; i++){
			if(column_tag == "MUc" || column_tag == "SIGMAc"){
				column_header = column_header + "\t"+column_names[i];
			}
			else if(column_tag == "MUf" || column_tag == "SIGMAf"){
				column_header = column_header + "\t"+"Factor-1";
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
	
	/**
	 * returns the index of the minimum double in the given list
	 * @param list
	 * @return
	 */
	public static int getMinindex(double[] list){
		double val = 100000.0;
		int ret=0;
		for(int i=0; i< list.length; i++){
			if(val > list[i]){
				ret = i;
				val = list[i];
			}
		}
		return ret;
	}
	
	public static Pair<Integer, Integer> getMinIndex(double[][] list){
		double val =  Double.MAX_VALUE;
		Pair<Integer,Integer> ret= new Pair<Integer, Integer>(0,0);
		int nrow = list.length;
		int ncol = list[0].length;
		for(int i=0; i<nrow; i++){
			for(int j=0; j<ncol; j++){
				if(list[i][j] < val){
					ret = new Pair<Integer, Integer>(i,j);
					val = list[i][j];
				}
			}
		}
		return ret;
	}
	
	/**
	 * returns the index of the maximum double in the given list
	 * @param list
	 * @return
	 */
	public static int getMaxindex(double[] list){
		double val = -100000.0;
		int ret=0;
		for(int i=0; i< list.length; i++){
			if(list[i] > val){
				ret =i;
				val = list[i];
			}
		}
		return ret;
	}
	
	public static boolean deleteDirectory(File path) {
	    if( path.exists() ) {
	      File[] files = path.listFiles();
	      for(int i=0; i<files.length; i++) {
	         if(files[i].isDirectory()) {
	           deleteDirectory(files[i]);
	         }
	         else {
	           files[i].delete();
	         }
	      }
	    }
	    return( path.delete() );
	}

}
