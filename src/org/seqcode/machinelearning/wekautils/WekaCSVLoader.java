package org.seqcode.machinelearning.wekautils;

import weka.core.converters.C45Loader;
import weka.core.converters.CSVLoader;


/**
 * Directly implements Weka's CVS loader. Use the same options in Weka's documentation
 * @author akshaykakumanu
 *
 */
public class WekaCSVLoader {
	
	public static void main(String[] args){
		CSVLoader loader = new CSVLoader();
		CSVLoader.runFileLoader(loader, args);
	}
	

}
