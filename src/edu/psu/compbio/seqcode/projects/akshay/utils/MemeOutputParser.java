package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class MemeOutputParser {
	
	protected static Float pseudo = (float)0.001; 
	
	public static void main(String[] args){
		try {
			ArgParser ap = new ArgParser(args);
			String meme_txt = Args.parseString(args, "memefile", "meme.txt");
			File memeOutFile = new File(meme_txt);
			if (!memeOutFile.exists()) {
				throw new FileNotFoundException("Can't find file " + memeOutFile.getName());
			}else{
				BufferedReader memeReader2 = new BufferedReader(new FileReader(memeOutFile));
	            List<Pair<WeightMatrix,Double>> currFM = parseMEMEResultsForFreqMatries(memeReader2);
	            memeReader2.close();
	            
	            for(int i=0; i< currFM.size(); i++){
	            	String out = WeightMatrix.printTransfacMatrix(currFM.get(i).car(),"Motif_"+Integer.toString(i));
	            	System.out.println(out);
	            }
	            
			}
			
		}catch (FileNotFoundException e) {
			e.printStackTrace();
		}catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static List<Pair<WeightMatrix,Double>> parseMEMEResultsForFreqMatries(BufferedReader memeOut){
		List<Pair<WeightMatrix,Double>> parsed = new ArrayList<Pair<WeightMatrix,Double>>();
		String line;
		int lineno = 1; 
		int motifCount=0;
		try {
			while((line = memeOut.readLine()) != null){
				while (line!=null && !line.matches(".*letter-probability matrix.*")) {
			      line = memeOut.readLine();
			      lineno++;
			    }
				if(line!=null){
					motifCount++;
				    String lenStr = line.replaceFirst("^.*w=\\s*", "");
				    lenStr = lenStr.replaceFirst("\\s*nsites=.*", "");
				    int length = Integer.parseInt(lenStr);
				    String EStr = line.replaceFirst("^.*E=\\s*", "");
				    double Eval = Double.parseDouble(EStr);
				    WeightMatrix matrix = new WeightMatrix(length);
				    matrix.setNameVerType("Motif"+motifCount, "freq", "MEME");
				    for (int i = 0; i < length; i++) {
				    	line = memeOut.readLine().replaceFirst("^\\s*", "");
				    	lineno++;
				    	try {
				    		String[] pieces = line.split("\\s+");
				    		float A = Float.parseFloat(pieces[0])+pseudo;
				    		float C = Float.parseFloat(pieces[1])+pseudo;
				    		float G = Float.parseFloat(pieces[2])+pseudo;
				    		float T = Float.parseFloat(pieces[3])+pseudo;
				    		float total = A+C+G+T;
				    		matrix.matrix[i]['A'] = A/total;
				    		matrix.matrix[i]['C'] = C/total;
				    		matrix.matrix[i]['G'] = G/total;
				    		matrix.matrix[i]['T'] = T/total;
				    	}
				    	catch (NumberFormatException ex) {
				    		System.err.println("At line " + lineno + ": " + line);
				    		ex.printStackTrace();
				    		throw ex;
				    	}
				    	catch (ArrayIndexOutOfBoundsException ex) {
				    		System.err.println("At line " + lineno + ": " + line);
				    		ex.printStackTrace();
				    		throw ex;
				    	}
				    }
				    matrix.setLogOdds();
				    parsed.add(new Pair<WeightMatrix,Double>(matrix, Eval));
				}
			}
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return parsed;
	}
}
