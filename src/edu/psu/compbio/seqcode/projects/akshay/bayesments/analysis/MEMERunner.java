package edu.psu.compbio.seqcode.projects.akshay.bayesments.analysis;

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

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.BayesmentsConfig;
import edu.psu.compbio.seqcode.projects.multigps.stats.StreamGobbler;



public class MEMERunner {
	
	protected BayesmentsConfig config;
	protected Genome gen;
	protected String MEMEpath;
	protected String MEMEargs;
	protected Float pseudo = (float)0.001; 
	
	//Constructor
	public MEMERunner(BayesmentsConfig c){
		config=c;
		MEMEpath = config.getMEMEpath();
		MEMEargs = config.getMEMEargs();
	}
	
	/**
	 * Run MEME on the sequences, and store the best motifs (aligned with each other) in the ExperimentConditions
	 * @param sequences: List of Lists of sequences. Indexed by condition
	 */
	public Pair<List<WeightMatrix>,List<WeightMatrix>> execute(List<String> sequences, String memeOutDirName, boolean bestOnly){
		List<WeightMatrix> wm = new ArrayList<WeightMatrix>();
		List<WeightMatrix> fm = new ArrayList<WeightMatrix>();
		File workingDir = config.getOutputInterDir();
		String memeOutDir = workingDir+"/meme_"+memeOutDirName;
		try {
			//Set up the input file
			File seqFile= File.createTempFile("seq", ".fa", workingDir);
			String seqFilename = seqFile.getCanonicalPath();
			FileWriter fout = new FileWriter(seqFile);
			int sCount=1;
			for(String seq : sequences){
				fout.write(">Seq"+sCount+"\n"+seq+"\n");
				sCount++;
			}
			fout.close();
			
			//Test if meme directory exists. If it does, recursively delete contents
			File memeOutPath = new File(memeOutDir);
			if(memeOutPath.exists())
				config.deleteDirectory(memeOutPath);
			
			//Call MEME
			String MEMEcmd = MEMEpath+"/meme ";
			Process proc = Runtime.getRuntime().exec(MEMEcmd+" "+seqFilename+" "+MEMEargs +" -o "+memeOutDir);
			// any error message? 
			StreamGobbler errorGobbler = new 
			StreamGobbler(proc.getErrorStream(), "MEME_ERR", true);
			// any output? 
			StreamGobbler outputGobbler = new 
			StreamGobbler(proc.getInputStream(), "MEME_OUT", true);
			// kick them off 
			errorGobbler.start(); 
			outputGobbler.start(); 
			// any error??? 
			int exitVal = proc.waitFor(); 
			System.err.println("MEME ExitValue: " + exitVal);

            File memeOutFile = new File(memeOutDir+"/meme.txt");
            if (!memeOutFile.exists()) {
            	//Clean up intermediate files (fasta, etc)
                if(seqFile.exists())
                	seqFile.delete();
            	throw new FileNotFoundException("Can't find file " + memeOutFile.getName());
            }else{
	            BufferedReader memeReader = new BufferedReader(new FileReader(memeOutFile));
	            Map<String,Double> back = parseMEMEResultsForBackground(memeReader);
	            memeReader.close();
	            BufferedReader memeReader2 = new BufferedReader(new FileReader(memeOutFile));
	            List<Pair<WeightMatrix,Double>> currFM = parseMEMEResultsForFreqMatries(memeReader2);
	            memeReader2.close();
	            
	            //Extract best motif?
				if(bestOnly){
					double minScore = Double.MAX_VALUE; 
					WeightMatrix bestMotif=null; 
					for(Pair<WeightMatrix,Double> m : currFM)
						if(m.cdr()<minScore){
							minScore = m.cdr();
							bestMotif = m.car();
						}
					fm.add(bestMotif);
					WeightMatrix wMatrix = WeightMatrix.getLogOddsVersion(bestMotif, back);
					wm.add(wMatrix);
					//System.out.println(bestMotif.getName()+"\t"+minScore);
					//System.out.println(WeightMatrix.printMatrix(bestMotif));
				}else{
					for(Pair<WeightMatrix,Double> m : currFM){
						fm.add(m.car());
						WeightMatrix wMatrix = WeightMatrix.getLogOddsVersion(m.car(), back);
						wm.add(wMatrix);
					}
				}
			}
         
            //Clean up intermediate files (fasta, etc)
            if(seqFile.exists())
            	seqFile.delete();
            proc.destroy();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		return new Pair<List<WeightMatrix>,List<WeightMatrix>>(wm,fm);
	}
	
	protected Map<String,Double> parseMEMEResultsForBackground(BufferedReader memeOut){
		Map<String,Double> back = new HashMap<String,Double>();
		try {
			String line=memeOut.readLine();
			int lineno = 1; 
			while (line!=null && !line.matches(".*Background letter frequencies.*")) {
		      line = memeOut.readLine();
		      lineno++;
		    }
			line = memeOut.readLine();
			if(line!=null){
		    	try {
		    		String[] pieces = line.split("\\s+");
		    		double A = Double.parseDouble(pieces[1]);
		    		double C = Double.parseDouble(pieces[3]);
		    		double G = Double.parseDouble(pieces[5]);
		    		double T = Double.parseDouble(pieces[7]);
		    		back.put("A", A);
		    		back.put("C", C);
		    		back.put("G", G);
		    		back.put("T", T);
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
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return back;
	}
	
	protected List<Pair<WeightMatrix,Double>> parseMEMEResultsForFreqMatries(BufferedReader memeOut){
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
