package edu.psu.compbio.seqcode.projects.shaun.dataexplorer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeSet;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.projects.gps.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.utils.io.BackgroundModelIO;
import edu.psu.compbio.seqcode.projects.shaun.FreqMatrixImport;

/** DataSourceLoader
 *  Loads information about each DataSource from a file
 *  File format:
 *  
 *  MNAME	motif	Thres	MotifFile MotifBack
 *  MNAME	secmotif	Thres	MotifFile MotifBack
 *	EXPT	fileexpt	ExptType	Ctrl	Thres	Weight	FileFormat ExptRep1File...
 *	EXPT	dbexpt	ExptType	Ctrl	Thres	Weight	DBExptRep1ID...
 *	EXPT	rdbexpt	ExptType	Ctrl	Thres	Weight	RDBExptRep1ID...
 *  TNAME	cons	Thres Weight
 *  TNAME	uniqueness	Thres Weight
 */
		
public class DataSourceLoader {
	private Genome gen;
	private boolean usingNonUniqueReads=true;
	private int readLen=32;
	private int readExt = 0;
	
	public DataSourceLoader(Genome g, int readExt){gen=g; this.readExt = readExt;}

	//Conservation DataSources
	public Collection<DataSource> getConservationDataSource(String dataInfoFile){
		ArrayList<DataSource> consSources= new ArrayList<DataSource>();
		
		//Open the file and read for conservation table information
		for(String[] words : getLinesFromFile(dataInfoFile, "cons")){
			if(words.length>=4){
				String name = words[0];
				Double thres = new Double(words[2]);
				Double weight = new Double(words[3]);
				consSources.add(new ConsDataSource(gen, name, thres, weight));
			}
		}
		return consSources;
	}
	
	//Uniqueness DataSources
	public Collection<DataSource> getUniquenessDataSource(String dataInfoFile){
		ArrayList<DataSource> uniqSources= new ArrayList<DataSource>();
		
		//Open the file and read for conservation table information
		for(String[] words : getLinesFromFile(dataInfoFile, "uniqueness")){
			if(words.length>=4){
				String name = words[0];
				Double thres = new Double(words[2]);
				Double weight = new Double(words[3]);
				uniqSources.add(new UniquenessDataSource(gen, name, thres, weight));
			}
		}
		return uniqSources;
	}
	
	//DeepSeq DataSources
	public Collection<DataSource> getDeepSeqDataSource(String dataInfoFile){
		ArrayList<DataSource> seqSources= new ArrayList<DataSource>();
		
		HashMap<String, DeepSeqExpt> experiments = new HashMap<String, DeepSeqExpt>(); 
		ArrayList<String> IPList = new ArrayList<String>();
		HashMap<String,String> ctrlMap = new HashMap<String, String>();
		HashMap<String,Double> thresMap = new HashMap<String, Double>();
		HashMap<String,Double> weightMap = new HashMap<String, Double>();
		ArrayList<String[]> lines = new ArrayList<String[]>();
		lines.addAll(getLinesFromFile(dataInfoFile, "fileexpt"));
		lines.addAll(getLinesFromFile(dataInfoFile, "dbexpt"));
		lines.addAll(getLinesFromFile(dataInfoFile, "rdbexpt"));
		
		
		//Read the experiments
		for(String[] words : lines){
			if(words.length>=7){
				String eName = words[0];
				String sourceType=words[1];
				String eType=words[2];
				String eCtrl = words[3];
				Double thres = new Double(words[4]);
				Double weight = new Double(words[5]);
				
				if(eType.equals("IP")){
					IPList.add(eName);
					if(eCtrl.equals("NONE"))
						ctrlMap.put(eName,null);
					else
						ctrlMap.put(eName, eCtrl);
				}
				thresMap.put(eName, thres);
				weightMap.put(eName, weight);
				
				if(sourceType.equals("fileexpt")){
					String fileFormat = words[6]; 
					ArrayList<File> files = new ArrayList<File>();
					for(int i=7; i<words.length; i++)
						files.add(new File(words[i]));					
					DeepSeqExpt dse = new DeepSeqExpt(gen, files, usingNonUniqueReads, fileFormat, readLen);
					dse.setThreePrimeExt(readExt);
					experiments.put(eName, dse);
				}else{
					ArrayList<SeqLocator> locs = new ArrayList<SeqLocator>();
					for(int i=6; i<words.length; i++){
						String[] pieces = words[i].split(";");
	                	Set<String> reps = new TreeSet<String>();
		            	if (pieces.length == 2) {
		            		locs.add(new SeqLocator(pieces[0], reps, pieces[1]));					
			            } else if (pieces.length == 3) {
			            	reps.add(pieces[1]);
			            	locs.add(new SeqLocator(pieces[0], reps, pieces[2]));
			            } else {
			                throw new RuntimeException("Couldn't parse a SeqLocator from " + words[i]);
			            }
					}
					if(sourceType.equals("dbexpt")){
						DeepSeqExpt dse = new DeepSeqExpt(gen, locs, "db", readLen);
						dse.setThreePrimeExt(readExt);
						experiments.put(eName, dse);
					}else if(sourceType.equals("rdbexpt")){
						DeepSeqExpt dse = new DeepSeqExpt(gen, locs, "readdb", readLen);
						dse.setThreePrimeExt(readExt);
						experiments.put(eName, dse);
					}
				}
			}
		}
		
		//Make the DataSources
		for(String e : IPList){
			seqSources.add(new DeepSeqDataSource(experiments.get(e),experiments.get(ctrlMap.get(e)), e, thresMap.get(e), weightMap.get(e)));
		}
		return(seqSources);
	}
	
	//Motif DataSources
	public Collection<DataSource> getMotifDataSources(String dataInfoFile) throws IOException, ParseException {return getMotifDataSources(dataInfoFile, false);}
	public Collection<DataSource> getMotifDataSources(String dataInfoFile, boolean secondary) throws IOException, ParseException {
		ArrayList<DataSource> motifSources = new ArrayList<DataSource>();
		
		FreqMatrixImport motifImport;
		HashMap<String,MarkovBackgroundModel> backs = new HashMap<String,MarkovBackgroundModel>();
		HashMap<String, HashMap<String,WeightMatrix>> motifSets = new HashMap<String, HashMap<String,WeightMatrix>>();
		
		//Open the file and read for motif information
		ArrayList<String[]> lines = secondary ? getLinesFromFile(dataInfoFile, "secmotif") : getLinesFromFile(dataInfoFile, "motif");
		for(String[] words : lines){
			if(words.length>=5){
				//Get the background model
				if(!backs.containsKey(words[5])){
					MarkovBackgroundModel back = BackgroundModelIO.parseMarkovBackgroundModel(words[5], gen);
			    	backs.put(words[5], back);
				}
				//Get the motifs in the designated file
				if(!motifSets.containsKey(words[4])){
					HashMap<String,WeightMatrix> set = new HashMap<String, WeightMatrix>();
					motifImport = new FreqMatrixImport();
					motifImport.setBackground(backs.get(words[5]));
					for(WeightMatrix wm : motifImport.readTransfacMatrices(words[4])){
						set.put(wm.getName(), wm);
					}
					motifSets.put(words[4], set);
				}
				
				double threshold = new Double(words[2]).doubleValue();
				double weight = new Double(words[3]).doubleValue();
				String currMotifName = words[0];
				WeightMatrix motif = motifSets.get(words[4]).containsKey(currMotifName) ?
						motifSets.get(words[4]).get(currMotifName) : null;
							
				if(motif!=null){
					//Finally, load the MotifDataSource
					MotifDataSource mds = new MotifDataSource(motif,backs.get(words[5]), threshold, weight);
					motifSources.add(mds);
				}else{
					System.err.println("Cannot find "+currMotifName+" in file "+words[3]);
				}
			}
		}
		return(motifSources);
	}
	
	private ArrayList<String[]> getLinesFromFile(String fileName, String tag){
		ArrayList<String[]> lines = new ArrayList<String[]>();	
		//Open the file and read for motif information
		try {
			File eFile = new File(fileName);
			if(!eFile.isFile()){System.err.println("Invalid file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(eFile));
			String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\t");
	         
	            if(words.length>1 && words[1].equals(tag)){
	            	lines.add(words);
	            }
	        }
        } catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return(lines);
	}
}
