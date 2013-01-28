package edu.psu.compbio.seqcode.projects.multigps.features;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.ewok.verbs.PointParser;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;

/**
 * Read a multi-condition GPS file
 * 
 * @author mahony
 *
 */
public class BindingEventFileReader {
	protected ExperimentSet experiments=null;
	protected Config config=null;
	ArrayList<BindingEvent> features = new ArrayList<BindingEvent>();
	
	public BindingEventFileReader(String fileName,  ExperimentSet ex, Config con){
		experiments = ex;
		config = con;
		BindingEvent.setExperimentSet(experiments);
		BindingEvent.setConfig(config);
		execute(fileName);
	}
	
	//Accessor
	public ArrayList<BindingEvent> getFeatures(){return features;}
	
	/**
	 * Read a multi-condition GPS file
	 * @param g
	 * @param fileName
	 * @return
	 */
	public ArrayList<BindingEvent> execute(String fileName){
		File f = new File(fileName);
		return(execute(f));		
	}
	public ArrayList<BindingEvent> execute(File f){
		features = new ArrayList<BindingEvent>();
		try {
			int numC = experiments.getConditions().size();
			HashMap<Integer, String> fileCondIndex = new HashMap<Integer, String>();
			
			BufferedReader reader = new BufferedReader(new FileReader(f));
			String line;
			while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	            String[] words = line.split("\\t");
	            
	            if(words[0].startsWith("#")){  //Header
	            	//Read Condition indices
	            	if(words[0].startsWith("#Condition") && !words[1].equals("Name")){
	            		fileCondIndex.put(new Integer(words[2]), words[1]);
	            	}
	            }else{ //Events
	            	PointParser pparser = new PointParser(config.getGenome());
	            	Point p = pparser.execute(words[0]);
	            	BindingEvent ev = new BindingEvent(p,null);
	            	
	            	int baseline = 1;
	            	
	            	//Load single condition stuff
	            	for(int currIndex=0; currIndex<numC; currIndex++){
	            		ExperimentCondition currCond = experiments.getNamedCondition(fileCondIndex.get(currIndex));
	            		
	            		ev.setCondSigHits(currCond, new Double(words[baseline+0]));
	            		ev.setCondCtrlHits(currCond, new Double(words[baseline+1]));
	            		ev.setCondSigVCtrlFold(currCond, new Double(words[baseline+2]));
	            		ev.setCondSigVCtrlP(currCond, Math.pow(10,-1*new Double(words[baseline+3])));
	            		baseline += BindingEvent.getNumSingleCondCols();
	            	}
	            	
	            	//Load inter-condition stuff
	            	for(int currIndexA=0; currIndexA<numC; currIndexA++){
	            		for(int currIndexB=0; currIndexB<numC; currIndexB++){
	            			if(currIndexA!=currIndexB){
	            				ExperimentCondition currCondA = experiments.getNamedCondition(fileCondIndex.get(currIndexA));
	            				ExperimentCondition currCondB = experiments.getNamedCondition(fileCondIndex.get(currIndexB));
	            				
	            				ev.setInterCondScMean(currCondA, currCondB, new Double(words[baseline+0]));
	            				ev.setInterCondFold(currCondA, currCondB, new Double(words[baseline+1]));
	            				ev.setInterCondP(currCondA, currCondB, Math.pow(10,-1*new Double(words[baseline+2])));
	            				baseline += BindingEvent.getNumInterCondCols();		
	            			}
	            		}
	            	}
	            	
	            	features.add(ev);
	            }	            
	            
			}reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return(features);
	}
}
