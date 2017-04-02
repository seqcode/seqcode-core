package org.seqcode.deepseq.events;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;


/**
 * BindingManager stores lists of binding events and motifs associated with experiment conditions, 
 * and binding models associated with replicates. 
 * These data structures used to be under the relevant experiment components, but we moved them here to 
 * allow experiment loading to be independent.  
 * 
 * @author mahony
 *
 */
public class BindingManager {

	protected EventsConfig config;
	protected ExperimentManager manager;
	protected List<BindingEvent> events;
	protected Map<ExperimentCondition, List <BindingEvent>> conditionEvents;
	protected Map<ExperimentCondition, WeightMatrix> motifs;
	protected Map<ExperimentCondition, WeightMatrix> freqMatrices;
	protected Map<ExperimentCondition, Integer> motifOffsets;
	protected Map<ControlledExperiment, BindingModel> models;
	protected Map<ExperimentCondition, Integer> maxInfluenceRange;
	
	
	public BindingManager(EventsConfig con, ExperimentManager exptman){
		config = con;
		manager = exptman;
		events  = new ArrayList<BindingEvent>();
		conditionEvents = new HashMap<ExperimentCondition, List <BindingEvent>>();
		motifs = new HashMap<ExperimentCondition, WeightMatrix>();
		freqMatrices = new HashMap<ExperimentCondition, WeightMatrix>();
		motifOffsets = new HashMap<ExperimentCondition, Integer>();
		models = new HashMap<ControlledExperiment, BindingModel>();
		maxInfluenceRange = new HashMap<ExperimentCondition, Integer>();
		for(ExperimentCondition cond : manager.getConditions()){
			conditionEvents.put(cond, new ArrayList<BindingEvent>());
			motifOffsets.put(cond,0);
			maxInfluenceRange.put(cond,0);
		}
	}
	
	public List<BindingEvent> getBindingEvents(){return events;}
	public List<BindingEvent> getConditionBindingEvents(ExperimentCondition ec){return conditionEvents.get(ec);}
	public WeightMatrix getMotif(ExperimentCondition ec){return motifs.get(ec);}
	public WeightMatrix getFreqMatrix(ExperimentCondition ec){return freqMatrices.get(ec);}
	public Integer getMotifOffset(ExperimentCondition ec){return motifOffsets.get(ec);}
	public BindingModel getBindingModel(ControlledExperiment ce){return models.get(ce);}
	public Integer getMaxInfluenceRange(ExperimentCondition ec){return maxInfluenceRange.get(ec);}

	public void setBindingEvents(List<BindingEvent> e){events =e;}
	public void setConditionBindingEvents(ExperimentCondition ec, List<BindingEvent> e){conditionEvents.put(ec, e);}
	public void setMotif(ExperimentCondition ec, WeightMatrix m){motifs.put(ec, m);}
	public void setFreqMatrix(ExperimentCondition ec, WeightMatrix m){freqMatrices.put(ec, m);}
	public void setMotifOffset(ExperimentCondition ec, Integer i){motifOffsets.put(ec, i);}
	public void setBindingModel(ControlledExperiment ce, BindingModel mod){models.put(ce, mod);}
	public void updateMaxInfluenceRange(ExperimentCondition ec){
		int max=0; 
		for(ControlledExperiment rep : ec.getReplicates()){
			if(getBindingModel(rep).getInfluenceRange()>max)
				max=getBindingModel(rep).getInfluenceRange();
		}maxInfluenceRange.put(ec, max);
	}
	
	
	
	/**
	 * Print the events
	 * @param outRoot
	 */
	public void printConditionEvents(ExperimentCondition ec, String outRoot){
		try{
			String outName = outRoot+"."+ec.getName()+".events";
			if(config.getEventsFileTXTExtension())
				outName = outName+".txt";
			FileWriter fw = new FileWriter(outName);
			
			fw.write(BindingEvent.fullHeadString()+"\n");
			for(BindingEvent e : getConditionBindingEvents(ec)){
				fw.write(e.toString()+"\n");
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Print the replicate counts at events
	 * @param outRoot
	 */
	public void printReplicateCounts(ExperimentCondition ec, String outRoot){
		try{
			String outName = outRoot+"."+ec.getName()+".repcounts";
			FileWriter fw = new FileWriter(outName);
			
			fw.write(BindingEvent.repCountHeadString()+"\n");
			for(BindingEvent e : getConditionBindingEvents(ec)){
				fw.write(e.getRepCountString()+"\n");
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * For each controlled experiment, simply calculate the proportion of reads in the provided 
	 * list of binding events to everything else. 
	 * @param regs
	 */
	public void estimateSignalVsNoiseFractions(List<BindingEvent> signalEvents){
		for(ExperimentCondition c : manager.getConditions()){
			for(ControlledExperiment r : c.getReplicates()){
				double repSigCount =0, repNoiseCount=0;
				for(BindingEvent event : signalEvents){
					if(event.isFoundInCondition(c))
						repSigCount += event.getRepSigHits(r);
				}
				repNoiseCount = r.getSignal().getHitCount() - repSigCount;
				r.setSignalVsNoiseFraction(repSigCount/r.getSignal().getHitCount());
				System.err.println(r.getName()+"\t"+r.getIndex()+"\tsignal-noise ratio:\t"+String.format("%.4f",r.getSignalVsNoiseFraction()));
			}
		}
	}
	/**
	 * Count the binding events present in a given condition
	 * @param cond
	 * @return
	 */
	public int countEventsInCondition(ExperimentCondition cond, double qMinThres){
		int count=0;
		for(BindingEvent e : events){
			if(e.isFoundInCondition(cond) && e.getCondSigVCtrlP(cond) <=qMinThres)
				count++;
		}
		return count;
	}
	/**
	 * Count the differential binding events present in a given pair of conditions
	 * @param cond
	 * @return
	 */
	public int countDiffEventsBetweenConditions(ExperimentCondition cond, ExperimentCondition othercond, double qMinThres, double diffPMinThres){
		int count=0;
		for(BindingEvent e : events){
			if(e.isFoundInCondition(cond) && e.getCondSigVCtrlP(cond) <=qMinThres)
    			if(e.getInterCondP(cond, othercond)<=diffPMinThres && e.getInterCondFold(cond, othercond)>0)
    				count++;
		}
		return count;
	}
    /**
     * Print all binding events to files
     */
    public void writeBindingEventFiles(String filePrefix, double qMinThres, boolean runDiffTests, double diffPMinThres){
    	if(events.size()>0){
    		
	    	try {
	    		//Full output table (all non-zero components)
	    		String filename = filePrefix+".all.events.table";
	    		if(config.getEventsFileTXTExtension())
    				filename = filename+".txt";
	    		FileWriter fout = new FileWriter(filename);
	    		fout.write(BindingEvent.fullHeadString()+"\n");
	    		for(BindingEvent e : events)
	    			fout.write(e.toString()+"\n");
				fout.close();
	    		
	    		//Per-condition event files
	    		for(ExperimentCondition cond : manager.getConditions()){
	    			//Sort on the current condition
	    			BindingEvent.setSortingCond(cond);
	    			Collections.sort(events, new Comparator<BindingEvent>(){
	    	            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareBySigCtrlPvalue(o2);}
	    	        });
	    			String condName = cond.getName(); 
	    			condName = condName.replaceAll("/", "-");
	    			//Print events in MultiGPS format
	    			filename = filePrefix+"_"+condName+".events";
	    			if(config.getEventsFileTXTExtension())
	    				filename = filename+".txt";
					fout = new FileWriter(filename);
					fout.write(BindingEvent.conditionHeadString(cond)+"\n");
			    	for(BindingEvent e : events){
			    		double Q = e.getCondSigVCtrlP(cond);
			    		//Because of the ML step and component sharing, I think that an event could be assigned a significant number of reads without being "present" in the condition's EM model.
			    		if(e.isFoundInCondition(cond) && Q <=qMinThres)
			    			fout.write(e.getConditionString(cond)+"\n");
			    	}
					fout.close();
					//Print events in BED
					if(config.getPrintBED()){
						String bedfilename = filePrefix+"_"+condName+".bed";
		    			fout = new FileWriter(bedfilename);
						fout.write(BindingEvent.conditionBEDHeadString(cond)+"\n");
				    	for(BindingEvent e : events){
				    		double Q = e.getCondSigVCtrlP(cond);
				    		//Because of the ML step and component sharing, I think that an event could be assigned a significant number of reads without being "present" in the condition's EM model.
				    		if(e.isFoundInCondition(cond) && Q <=qMinThres)
				    			fout.write(e.getConditionBED(cond)+"\n");
				    	}
						fout.close();
					}
	    		}
	    		
	    		//Differential event files
	    		if(manager.getNumConditions()>1 && runDiffTests){
	    			for(ExperimentCondition cond : manager.getConditions()){
		    			//Sort on the current condition
		    			BindingEvent.setSortingCond(cond);
		    			Collections.sort(events, new Comparator<BindingEvent>(){
		    	            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareBySigCtrlPvalue(o2);}
		    	        });
		    			
		    			for(ExperimentCondition othercond : manager.getConditions()){
		    				if(!cond.equals(othercond)){
				    			//Print diff events
				    			String condName = cond.getName(); 
				    			String othercondName = othercond.getName(); 
				    			condName = condName.replaceAll("/", "-");
				    			othercondName = othercondName.replaceAll("/", "-");
				    			
				    			//MultiGPS format
				    			filename = filePrefix+"_"+condName+"_gt_"+othercondName+".diff.events";
				    			if(config.getEventsFileTXTExtension())
				    				filename = filename+".txt";
								fout = new FileWriter(filename);
								fout.write(BindingEvent.conditionShortHeadString(cond)+"\n");
						    	for(BindingEvent e : events){
						    		double Q = e.getCondSigVCtrlP(cond);
						    		//Because of the ML step and component sharing, I think that an event could be assigned a significant number of reads without being "present" in the condition's EM model.
						    		if(e.isFoundInCondition(cond) && Q <=qMinThres){
						    			if(e.getInterCondP(cond, othercond)<=diffPMinThres && e.getInterCondFold(cond, othercond)>0){
						    				fout.write(e.getConditionString(cond)+"\n");
						    			}
						    		}
						    	}
								fout.close();
								
								//Print events in BED
								if(config.getPrintBED()){
									filename = filePrefix+"_"+condName+"_gt_"+othercondName+".diff.bed";
									fout = new FileWriter(filename);
									fout.write(BindingEvent.diffConditionBEDHeadString(cond, othercond)+"\n");
							    	for(BindingEvent e : events){
							    		double Q = e.getCondSigVCtrlP(cond);
							    		//Because of the ML step and component sharing, I think that an event could be assigned a significant number of reads without being "present" in the condition's EM model.
							    		if(e.isFoundInCondition(cond) && Q <=qMinThres){
							    			if(e.getInterCondP(cond, othercond)<=diffPMinThres && e.getInterCondFold(cond, othercond)>0){
							    				fout.write(e.getConditionBED(cond)+"\n");
							    			}
							    		}
							    	}
									fout.close();
								}
		    				}
		    			}
		    		}
	    		}
	    		
	    		//If necessary, print per-replicate event base composition files
	    		if(config.getCalcEventBaseCompositions()){
	    			for(ExperimentCondition cond : manager.getConditions()){
	    				for(ControlledExperiment rep : cond.getReplicates()){
		    				//Print events
			    			String repName = rep.getName(); 
			    			repName = repName.replaceAll("/", "-");
			    			repName = repName.replaceAll(":", "-");
			    			filename = filePrefix+"_"+repName+".eventsbasecomps.txt";
							fout = new FileWriter(filename);
							fout.write("#Position\tEventRegA\tEventRegC\tEventRegG\tEventRegT\tEventTagA\tEventTagC\tEventTagG\tEventTagT\tBubbleRegA\tBubbleRegC\tBubbleRegG\tBubbleRegT\tBubbleTagA\tBubbleTagC\tBubbleTagG\tBubbleTagT\tBubbleIndex\n");
					    	for(BindingEvent e : events){
					    		double Q = e.getCondSigVCtrlP(cond);
					    		if(e.isFoundInCondition(cond) && Q <=qMinThres){
					    			float[] eventTagBases = e.getRepEventTagBases(rep);
					    			float[] bubbleTagBases = e.getRepBubbleTagBases(rep);
					    			float[] eventBases = e.getEventBases();
					    			float[] bubbleBases = e.getBubbleBases();
					    			float bubbleIndex = -1;
					    			float expectedT=0, expectedACG=0;
					    			if(bubbleBases[3]>0){
					    				expectedT = bubbleTagBases[3] / bubbleBases[3];
					    				float baseACG = bubbleBases[0]+bubbleBases[1]+bubbleBases[2];
					    				if(baseACG>0)
					    					expectedACG = (bubbleTagBases[0]+bubbleTagBases[1]+bubbleTagBases[2])/baseACG;
					    				if(expectedACG>0)
					    					bubbleIndex = expectedT/expectedACG;
					    			}
					    			fout.write(e.getPoint()+"\t"+
					    					eventBases[0]+"\t"+eventBases[1]+"\t"+eventBases[2]+"\t"+eventBases[3]+"\t"+
					    					eventTagBases[0]+"\t"+eventTagBases[1]+"\t"+eventTagBases[2]+"\t"+eventTagBases[3]+"\t"+
					    					bubbleBases[0]+"\t"+bubbleBases[1]+"\t"+bubbleBases[2]+"\t"+bubbleBases[3]+"\t"+
					    					bubbleTagBases[0]+"\t"+bubbleTagBases[1]+"\t"+bubbleTagBases[2]+"\t"+bubbleTagBases[3]+"\t"+
					    					bubbleIndex+
					    					"\n");
					    			
					    		}
					    	}
							fout.close();
	    				}
	    			}
	    		}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
    }
 
    /**
     * Print all binding events to files
     */
    public void writeFullEventFile(String filename){
    	if(events.size()>0){
    		try {
    			if(config.getEventsFileTXTExtension())
    				filename = filename+".txt";
	    		//Full dataset table
		    	FileWriter fout = new FileWriter(filename);
		    	fout.write(BindingEvent.fullHeadString()+"\n");
		    	for(BindingEvent e : events)
		    		fout.write(e.toString()+"\n");
				fout.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
    }
    
    /**
     * Print all motifs to file
     */
    public void writeMotifFile(String filename){
		try {
			if(config.getEventsFileTXTExtension())
				filename = filename+".txt";
    		//Full dataset table
	    	FileWriter fout = new FileWriter(filename);
	    	for(ExperimentCondition cond : manager.getConditions()){
	    		if(getFreqMatrix(cond)!=null)
	    			fout.write(WeightMatrix.printTransfacMatrix(getFreqMatrix(cond), cond.getName()));
	    	}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
 
    
    /**
     * Print replicate data counts to a file
     */
    public void writeReplicateCounts(String filename){
    	if(events.size()>0){
    		try {
	    		//Full dataset table
	    		FileWriter fout = new FileWriter(filename);
	    		fout.write(BindingEvent.repCountHeadString()+"\n");
	    		for(BindingEvent e : events)
	    			fout.write(e.getRepCountString()+"\n");
				fout.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
    }
    /**
     * Print all binding events to screen
     * TESTING ONLY
     */
    public void printBindingEvents(){
    	System.err.println(events.size()+" events found");
    	for(BindingEvent e : events){
    		System.err.println(e.toString());
    	}
    }
}
