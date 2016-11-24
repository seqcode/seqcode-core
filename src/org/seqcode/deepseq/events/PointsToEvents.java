package org.seqcode.deepseq.events;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;


/**
 * PointsToEvents: this class converts a list of points into pseudo-events with tag counts.
 * 
 * @author mahony
 *
 */
public class PointsToEvents {
	private EventsConfig config;
	private ExperimentManager manager;
	private BindingManager bindingManager;
	private List<Point> points;
	private int regionWin;
	private boolean assignReadsWithModel=true;
	
	/**
	 * Constructor
	 * @param c
	 * @param man
	 * @param p
	 * @param win
	 */
	public PointsToEvents(EventsConfig c, ExperimentManager man, BindingManager bindMan, List<Point> p, int win, boolean assignWithModel){
		config = c;
		manager = man;
		bindingManager = bindMan;
		points = p;
		regionWin = win;
		assignReadsWithModel=assignWithModel;
	}
	
	/**
	 * Convert points to pseudo-events
	 * @return
	 */
	public List<BindingEvent> execute(){
		List<BindingEvent> events = new ArrayList<BindingEvent>();
		BindingEvent.setExperimentManager(manager);
		BindingEvent.setConfig(config);
		BindingEvent.setSortingCond(manager.getConditions().get(0));
		Collections.sort(points); //Sort for efficient experiment file cache loading
		//For each point
		for(Point p : points){
			//Expand point into region
			Region potentialReg = p.expand(regionWin/2);
			
			//Initialize binding event
			BindingEvent e = new BindingEvent(p,potentialReg);
			
			//for each condition
			for(ExperimentCondition c : manager.getConditions()){
				boolean[] addedToCond = new boolean[manager.getSamples().size()];
				for(int a=0; a<manager.getSamples().size(); a++)
					addedToCond[a]=false;

				//Signal hits for each replicate
				for(ControlledExperiment r : c.getReplicates()){
					//Get the hits for this replicate in this region
					List<StrandedBaseCount> sigHits = r.getSignal().getBases(potentialReg);
					List<StrandedBaseCount> ctrlHits=null;
					if(r.getControl()!=null)
						ctrlHits = r.getControl().getBases(potentialReg);
					
					double sigResp=0.0, ctrlResp=0.0;
					
					if(assignReadsWithModel){
						//Scan region with binding distribution to find ML position
						Point maxSigPoint = findMaxWithBindingModel(sigHits, potentialReg, bindingManager.getBindingModel(r));
						Point maxCtrlPoint = ctrlHits==null ? null : findMaxWithBindingModel(ctrlHits, potentialReg, bindingManager.getBindingModel(r));
						
						//ML assign reads (single binding event)
						sigResp = assignReadsSingleEvent(sigHits, maxSigPoint, bindingManager.getBindingModel(r));
						ctrlResp = ctrlHits==null ? 0 : assignReadsSingleEvent(ctrlHits, maxCtrlPoint, bindingManager.getBindingModel(r)); 
					}else{
						sigResp = simpleAssignReadsToEvent(sigHits, e.getPoint(), potentialReg);
						ctrlResp = ctrlHits==null ? 0 : simpleAssignReadsToEvent(ctrlHits, e.getPoint(), potentialReg);
					}
					
					//Set the replicate responsibilities
					e.setRepSigHits(r, sigResp);
					e.setRepCtrlHits(r, ctrlResp);
					
					//Increment the condition responsibilities
					//TODO: Is this the right way to integrate replicates?
					if(!addedToCond[r.getSignal().getIndex()])
						e.setCondSigHits(c, e.getCondSigHits(c)+sigResp);
					addedToCond[r.getSignal().getIndex()]=true;
					
					if(r.getControl()!=null){
						if(!addedToCond[r.getControl().getIndex()])
							e.setCondCtrlHits(c, e.getCondCtrlHits(c)+ctrlResp);
						addedToCond[r.getControl().getIndex()]=true;
					}else{
						e.setCondCtrlHits(c, 0.0);
					}
				}
				e.setIsFoundInCondition(c,true);
			}
			if(config.isAddingAnnotations())
				e.addClosestGenes();
			
			events.add(e);
		}
		return(events);
	}
	
	
	/* Find exact peak using a BindingModel */
	protected Point findMaxWithBindingModel(List<StrandedBaseCount> hits, Region coords, BindingModel model){
		int maxPos=0; double maxScore=0;
		double [] p = new double[coords.getWidth()+1];
		for(int k=0; k<=coords.getWidth(); k++){p[k]=0;}
		for(StrandedBaseCount x : hits){
			int readStart = x.getCoordinate();
			if(readStart>=coords.getStart()-model.getMax() && readStart<=coords.getEnd()+model.getMax()){
				int offset = readStart-coords.getStart();
				if(x.getStrand()=='+')
					for(int i=Math.max(model.getMin()+offset, 0); i<=Math.min(coords.getWidth(), offset+model.getMax()); i++)
						p[i]+=model.probability(i-offset) * x.getCount();
				else
					for(int i=Math.max(offset-model.getMax(), 0); i<=Math.min(coords.getWidth(), offset-model.getMin()); i++)
						p[i]+=model.probability(offset-i) * x.getCount();
			}
		}
		for(int k=0; k<=coords.getWidth(); k++)
			if(p[k]>maxScore){maxScore=p[k]; maxPos=k;}
			
		Point pt = new Point(coords.getGenome(), coords.getChrom(), coords.getStart()+maxPos);
		return(pt);
	}
	
	/**
	 * Assign all reads in the list that are within the binding model boundaries to the event at point
	 * @param hits
	 * @param point
	 * @param model
	 * @return
	 */
	protected double assignReadsSingleEvent(List<StrandedBaseCount> hits, Point point, BindingModel model){
		double count =0;
		for(StrandedBaseCount x : hits){
			int readStart = x.getCoordinate();
			if(readStart>=point.getLocation()-model.getMax() && readStart<=point.getLocation()+model.getMax()){
				count += x.getCount();
			}
		}
		return count;
	}
	
	/**
	 * Assign all reads in the list that are within the region to the event at point
	 * @param hits
	 * @param point
	 * @param region
	 * @return
	 */
	protected double simpleAssignReadsToEvent(List<StrandedBaseCount> hits, Point point, Region reg){
		double count =0;
		for(StrandedBaseCount x : hits){
			int readStart = x.getCoordinate();
			if(readStart>=reg.getStart() && readStart<=reg.getEnd()){
				count += x.getCount();
			}
		}
		return count;
	}
}
