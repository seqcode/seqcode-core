package org.seqcode.projects.sequtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.projects.multigps.features.BindingEvent;
import org.seqcode.projects.multigps.framework.BindingManager;
import org.seqcode.projects.multigps.framework.MultiGPSConfig;
import org.seqcode.projects.multigps.utilities.Utils;


public class ReplicateDataExporter {
	protected MultiGPSConfig config;
	protected ExperimentManager manager;
	protected BindingManager bindingManager;
	protected List<Point> potentialSites;
	protected int siteJoinWin = 200;
	protected int searchRegionWin = 500;
	
	//Constructor
	public ReplicateDataExporter(MultiGPSConfig config, ExperimentManager manager, BindingManager bindMan, List<Point> potentialSites) {
		this.config = config;
		this.manager = manager;
		this.bindingManager= bindMan;
		this.potentialSites = cleanPotentialSites(potentialSites, siteJoinWin);
	}

	
	
	/**
	 * Execute the enrichment tester
	 */
	public void execute(){
		//Convert our points to events
		PointsToEvents p2e = new PointsToEvents(config, manager, bindingManager, potentialSites, searchRegionWin, true);
		List<BindingEvent> events = p2e.execute();
		
		//Estimate signal-noise ratios
		bindingManager.estimateSignalVsNoiseFractions(events);
		
		bindingManager.setBindingEvents(events);
		
		bindingManager.writeReplicateCounts(config.getOutputParentDir()+File.separator+config.getOutBase()+".replicates.counts");
	}
	
	
	/**
	 * Sorts sites and combines nearby ones
	 * @param sites
	 * @return
	 */
	protected List<Point> cleanPotentialSites(List<Point> sites, int joinWin){
		List<Point> cleanSites = new ArrayList<Point>();
		System.err.println(sites.size()+" sites before cleaning.");
		
		Collections.sort(sites);
		Point lastSite = null;
		Point currSite = null;
		for(Point p : sites){
			if(lastSite == null){
				currSite = p;
			}else if(lastSite.getChrom().equals(p.getChrom()) && lastSite.distance(p)<=siteJoinWin){
				currSite = new Point(p.getGenome(), p.getChrom(), (p.getLocation()+lastSite.getLocation())/2);
			}else{
				cleanSites.add(lastSite);
				currSite=p;
			}
			lastSite = currSite;
		}
		if(lastSite!=null)
			cleanSites.add(lastSite);
		
		System.err.println(cleanSites.size()+" sites after cleaning.");
		return(cleanSites);
	}
	
	//Main
	public static void main(String[] args){
		List<Point> potentialSites = new ArrayList<Point>();
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		MultiGPSConfig config = new MultiGPSConfig(gcon, args, false);
		if(config.helpWanted()){
			System.err.println("ReplicateDataExporter:");
			System.err.println("\t--sites <potential site coords>\n" +
					"\t--scalingset <optional set in which to est. scaling ratio>");
			System.err.println(config.getArgsList());			
		}else{
			ExperimentManager manager = new ExperimentManager(econ);
			BindingManager bMan = new BindingManager(config, manager);
			//Just a test to see if we've loaded all conditions
			System.err.println("Conditions:\t"+manager.getConditions().size());
			for(ExperimentCondition c : manager.getConditions()){
				System.err.println("Condition "+c.getName()+":\t#Replicates:\t"+c.getReplicates().size());
			}
			for(ExperimentCondition c : manager.getConditions()){
				for(ControlledExperiment r : c.getReplicates()){
					System.err.println("Condition "+c.getName()+":\tRep "+r.getName());
					if(r.getControl()==null)
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount());
					else
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount()+"\tControl:\t"+r.getControl().getHitCount());
				}
			}
			
			//Load potential sites
			Collection<String> siteFiles = Args.parseStrings(args, "sites");
			for(String sf : siteFiles)
				potentialSites.addAll(Utils.loadPointsFromFile(sf, config.getGenome()));
			
			ReplicateDataExporter tester = new ReplicateDataExporter(config, manager, bMan, potentialSites);
			tester.execute();
			
			manager.close();
		}
	}
}
