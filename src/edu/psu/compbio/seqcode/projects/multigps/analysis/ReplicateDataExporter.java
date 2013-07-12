package edu.psu.compbio.seqcode.projects.multigps.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.utilities.Utils;

public class ReplicateDataExporter {
	protected Config config;
	protected ExperimentManager manager;
	protected List<Point> potentialSites;
	protected int siteJoinWin = 200;
	protected int searchRegionWin = 500;
	
	//Constructor
	public ReplicateDataExporter(Config config, ExperimentManager manager, List<Point> potentialSites) {
		this.config = config;
		this.manager = manager;
		this.potentialSites = cleanPotentialSites(potentialSites, siteJoinWin);
	}

	
	
	/**
	 * Execute the enrichment tester
	 */
	public void execute(){
		//Convert our points to events
		PointsToEvents p2e = new PointsToEvents(config, manager, potentialSites, searchRegionWin, true);
		List<BindingEvent> events = p2e.execute();
		
		//Estimate signal-noise ratios
		manager.estimateSignalProportion(events);
		
		manager.setEvents(events);
		
		manager.writeReplicateCounts(config.getOutputParentDir()+File.separator+config.getOutBase()+".replicates.counts");
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
		Config config = new Config(args);
		if(config.helpWanted()){
			System.err.println("ReplicateDataExporter:");
			System.err.println("\t--sites <potential site coords>\n" +
					"\t--scalingset <optional set in which to est. scaling ratio>");
			System.err.println(config.getArgsList());			
		}else{
			ExperimentManager manager = new ExperimentManager(config);
			
			//Just a test to see if we've loaded all conditions
			ExperimentSet eset = manager.getExperimentSet();
			System.err.println("Conditions:\t"+eset.getConditions().size());
			for(ExperimentCondition c : eset.getConditions()){
				System.err.println("Condition "+c.getName()+":\t#Replicates:\t"+c.getReplicates().size());
			}
			for(ExperimentCondition c : eset.getConditions()){
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
			
			ReplicateDataExporter tester = new ReplicateDataExporter(config, manager, potentialSites);
			tester.execute();
			
			manager.close();
		}
	}
}
