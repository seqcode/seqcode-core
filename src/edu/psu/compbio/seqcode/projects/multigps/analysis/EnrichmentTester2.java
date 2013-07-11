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
import edu.psu.compbio.seqcode.projects.multigps.framework.EnrichmentSignificance;
import edu.psu.compbio.seqcode.projects.multigps.stats.CountsDataset;
import edu.psu.compbio.seqcode.projects.multigps.stats.MedianRatiosNormalization;
import edu.psu.compbio.seqcode.projects.multigps.stats.Normalization;
import edu.psu.compbio.seqcode.projects.multigps.utilities.Utils;

public class EnrichmentTester2 {

	protected Config config;
	protected ExperimentManager manager;
	protected Normalization normalizer;
	protected List<Point> potentialSites;
	protected List<Point> scalingSites;
	protected int siteJoinWin = 200;
	protected int searchRegionWin = 500;
	protected boolean simpleReadAssignment=false;
	
	//Constructor
	public EnrichmentTester2(Config config, ExperimentManager manager, List<Point> potentialSites, List<Point> scalingSites) {
		this.config = config;
		this.manager = manager;
		this.potentialSites = potentialSites;
		this.scalingSites = scalingSites;
		config.makeGPSOutputDirs(false);
	}
	
	
	/**
	 * Execute the enrichment tester
	 */
	public void execute(){
		potentialSites = cleanPotentialSites(potentialSites, siteJoinWin);
		scalingSites = cleanPotentialSites(scalingSites, siteJoinWin);
		//normalizer = new MedianRatiosNormalization(manager.getExperimentSet().getReplicates().size());
		
		
		//Convert our points to events
		PointsToEvents p2e = new PointsToEvents(config, manager, potentialSites, searchRegionWin,simpleReadAssignment);
		List<BindingEvent> events = p2e.execute();
		
		//Estimate signal fraction
		manager.estimateSignalProportion(events);
		
		/*//Get the scaling ratio from the scaling events if appropriate
		if(scalingSites.size()>0){
			PointsToEvents p2e_scaling = new PointsToEvents(config, manager, scalingSites, searchRegionWin, simpleReadAssignment);
			List<BindingEvent> scalingEvents = p2e_scaling.execute();
			//Convert events to a CountsDataset for the normalizer
			CountsDataset data = new CountsDataset(manager, scalingEvents, 0);
			normalizer.normalize(data);
		}else{
			//Convert events to a CountsDataset for the normalizer
			CountsDataset data = new CountsDataset(manager, events, 0);
			normalizer.normalize(data);
		}*/
		
		EnrichmentSignificance tester = new EnrichmentSignificance(config, manager.getExperimentSet(), events, config.getMinEventFoldChange(), config.getMappableGenomeLength());
		tester.execute();
		
		manager.setEvents(events);
		
		manager.writeFullEventFile(config.getOutputParentDir()+File.separator+config.getOutBase()+".all.events.table");
		manager.writeReplicateCounts(config.getOutputParentDir()+File.separator+config.getOutBase()+".replicates.counts");
		manager.writeBindingEventFiles(config.getOutputParentDir()+File.separator+config.getOutBase());
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
	
	//Accessors
	void setSiteJoinWin(int sjw){siteJoinWin = sjw;}
	void setSearchRegionWin(int srw){searchRegionWin = srw;}
	void setSimpleReadAssignment(boolean sra){simpleReadAssignment = sra;}
	
	//Main
	public static void main(String[] args){
		List<Point> potentialSites = new ArrayList<Point>();
		List<Point> scalingSites = new ArrayList<Point>();
		Config config = new Config(args, false);
		if(config.helpWanted()){
			System.err.println("EnrichmentTester:");
			System.err.println("\t--sites <potential site coords>\n" +
					"\t--scalingset <optional set in which to est. scaling ratio>\n" +
					"\t--win <window around events>\n" +
					"\t--joinwin <window to join sites>\n" +
					"\t--simple [assign all reads in window to event]\n");
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
			
			//Load scaling sites 
			Collection<String> scalingFiles = Args.parseStrings(args, "scalingset");
			for(String sf : scalingFiles)
				scalingSites.addAll(Utils.loadPointsFromFile(sf, config.getGenome()));
			
			int win = Args.parseInteger(args, "win", 200);
			int joinWin = Args.parseInteger(args, "joinwin", 500);
			boolean useSimpleAssignment = Args.parseFlags(args).contains("simple");
			
			EnrichmentTester2 tester = new EnrichmentTester2(config, manager, potentialSites, scalingSites);
			tester.setSiteJoinWin(joinWin);
			tester.setSearchRegionWin(win);
			tester.setSimpleReadAssignment(useSimpleAssignment);
			
			tester.execute();
			
			manager.close();
		}
	}
}
