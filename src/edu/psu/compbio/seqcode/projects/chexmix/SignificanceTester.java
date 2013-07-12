package edu.psu.compbio.seqcode.projects.chexmix;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.io.parsing.GFFEntry;
import edu.psu.compbio.seqcode.gse.utils.io.parsing.ParseGFF;
import edu.psu.compbio.seqcode.projects.multigps.analysis.PointsToEvents;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.framework.EnrichmentSignificance;
import edu.psu.compbio.seqcode.projects.multigps.stats.Normalization;

/**
 * Utility to quantify enrichment and assess significance of a set of peaks (from GFF file). 
 * Meant to be used by Pugh lab pipeline.
 * 
 * 	Input: 
 *		- Genome
 * 		- GFF file of potential peak locations
 * 		- Window around peaks in which to count reads
 * 	Output:
 * 		- Annotated GFF file (all peaks)
 * 		- Annotated GFF file (peaks passing threshold)
 * 		- Replicate count file
 * 
 * SignificanceTester performs no cleaning, merging, or shifting of input locations.
 *  
 * @author mahony
 *
 */
public class SignificanceTester {

	protected Config config;
	protected ExperimentManager manager;
	protected Normalization normalizer;
	protected HashMap<Point,GFFEntry> siteToGFFEntry;
	protected List<Point> potentialSites;
	protected int searchRegionWin = 50;
	protected double qThres = 0.01;
	protected double minFold = 2;
	protected boolean simpleReadAssignment=true;
	protected String outDirName, outFileBase;
	protected File outDir;
	
	//Constructor
	public SignificanceTester(Config config, ExperimentManager manager, String gffFileName, int win, double q, double minfold) {
		this.config = config;
		this.manager = manager;
		potentialSites = new ArrayList<Point>();
		siteToGFFEntry = new HashMap<Point, GFFEntry>();
		searchRegionWin = win;
		qThres = q;
		this.minFold = minfold;
		
		//Load GFF file
		String gffDir="";
		try {
			File gffFile = new File(gffFileName);
			ParseGFF parser = new ParseGFF(gffFile);
			while(parser.hasNext()){
				GFFEntry site = parser.next();
				Point currPt = new Point(config.getGenome(), site.getChr(), site.getMidPoint());
				potentialSites.add(currPt);
				siteToGFFEntry.put(currPt, site);
			}
			if(gffFile.getParent() != null)
				gffDir = gffFile.getParent()+File.separator;
			outFileBase = gffFile.getName().replaceAll(".gff", "");
		} catch (IOException e) {
			//Silent exceptions
		}
		
		outDirName = gffDir+"signif_w"+searchRegionWin+"_q"+String.format("%e", qThres)+"_minfold"+String.format("%.1f", minFold);
		outDir = new File(outDirName);
		outDir.mkdir();
	}
	
	
	/**
	 * Execute the enrichment tester
	 */
	public void execute(){
		//normalizer = new MedianRatiosNormalization(manager.getExperimentSet().getReplicates().size());
		
		
		//Convert our points to events
		PointsToEvents p2e = new PointsToEvents(config, manager, potentialSites, searchRegionWin,simpleReadAssignment);
		List<BindingEvent> events = p2e.execute();
				
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
		
		manager.writeReplicateCounts(outDirName+File.separator+outFileBase+"_replicatecounts.txt");
	}
	
		
	//Main
	public static void main(String[] args){
		Config config = new Config(args, false);
		config.setMedianScaling(true);

		if(config.helpWanted()){
			printHelp();	
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
			
			//Arguments
			String siteFile = Args.parseString(args, "gff", null);
			if(siteFile==null){
				printHelp(); System.exit(1);
			}
			int win = Args.parseInteger(args, "win", 50);
			double qThres = Args.parseDouble(args, "q", 0.01);
			double minFold = Args.parseDouble(args, "minfold", 2);
			
			SignificanceTester tester = new SignificanceTester(config, manager, siteFile, win, qThres, minFold); 
			
			tester.execute();
			
			manager.close();
		}
	}
	
	public static void printHelp(){
		System.err.println("SignificanceTester:\n" +
				"\tCounts signal and control tags in windows around specified binding events, \n" +
				"\tnormalizes signal vs control using median scaling, and assesses the significance\n" +
				"\tof signal vs control enrichment using a Binomial test. \n" +
				"\tNote that median scaling is probably not appropriate in yeast ChIP-exo datasets.\n" +
				"\tAlternatives will be added soon. \n");
		System.err.println("\t--gen <genome version>\n" +
				"\t--format <format of seq data files (SAM/IDX)>  default=SAM\n" +
				"\t--exptNAME-REP <file name of signal experiment, where NAME and REP are experiment name and replicate number>\n" +
				"\t--ctrlNAME-REP <file name of control experiment, where NAME and REP are experiment name and replicate number>\n" +
				"\t--gff <potential site coords>\n" +
				"\t--win <window around events>  default=50bp\n" +
				"\t--q <Q-value minimum (corrected p-value)>  default=0.01\n" +
				"\t--minfold <min event fold-change>  default=2\n" +
				"\t--fixedpb <fixed per base limit>\n" +
				"\n" +
				"\t--design <experiment design file>  optional: can use design file instead of --expt, --ctrl and --format\n");
	}

}
