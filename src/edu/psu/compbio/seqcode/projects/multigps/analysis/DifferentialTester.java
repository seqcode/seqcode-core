package edu.psu.compbio.seqcode.projects.multigps.analysis;

import java.io.File;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.projects.multigps.analysis.PointsToEvents;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.framework.EnrichmentSignificance;
import edu.psu.compbio.seqcode.projects.multigps.stats.CountsDataset;
import edu.psu.compbio.seqcode.projects.multigps.stats.DifferentialEnrichment;
import edu.psu.compbio.seqcode.projects.multigps.stats.EdgeRDifferentialEnrichment;
import edu.psu.compbio.seqcode.projects.multigps.stats.Normalization;
import edu.psu.compbio.seqcode.projects.shaun.Utilities;

/**
 * Utility to quantify differential enrichment of a defined set of peaks across conditions. 
 * 
 * 
 * 	Input: 
 *		- Genome
 * 		- Peak file of potential peak locations
 * 		- Window around peaks in which to count reads
 * 	Output:
 * 		- Replicate count file
 * 		- Differential analysis
 *  
 * @author mahony
 *
 */
public class DifferentialTester {

	protected Config config;
	protected ExperimentManager manager;
	protected Normalization normalizer;
	protected List<Point> potentialSites;
	protected int searchRegionWin = 50;
	protected double qThres = 0.01;
	protected double minFold = 2;
	protected boolean simpleReadAssignment=true;
	protected String outDirName, outFileBase;
	protected File outDir;
	
	//Constructor
	public DifferentialTester(Config config, ExperimentManager manager, String peaksFileName, int win, double q, double minfold, String outRoot) {
		this.config = config;
		this.manager = manager;
		searchRegionWin = win;
		qThres = q;
		this.minFold = minfold;

		potentialSites = Utilities.loadPeaksFromPeakFile(config.getGenome(), peaksFileName, win);
		outFileBase=outRoot;		
		outDirName = outRoot;
		outDir = new File(outDirName);
		outDir.mkdir();
	}
	
	
	/**
	 * Execute the enrichment tester
	 */
	public void execute(){
		
		//Convert our points to events
		PointsToEvents p2e = new PointsToEvents(config, manager, potentialSites, searchRegionWin,simpleReadAssignment);
		List<BindingEvent> events = p2e.execute();
		
		//Estimate signal fraction (necessary for calculating statistics)
		manager.estimateSignalProportion(events);
				
		EnrichmentSignificance tester = new EnrichmentSignificance(config, manager.getExperimentSet(), events, config.getMinEventFoldChange(), config.getMappableGenomeLength());
		tester.execute();
		
		manager.setEvents(events);
		manager.writeReplicateCounts(outDirName+File.separator+outFileBase+".replicates.counts");
		
		//Statistical analysis: inter-condition differences
		if(manager.getNumConditions()>1 && config.getRunDiffTests()){
			DifferentialEnrichment edgeR = new EdgeRDifferentialEnrichment(config);
			CountsDataset data;
			File outImagesDir = new File(outDirName+File.separator+"images");
			outImagesDir.mkdir();
			
			for(int ref=0; ref<manager.getNumConditions(); ref++){
				data = new CountsDataset(manager, manager.getEvents(), ref);
				data = edgeR.execute(data);
				data.updateEvents(manager.getEvents(), manager);
				
				//Print MA scatters (inter-sample & inter-condition)
				data.savePairwiseConditionMAPlots(config.getDiffPMinThres(), outImagesDir.getAbsolutePath()+File.separator, true);
				
				//Print XY scatters (inter-sample & inter-condition)
				data.savePairwiseFocalSampleXYPlots(outImagesDir.getAbsolutePath()+File.separator, true);
				data.savePairwiseConditionXYPlots(manager, config.getDiffPMinThres(), outImagesDir.getAbsolutePath()+File.separator, true);
			}
		}
		manager.writeBindingEventFiles(outDirName+File.separator+outFileBase);
		
		System.out.println("Output files written to: "+outDirName);
	}
	
		
	//Main
	public static void main(String[] args){
		//Hack to set a default fixedpb limit. Need a better way to set per-application defaults
		String [] newargs = new String[args.length+2];
		for(int a=0; a<args.length; a++)
			newargs[a] = args[a];
		newargs[newargs.length-2] = "--fixedpb";
		newargs[newargs.length-1] = "1000";
		
		//Initialize Config
		Config config = new Config(newargs, false);
		config.setMedianScaling(true);
		config.setScalingSlidingWindow(50000);
		
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
			String siteFile = Args.parseString(args, "peaks", null);
			if(siteFile==null){
				printHelp(); System.exit(1);
			}
			int win = Args.parseInteger(args, "win", 50);
			double qThres = Args.parseDouble(args, "q", 0.01);
			double minFold = Args.parseDouble(args, "minfold", 2);
			DifferentialTester tester = new DifferentialTester(config, manager, siteFile, win, qThres, minFold, config.getOutName()); 
			
			tester.execute();
			
			manager.close();
		}
	}
	
	public static void printHelp(){
		System.err.println("DifferentialTester:\n" +
				"\tDifferential analysis of read counts in windows around a fixed set of peak positions using edgeR. \n");
		System.err.println("\t--gen <genome version>\n" +
				"\t--format <format of seq data files (SAM/IDX)>  default=SAM\n" +
				"\t--exptNAME-REP <file name of signal experiment, where NAME and REP are experiment name and replicate number>\n" +
				"\t--ctrlNAME-REP <file name of control experiment, where NAME and REP are experiment name and replicate number>\n" +
				"\t--peaks <potential site coords>\n" +
				"\t--win <window around events>  default=50bp\n" +
				"\t--q <Q-value minimum (corrected p-value)>  default=0.01\n" +
				"\t--minfold <min event fold-change>  default=2\n" +
				"\n" +
				"\t--design <experiment design file>  optional: can use design file instead of --expt, --ctrl and --format\n");
		System.err.println("\tExample Usage:\n" +
				"\tjava -Xmx2G edu.psu.compbio.seqcode.projects.multigps.analysis.DifferentialTester --gen mm10 --format IDX --exptFoxa2-r1 FoxA2_07-633_liver_-_-_-_XO111_kaz1-S001_Pugh40203mm10.tab --exptFoxa2-r2 FoxA2_07-633_liver_-_-_-_XO211_kaz1-S001_Pugh40205mm10.tab --ctrlFoxa2 IgG_12-370_liver_-_-_-_XO_kaz1-S001_Pugh4020mm10.idx --peaks foxa2.peaks --win 50 --q 0.01 --minfold 2\n" +
				"");
	}

}
