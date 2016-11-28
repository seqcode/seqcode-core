package org.seqcode.projects.galaxyexo;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import org.seqcode.data.io.parsing.GFFEntry;
import org.seqcode.data.io.parsing.ParseGFF;
import org.seqcode.deepseq.events.BindingEvent;
import org.seqcode.deepseq.events.BindingManager;
import org.seqcode.deepseq.events.EnrichmentSignificance;
import org.seqcode.deepseq.events.EventsConfig;
import org.seqcode.deepseq.events.PointsToEvents;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.gseutils.Args;
import org.seqcode.math.diff.CountsDataset;
import org.seqcode.math.diff.DifferentialEnrichment;
import org.seqcode.math.diff.EdgeRDifferentialEnrichment;
import org.seqcode.math.diff.Normalization;


/**
 * Utility to quantify enrichment and assess significance of a set of peaks (from GFF file). 
 * If multiple conditions are specified, differential binding analysis is performed via edgeR.
 * Meant to be used in the galaxy-exo pipeline.
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

	protected ExptConfig econfig;
	protected EventsConfig config;
	protected ExperimentManager manager=null;
	protected BindingManager bindingManager;
	protected Normalization normalizer;
	protected HashMap<Point,GFFEntry> siteToGFFEntry;
	protected List<Point> potentialSites;
	protected int searchRegionWin = 50;
	protected double qThres = 0.01;
	protected double minFold = 2;
	protected boolean simpleReadAssignment=true;
	protected String outDirName, outFileBase;
	protected File outDir;
	protected boolean rankByQ;
	
	//Constructor
	public SignificanceTester(ExptConfig econ, EventsConfig config, ExperimentManager manager, BindingManager bman, String gffFileName, int win, double q, double minfold, boolean rankByQ, String outDirName) {
		this.econfig = econ;
		this.config = config;
		this.manager = manager;
		this.bindingManager = bman;
		potentialSites = new ArrayList<Point>();
		siteToGFFEntry = new HashMap<Point, GFFEntry>();
		searchRegionWin = win;
		qThres = q;
		this.minFold = minfold;
		this.rankByQ = rankByQ;
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
		
		this.outDirName = outDirName==null?
				gffDir+"signif_w"+searchRegionWin+"_q"+String.format("%.2e", qThres)+"_minfold"+String.format("%.1f", minFold):
				outDirName+File.separator+"signif_w"+searchRegionWin+"_q"+String.format("%.2e", qThres)+"_minfold"+String.format("%.1f", minFold);
		outDir = new File(outDirName);
		outDir.mkdirs();
	}
	
	
	/**
	 * Execute the enrichment tester
	 */
	public void execute(){
		
		//Convert our points to events
		PointsToEvents p2e = new PointsToEvents(config, manager, bindingManager, potentialSites, searchRegionWin,!simpleReadAssignment);
		List<BindingEvent> events = p2e.execute();
		bindingManager.setBindingEvents(events);
		
		//Estimate signal fraction (necessary for calculating statistics)
		bindingManager.estimateSignalVsNoiseFractions(events);
		
		EnrichmentSignificance tester = new EnrichmentSignificance(config, manager, bindingManager, config.getMinEventFoldChange(), econfig.getMappableGenomeLength());
		tester.execute(searchRegionWin);
		
		bindingManager.writeReplicateCounts(outDirName+File.separator+outFileBase+"_replicatecounts.txt");
		writeBindingEventGFFFiles(outDirName+File.separator+outFileBase, events);
		
		//Statistical analysis: inter-condition differences
		if(manager.getNumConditions()>1 && config.getRunDiffTests()){
			DifferentialEnrichment edgeR = new EdgeRDifferentialEnrichment(config, outDir, outFileBase);
			CountsDataset data;
			File outImagesDir = new File(outDirName+File.separator+"images");
			outImagesDir.mkdir();
			
			for(int ref=0; ref<manager.getNumConditions(); ref++){
				data = new CountsDataset(manager, bindingManager.getBindingEvents(), ref);
				data = edgeR.execute(data);
				data.updateEvents(bindingManager.getBindingEvents(), manager);
				
				//Print MA scatters (inter-sample & inter-condition)
				data.savePairwiseConditionMAPlots(config.getDiffPMinThres(), outImagesDir.getAbsolutePath()+File.separator, true);
				
				//Print XY scatters (inter-sample & inter-condition)
				data.savePairwiseFocalSampleXYPlots(outImagesDir.getAbsolutePath()+File.separator, true);
				data.savePairwiseConditionXYPlots(manager, bindingManager, config.getDiffPMinThres(), outImagesDir.getAbsolutePath()+File.separator, true);
			}
		}
		
		System.out.println("Output files written to: "+outDirName);
	}
	
	/**
     * Print all binding events to files
     */
    public void writeBindingEventGFFFiles(String filePrefix, List<BindingEvent> events){
    	if(events.size()>0){
    		try {
	    		//Per-condition event files
	    		for(ExperimentCondition cond : manager.getConditions()){
	    			//Sort on the current condition
	    			BindingEvent.setSortingCond(cond);
	    			Collections.sort(events, new Comparator<BindingEvent>(){
	    	            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareBySigCtrlPvalue(o2);}
	    	        });
	    			//Print
	    			String condName = cond.getName(); 
	    			condName = condName.replaceAll("/", "-");
	    			String allFilename = filePrefix+"_all_"+condName+".gff";
	    			FileWriter allFout = new FileWriter(allFilename);
	    			String signifFilename = filePrefix+"_signif_"+condName+".gff";
	    			FileWriter signifFout = new FileWriter(signifFilename);
					for(BindingEvent e : events){
						Point currPt = e.getPoint();
		    			GFFEntry origGFF = siteToGFFEntry.get(currPt);
			    		double P = e.getCondSigVCtrlP(cond);
						double logP = Math.log(e.getCondSigVCtrlP(cond))/config.LOG2;
						double logF = Math.log(e.getCondSigVCtrlFold(cond))/config.LOG2;
						double sigHits = e.getCondSigHits(cond);
						double ctrlHits = e.getCondCtrlHits(cond);
						String attrib = origGFF.getAttribString()+String.format(";sig_tags=%.1f;ctrl_tags=%.1f;log2_fold_sigctrl=%.3f;log2_qval_sigctrl=%.3f", sigHits, ctrlHits, logF, logP);
			    		
						double score = sigHits;
						if(rankByQ)
							score = -logP;
			    		if(e.isFoundInCondition(cond) && P <=config.getQMinThres()){
			    			signifFout.write("chr"+origGFF.getChr()+"\t"+origGFF.getSource()+"\t"+origGFF.getFeature()+"\t"+origGFF.getStart()+"\t"
			    	    			+origGFF.getEnd()+"\t"+score+"\t"+origGFF.getStrand()+"\t"+origGFF.getFrame()
			    	    			+"\t"+attrib+"\n");
			    		}
			    		allFout.write("chr"+origGFF.getChr()+"\t"+origGFF.getSource()+"\t"+origGFF.getFeature()+"\t"+origGFF.getStart()+"\t"
		    	    			+origGFF.getEnd()+"\t"+score+"\t"+origGFF.getStrand()+"\t"+origGFF.getFrame()
		    	    			+"\t"+attrib+"\n");
			    	}
					signifFout.close();
					allFout.close();
	    		}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
    }
		
	//Main
	public static void main(String[] args){
		//Hack to set a default fixedpb limit. Need a better way to set per-application defaults
		String [] newargs = new String[args.length+2];
		for(int a=0; a<args.length; a++)
			newargs[a] = args[a];
		newargs[newargs.length-2] = "--fixedpb";
		newargs[newargs.length-1] = "1000";
		
		//Initialize MultiGPSConfig
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		EventsConfig config = new EventsConfig(gcon, args);
		econ.setScalingSlidingWindow(50000);
		
		if(config.helpWanted()){
			printHelp();	
		}else{
			ExperimentManager manager = new ExperimentManager(econ);
			BindingManager bman = new BindingManager(config, manager);

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
			
			//Arguments
			String siteFile = Args.parseString(args, "gff", null);
			if(siteFile==null)
				printHelp(); System.exit(1);
			String outDirName = Args.parseString(args, "outdir", null);
			int win = Args.parseInteger(args, "win", 50);
			double qThres = Args.parseDouble(args, "q", 0.01);
			double minFold = Args.parseDouble(args, "minfold", 2);
			boolean rankbyq = Args.parseFlags(args).contains("rankbyq");
			SignificanceTester tester = new SignificanceTester(econ, config, manager, bman, siteFile, win, qThres, minFold, rankbyq, outDirName); 
			
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
				"\t--outdir <output directory, default = based on input gff file directory name>\n" +
				"\t--win <window around events>  default=50bp\n" +
				"\t--q <Q-value minimum (corrected p-value)>  default=0.01\n" +
				"\t--minfold <min event fold-change>  default=2\n" +
				"\t--rankbyq [flag to rank by Q-value] default=rank by signal tags\n" +
				"\n" +
				"\t--design <experiment design file>  optional: can use design file instead of --expt, --ctrl and --format\n");
		System.err.println("\tOutput:\n" +
				"\t- GFF file containing all significant peak-pairs (i.e. passing Q-value and fold thresholds) with annotation.\n" +
				"\t- GFF file containing all input peak-pairs with annotation.\n" +
				"\t- Text file containing per-replicate signal tag counts for all peak-pairs.\n");
		System.err.println("\tExample Usage:\n" +
				"\tjava -Xmx2G org.seqcode.projects.galaxyexo.SignificanceTester --gen mm10 --format IDX --exptFoxa2-r1 FoxA2_07-633_liver_-_-_-_XO111_kaz1-S001_Pugh40203mm10.tab --exptFoxa2-r2 FoxA2_07-633_liver_-_-_-_XO211_kaz1-S001_Pugh40205mm10.tab --ctrlFoxa2 IgG_12-370_liver_-_-_-_XO_kaz1-S001_Pugh4020mm10.idx --gff genetrack_s5e10F1/cwpair_output_mode_f0u5d25b1/S_FoxA2_07-633_liver_-_-_-_XO_kaz1-S001_Pugh4020mm10_s5e10F1.gff --win 50 --q 0.01 --minfold 2\n" +
				"");
	}

}
