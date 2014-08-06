package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;
import edu.psu.compbio.seqcode.projects.multigps.framework.BackgroundDetector;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingManager;
import edu.psu.compbio.seqcode.projects.multigps.framework.MultiGPSConfig;
import edu.psu.compbio.seqcode.projects.multigps.framework.EnrichmentSignificance;
import edu.psu.compbio.seqcode.projects.sequtils.PointsToEvents;

public class BackgroundDetectingEnrichmentTester {

	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected MultiGPSConfig mgpsconfig;
	private List<BindingEvent> events;
	private double backgroundProportion;
	private int binWidth;
	private ExperimentManager manager;
	private BindingManager bindingManager;
	
	public BackgroundDetectingEnrichmentTester(GenomeConfig gcon, ExptConfig econ, MultiGPSConfig mgpscon, double backProp, int binW){
		backgroundProportion=backProp;
		binWidth = binW;
		this.gconfig = gcon;
		this.econfig = econ;
		this.mgpsconfig = mgpscon;
		manager = new ExperimentManager(econfig);
		bindingManager = new BindingManager(manager);
		BindingEvent.setExperimentManager(manager);
		BindingEvent.setConfig(mgpsconfig);
	}
	
	public void execute(String eventsFile){
		//Read the points from a file
		List<Point> pts = getPoints(eventsFile);
		
		//Calculate the expected number of background reads in a binWidth window
		double expectedBack = manager.getIndexedCondition(0).getTotalSignalCount() * 
								backgroundProportion / gconfig.getGenome().getGenomeLength() * binWidth;
		System.err.println("Expected number of background reads in "+binWidth+"bp window = "+expectedBack);
		
		//Quantify each of the potential events
		PointsToEvents eventConvertor = new PointsToEvents(mgpsconfig, manager, bindingManager, pts, binWidth, true);
		events = eventConvertor.execute();
		
		//Set background rates
		for(BindingEvent ev : events){
			ev.setIsFoundInCondition(manager.getIndexedCondition(0), true);
			ev.setRepCtrlHits(manager.getIndexedCondition(0).getIndexedReplicate(0), expectedBack);
			ev.setRepSigVCtrlFold(manager.getIndexedCondition(0).getIndexedReplicate(0), ev.getRepSigHits(manager.getIndexedCondition(0).getIndexedReplicate(0))/expectedBack);
			ev.setCondCtrlHits(manager.getIndexedCondition(0), expectedBack);
			ev.setCondSigVCtrlFold(manager.getIndexedCondition(0), ev.getCondSigHits(manager.getIndexedCondition(0))/expectedBack);
		}
		
		//EnrichmentSignificance
		EnrichmentSignificanceTesting signifTester = new EnrichmentSignificanceTesting(gconfig, econfig, events, 1, gconfig.getGenome().getGenomeLength());
		signifTester.execute();
		
		//Write files
		writeBindingEventFiles();
	}

	/**
     * Print all binding events to files
     */
    public void writeBindingEventFiles(){
    	if(events.size()>0){
    		
	    	try {
	    		//Full output table (all non-zero components)
	    		String filename = mgpsconfig.getOutBase()+".all.events.table";
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
	    			//Print
	    			String condName = cond.getName(); 
	    			condName = condName.replaceAll("/", "-");
	    			filename = mgpsconfig.getOutBase()+"_"+condName+".events";
					fout = new FileWriter(filename);
					fout.write(BindingEvent.conditionHeadString(cond)+"\n");
			    	for(BindingEvent e : events){
			    		double Q = e.getCondSigVCtrlP(cond);
			    		//Because of the ML step and component sharing, I think that an event could be assigned a significant number of reads without being "present" in the condition's EM model.
			    		if( Q <=mgpsconfig.getQMinThres())
			    			fout.write(e.getConditionString(cond)+"\n");
			    	}
					fout.close();
	    		}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
    }
	public List<Point> getPoints(String fileName){
		File f = new File(fileName);
		List<Point> points = new ArrayList<Point>();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(f));
			String line;
			while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	            String[] words = line.split("\\t");
	            
	            if(words[0].startsWith("#") || words[0].equals("chrom")){  //Header
	            	
	            }else{ //Events
	            	String chrom = words[0];
	            	chrom = chrom.replaceFirst("chr", "");
	            	int pos = new Integer(words[10]);
	            	points.add(new Point(gconfig.getGenome(), chrom, pos));
	            }	            
	            
			}reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return(points);
	}
	/**
	 * This main method is for testing the BackgroundDetectingEnrichmentTester
	 * @param args
	 */
	public static void main(String[] args){
		
		GenomeConfig gconfig = new GenomeConfig(args);
		ExptConfig econfig = new ExptConfig(gconfig.getGenome(), args);
		MultiGPSConfig mgpsconfig = new MultiGPSConfig(gconfig, args, false);
		if(econfig.helpWanted()){
			System.err.println("BackgroundDetectingEnrichmentTester:");
			System.err.println("Genome:" +
					"\t--species <Organism;Genome>\n" +
					"\tOR\n" +
					"\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" +
					"Experiments:\n" +
					"\t--expt <read file name> AND --format <SAM/BED/IDX/BOWTIE/NOVO/ELAND>\n" +
					"AND/OR" +
					"\t--rdbexpt <ReadDB experiment identifier>\n" +
					"Sites to test:\n" +
					"\t--sites <file name>\n" +
					"Other:\n" +
					"\t--threads <number of threads to use>\n" +
					"\t--binwidth <bin width>\n" +
					"\t--binstep <bin step>\n" +
					"\t--fixedpb <fixed per base limit>\n" +
					"\t--poissongausspb <filter per base using Poisson Gaussian sliding window>\n" +
					"");
		}else{
			int binW = Args.parseInteger(args,"binwidth", 50);
			int binS = Args.parseInteger(args,"binstep", 25);
			ExperimentManager manager = new ExperimentManager(econfig);
			
			BackgroundDetector detector = new BackgroundDetector(econfig, mgpsconfig, manager, binW, binS);
			HashMap<Sample, Double> backProps = detector.execute();
			double backProp = backProps.get(backProps.keySet().toArray()[0]);
			String pFile = Args.parseString(args, "sites", null);
			
			BackgroundDetectingEnrichmentTester bdet = new BackgroundDetectingEnrichmentTester(gconfig, econfig, mgpsconfig, backProp, binW);
			bdet.execute(pFile);
			manager.close();
		}
	}
}
