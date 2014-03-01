package edu.psu.compbio.seqcode.projects.multigps.utilities;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.utils.RealValuedHistogram;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.motifs.MotifPlatform;

public class EventsPostAnalysis {

	protected Config config;
	protected ExperimentManager manager;
	protected MotifPlatform motifFinder = null;
	protected List<BindingEvent> events;
	protected double motifThres = 0.6;  // fraction of max score threshold
	
	public EventsPostAnalysis(Config c, ExperimentManager man, List<BindingEvent> ev, MotifPlatform mp){
		config = c;
		manager = man;
		events = ev;
		motifFinder = mp;
	}
	
	/**
	 * Run post-analysis of peaks.
	 * 	1) Histograms of peak-closestMotif distances
	 * 	2) Histograms of peak-peak distances (same condition)
	 * 	3) Histograms of peak-peak distances (inter-condition) 
	 */
	public void execute(int histoWin){
		System.err.println("Events post-analysis");
		String pcmfilename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+".peaks2motifs.histo.txt";
		String ppdscfilename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+".intraCondPeakDistances.histo.txt";
		String ppdicfilename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+".interCondPeakDistances.histo.txt";
		String htmlfilename = config.getOutputParentDir()+File.separator+"multiGPS_"+config.getOutBase()+"_results.html";
		
		//0) Set up hash map structure for events by chromosome
		List<HashMap<String,List<Integer>>> eventStruct = new ArrayList<HashMap<String,List<Integer>>>();
		for(int c=0; c<manager.getNumConditions(); c++){
			ExperimentCondition cond = manager.getExperimentSet().getIndexedCondition(c);
			eventStruct.add(new HashMap<String, List<Integer>>());
			for(String chr : config.getGenome().getChromList())
				eventStruct.get(c).put(chr, new ArrayList<Integer>());
			for(BindingEvent ev : events){
				double Q = ev.getCondSigVCtrlP(cond);
	    		if(ev.isFoundInCondition(cond) && Q <=config.getQMinThres()){
					String chr = ev.getPoint().getChrom();
					int loc = ev.getPoint().getLocation();
					eventStruct.get(c).get(chr).add(loc);
				}
			}
		}
		
		//1) Histograms of peak-closestMotif distances
		try {
			if(config.getFindingMotifs()){
				System.err.println("\tPeak-motif distance histograms");	    		
	    		FileWriter fout = new FileWriter(pcmfilename);
	    		fout.write("#Peaks to closest motifs distance histograms\n\n");
				for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
					if(cond.getMotif()!=null){
						fout.write("#Condition:"+cond.getName()+"\n");
						RealValuedHistogram peakMotifHisto = new RealValuedHistogram(0, histoWin, histoWin/5);
						double currThreshold = cond.getMotif().getMaxScore() * motifThres;
						for(BindingEvent ev : events){
							double Q = ev.getCondSigVCtrlP(cond);
				    		if(ev.isFoundInCondition(cond) && Q <=config.getQMinThres()){
								int loc = ev.getPoint().getLocation();
								if(ev.getContainingRegion()!=null){
									if(loc - ev.getContainingRegion().getStart() > histoWin && ev.getContainingRegion().getEnd()-loc >histoWin){
										double[] scores = motifFinder.scanRegionWithMotif(ev.getContainingRegion(), cond);
										int index = loc - ev.getContainingRegion().getStart();
										int closestMatch = Integer.MAX_VALUE;
										for(int x=0; x<scores.length; x++){
											if(scores[x]>=currThreshold && Math.abs(x-index)<closestMatch){
												closestMatch = Math.abs(x-index);
											}
										}
										peakMotifHisto.addValue(closestMatch);
									}
								}
							}
						}
						fout.write(peakMotifHisto.contentsToString()+"\n");
					}
				}
				fout.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//2) Histograms of peak-peak distances (same condition)
		try {
			System.err.println("\tPeak-peak distance histograms (same condition)");    		
    		FileWriter fout = new FileWriter(ppdscfilename);
    		fout.write("#Peaks to other peaks in same condition distance histograms\n\n");
			for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
				RealValuedHistogram peakPeakHisto = new RealValuedHistogram(0, histoWin, histoWin/5);
				fout.write("#Condition: "+cond.getName()+"\n");
				for(String chr : config.getGenome().getChromList()){
					List<Integer> currCondChrLocs = eventStruct.get(cond.getIndex()).get(chr);
					for(int x=0; x<currCondChrLocs.size(); x++){
						int xLoc = currCondChrLocs.get(x);
						int closestPeak = Integer.MAX_VALUE;
						for(int y=0; y<currCondChrLocs.size(); y++){ if(x!=y){
							int yLoc = currCondChrLocs.get(y);
							int dist = Math.abs(xLoc-yLoc);
							if(dist<closestPeak)
								closestPeak = dist;
						}}
						peakPeakHisto.addValue(closestPeak);	
					}
				}
				fout.write(peakPeakHisto.contentsToString()+"\n");
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//3) Histograms of peak-peak distances (inter-condition)
		try {
			if(manager.getNumConditions()>1){
				System.err.println("\tPeak-peak distance histograms (different conditions)");	    		
	    		FileWriter fout = new FileWriter(ppdicfilename);
	    		fout.write("#Peaks to peaks in other conditions distance histograms\n\n");
				for(ExperimentCondition condA : manager.getExperimentSet().getConditions()){
					for(ExperimentCondition condB : manager.getExperimentSet().getConditions()){if(condA != condB){
						RealValuedHistogram peakPeakHisto = new RealValuedHistogram(0, histoWin, histoWin/5);
						fout.write("#Condition: "+condA.getName()+" vs "+condB.getName()+"\n");
						for(String chr : config.getGenome().getChromList()){
							List<Integer> currCondChrLocsA = eventStruct.get(condA.getIndex()).get(chr);
							List<Integer> currCondChrLocsB = eventStruct.get(condB.getIndex()).get(chr);
							for(int x=0; x<currCondChrLocsA.size(); x++){
								int xLoc = currCondChrLocsA.get(x);
								int closestPeak = Integer.MAX_VALUE;
								for(int y=0; y<currCondChrLocsB.size(); y++){ if(x!=y){
									int yLoc = currCondChrLocsB.get(y);
									int dist = Math.abs(xLoc-yLoc);
									if(dist<closestPeak)
										closestPeak = dist;
								}}
								peakPeakHisto.addValue(closestPeak);	
							}
						}
						fout.write(peakPeakHisto.contentsToString()+"\n");
					}}
				}
				fout.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//4) HTML report
		try {
			System.err.println("Writing results report to: "+htmlfilename);
			
			//Write motif images
			HashMap<ExperimentCondition, String> motifImageNames = new HashMap<ExperimentCondition, String>();
			HashMap<ExperimentCondition, String> motifRCImageNames = new HashMap<ExperimentCondition, String>();
			if(config.getFindingMotifs()){
				for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
					if(cond.getMotif()!=null){
						String imName = config.getOutputImagesDir()+File.separator+config.getOutBase()+"_"+cond.getName()+"_motif.png";
						String imName2 = "images/"+config.getOutBase()+"_"+cond.getName()+"_motif.png";
						String motifLabel = cond.getName()+" motif, MEME";
						Utils.printMotifLogo(cond.getMotif(), new File(imName), 75, motifLabel);
						motifImageNames.put(cond,  imName2);
						WeightMatrix wm_rc = WeightMatrix.reverseComplement(cond.getMotif());
						imName = config.getOutputImagesDir()+File.separator+config.getOutBase()+"_"+cond.getName()+"_motif_rc.png";
						imName2 = "images/"+config.getOutBase()+"_"+cond.getName()+"_motif_rc.png";
						motifLabel = cond.getName()+" revcomp motif, MEME";
						Utils.printMotifLogo(wm_rc, new File(imName), 75, motifLabel);
						motifRCImageNames.put(cond,  imName2);
					}else{
						motifImageNames.put(cond,  null);
						motifRCImageNames.put(cond,  null);
					}
				}
			}
			
			//Build up the HTML file
			
			//Header and run information 
	    	FileWriter fout = new FileWriter(htmlfilename);
	    	fout.write("<html>\n" +
	    			"\t<head><title>MultiGPS results ("+config.getOutBase()+")</title></head>\n" +
	    			"\t<style type='text/css'>/* <![CDATA[ */ table, th{border-color: #600;border-style: solid;} td{border-color: #600;border-style: solid;} table{border-width: 0 0 1px 1px; border-spacing: 0;border-collapse: collapse;} th{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} td{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} /* ]]> */</style>\n" +
	    			"\t<script language='javascript' type='text/javascript'><!--\nfunction motifpopitup(url) {	newwindow=window.open(url,'name','height=75');	if (window.focus) {newwindow.focus()}	return false;}// --></script>\n" +
	    			"\t<script language='javascript' type='text/javascript'><!--\nfunction fullpopitup(url) {	newwindow=window.open(url,'name');	if (window.focus) {newwindow.focus()}	return false;}// --></script>\n" +
	    			"\t<body>\n" +
	    			"\t<h1>MultiGPS results ("+config.getOutBase()+")</h1>\n" +
	    			"");
	    	DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
	    	Date date = new Date();
	    	fout.write("\t<p>MultiGPS version "+config.version+" run completed on: "+dateFormat.format(date));
	    	fout.write(" with arguments:\n "+config.getArgs()+"\n</p>\n");
	    	
	    	
	    	//Binding event information (per condition)
	    	fout.write("\t<h2>Binding events</h2>\n" +
	    			"\t<table>\n");
	    	fout.write("\t\t<tr>" +
	    			"\t\t<th>Condition</th>\n" +
	    			"\t\t<th>Events</th>\n" +
	    			"\t\t<th>File</th>\n");
	    	if(config.getFindingMotifs())
	    		fout.write("\t\t<th>Positional Prior Motif</th>\n" +
	    				"\t\t<th>Motif Relative Offset</th>\n");
	    	fout.write("\t\t</tr>\n");
	    	for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
	    		String eventFileName=config.getOutBase()+"_"+cond.getName()+".events";
	    		fout.write("\t\t<tr>" +
		    			"\t\t<td>"+cond.getName()+"</td>\n" +
	    				"\t\t<td>"+manager.countEventsInCondition(cond)+"</td>\n" +
		    			"\t\t<td><a href='"+eventFileName+"'>"+eventFileName+"</a></td>\n");
		    	if(config.getFindingMotifs()){
		    		if(motifImageNames.get(cond)!=null)
		    			fout.write("\t\t<td><img src='"+motifImageNames.get(cond)+"'><a href='#' onclick='return motifpopitup(\""+motifRCImageNames.get(cond)+"\")'>rc</a></td>\n" +
		    					"\t\t<td>"+cond.getMotifOffset()+"</td>\n");
		    		else
		    			fout.write("\t\t<td>No motif found</td>\n" +
		    					"\t\t<td>NA</td>\n");
		    	}
		    	fout.write("\t\t</tr>\n");
			}fout.write("\t</table>\n");
			
	    	
	    	//Input data read counts and read distribs (per replicate)
	    	fout.write("\t<h2>Input data</h2>\n" +
	    			"\t<table>\n");
	    	fout.write("\t\t<tr>" +
	    			"\t\t<th>Replicate</th>\n" +
	    			"\t\t<th>ReadCount</th>\n" +
	    			"\t\t<th>CtrlScaling</th>\n" +
	    			"\t\t<th>SignalFraction</th>\n" +
	    			"\t\t<th>ReadDistributionModel</th>\n");
	    	fout.write("\t\t</tr>\n");
	    	for(ControlledExperiment rep : manager.getExperimentSet().getReplicates()){
	    		String replicateName = rep.getCondName()+"-"+rep.getRepName();
				String distribFilename = "images/"+config.getOutBase()+"_"+replicateName+"_Read_Distributions.png";
				String tmpscale = rep.hasControl()?String.format("%.3f",rep.getControlScaling()):"NA";
	    		fout.write("\t\t<tr>" +
		    			"\t\t<td>"+rep.getCondName()+" "+rep.getRepName()+"</td>\n" +
	    				"\t\t<td>"+rep.getSignal().getHitCount()+"</td>\n" +
	    				"\t\t<td>"+tmpscale+"</td>\n" +
	    				"\t\t<td>"+String.format("%.3f",rep.getSigProp())+"</td>\n");
	    		fout.write("\t\t<td><a href='#' onclick='return fullpopitup(\""+distribFilename+"\")'><img src='"+distribFilename+"' height='300'></a></td>\n");
	    		fout.write("\t\t</tr>\n");
			}fout.write("\t</table>\n");
	    	
			
	    	if(manager.getNumConditions()>1 && config.getRunDiffTests()){
	    		
	    		//Differential site count matrix
				fout.write("\t<h2>Differentially enriched events</h2>\n" +
						"\t<p>The table displays counts of events that are significantly enriched in the row condition with respect to the column condition. Click on the counts to get the differential event data file.</p>\n" +
						"\t<table>\n" +
						"\t\t<tr>\n" +
						"\t\t<th>Diff</th>\n");
				for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
					fout.write("\t\t<th>"+cond.getName()+"</th>\n");
				}fout.write("\t\t</tr>\n");
				for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
					fout.write("\t\t<tr>\n" +
							"\t\t<td>"+cond.getName()+"</td>\n");
					for(ExperimentCondition othercond : manager.getExperimentSet().getConditions()){
						if(cond.equals(othercond)){
							fout.write("\t\t<td>-</td>\n");
						}else{
							String filename = config.getOutBase()+"_"+cond.getName()+"_gt_"+othercond.getName()+".diff.events";
							fout.write("\t\t<td><a href='"+filename+"'>"+manager.countDiffEventsBetweenConditions(cond, othercond)+"</a></td>\n");
						}
					}fout.write("\t\t</tr>\n");
				}fout.write("\t</table>\n");
				
				//Differential scatterplots matrix
				fout.write("\t<h2>Differential enrichment scatterplots</h2>\n" +
						"\t<table>\n" +
						"\t\t<tr>\n" +
						"\t\t<th>Diff</th>\n");
				for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
					fout.write("\t\t<th>"+cond.getName()+"</th>\n");
				}fout.write("\t\t</tr>\n");
				for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
					fout.write("\t\t<tr>\n" +
							"\t\t<td>"+cond.getName()+"</td>\n");
					for(ExperimentCondition othercond : manager.getExperimentSet().getConditions()){
						if(cond.equals(othercond)){
							fout.write("\t\t<td> </td>\n");
						}else{
							String filename = "images/"+cond.getName()+"_vs_"+othercond.getName()+".XY.png";
							fout.write("\t\t<td><a href='#' onclick='return fullpopitup(\""+filename+"\")'><img src='"+filename+"' height='200'></a></td>\n");
						}
					}fout.write("\t\t</tr>\n");
				}fout.write("\t</table>\n");
				
				//Differential MA plots matrix
				fout.write("\t<h2>Differential enrichment MA plots</h2>\n" +
						"\t<table>\n" +
						"\t\t<tr>\n" +
						"\t\t<th>Diff</th>\n");
				for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
					fout.write("\t\t<th>"+cond.getName()+"</th>\n");
				}fout.write("\t\t</tr>\n");
				for(ExperimentCondition cond : manager.getExperimentSet().getConditions()){
					fout.write("\t\t<tr>\n" +
							"\t\t<td>"+cond.getName()+"</td>\n");
					for(ExperimentCondition othercond : manager.getExperimentSet().getConditions()){
						if(cond.equals(othercond)){
							fout.write("\t\t<td> </td>\n");
						}else{
							String filename = "images/"+cond.getName()+"_vs_"+othercond.getName()+".MA.png";
							fout.write("\t\t<td><a href='#' onclick='return fullpopitup(\""+filename+"\")'><img src='"+filename+"' height='200'></a></td>\n");
						}
					}fout.write("\t\t</tr>\n");
				}fout.write("\t</table>\n");
				
			}
			
			
			//File list of extras (histograms, etc)
			fout.write("\t<h2>Miscellaneous files</h2>\n");
			if(config.getFindingMotifs())
				fout.write("\t<p><a href='"+config.getOutBase()+".motifs'>Positional prior motifs.</a> Try inputting these motifs into <a href='http://www.benoslab.pitt.edu/stamp/'>STAMP</a> for validation.</p>\n");
			fout.write("\t<p><a href='intermediate-results/"+config.getOutBase()+".intraCondPeakDistances.histo.txt'>Peak-peak distance histograms (same condition)</a></p>\n");
			if(manager.getNumConditions()>1)
				fout.write("\t<p><a href='intermediate-results/"+config.getOutBase()+".interCondPeakDistances.histo.txt'>Peak-peak distance histograms (between conditions)</a></p>\n");
			if(config.getFindingMotifs())
				fout.write("\t<p><a href='intermediate-results/"+config.getOutBase()+".peaks2motifs.histo.txt'>Peak-motif distance histograms</a></p>\n");
	    	
	    	
	    	fout.write("\t</body>\n</html>\n");
	    	fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
}
