package edu.psu.compbio.seqcode.projects.multigps.utilities;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.gse.utils.RealValuedHistogram;
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
				System.err.println("Peak-motif distance histograms");
	    		String filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+".peaks2motifs.histo.txt";
	    		FileWriter fout = new FileWriter(filename);
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
			System.err.println("Peak-peak distance histograms (same condition)");
    		String filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+".intraCondPeakDistances.histo.txt";
    		FileWriter fout = new FileWriter(filename);
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
				System.err.println("Peak-peak distance histograms (different conditions)");
	    		String filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+".interCondPeakDistances.histo.txt";
	    		FileWriter fout = new FileWriter(filename);
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
	}
}
