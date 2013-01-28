package edu.psu.compbio.seqcode.projects.multigps.utilities;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEvent;
import edu.psu.compbio.seqcode.projects.multigps.features.BindingEventFileReader;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;

public class BindingEventFilter {

	protected ExperimentSet experiments;
	protected Config config;
	protected List<BindingEvent> events;
	protected double sigThres=0.01;
	protected double maxLogFold=8.0, minLogFold=-8.0;
	
	//Constructor
	public BindingEventFilter(ExperimentSet e, Config c, List<BindingEvent> ev, double st){
		experiments = e;
		config = c;
		events = ev;
		sigThres = st;
	}
	
	/**
	 * Filter the events a number of different ways
	 */
	public void execute(){
		/*//TEST
		try{
			String testName = "TEST.events";
			FileWriter fw1 = new FileWriter(testName);
			fw1.write(BindingEvent.headString()+"\n");
			for(BindingEvent e : events){
				fw1.write(e.toString()+"\n");
			}
			fw1.close();
		} catch (IOException e) {
			e.printStackTrace();
		}*/
		
		for(ExperimentCondition c1 : experiments.getConditions()){
			for(ExperimentCondition c2 : experiments.getConditions()){
				if(c1!=c2){
					try{
						//1) List of sites enriched in one condition more than another
						String outName = c1.getName()+"_gt_"+c2.getName()+".events";
						FileWriter fw = new FileWriter(outName);
						fw.write(BindingEvent.fullHeadString()+"\n");
						for(BindingEvent e : events){
							if(e.getCondSigVCtrlP(c1)<=sigThres && e.getInterCondP(c1, c2)<=sigThres)
								fw.write(e.toString()+"\n");
						}
						fw.close();
		
						//2) List of sites enriched in either condition
						outName = c1.getName()+"_or_"+c2.getName()+".events";
						fw = new FileWriter(outName);
						fw.write(BindingEvent.fullHeadString()+"\n");
						for(BindingEvent e : events){
							if(e.getCondSigVCtrlP(c1)<=sigThres || e.getCondSigVCtrlP(c2)<=sigThres)
								fw.write(e.toString()+"\n");
						}
						fw.close();
						
						//3) List of sites enriched similarly in both conditions
						outName = c1.getName()+"_sim_"+c2.getName()+".events";
						fw = new FileWriter(outName);
						fw.write(BindingEvent.fullHeadString()+"\n");
						for(BindingEvent e : events){
							if(e.getCondSigVCtrlP(c1)<=sigThres && e.getCondSigVCtrlP(c2)<=sigThres && e.getInterCondP(c1, c2)>sigThres && e.getInterCondFold(c1, c2)<2)
								fw.write(e.toString()+"\n");
						}
						fw.close();
						
						//4) Scatter plot data
						String outNameA = c1.getName()+"_vs_"+c2.getName()+".unscaled.scatter";
						FileWriter fwAll = new FileWriter(outNameA);
						String outNameB = c1.getName()+"_vs_"+c2.getName()+".outliers.unscaled.scatter";
						FileWriter fwOut = new FileWriter(outNameB);
						for(BindingEvent e : events){
							if(e.getCondSigVCtrlP(c1)<=sigThres || e.getCondSigVCtrlP(c2)<=sigThres)
								fwAll.write(e.getPoint().getLocationString()+"\t"+e.getCondSigHits(c1)+"\t"+e.getCondSigHits(c2)+"\n");
							if(e.getCondSigVCtrlP(c1)<=sigThres && e.getInterCondP(c1, c2)<=sigThres)
								fwOut.write(e.getPoint().getLocationString()+"\t"+e.getCondSigHits(c1)+"\t"+e.getCondSigHits(c2)+"\n");
							if(e.getCondSigVCtrlP(c2)<=sigThres && e.getInterCondP(c2, c1)<=sigThres)
								fwOut.write(e.getPoint().getLocationString()+"\t"+e.getCondSigHits(c1)+"\t"+e.getCondSigHits(c2)+"\n");
						}
						fwAll.close();
						fwOut.close();
						
						//5) Mean/fold data
						outNameA = c1.getName()+"_vs_"+c2.getName()+".scaledmean.meanfold";
						fwAll = new FileWriter(outNameA);
						outNameB = c1.getName()+"_vs_"+c2.getName()+".outliers.scaledmean.meanfold";
						fwOut = new FileWriter(outNameB);
						for(BindingEvent e : events){
							double logFold = e.getInterCondFold(c1, c2)==0? minLogFold : Math.log(e.getInterCondFold(c1, c2))/Math.log(2);
							if(logFold<minLogFold)								
								logFold=minLogFold;
							if(logFold>maxLogFold)
								logFold=maxLogFold;
							if(e.getCondSigVCtrlP(c1)<=sigThres || e.getCondSigVCtrlP(c2)<=sigThres)
								fwAll.write(e.getPoint().getLocationString()+"\t"+e.getInterCondScMean(c1, c2)+"\t"+logFold+"\n");
							if(e.getCondSigVCtrlP(c1)<=sigThres && e.getInterCondP(c1, c2)<=sigThres)
								fwOut.write(e.getPoint().getLocationString()+"\t"+e.getInterCondScMean(c1, c2)+"\t"+logFold+"\n");
							if(e.getCondSigVCtrlP(c2)<=sigThres && e.getInterCondP(c2, c1)<=sigThres)
								fwOut.write(e.getPoint().getLocationString()+"\t"+e.getInterCondScMean(c1, c2)+"\t"+logFold+"\n");
						}
						fwAll.close();
						fwOut.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}
		}
	}
	
	//Main
	public static void main(String[] args){
		Config config = new Config(args);
		if(config.helpWanted()){
			System.err.println("BindingEventFilter:");
			System.err.println("\t--events <event file>");
			System.err.println("\t--sigthres <significance threshold (log)>");
			System.err.println(config.getArgsList());			
		}else{
			ExperimentManager manager = new ExperimentManager(config, false);
			ExperimentSet eset = manager.getExperimentSet();
			
			String eFile = Args.parseString(args, "events", null);
			BindingEventFileReader reader = new BindingEventFileReader(eFile, eset, config);
			List<BindingEvent> events = reader.execute(eFile);
			
			double sigThres = Args.parseDouble(args, "sigthres", 0.01);
			
			BindingEventFilter filter = new BindingEventFilter(eset, config, events, sigThres);
			filter.execute();
			
			manager.close();
		}
	}
}
