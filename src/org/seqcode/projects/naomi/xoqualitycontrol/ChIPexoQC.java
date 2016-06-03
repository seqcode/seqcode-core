package org.seqcode.projects.naomi.xoqualitycontrol;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExperimentScaler;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.ArgParser;


public class ChIPexoQC {
	
	public boolean printBinCounts = false;
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ExperimentManager manager;
	
	public ChIPexoQC(GenomeConfig gcon, ExptConfig econ, ExperimentManager man){	
		gconfig = gcon;
		econfig = econ;
		manager = man;
	}
	
	public void setBinCounts(boolean printBin){printBinCounts = true;}
	 
	public void printQCMetrics(){
		
		double ncis, signalHits, controlHits;
		double IPstrength=0;
				
		for(ExperimentCondition exptCond: manager.getConditions()){
			for(ControlledExperiment rep : exptCond.getReplicates()){
				ncis = rep.getControlScaling();
				signalHits = rep.getSignal().getHitCount();
				controlHits = rep.getControl().getHitCount();
				IPstrength = 1-(ncis/(signalHits/controlHits));
				if (IPstrength<0)
					IPstrength=0;
				System.out.println("Condition:"+rep.getCondName()+"\tSignal:"+signalHits+"\tControl:"+controlHits+"\tScalingFactor:"+ncis+"\tIPstrength: "+IPstrength);

			}
		}
		manager.close();
	}
	
	public void printPairedBinCounts(int scalingWindowSize) throws FileNotFoundException, UnsupportedEncodingException{

		
		Genome genome = econfig.getGenome();
		Map<Sample, List<Float>> sampleWindowCounts = new HashMap<Sample, List<Float>>();
		List<Sample> allSamples = new ArrayList<Sample>();
		List<Sample> signalSamples = new ArrayList<Sample>();
		List<Sample> controlSamples = new ArrayList<Sample>();
			
		for(ExperimentCondition exptCond: manager.getConditions()){			
			for(ControlledExperiment rep : exptCond.getReplicates()){
				signalSamples.add(rep.getSignal());
				controlSamples.add(rep.getControl());
			}
			
		allSamples.addAll(signalSamples);
		allSamples.addAll(controlSamples);
		int listSize=0;
		for(Sample samp : allSamples){
			List<Float> currSampCounts = new ArrayList<Float>();
			for(String chrom:genome.getChromList()) {
				int chrlen = genome.getChromLength(chrom);
				for (int start = 1; start  < chrlen - scalingWindowSize; start += scalingWindowSize) {
					Region r = new Region(genome, chrom, start, start + scalingWindowSize);
					currSampCounts.add(samp.countHits(r));
		        }
		    }
			sampleWindowCounts.put(samp, currSampCounts);
			listSize = currSampCounts.size();
			}
		}
		for(ExperimentCondition exptCond: manager.getConditions()){
			
			for(ControlledExperiment rep : exptCond.getReplicates()){
				
				PrintWriter writer = new PrintWriter(rep.getName()+rep.getCondName()+".counts.txt","UTF-8");
				writer.println("#signalBinCounts : controlBinCounts"+"\t"+"#"+rep.getCondName());

				
				List<Float> signalSampCounts = new ArrayList<Float>();
				List<Float> controlSampCounts = new ArrayList<Float>();
				signalSampCounts = sampleWindowCounts.get(rep.getSignal());
				controlSampCounts = sampleWindowCounts.get(rep.getControl());
				for (int i = 0; i < signalSampCounts.size(); i ++){
					if (signalSampCounts.get(i)+controlSampCounts.get(i)>0)
						writer.println(signalSampCounts.get(i)+":"+controlSampCounts.get(i));				
				}
			}
		}
		
		manager.close();
	}
		
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
		
		
		/***
		 * example: 
		 * java org.seqcode.projects.naomi.xoqualitycontrol.ChIPexoQC --species "Homo sapiens;hg19" --design Ikaros20160328.design --scalewin 5000 --plotscaling >NCIS.out
		 *  
		 * if you want to print binCounds --binCounts 
		 * if you want to plot scaling --plotscaling
		 ***/
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
		
		econf.setPerBaseReadFiltering(false);
		
		ExperimentManager manager = new ExperimentManager(econf);
		ChIPexoQC exoQC = new ChIPexoQC(gconf, econf, manager); 
		
		exoQC.printQCMetrics();
		
		ArgParser ap = new ArgParser(args);
		
		if (ap.hasKey("binCounts")){
			exoQC.printPairedBinCounts(econf.getScalingSlidingWindow());
		}
		
		manager.close();
	}
}
