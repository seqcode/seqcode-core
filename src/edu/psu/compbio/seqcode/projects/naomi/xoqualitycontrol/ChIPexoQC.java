package edu.psu.compbio.seqcode.projects.naomi.xoqualitycontrol;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentScaler;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;

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
				
				if (printBinCounts ==true){
					
				}
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
				
				PrintWriter writer = new PrintWriter(rep.getName()+rep.getCondName()+rep.getRepName(),"UTF-8");
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
		 * You need to specify --fixedpb otherwise upper count limit would be set to a random number
		 * if you want to print binCounds --binCounts 
		 ***/
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
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
