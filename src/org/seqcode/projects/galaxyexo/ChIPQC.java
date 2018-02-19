package org.seqcode.projects.galaxyexo;

import java.io.File;
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
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.gseutils.ArgParser;

/**
 * Utility to output a signal fraction calculated based on a normalization. It will also output a scaling factor, number of reads etc.
 * It includes an option to output hit counts of signal and control experiments in a given bin size in genome wide.
 * 
 * Input:
 * 		- Genome
 * 		- Signal experiment
 * 		- Control experiment
 * Output:
 * 		- A text file listing condition, signal hits, control hits, scaling factor, and a signal fraction.
 * 
 * @author naomi yamada
 */
public class ChIPQC {	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ExperimentManager manager;
	protected Genome genome;
	
	public ChIPQC(GenomeConfig gcon, ExptConfig econ, ExperimentManager man){	
		gconfig = gcon;
		econfig = econ;
		manager = man;
		genome = econfig.getGenome();
	}	
	 
	public void printQCMetrics(){		
		double ncis, signalHits, controlHits;
		double IPstrength=0;
				
		for(ExperimentCondition exptCond: manager.getConditions()){
			for(ControlledExperiment rep : exptCond.getReplicates()){
				if (!rep.hasControl()){
					System.err.println("Please provide a control experiment");
					System.exit(1);
				}				
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
	
	public void printGenomeBins(int scalingWindowSize) throws FileNotFoundException{
		File outFile = new File(System.getProperty("user.dir")+File.separator+"genome_windows.bed");
		PrintWriter writer = new PrintWriter(outFile);
		for(String chrom:genome.getChromList()) {
			int chrlen = genome.getChromLength(chrom);
			for (int start = 1; start  < chrlen - scalingWindowSize; start += scalingWindowSize) {
				Region r = new Region(genome, chrom, start, start + scalingWindowSize);
				writer.write("chrm"+chrom.toString()+"\t"+start+"\t"+start + scalingWindowSize+"\n");
			}
		}
		writer.close();
	}
		
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {	
		/***
		 * example: 
		 * java org.seqcode.projects.galaxyexo.ChIPQC --species "Homo sapiens;hg19" --expt expt.bam --ctrl ctrl.bam --format BAM --scalewin 5000 --plotscaling --binCounts > NCIS.out
		 ***/
		ArgParser ap = new ArgParser(args);	
        if((!ap.hasKey("species") && !ap.hasKey("geninfo"))) { 
            System.err.println("Usage:\n" +
                               "RegionCountSorter\n" +
                               "\t--species <organism;genome> OR\n" +
                               "\t--geninfo <genome info file> AND --seq <path to seqs>\n" +
                               "\t--expt <experiments> \n" +
                               "\nOPTIONS:\n" +
                               "\t--scalewin <window size for scaling procedure (default=10000)>\n" +
                               "\t--binCounts [flag to print bin counts] \n" +
            				   "\t--plotscaling [flag to plot diagnostic information for the chosen scaling method]\n" +
            				   "\t--printBins [flag to print genomic bin coordinates in bed file]\n" +
                               "");
            System.exit(0);
        }
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);		
		econf.setPerBaseReadFiltering(false);		
		ExperimentManager manager = new ExperimentManager(econf);
		ChIPQC exoQC = new ChIPQC(gconf, econf, manager); 		
		exoQC.printQCMetrics();			
		if (ap.hasKey("binCounts"))
			exoQC.printPairedBinCounts(econf.getScalingSlidingWindow());
		if (ap.hasKey("printBins"))
			exoQC.printGenomeBins(econf.getScalingSlidingWindow());
		manager.close();
	}
}
