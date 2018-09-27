package org.seqcode.projects.galaxyexo;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExperimentScaler;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;

/**
 * Utility to output a signal fraction using a normalization factor between signal and control experiments. 
 * It includes an option to output binned counts of signal and control experiments.
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
	protected List<Region> regionsToIgnore = new ArrayList<Region>(); 
	protected List<Region> regionsToCount = new ArrayList<Region>();
	protected Map<Sample, List<Float>> sampleWindowCounts = new HashMap<Sample, List<Float>>();
	
	public ChIPQC(GenomeConfig gcon, ExptConfig econ, ExperimentManager man){	
		gconfig = gcon;
		econfig = econ;
		manager = man;
		genome = econfig.getGenome();
	}	
	
	// setters
	public void setExcludeRegions(List<Region> regionsToIgnore){this.regionsToIgnore = regionsToIgnore;}
	 
	public void execute(int scalingWindowSize) throws FileNotFoundException, UnsupportedEncodingException{
		
		ExperimentScaler scaler = new ExperimentScaler();
		Genome genome = econfig.getGenome();
		
		List<Region> regs = new ArrayList<Region>();
		for(String chrom:genome.getChromList()) {
			int chrlen = genome.getChromLength(chrom);
			for (int start = 1; start  < chrlen - scalingWindowSize; start += scalingWindowSize) {
				Region r = new Region(genome, chrom, start, start + scalingWindowSize);
				regs.add(r);
			}
		}
		regionsToCount = filterExcluded(regs);
			
		for(ExperimentCondition exptCond: manager.getConditions()){	
			Map<Sample, List<Float>> sampleWindowCounts = new HashMap<Sample, List<Float>>();
			List<Sample> allSamples = new ArrayList<Sample>();
			List<Sample> signalSamples = new ArrayList<Sample>();
			List<Sample> controlSamples = new ArrayList<Sample>();
			
			for(ControlledExperiment rep : exptCond.getReplicates()){
				signalSamples.add(rep.getSignal());
				controlSamples.add(rep.getControl());
			}
			
			allSamples.addAll(signalSamples);
			allSamples.addAll(controlSamples);
			int listSize=0;
			for(Sample samp : allSamples){
				List<Float> currSampCounts = new ArrayList<Float>();
				for (Region r : regionsToCount)
					currSampCounts.add(samp.countHits(r));				
				sampleWindowCounts.put(samp, currSampCounts);
				listSize = currSampCounts.size();		
			}
		
			for(ControlledExperiment expt : exptCond.getReplicates()){
				if(expt.hasControl()){
					String scalingPlotFilename =null;
					if(econfig.getPlotScaling())
						scalingPlotFilename = expt.getName();
					scaler.scalingRatioByNCIS(sampleWindowCounts.get(expt.getSignal()), sampleWindowCounts.get(expt.getControl()), scalingPlotFilename,econfig.getNCISMinBinFrac());
				}
			}
		
			//Calculate scaling factor for pooled signal vs pooled control for this condition
			List<Float> pooledSignal = new ArrayList<Float>();
			List<Float> pooledControl = new ArrayList<Float>();
			for(int x=0; x<listSize; x++){
				float sumSig = 0;
				for(Sample s : signalSamples)
					sumSig+=sampleWindowCounts.get(s).get(x);
				pooledSignal.add(sumSig);
			
				float sumCtrl = 0;
				for(Sample s : controlSamples)
					sumCtrl+=sampleWindowCounts.get(s).get(x);
				pooledControl.add(sumCtrl);
			}
			
			double pooledSampleControlScaling = scaler.scalingRatioByNCIS(pooledSignal, pooledControl, null, econfig.getNCISMinBinFrac());
			double pooledsignal = exptCond.getTotalSignalCount();
			double pooledcontrl = exptCond.getTotalControlCount();
			double pooledIPstrength =1-(pooledSampleControlScaling/(pooledsignal/pooledcontrl));
			if (pooledIPstrength<0)
				pooledIPstrength=0;
			System.out.println("Condition:"+exptCond.getName()+"\tSignal:"+pooledsignal+"\tControl:"+pooledcontrl+"\tScalingFactor:"+pooledSampleControlScaling+"\tIPstrength: "+pooledIPstrength);
		}
	}
	
	public void printPairedBinCounts() throws FileNotFoundException, UnsupportedEncodingException{
		
		//Calculate scaling factors for each replicate's signal vs control		
		for(ExperimentCondition exptCond: manager.getConditions()){			
			for(ControlledExperiment rep : exptCond.getReplicates()){	
				if (rep.hasControl()){			
					PrintWriter writer = new PrintWriter(rep.getName()+rep.getCondName()+".counts.txt","UTF-8");
					writer.println("#signalBinCounts : controlBinCounts"+"\t"+"#"+rep.getCondName());
			
					List<Float> signalSampCounts = new ArrayList<Float>();
					List<Float> controlSampCounts = new ArrayList<Float>();
					signalSampCounts = sampleWindowCounts.get(rep.getSignal());
					controlSampCounts = sampleWindowCounts.get(rep.getControl());
					for (int i = 0; i < signalSampCounts.size(); i ++){
						if (signalSampCounts.get(i)+controlSampCounts.get(i)>0)
							writer.println(signalSampCounts.get(i)+":"+controlSampCounts.get(i));		
					}}}}}					
	
	public void printGenomeBins() throws FileNotFoundException{
		File outFile = new File(System.getProperty("user.dir")+File.separator+"genome_windows.bed");
		PrintWriter writer = new PrintWriter(outFile);
		for (Region r: regionsToCount)
			writer.write("chr"+r.getChrom()+"\t"+r.getStart()+"\t"+r.getEnd()+"\n");
		writer.close();
	}
	
	//Filter out pre-defined regions to ignore (e.g. tower regions)
    protected List<Region> filterExcluded(List<Region> testRegions) {
		List<Region> filtered = new ArrayList<Region>();
		if(regionsToIgnore.size()==0)
			return testRegions;
		
		for(Region t : testRegions){
			boolean ignore = false;
			for(Region i : regionsToIgnore){
				if(t.overlaps(i)){
					ignore = true; break;
				}
			}
			if(!ignore)
				filtered.add(t);
		}
		return filtered;
	}
		
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {	
		/***
		 * example: 
		 * java org.seqcode.projects.galaxyexo.ChIPQC --species "Homo sapiens;hg19" --expt expt.bam --ctrl ctrl.bam --format BAM --scalewin 5000 --exclude path-to-exclude-file --plotscaling --binCounts > NCIS.out
		 ***/
		ArgParser ap = new ArgParser(args);	
        if((!ap.hasKey("species") && !ap.hasKey("geninfo"))) { 
            System.err.println("Usage:\n" +
                               "ChIPQC\n" +
                               "\t--species <organism;genome> OR\n" +
                               "\t--geninfo <genome info file> AND --seq <path to seqs>\n" +
                               "\t--expt <signal experiment> \n" +
                               "\t--expt <control experiment> \n" +
	                       "\t--format <BAM/IDX/BED/etc> \n" +
                               "\nOPTIONS:\n" +
                               "\t--scalewin <window size for scaling procedure (default=10000)>\n" +
                               "\t--binCounts [flag to print bin counts] \n" +
            				   "\t--plotscaling [flag to plot diagnostic information for the chosen scaling method]\n" +
            				   "\t--printBins [flag to print genomic bin coordinates in bed file]\n" +
            				   "\t--exclude <file of regions to ignore>\n" +
                               "");
            System.exit(0);
        }
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);		
		econf.setPerBaseReadFiltering(false);		
		ExperimentManager manager = new ExperimentManager(econf);
		ChIPQC exoQC = new ChIPQC(gconf, econf, manager); 
		//Regions to ignore during scaling factor calculations
		if(ap.hasKey("exclude")){
			List<Region> regs = RegionFileUtilities.loadRegionsFromFile(Args.parseString(args, "exclude", null), gconf.getGenome(), -1);
			exoQC.setExcludeRegions(regs);
		}
		exoQC.execute(econf.getScalingSlidingWindow());			
		if (ap.hasKey("binCounts"))
			exoQC.printPairedBinCounts();
		if (ap.hasKey("printBins"))
			exoQC.printGenomeBins();
		manager.close();
	}
}
