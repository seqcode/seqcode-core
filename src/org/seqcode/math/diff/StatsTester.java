package org.seqcode.math.diff;

import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.TimeZone;

import org.seqcode.deepseq.events.EventsConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;

/**
 * StatsTester: Tests the statistics classes
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class StatsTester {

	protected EventsConfig config;
	private CountsDataset data;
	private Normalization normalizer;
	private int focalCondition=0;
	protected String outName="test", outBase="test";
	protected File outDir=null;
	
	
	/**
	 * Main: testing class for statistical methods
	 * @param args
	 */
	public static void main(String[] args) {
		String dataFile = "";
		String norm = "TMM";
		String diff = "EDGER";
		int focal=0;
		ArgParser ap = new ArgParser(args);
		if(args.length==0 || ap.hasKey("h") || !ap.hasKey("data")){
			System.err.println("StatsTester:\n" +
					"\t--out <output dir/name>\n"+
					"\t--rpath <path to the R bin dir (default: R is in $PATH). Note that you need to install edgeR separately>\n" +
					"\t--edgerod <EdgeR overdispersion (default=0.15)>\n" +
					"\t--data <counts data file>\n" +
					"\t--norm <TMM/MR>\n" +
					"\t--focal <focal condition>\n" +
					"\t--diff <DE/EDGER>\n" +
					"\t--out <out name>\n");
		}else{
			GenomeConfig gcon = new GenomeConfig(args);
			EventsConfig config = new EventsConfig(gcon, args);
			//Output path
			DateFormat df = new SimpleDateFormat("yyyy-MM-dd-hh-mm-ss");  
		    df.setTimeZone(TimeZone.getTimeZone("EST"));
			String outName = Args.parseString(args, "out", "test_"+df.format(new Date()));
			File outDir =  new File(outName); //Output directory
			String outBase = outDir.getName(); //Last part of name

			if(ap.hasKey("data"))
				dataFile = ap.getKeyValue("data");
			if(ap.hasKey("norm"))
				norm = ap.getKeyValue("norm");
			if(ap.hasKey("focal"))
				focal = new Integer(ap.getKeyValue("focal"));
			StatsTester tester = new StatsTester(config, dataFile, norm, diff, focal, outDir, outBase);
		}
	}
	
	/**
	 * Constructor: initialize a statistics tester. 
	 * @param dataFile
	 * @param normMethod
	 * @param focalCond
	 */
	public StatsTester(EventsConfig econ, String dataFile, String normMethod, String diffMethod, int focalCond, File outDir, String outBase){
		config = econ;
		focalCondition = focalCond;
		data = CountsDatasetLoader.loadCountsDataFile(dataFile);
		data.setFocalCondition(focalCondition);
		this.outDir = outDir;
		this.outBase = outBase;
		
		if(diffMethod.equals("EDGER")){
			DifferentialEnrichment EdgeR = new EdgeRDifferentialEnrichment(
					config, outDir, outBase);
			
			EdgeR.execute(data);
		}else if(diffMethod.equals("DE")){
			
			if(normMethod.equals("MR"))
				normalizer = new MedianRatiosNormalization(data.getNumSamples());
			else
				normalizer = new TMMNormalization(data.getNumSamples(), 0.3, 0.05);
			//normalizer.normalize(data);
			//normalizer.printPairwiseMAData(data);
			//normalizer.savePairwiseMAPlots(data, true);
			DifferentialEnrichment DESeq = new DESeqDifferentialEnrichment(normalizer);
			DESeq.execute(data);
		}
		
	}

}
