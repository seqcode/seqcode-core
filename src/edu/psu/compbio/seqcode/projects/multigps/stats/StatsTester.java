package edu.psu.compbio.seqcode.projects.multigps.stats;

import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;

/**
 * StatsTester: Tests the statistics classes
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class StatsTester {

	private CountsDataset data;
	private Normalization normalizer;
	private int focalCondition=0;
	private Config config;
	
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
					"\t--data <counts data file>\n" +
					"\t--norm <TMM/MR>\n" +
					"\t--focal <focal condition>\n" +
					"\t--diff <DE/EDGER>\n" +
					"\t--out <out name>\n");
		}else{
			if(ap.hasKey("data"))
				dataFile = ap.getKeyValue("data");
			if(ap.hasKey("norm"))
				norm = ap.getKeyValue("norm");
			if(ap.hasKey("focal"))
				focal = new Integer(ap.getKeyValue("focal"));
				
			Config conf = new Config(args, false);
			StatsTester tester = new StatsTester(conf, dataFile, norm, diff, focal);
		}
	}
	
	/**
	 * Constructor: initialize a statistics tester. 
	 * @param dataFile
	 * @param normMethod
	 * @param focalCond
	 */
	public StatsTester(Config conf, String dataFile, String normMethod, String diffMethod, int focalCond){
		config = conf;
		focalCondition = focalCond;
		data = CountsDatasetLoader.loadCountsDataFile(dataFile);
		data.setFocalCondition(focalCondition);
	
		if(diffMethod.equals("EDGER")){
			DifferentialEnrichment EdgeR = new EdgeRDifferentialEnrichment(config);
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
