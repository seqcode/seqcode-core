package edu.psu.compbio.seqcode.projects.seed;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.TimeZone;

import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.projects.multigps.utilities.Utils;

public class SEEDConfig {

	protected GenomeConfig gconfig;
	//Settable variables
	protected String outName="seed", outBase="seed"; //output names
	protected File outDir=null;
	protected int tagShift=0;			//Shift tags this length in the 3' direction
	protected int tag3PrimeExtension=0; //Extend tags this number of bp in the 3' direction (starting from 5' position)
	protected int tag5PrimeExtension=0; //Extend tags this number of bp in the 5' direction (starting from 5' position)
	protected float tagGaussSigma=0;    //Sigma for Gaussian tag smoothing (std dev)
	protected float tagGaussWidthFactor=5; //If smoothing with a Gaussian, each tag's prob density 'extends' over this factor times sigma 
	protected int binWidth=50;		//Bin size for feature scanning
	protected int binStep=25;		//Bin step for feature scanning
	protected int featureMergeWindow=100; //Merge neighboring domains that are this close to one another (in bp)
	protected List<Integer> localBackWins=new ArrayList<Integer>(); //MACS-style local background windows
	protected double perBinPoissonLogPThres=-7; //Log base 10 confidence threshold for use with Poisson background models when scoring bins
	protected double perBinBinomialPThres = 0.01; //threshold to use with per-bin Binomial test (we can afford to be lax here, as the q-value threshold filters for final results)
	protected double perFeatureBinomialQThres = 0.01; //threshold to use with per-event Binomial tests 
	protected double minSigCtrlFoldDifference = 1; //minimum signal/control difference to use in Binomial testing (per bin & per event)
	protected PeakFindingMethod peakFinding = PeakFindingMethod.MAXDENSITY; //Approach for finding peak point(s) in enriched domains 
	protected int maxThreads=1;		//Number of threads to use. Default is 1 for single processor machines.
	protected List<Region> regionsToIgnore = new ArrayList<Region>(); //List of regions that will be ignored during EM training (i.e. known towers, etc)
	protected boolean printHelp=false;
	protected boolean outputGFF=false; //Print the output files in GFF format
	protected int superStitchWindow = 12500; // Window to stitch typical enhancers to super enhancers
	protected int minDistalDistance = 2000; // Minmum distance from TSSs to call an element Distal
	protected List<StrandedPoint> refTSSs = new ArrayList<StrandedPoint>(); // refTSSs; currently loading from a stranded peak file 
	
	
	//Constants
	public final int MAXSECTION = 50000000; //Size of genomic sections that each thread analyzes

	
	protected String[] args;
	public String getArgs(){
		String a="";
		for(int i=0; i<args.length; i++)
			a = a+" "+args[i];
		return a;
	}
	
	public SEEDConfig(GenomeConfig gcon, String[] args){
		this.args = args;
		gconfig = gcon;
		ArgParser ap = new ArgParser(args);
		if(args.length==0 || ap.hasKey("h")){
			printHelp=true;			
		}else{
			try{
				//Test for a config file... if there is concatenate the contents into the args
				if(ap.hasKey("config")){
					ArrayList<String> confArgs = new ArrayList<String>();
					String confName = ap.getKeyValue("config");
					File confFile = new File(confName);
					if(!confFile.isFile())
						System.err.println("\nCannot find configuration file: "+confName);
					BufferedReader reader = new BufferedReader(new FileReader(confFile));
				    String line;
			        while ((line = reader.readLine()) != null) {
			        	line = line.trim();
			        	String[] words = line.split("\\s+");
			        	if(!words[0].startsWith("--"))
			        		words[0] = new String("--"+words[0]);
			        	confArgs.add(words[0]); 
			        	if(words.length>1){
				        	String rest=words[1];
				        	for(int w=2; w<words.length; w++)
				        		rest = rest+" "+words[w];
				        	confArgs.add(rest);
			        	}
			        }
			        String [] confArgsArr = confArgs.toArray(new String[confArgs.size()]);
			        String [] newargs =new String[args.length + confArgsArr.length];
			        System.arraycopy(args, 0, newargs, 0, args.length);
			        System.arraycopy(confArgsArr, 0, newargs, args.length, confArgsArr.length);
			        args = newargs;
			        ap = new ArgParser(args);
				}
				
				/******* General options *******/
				//Output path
				DateFormat df = new SimpleDateFormat("yyyy-MM-dd-hh-mm-ss");  
			    df.setTimeZone(TimeZone.getTimeZone("EST"));
				outName = Args.parseString(args, "out", outName+"_"+df.format(new Date()));
				outDir =  new File(outName); //Output directory
				outBase = outDir.getName(); //Last part of name
				outputGFF = ap.hasKey("gffout");
				
				//Threads
				maxThreads = Args.parseInteger(args,"threads",maxThreads);
				maxThreads = Math.min(maxThreads, java.lang.Runtime.getRuntime().availableProcessors());
				
				//Feature detection 
				binWidth = Args.parseInteger(args,"binwidth",binWidth);
				binStep = Args.parseInteger(args,"binstep",binStep);
				featureMergeWindow = Args.parseInteger(args,"mergewin",featureMergeWindow);
				
				//Regions to ignore during analysis
				if(ap.hasKey("exclude"))
					regionsToIgnore = Utils.loadRegionsFromFile(Args.parseString(args, "exclude", null), gconfig.getGenome(), -1);
				
				//Tag manipulation
				tagShift = Args.parseInteger(args,"tagshift",tagShift);
				tag3PrimeExtension = Args.parseInteger(args,"tag3ext",tag3PrimeExtension);
				tag5PrimeExtension = Args.parseInteger(args,"tag5ext",tag5PrimeExtension);
				tagGaussSigma = Args.parseFloat(args,"tagsigma",tag5PrimeExtension);
				
				//Statistical thresholds
				perBinPoissonLogPThres = Args.parseDouble(args,"poisslogpthres",perBinPoissonLogPThres);
				perBinBinomialPThres = Args.parseDouble(args,"binpthres",perBinBinomialPThres);
				perFeatureBinomialQThres = Args.parseDouble(args,"featureqthres",perFeatureBinomialQThres);
				minSigCtrlFoldDifference = Args.parseDouble(args,"minfolddiff",minSigCtrlFoldDifference);
				
				//Local background models
				localBackWins= (List<Integer>) Args.parseIntegers(args, "localbackwin");
				
				//Peak-finding method
				String pfmethod = Args.parseString(args,"peakfinding","max").toLowerCase();
				if(pfmethod.equals("max")){peakFinding = PeakFindingMethod.MAXDENSITY;}
				else if(pfmethod.equals("lrbal")){peakFinding = PeakFindingMethod.LRBALANCE;}
				else if(pfmethod.equals("distrib")){peakFinding = PeakFindingMethod.TAGDISTRIBUTION;}
				else{peakFinding = PeakFindingMethod.MAXDENSITY;}
				
				// Super Feature detection method optioins
				superStitchWindow = Args.parseInteger(args, "supStitchWin", 12500);
				minDistalDistance  = Args.parseInteger(args, "distalDistance", 2000);
				if(ap.hasKey("refTSSs")){
					refTSSs = Utils.loadStrandedPointsFromFile(gconfig.getGenome(), Args.parseString(args, "refTSSs",""));
				}
				
						
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	//Accessors
	public boolean helpWanted(){return printHelp;}
	public String getOutName(){return outName;}
	public String getOutBase(){return outBase;}
	public File getOutputParentDir(){return outDir;}
	public int getMaxThreads(){return maxThreads;}
	public List<Region> getRegionsToIgnore(){return regionsToIgnore;}
	public int getTagShift(){return tagShift;}
	public int getTag3PrimeExtension(){return tag3PrimeExtension;}
	public int getTag5PrimeExtension(){return tag5PrimeExtension;}
	public float getTagGaussSigma(){return tagGaussSigma;}
	public int getTagGaussWidth(){return (int)(tagGaussSigma * tagGaussWidthFactor);}
	public int getBinWidth(){return binWidth;}
	public int getBinStep(){return binStep;}
	public int getFeatureMergeWindow(){return featureMergeWindow;}
	public List<Integer> getLocalBackWins(){return localBackWins;}
	public double getPerBinPoissonLogPThres(){return perBinPoissonLogPThres;}
	public double getPerBinBinomialPThres(){return perBinBinomialPThres;}
	public double getPerFeatureBinomialQThres(){return perFeatureBinomialQThres;}
	public double getMinSigCtrlFoldDifference(){return minSigCtrlFoldDifference;}
	public PeakFindingMethod getPeakFindingApproach(){return peakFinding;}
	public boolean outputGFF(){return outputGFF;}
	public int getSuperStitchWin(){return superStitchWindow;}
	public int getMinDistalDistance(){return minDistalDistance;}
	public List<StrandedPoint> getRefTSSs(){return refTSSs;}

	//Settors
	public void setTagShift(int s){tagShift=s;}
	public void setTag3PrimeExtension(int e){tag3PrimeExtension = e;}
	public void setTag5PrimeExtension(int e){tag5PrimeExtension = e;}
	public void setTagGaussSigma(float s){tagGaussSigma = s;}
	public void setPeakFindingApproach(PeakFindingMethod p){peakFinding=p;}
	
	
	/**
	 * Make some output directories used by SEED implementations
	 */
	public void makeSEEDOutputDirs(){
		//Test if output directory already exists. If it does,  recursively delete contents
		outDir =  new File(outName);
		if(outDir.exists())
			deleteDirectory(outDir);
		outBase = outDir.getName();
		//(re)make the output directory
		outDir.mkdirs();
	}
	/**
	 * Delete a direcctory
	 */
	public boolean deleteDirectory(File path) {
	    if( path.exists() ) {
	      File[] files = path.listFiles();
	      for(int i=0; i<files.length; i++) {
	         if(files[i].isDirectory()) {
	           deleteDirectory(files[i]);
	         }
	         else {
	           files[i].delete();
	         }
	      }
	    }
	    return( path.delete() );
	}
	
	/**
	 * returns a string describing the arguments handled by this parser. 
	 * @return String
	 */
	public String getArgsList(){
		return(new String("" +
				"SEED arguments:\n"+
				"\t--out <output file prefix>\n" +
				"\t--tagshift <shift tags this number of bases in 3' direction (default: shift="+tagShift+")>\n" +
				"\t--tag3ext <extend tags this number of bp from (shifted) 5' position in the 3' direction (default="+tag3PrimeExtension+")>\n" +
				"\t--tag5ext <extend tags this number of bp from (shifted) 5' position in the 5' direction (default="+tag5PrimeExtension+")>\n" +
				"\t--tagsigma <smooth tags using Gaussian with this sigma (default="+tagGaussSigma+")>\n" +
				"\t--binwidth <width of bin in bp for feature scanning (default="+binWidth+")>\n" +
				"\t--binstep <bin step size in bp for feature scanning (default="+binStep+")>\n" +
				"\t--mergewin <size of window in which to merge features/domains (default="+featureMergeWindow+")>\n"+
				"\t--localbackwin <size of local window(s) for MACS-style calculation of expected tag count (default=none)>\n" +
				"\t--poisslogpthres <log base 10 Poisson confidence threshold for use in finding enriched bins (default="+perBinPoissonLogPThres+")>\n" +
				"\t--binpthres <uncorrected p-value threshold for use in Binomial test for enriched bins (default="+perBinBinomialPThres+")>\n" +
				"\t--featureqthres <q-value threshold for use in Binomial test for enriched events (default="+perFeatureBinomialQThres+")>\n" +
				"\t--minfolddiff <minimum signal/control fold difference for use in Binomial test for enriched domains (default="+minSigCtrlFoldDifference+")>\n" +
				"\t--peakfinding <max/lrbal/distrib: approach for finding the peak in a domain; max. density, left-right balance point, or tag density scan, resp. (default=max)>\n" +
				"\t--threads <number of threads to use>\n" +
				"\t--config <config file: all options can be specified in a name<space>value text file, over-ridden by command-line args>\n" +
				"\t--exclude <file of regions to ignore>\n" +
				"\t--gffout <output GFF format files>\n" +
				""));
	}

	/**
	 * Enumerated type for peak-finding approach
	 * @author mahony
	 */
	public enum PeakFindingMethod{ 
		MAXDENSITY, LRBALANCE, TAGDISTRIBUTION
	}
}

