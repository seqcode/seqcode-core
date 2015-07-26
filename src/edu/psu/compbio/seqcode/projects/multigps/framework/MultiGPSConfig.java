package edu.psu.compbio.seqcode.projects.multigps.framework;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TimeZone;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.projects.multigps.utilities.AnnotationLoader;
import edu.psu.compbio.seqcode.projects.multigps.utilities.Utils;

/**
 * MultiGPSConfig: 
 * 		Maintains all constants needed by MultiGPS. 
 *     
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class MultiGPSConfig {
	public static String version = "0.7";
	public boolean isGPS=true;
	protected GenomeConfig gconfig;
	protected Genome gen=null;
	protected String outName="multigps", outBase="multigps";
	protected File outDir=null, interDir=null, imagesDir=null;
	protected boolean printHelp=false;
	protected BindingModel defaultModel=null;
	protected double sigLogConf=-7; 
	protected double prLogConf=-6; 
	protected double qMinThres=0.001;		//Minimum  Q-value for reported binding events
	protected double differentialSignificanceP = 0.01;
	protected int maxModelUpdateRounds=3;
	protected int maxThreads=1;				//Number of threads to use. Default is 1 for single processor machines. 
	protected double alphaScalingFactor = 1.0; //Scale the condition-specific alpha value by this factor
	protected boolean multicondition_posprior=true; //Multiple condition positional prior
	protected double prob_shared_binding=0.9; //Prior probability that binding sites are shared between conditions (Constant used to build positional priors between conditions)
	protected int bmAnalysisWindowMax=10000;
	protected int minComponentsForBMUpdate = 500;
	protected double minComponentReadFactorForBM = 3; //Components must have (this factor times the condition alpha) number of reads assigned before being included in BM update
	protected boolean smoothingBMDuringUpdate=true;
	protected boolean gaussianSmoothingBMDuringUpdate=false;
	protected boolean updateBM=true; //Set to false to turn off binding model update
	protected boolean includeJointEventsInBMUpdate=false; 
	protected double bindingmodel_spline_smooth = 30; //Smoothing step for cubic spline in binding model reestimation
    protected double bindingmodel_gauss_smooth = 2; //Variance for Gaussian smoothing in binding model reestimation
	protected int addFlankingComponentSpacing=20; //In non-first rounds of EM, the components are initialized using the positions from the last round with additional flanking components added at this spacing
	protected double minEventFoldChange=1.5;
	protected boolean addAnnotations=false;
	protected boolean addSequences=true;
	protected List<AnnotationLoader> geneAnnotations = new ArrayList<AnnotationLoader>();
	protected int maxAnnotDistance=50000;
	protected boolean annotOverlapOnly=false;
	protected List<Region> regionsToPlot = new ArrayList<Region>(); //List of regions that will be printed during EM training (for debugging/demonstration)
	protected List<Region> regionsToIgnore = new ArrayList<Region>(); //List of regions that will be ignored during EM training (i.e. known towers, etc)
	protected boolean fixedModelRange = false;
	protected boolean MLSharedComponentConfiguration = true; //For ML assignment: use a component configuration shared across all conditions or have condition-specific configs.
	protected boolean findMotifs = true; //Run motif-finding for motif prior
	protected boolean motif_posprior=true; //You can have motif-finding without using the motif-prior
	protected String MEMEpath="";
	protected String MEMEargs=" -dna -mod zoops -revcomp -nostatus ";    //Do not try using -p here; it leads to MEME runtime errors
	public int MEMEminw=6;
	public int MEMEmaxw=18;
	protected boolean runDiffTests = true; //Run differential enrichment testing
	protected String Rpath="";
	protected double edger_overdispersion = 0.15; //Overdispersion used by EdgeR differential enrichment tests
	protected boolean verbose = false; //Print extra output
	 
    
	//Constants
	public final double LOG2 = Math.log(2);
	public final int POTREG_BIN_STEP = 100; //Sliding window step in potential region scanner
	public final int MAXSECTION = 50000000;
    public final int INIT_COMPONENT_SPACING=30;  //Initial component spacing
    public final int MAX_EM_ITER=2000;
    public final int EM_ML_ITER=100;     				//Run EM up until <tt>ML_ITER</tt> without using sparse prior
    public final int ML_ML_ITER=100;     				//Run ML up until <tt>ML_ITER</tt> without using sparse prior
    public final int ALPHA_ANNEALING_ITER=100;     //Run EM up until <tt>ALPHA_ANNEALING_ITER</tt> with smaller alpha based on the current iteration
    public final int POSPRIOR_ITER=150;     //Run EM up until <tt>ALPHA_ANNEALING_ITER</tt> with uniform positional prior and then up until at least <tt>POSPRIOR_ANNEALING_ITER</tt> with activated positional prior
    public final int EM_MU_UPDATE_WIN=100; //Half the window size in which to look for mu maximization (i.e. component position) during EM.
    public final double EM_CONVERGENCE = 1e-10; //EM convergence between the likelihood of the current and the previous step
    public final double EM_STATE_EQUIV_THRES = 1e-10; //EM state equivalence threshold 
    public final int EM_STATE_EQUIV_ROUNDS = 3; //Number of training rounds where the EM states have to be equivalent
    public final double NOISE_EMISSION_MIN = 0.01; //Arbitrary floor on the emission probability of noise (must be non-zero to mop up noise reads)
    public final double NOISE_EMISSION_MAX = 0.95; //Arbitrary ceiling on the emission probability of noise
    public final int NOISE_DISTRIB_SMOOTHING_WIN = 50; //Smoothing window for the noise distribution used in the BindingMixture
    public final int MAX_BINDINGMODEL_WIDTH=800; //Maximum width for binding models (affects how large the read profiles are for binding components    
    public final int MOTIF_FINDING_SEQWINDOW=80; //Bases to extract around components for motif-finding
    public final int MOTIF_FINDING_TOPSEQS=500; //Number of top components to analyze
    public final double MOTIF_FINDING_ALLOWED_REPETITIVE = 0.2; //Percentage of the examined sequence window allowed to be lowercase or N
    public final int MOTIF_FINDING_NEGSEQ=5000; //Number of negative sequences for motif significance tests
    public final double MOTIF_MIN_ROC = 0.7; //Motif prior is used only if the ROC is greater than this .
    public final double LOG_FC_LIMIT = 10; //Maximum absolute log fold-change reported
	public final boolean CALC_LL=false; //Calculate the log-likelihood during EM.
	public final boolean CALC_COMP_LL=false; //Calculate component-wise log-likelihoods during ML 
    
	protected String[] args;
	public String getArgs(){
		String a="";
		for(int i=0; i<args.length; i++)
			a = a+" "+args[i];
		return a;
	}
	
	public MultiGPSConfig(GenomeConfig gcon, String[] arguments){this(gcon, arguments, true);}
	public MultiGPSConfig(GenomeConfig gcon, String[] arguments, boolean isGPS){
		System.setProperty("java.awt.headless", "true");
		gconfig = gcon;
		gen = gconfig.getGenome();
		this.args=arguments; 
		this.isGPS=isGPS;
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
				
				//Read distribution file
				String modelFile = Args.parseString(args, "d", null);	// read distribution file
				if (modelFile != null){
					File pFile = new File(modelFile);
					if(!pFile.isFile()){
						System.err.println("\nCannot find read distribution file: "+modelFile);
						System.exit(1);
					}
					defaultModel = new BindingModel(pFile);
				}
				
				
				/****Gene Annotation****/
				//Gene Annotations
				Collection<String> tfiles = Args.parseStrings(args,"transcripts");
				Collection<String> dbgenes = Args.parseStrings(args,"dbgenes");
				for(String s:dbgenes)
		        	geneAnnotations.add(new AnnotationLoader(gen, s, "refGene", maxAnnotDistance, annotOverlapOnly));
		        for(String s:tfiles)
		        	geneAnnotations.add(new AnnotationLoader(gen, s, "file", maxAnnotDistance, annotOverlapOnly));
				if(geneAnnotations.size()>0)
					addAnnotations=true;
				
				/****Miscellaneous arguments****/
				
				//Maximum number of model update rounds
				maxModelUpdateRounds = Args.parseInteger(args,"r", 3);
				//Turn off binding model updates
				updateBM = Args.parseFlags(args).contains("nomodelupdate") ? false : true;
				//Minimum number of components to support a binding model update
				minComponentsForBMUpdate = Args.parseInteger(args,"minmodelupdateevents",minComponentsForBMUpdate);
				//Turn off smoothing during binding model updates 
				smoothingBMDuringUpdate = Args.parseFlags(args).contains("nomodelsmoothing") ? false : true;
				//Parameter for spline smoothing
				bindingmodel_spline_smooth = Args.parseDouble(args,"splinesmoothparam",bindingmodel_spline_smooth); 
				//Turn on Gaussian smoothing during binding model updates
				gaussianSmoothingBMDuringUpdate = Args.parseFlags(args).contains("gaussmodelsmoothing") ? true : false;
				//Parameter for Gaussian smoothing (std. dev.)
				bindingmodel_gauss_smooth = Args.parseDouble(args,"gausssmoothparam",bindingmodel_gauss_smooth);
				//Allow joint events in binding model updates
				includeJointEventsInBMUpdate = Args.parseFlags(args).contains("jointinmodel");
				//Fixed binding model range
				fixedModelRange = Args.parseFlags(args).contains("fixedmodelrange");
				//Output path
				DateFormat df = new SimpleDateFormat("yyyy-MM-dd-hh-mm-ss");  
			    df.setTimeZone(TimeZone.getTimeZone("EST"));
				outName = Args.parseString(args, "out", outName+"_"+df.format(new Date()));
				outDir =  new File(outName); //Output directory
				outBase = outDir.getName(); //Last part of name
				//Background model parameters
				sigLogConf = Args.parseDouble(args,"highlogconf",sigLogConf);
				prLogConf = Args.parseDouble(args,"prlogconf",prLogConf);
				//Q-value threshold
				qMinThres = Args.parseDouble(args,"q",qMinThres);
				//differential p-value threshold
				differentialSignificanceP = Args.parseDouble(args,"diffp",differentialSignificanceP);
				//Threads
				maxThreads = Args.parseInteger(args,"threads",maxThreads);
				maxThreads = Math.min(maxThreads, java.lang.Runtime.getRuntime().availableProcessors());
				//Alpha scaling factor
				alphaScalingFactor = Args.parseDouble(args,"alphascale",alphaScalingFactor);
				//Event Fold-change minimum
				minEventFoldChange = Args.parseDouble(args,"minfold",minEventFoldChange);
				//Regions to print during EM training
				if(ap.hasKey("plotregions"))
					regionsToPlot = Utils.loadRegionsFromFile(Args.parseString(args, "plotregions", null), gen, -1);
				//Regions to ignore during EM training
				if(ap.hasKey("exclude"))
					regionsToIgnore = Utils.loadRegionsFromFile(Args.parseString(args, "exclude", null), gen, -1);
				//Turn off multi-condition positional prior
				multicondition_posprior = Args.parseFlags(args).contains("noposprior") ? false : true;
				//Set a value for the multi-condition positional prior
				prob_shared_binding = Args.parseDouble(args,"probshared",prob_shared_binding);
				//Turn off motif-finding 
				findMotifs = Args.parseFlags(args).contains("nomotifs") ? false : true;
				//Turn off motif prior only
				motif_posprior = (findMotifs && Args.parseFlags(args).contains("nomotifprior")) ? false : true;				
				//Check whether sequence is available (affects motif-finding behavior)
				if(isGPS && !gconfig.getSequenceGenerator().usingLocalFiles()){
					findMotifs=false;
					motif_posprior=false;
					System.err.println("No genome sequence data was provided with --seq, so motif-finding and the motif prior are switched off.");
				}
				
				//MEME path
				MEMEpath = Args.parseString(args, "memepath", MEMEpath);
				if(!MEMEpath.equals("") && !MEMEpath.endsWith("/")){ MEMEpath= MEMEpath+"/";}
				//MEME args
				MEMEargs = Args.parseString(args, "memeargs", MEMEargs);
				//MEME minw
				MEMEminw = Args.parseInteger(args, "mememinw", MEMEminw);
				//MEME maxw
				MEMEmaxw = Args.parseInteger(args, "mememaxw", MEMEmaxw);
				//MEME nmotifs option
				int MEMEnmotifs = Args.parseInteger(args,"memenmotifs", 3);
				MEMEargs = MEMEargs + " -nmotifs "+MEMEnmotifs + " -minw "+MEMEminw+" -maxw "+MEMEmaxw;
				
				//Turn off DE testing
				runDiffTests = Args.parseFlags(args).contains("nodifftests") ? false : true;
				//R path
				Rpath = Args.parseString(args, "rpath", Rpath);
				if(!Rpath.equals("") && !Rpath.endsWith("/")){ Rpath= Rpath+"/";}
				//EdgeR overdispersion parameter
				edger_overdispersion = Args.parseDouble(args,"edgerod",edger_overdispersion);
				
				//Extra output
				verbose = Args.parseFlags(args).contains("verbose") ? true : false;
				//Shared component config in ML step
				//MLSharedComponentConfiguration = Args.parseFlags(args).contains("mlsharedconfig") ? true : false;
				MLSharedComponentConfiguration = Args.parseFlags(args).contains("mlconfignotshared") ? false : true;
				
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Merge a set of estimated genomes 
	 * @param estGenomes
	 * @return
	 */
	public Genome mergeGenomes(List<Genome> estGenomes){
		//Combine the chromosome information
		HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
		for(Genome e : estGenomes){
			Map<String, Integer> currMap = e.getChromLengthMap();
			for(String s: currMap.keySet()){
				if(!chrLenMap.containsKey(s) || chrLenMap.get(s)<currMap.get(s))
					chrLenMap.put(s, currMap.get(s));
			}
		}
		gen =new Genome("Genome", chrLenMap);
		return gen;		
	}
	
	//Accessors
	public Genome getGenome(){return gen;}
	public boolean helpWanted(){return printHelp;}
	public double getQMinThres(){return qMinThres;}
	public double getDiffPMinThres(){return differentialSignificanceP;}
	public double getMinEventFoldChange(){return minEventFoldChange;}
	public double getSigLogConf(){return sigLogConf;}
	public double getPRLogConf(){return prLogConf;}
	public int getMaxThreads(){return maxThreads;}
	public double getAlphaScalingFactor(){return alphaScalingFactor;}
	public boolean useMultiConditionPosPrior(){return multicondition_posprior;}
	public double getProbSharedBinding(){return prob_shared_binding;}
	public int getBMAnalysisWindowMax(){return bmAnalysisWindowMax;}
	public int getAddFlankingComponentSpacing(){return addFlankingComponentSpacing;}
	public boolean isAddingAnnotations(){return addAnnotations;}
	public boolean isAddingSequences(){return addSequences;}
	public List<AnnotationLoader> getGeneAnnotations(){return geneAnnotations;}
	public int getMaxAnnotDistance(){return maxAnnotDistance;}
	public boolean getAnnotOverlapOnly(){return annotOverlapOnly;}
	public List<Region> getRegionsToPlot(){return regionsToPlot;}
	public List<Region> getRegionsToIgnore(){return regionsToIgnore;}
	public BindingModel getDefaultBindingModel(){return defaultModel;}
	public boolean doBMUpdate(){return updateBM;}
	public int getMinComponentsForBMUpdate(){return minComponentsForBMUpdate;}
	public double getMinComponentReadFactorForBM(){return minComponentReadFactorForBM;}
	public boolean getIncludeJointEventsInBMUpdate(){return includeJointEventsInBMUpdate;}
	public boolean getSmoothingBMDuringUpdate(){return smoothingBMDuringUpdate;}
	public double getBindingModelSplineSmoothParam(){return bindingmodel_spline_smooth;}
	public boolean getGaussBMSmooth(){return gaussianSmoothingBMDuringUpdate;}
	public double getBindingModelGaussSmoothParam(){return bindingmodel_gauss_smooth;}
	public int getMaxModelUpdateRounds(){return maxModelUpdateRounds;}
	public boolean getFixedModelRange(){return fixedModelRange;}
	public boolean getMLSharedComponentConfiguration(){return MLSharedComponentConfiguration;}
	public boolean getFindingMotifs(){return findMotifs;}
	public boolean useMotifPrior(){return motif_posprior;}
	public String getMEMEpath(){return MEMEpath;}
	public String getMEMEargs(){return MEMEargs;}
	public boolean getRunDiffTests(){return runDiffTests;}
	public String getRpath(){return Rpath;}
	public double getEdgeROverDisp(){return edger_overdispersion;}
	public boolean isVerbose(){return verbose;}
	
	/**
	 * Make some output directories used by multiGPS
	 */
	public void makeGPSOutputDirs(boolean makeInterAndImageDirs){
		//Test if output directory already exists. If it does,  recursively delete contents
		outDir =  new File(outName);
		if(outDir.exists())
			deleteDirectory(outDir);
		outBase = outDir.getName();
		//(re)make the output directory
		outDir.mkdirs();
		if(makeInterAndImageDirs){
			//Make the gps intermediate results output directory
			interDir = new File(outDir.getAbsolutePath()+File.separator+"intermediate-results");
			interDir.mkdirs();
			//Make the image results output directory
			imagesDir = new File(outDir.getAbsolutePath()+File.separator+"images");
			imagesDir.mkdirs();
		}
	}
	public String getOutName(){return outName;}
	public String getOutBase(){return outBase;}
	public File getOutputParentDir(){return outDir;}
	public File getOutputIntermediateDir(){return interDir;}
	public File getOutputImagesDir(){return imagesDir;}
	
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
				"Genome:" +
				"\t--species <Species;Genome>\n" +
				"\tOR\n" +
				"\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" +
				"General:\n" +
				"\t--d <read distribution model file>\n" +
				"\t--r <max. model update rounds (default="+maxModelUpdateRounds+">\n" +
				"\t--out <out name (default="+outBase+">\n" +
				"\t--nonunique [flag to use non-unique reads]\n" +
				"\t--threads <number of threads to use (default="+maxThreads+")>\n" +
				"Annotations:\n" +
				"\t--transcripts <transcripts file>\n" +
				"\t--dbgenes refGene\n" +
				"Experiment Design File:\n" +
				"\t--design <file name>\n" +
				"Miscellaneous:\n" +
				"\t--q <Q-value minimum, i.e corrected p-value(default="+qMinThres+")>\n" +
				"\t--minfold <min event fold-change (default="+minEventFoldChange+")>\n" +
				"\t--prlogconf <Poisson log threshold for potential region scanning(default="+prLogConf+")>\n" +
				"\t--alphascale <alpha scaling factor(default="+alphaScalingFactor+">\n" +
				"\t--nomodelupdate [flag to turn off binding model updates]\n" +
				"\t--minmodelupdateevents <minimum number of events to support an update(default="+minComponentsForBMUpdate+")>\n" +
				"\t--nomodelsmoothing [flag to turn off binding model smoothing]\n" +
				"\t--splinesmoothparam <spline smoothing parameter (default="+bindingmodel_spline_smooth+">\n" +
				"\t--gaussmodelsmoothing [flag to turn o Gaussian model smoothing (default = cubic spline)]\n" +
				"\t--gausssmoothparam <Gaussian smoothing std dev (default="+bindingmodel_gauss_smooth+">\n" +
				"\t--jointinmodel [flag to allow joint events in model updates]\n" +
				"\t--fixedmodelrange [flag to keep binding model range constant]\n" +
				//"\t--mlsharedconfig [flag to share the component config in the ML step]\n" +
				"\t--mlconfignotshared [flag to not share component configs in the ML step]\n" +
				"\t--exclude <file of regions to ignore>\n" +
				"\t--plotregions <regions to print during EM training>\n" +
				"\t--noposprior [flag to turn off multi-cond pos prior]\n" +
				"\t--probshared <probability that events are shared across conditions, i.e. pos prior(default="+prob_shared_binding+")>\n" +
				"\t--nomotifs [flag to turn off motif-finding & motif priors]\n" +
				"\t--nomotifprior [flag to turn off motif priors only]\n" +
				"\t--memepath <path to the meme bin dir (default: meme is in $PATH)>\n" +
				"\t--memenmotifs <number of motifs MEME should find for each condition (default=3)>\n" +
				"\t--mememinw <minw arg for MEME (default="+MEMEminw+")>\n"+
				"\t--mememaxw <maxw arg for MEME (default="+MEMEmaxw+")>\n"+
				"\t--nodifftests [flag to turn off DE tests]\n" +
				"\t--rpath <path to the R bin dir (default: R is in $PATH). Note that you need to install edgeR separately>\n" +
				"\t--edgerod <EdgeR overdispersion (default="+edger_overdispersion+")>\n" +
				"\t--diffp <minimum p-value for differential enrichment (default="+differentialSignificanceP+")>\n" +
				"\t--verbose [flag to print intermediate files and extra output]\n" +
				"\t--config <config file: all options can be specified in a name<space>value text file, over-ridden by command-line args>\n" +
				""));
	}
	}
