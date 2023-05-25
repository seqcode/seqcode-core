package org.seqcode.deepseq.experiments;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.Pair;


/**
 * ExptConfig: 
 * A config parser that:		
 * 		Handles specification of experiments and read loaders from command-line, config files and design files.
 * 		Maintains constants needed by experiment loading. 
 *     
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExptConfig {
	protected Genome gen=null;
	protected List<ExptDescriptor> expts = new ArrayList<ExptDescriptor>();
	protected boolean nonUnique=false;
	protected boolean printHelp=false;
	protected boolean printLoadingProgress=true;
	protected List<Integer> localBackgroundWindows=new ArrayList<Integer>(); 
	protected double perBaseLogConf=-7;
	protected boolean poissonGaussWinPerBaseFilter = false; //Filter using a poisson Gaussian window
	protected float perBaseReadLimit = -1;
	protected boolean perBaseReadFiltering = true;
	protected double mappableGenome = 0.8;
	protected boolean estimateScaling=true;
	protected boolean scalingByRegression=false; //Default is to scale by NCIS
	protected boolean scalingBySES = false; //Default is to estimate scaling by NCIS
	protected boolean scalingByMedian = false; //Default is to estimate scaling by NCIS
	protected boolean scalingByHitRatioAndNCIS = false; //Default is to estimate scaling by NCIS
	protected float fixedScalingFactor = 1; //Default is to estimate scaling by NCIS
	protected int scalingSlidingWindow = 10000; 
	protected boolean plotScaling = false; //Make a scaling method plot
	protected boolean cacheAllHits=true; //Cache all hits
	protected String fileCacheDir = "hitcache";
	protected List<Region> initialCachedRegions=null;
	//Different loaders will have different behaviors in the following
	//For example, some file formats cannot store pairs. ReadDB ignores the difference between R1 & R2 in single-end, etc.
	protected boolean loadType1Reads = true; //Load Type1 reads
	protected boolean loadType2Reads = false; //Load Type2 reads (if exists and distinguishable)
	protected boolean loadRead2=true; //Load second in pair reads (only used by BAM loader for now)
	protected boolean loadReads = true;  //Load read information
	protected boolean loadPairs = false; //Load pair information (if exists)
	protected boolean sortMid = false; //Sort pairs according to midpoint
	protected boolean keepHDF5 = false; // Keep the hdf5 file for future use
	protected double NCISMinBinFrac = 0.75; //NCIS estimates begin using the lower fraction of the genome (based on total tags)
	
	    
	protected String[] args;
	public String getArgs(){
		String a="";
		for(int i=0; i<args.length; i++)
			a = a+" "+args[i];
		return a;
	}
	
	public ExptConfig(Genome g, String[] arguments){
		this.gen = g;
		this.args=arguments; 
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
			        reader.close();
				}
				//////////////////////
				//Which reads to load?
				////////////////////////
				loadPairs = Args.parseFlags(args).contains("loadpairs");
				loadType1Reads = !Args.parseFlags(args).contains("not1reads");
				loadType2Reads = Args.parseFlags(args).contains("loadt2reads");
				loadRead2 = !Args.parseFlags(args).contains("noread2");
				sortMid = Args.parseArgs(args).contains("sortMid");
				
				//////////////////////
				//keep the HDF5Cache file?
				////////////////////////
				keepHDF5 = Args.parseFlags(args).contains("keepHDF5Cache");
				
				////////////////////////
				//Read limit parameters
				////////////////////////
				//Fixed per base read limits
				perBaseReadLimit = Args.parseInteger(args,"fixedpb",-1);
				//Do per base filtering with a Poisson Gaussian window
				poissonGaussWinPerBaseFilter = Args.parseFlags(args).contains("poissongausspb");
				if(poissonGaussWinPerBaseFilter)
					perBaseReadLimit=0;
				//Switch off per-base read filtering?
				perBaseReadFiltering = !Args.parseFlags(args).contains("nopbfilter");
				//////////////////////
				//Scaling parameters
				//////////////////////
				//Turn off scaling estimation
				estimateScaling = Args.parseFlags(args).contains("noscaling") ? false : true;
				//Fixed scaling factor if not estimating scaling
				if(ap.hasKey("fixedscaling")){
					estimateScaling=false;
					fixedScalingFactor = Args.parseFloat(args,"fixedscaling",1);
				}
				//Scale by NCIS is default
				NCISMinBinFrac = Args.parseDouble(args, "ncisbinmin", NCISMinBinFrac);	//NCIS estimates begin using the lower fraction of the genome (based on total tags)
				//Scale by median 
				scalingByMedian = Args.parseFlags(args).contains("medianscale") ? true : false;
				//Scale by SES
				scalingBySES = Args.parseFlags(args).contains("sesscale") ? true : false;
				//Scale by regression
				scalingByRegression = Args.parseFlags(args).contains("regressionscale") ? true : false;
				//Scale by total tag followed by NCIS
				scalingByHitRatioAndNCIS = Args.parseFlags(args).contains("normalizedncisscale") ? true : false;
				//Scaling window
				scalingSlidingWindow = Args.parseInteger(args,"scalewin",scalingSlidingWindow);
				//Make a scaling method plot
				plotScaling = Args.parseFlags(args).contains("plotscaling") ? true : false;
				////////////////////////
				//Misc parameters
				////////////////////////
				//Mappability 
				mappableGenome =  Args.parseDouble(args,"mappability",0.8);
				//Use non-unique reads?
				nonUnique = Args.parseFlags(args).contains("nonunique") ? true : false;
				//Local dynamic backgrounds
				localBackgroundWindows= (List<Integer>) Args.parseIntegers(args, "dynback");
				if(localBackgroundWindows.size()==0){localBackgroundWindows.add(10000);}
				//Caching
				cacheAllHits = Args.parseFlags(args).contains("nocache") ? false : true;
				
				//Parse command-line experiments (optional experiment and replicate names can be specified within the argument name - e.g. --exptName-Rep )
				String fileFormat = Args.parseString(args, "format", "SAM").toUpperCase();
				ArrayList<String> exptTags = new ArrayList<String>();
				for(String s : ap.getKeys()){
					if(!exptTags.contains(s)){
			        	if(s.contains("expt") || s.contains("ctrl")){
			        		String datatype = "";
			        		String condrep = "";
			        		boolean signal = true;
			        		if(s.contains("ctrl"))
			        			signal = false;
			        		if(s.startsWith("rdb")){
			        			datatype = "READDB";
			        			condrep = s.replaceFirst("rdbexpt", "");
			        			condrep = condrep.replaceFirst("rdbctrl", "");
			        		}else{
			        			datatype = fileFormat;
			        			condrep = s.replaceFirst("expt", "");
			        			condrep = condrep.replaceFirst("ctrl", "");
			        		}
			        		//Parse the condition & replicate names
			        		String cond = condrep.split("-")[0];
			        		if(cond.length()==0 || cond.equals(""))
			        			cond = signal ? "experiment" : "DEFAULT";
			        		String rep = condrep.split("-").length>1 ? condrep.split("-")[1] : "";
			        		if(rep.length()==0 || rep.equals("")){
			        			rep = signal ? "rep1" : "DEFAULT";
			        		}
			        		//Parse the file/rdb names
			        		Collection<String> sourceNames = Args.parseStrings(args, s);
			        		List<Pair<String,String>> sources = new ArrayList<Pair<String,String>>();
			        		for(String n : sourceNames)
			        			sources.add(new Pair<String,String>(n,datatype));
			        		
			        		expts.add(new ExptDescriptor("", "", cond, rep, signal, sources, perBaseReadLimit));			        		
			        		exptTags.add(s);
			        	}
		        	}
				}
		     
				//Parse experiment design file
				// Format: (tab separated)
				// SrcName	Signal/Control   DataType   Condition   Replicate	[ExptType]	[per-base max]	[Target]
				if(ap.hasKey("design")){
					String dfile = ap.getKeyValue("design");
					File df = new File(dfile);
					BufferedReader reader = new BufferedReader(new FileReader(df));
				    String line;
			        while ((line = reader.readLine()) != null) {
			        	if(!line.startsWith("#")){
				            line = line.trim();
				            String[] words = line.split("\\t");
				            if(words.length >=3){
				            	String cond="", rep="", expttype="", expttarget="";
					            
				            	boolean validLine = true;
					            Pair<String, String> src=null;
					            boolean signal= true;
					            if(words[0].toUpperCase().equals("SIGNAL") || words[0].toUpperCase().equals("CONTROL")){
					            	signal = words[0].toUpperCase().equals("SIGNAL") ? true : false;
					            	src = new Pair<String, String>(words[1], words[2]);
					            }else if(words[1].toUpperCase().equals("SIGNAL") || words[1].toUpperCase().equals("CONTROL")){
					            	signal = words[1].toUpperCase().equals("SIGNAL") ? true : false;
					            	src = new Pair<String, String>(words[0], words[2]);
					            }else{
					            	System.err.println("Incorrectly formatted line in design file:\n\t"+line+"\n");
					            	validLine=false;
					            }
					            
					            if(validLine){
					            	//Condition in field 3, replicate in field 4
						            if(words.length>=4 && !words[3].equals("")){
						            	cond = words[3];
						            }else
						            	cond = signal ? "experiment" :"DEFAULT";
						            
						            if(words.length>=5 &&  !words[4].equals("")){
						            	rep = words[4];
						            }else{
						            	rep = signal ? "rep1" : "DEFAULT";
						            }
						            
						            //Expt type in field 5 (optional)
						            if(words.length>=6 && !words[5].equals("") && !words[5].equals("NULL"))
						            	expttype = words[5];

						            //Per-base read limit in field 6
						            float currCondPerBaseReads = perBaseReadLimit;
						            if(words.length>=7 && !words[6].equals("") && !words[6].equals("NULL")){
						            	if(words[6].equals("P"))
						            		currCondPerBaseReads=0;
						            	else
						            		currCondPerBaseReads = new Float(words[6]);
						            }
						            
						            //Expt target in field 7
						            if(words.length>=8 && !words[7].equals("") && !words[7].equals("NULL"))
						            	expttarget = words[7];
						            
						            //Check if we have other entries for this experiment
						            boolean found=false;
						            for(ExptDescriptor e : expts){
						            	if(e.signal==signal && e.expttype.equals(expttype) && e.target.equals(expttarget) && e.condition.equals(cond) && e.replicate.equals(rep)){
						            		found = true;
						            		e.sources.add(src);
						            	}
						            }
						            if(!found){
						            	expts.add(new ExptDescriptor(expttype, expttarget, cond, rep, signal, src, currCondPerBaseReads));
						            }
					            }
				            }else{
				            	System.err.println("Incorrectly formatted line in design file:\n\t"+line+"\n");
				            }
			            }
			    	}
			        reader.close();
				}
				
				
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	//Accessors
	public Genome getGenome(){return gen;}
	public boolean getNonUnique(){return nonUnique;}
	public List<ExptDescriptor> getExperimentDescriptors(){return expts;}
	public boolean helpWanted(){return printHelp;}
	public boolean getPrintLoadingProgress(){return printLoadingProgress;}
	public boolean doPoissonGaussWinPerBaseFiltering(){return poissonGaussWinPerBaseFilter;}
	public boolean doPerBaseFiltering(){return perBaseReadFiltering;}
	public double getPerBaseLogConf(){return perBaseLogConf;}
	public float getPerBaseMax(){return perBaseReadLimit;}
	public double getMappableGenomeProp(){return mappableGenome;}
	public double getMappableGenomeLength(){return mappableGenome*gen.getGenomeLength();}
	public List<Integer> getLocalBackgroundWindows(){return localBackgroundWindows;}
	public boolean getEstimateScaling(){return estimateScaling;}
	public float getFixedScalingFactor(){return fixedScalingFactor;}
	public boolean getScalingByMedian(){return scalingByMedian;}
	public boolean getScalingByRegression(){return scalingByRegression;}
	public boolean getScalingBySES(){return scalingBySES;}
	public boolean getScalingByHitRatioAndNCIS(){return scalingByHitRatioAndNCIS;}
	public int getScalingSlidingWindow(){return scalingSlidingWindow;}
	public boolean getPlotScaling(){return plotScaling;}
	public boolean getCacheAllData(){return cacheAllHits;}
	public String getFileCacheDirName(){return fileCacheDir;}
	public List<Region> getInitialCachedRegions(){return initialCachedRegions;}
	public boolean getLoadType1Reads(){return loadType1Reads;}
	public boolean getLoadType2Reads(){return loadType2Reads;}
	public boolean getLoadRead2(){return loadRead2;}
	public boolean getLoadReads() {return loadReads;}
	public boolean getLoadPairs(){return loadPairs;}
	public boolean getKeepHDF5() {return keepHDF5;}
	public double getNCISMinBinFrac(){return NCISMinBinFrac;}
	
	//Some accessors to allow modification of options after config .
	public void setPrintLoadingProgress(boolean plp){printLoadingProgress = plp;}
	public void setPerBaseReadFiltering(boolean pbrf){perBaseReadFiltering = pbrf;}
	public void setPerBaseReadLimit(float pbrl){perBaseReadLimit = pbrl; perBaseReadFiltering=true;}
	public void setMedianScaling(boolean ms){scalingByMedian = ms;}
	public void setRegressionScaling(boolean rs){scalingByRegression = rs;}
	public void setSESScaling(boolean ses){scalingBySES = ses;}
	public void setScalingSlidingWindow(int ssw){scalingSlidingWindow = ssw;}
	public void setFileCacheDirName(String d){fileCacheDir = d;}
	public void setLoadType1Reads(boolean l){loadType1Reads = l;}
	public void setLoadType2Reads(boolean l){loadType2Reads = l;}
	public void setLoadRead2(boolean l){loadRead2 = l;}
	public void setLoadReads(boolean l){loadReads = l;}
	public void setLoadPairs(boolean l){loadPairs = l;}
	public void setSortMid(boolean l) {sortMid = l;}
	
	
	/**
	 * Merge a set of genomes that have been estimated from individual Samples 
	 * @param estGenomes
	 * @return
	 */
	public Genome mergeEstGenomes(List<Genome> estGenomes){
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
	
	
	/**
	 * returns a string describing the arguments handled by this config parser. 
	 * @return String
	 */
	public static String getArgsList(){
		return(new String("" +
				"Experiments:\n" +
				"\t--design <design file name>\n" +
				"\tOR\n" +
				"\t--expt/--ctrl <signal/control experiment file name> AND --format <SAM/BED/SCIDX/BOWTIE/NOVO>\n" +
				"\tAND/OR" +
				"\t--rdbexpt/--rdbctrl <signal/control ReadDB experiment identifier>\n" +
				"\t\tNote that if you use --expt/--ctrl or --rdbexpt/--rdbctrl, you can specify the names of the experiment & replicate\n" +
				"\t\tdirectly in the argument. Here's an example: --exptConditionA-Rep1 somefile.bam\n" +
				"Scaling control vs signal counts:\n" +
				"\t--noscaling [flag to turn off auto estimation of signal vs control scaling factor]\n" +
				"\t--medianscale [flag to use scaling by median ratio (default = scaling by NCIS)]\n" +
				"\t--regressionscale [flag to use scaling by regression (default = scaling by NCIS)]\n" +
				"\t--sesscale [flag to use scaling by SES (default = scaling by NCIS)]\n" +
				"\t--normalizedncisscale [flag to use scaling by total tag followed by NCIS (default = scaling by NCIS)]\n" +
				"\t--fixedscaling <multiply control counts by total tag count ratio and then by this factor if not estimating scaling>\n" +
				"\t--scalewin <window size for scaling procedure (default=10000)>\n" +
				"\t--plotscaling [flag to plot diagnostic information for the chosen scaling method]\n" +
				"Miscellaneous Experiment Loading Args:\n" +
				"\t--fixedpb <fixed per base limit>\n" +
				"\t--poissongausspb <filter per base using a Poisson threshold parameterized by a local Gaussian sliding window>\n" +
				"\t--mappability <fraction of the genome that is mappable for these experiments>\n" +
				"\t--nocache [flag to turn off caching of the entire set of experiments (i.e. run slower with less memory)]\n" +
				"\t--not1reads / --loadt2reads [flags to use Type1 or Type2 reads] (Type1 loaded by default)\n" +
				"\t--noread2 [flag to ignore second reads in paired-end]\n" +
				"\t--sortMid [flag to decide if sort read pairs by midpoint or 5' end (default: 5' end)]\n" +
				"\t--loadpairs [flag to load pair-end reads]\n" + 
				"\t--keepHDF5Cache" +
				""));
	}
}
