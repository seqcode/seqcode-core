package edu.psu.compbio.seqcode.deepseq.experiments;

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

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.Pair;

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
	protected List<Integer> localBackgroundWindows=new ArrayList<Integer>(); 
	protected double perBaseLogConf=-7;
	protected boolean poissonGaussWinPerBaseFilter = false; //Filter using a poisson Gaussian window
	protected float perBaseReadLimit = -1;
	protected boolean perBaseReadFiltering = true;
	protected double mappableGenome = 0.8;
	protected boolean estimateScaling=true;
	protected boolean scalingByMedian = false; //Default is to estimate scaling by regression
	protected boolean scalingBySES = false; //Default is to estimate scaling by regression
	protected int scalingSlidingWindow = 10000; 
	protected boolean cacheAllHits=true; //Cache all hits
	protected String fileCacheDir = "hitcache";
	protected List<Region> initialCachedRegions=null;
	
	    
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
				//Scale by median or regression
				scalingByMedian = Args.parseFlags(args).contains("medianscale") ? true : false;
				//Scale by SES or regression
				scalingBySES = Args.parseFlags(args).contains("sesscale") ? true : false;
				//Scaling window
				scalingSlidingWindow = Args.parseInteger(args,"scalewin",scalingSlidingWindow);
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
	public boolean doPoissonGaussWinPerBaseFiltering(){return poissonGaussWinPerBaseFilter;}
	public boolean doPerBaseFiltering(){return perBaseReadFiltering;}
	public double getPerBaseLogConf(){return perBaseLogConf;}
	public double getMappableGenomeProp(){return mappableGenome;}
	public double getMappableGenomeLength(){return mappableGenome*gen.getGenomeLength();}
	public List<Integer> getLocalBackgroundWindows(){return localBackgroundWindows;}
	public boolean getEstimateScaling(){return estimateScaling;}
	public boolean getScalingByMedian(){return scalingByMedian;}
	public boolean getScalingBySES(){return scalingBySES;}
	public int getScalingSlidingWindow(){return scalingSlidingWindow;}
	public boolean getCacheAllData(){return cacheAllHits;}
	public String getFileCacheDirName(){return fileCacheDir;}
	public List<Region> getInitialCachedRegions(){return initialCachedRegions;}
	
	//Some accessors to allow modification of options after config .
	public void setPerBaseReadFiltering(boolean pbrf){perBaseReadFiltering = pbrf;}
	public void setMedianScaling(boolean ms){scalingByMedian = ms;}
	public void setSESScaling(boolean ses){scalingBySES = ses;}
	public void setScalingSlidingWindow(int ssw){scalingSlidingWindow = ssw;}
	public void setFileCacheDirName(String d){fileCacheDir = d;}
	
	
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
	public String getArgsList(){
		return(new String("" +
				"Experiments:\n" +
				"\t--design <design file name>\n" +
				"\tOR\n" +
				"\t--expt <read file name> AND --format <SAM/BED/IDX/BOWTIE/NOVO>\n" +
				"\tAND/OR" +
				"\t--rdbexpt <ReadDB experiment identifier>\n" +
				"Miscellaneous Experiment Loading Args:\n" +
				"\t--fixedpb <fixed per base limit>\n" +
				"\t--poissongausspb <filter per base using a Poisson threshold parameterized by a local Gaussian sliding window>\n" +
				"\t--mappability <fraction of the genome that is mappable for these experiments>\n" +
				"\t--nocache [flag to turn off caching of the entire set of experiments (i.e. run slower with less memory)]\n" +
				""));
	}
}
