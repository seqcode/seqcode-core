package edu.psu.compbio.seqcode.deepseq.experiments;

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

import edu.psu.compbio.seqcode.deepseq.experiments.ExptDescriptor;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

/**
 * ExptConfig: 
 * 		Handles specification of experiments and read loaders from command-line.
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
	protected double mappableGenome = 0.8;
	protected List<Integer> localBackgroundWindows=new ArrayList<Integer>(); 
	protected double perBaseLogConf=-7;
	protected boolean poissonGaussWinPerBaseFilter = false; //Filter using a poisson Gaussian window
	protected float perBaseReadLimit = -1;
	protected boolean perBaseReadFiltering = true;
	    
	protected String[] args;
	public String getArgs(){
		String a="";
		for(int i=0; i<args.length; i++)
			a = a+" "+args[i];
		return a;
	}
	
	public ExptConfig(String[] arguments){
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
				}
				
				//Load genome
				if(ap.hasKey("species") || ap.hasKey("genome") || ap.hasKey("gen")){
					Pair<Organism, Genome> pair = Args.parseGenome(args);
					if(pair != null){
						gen = pair.cdr();
					}
				}else{
					if(ap.hasKey("geninfo") || ap.hasKey("g")){
						//Make fake genome... chr lengths provided
						String fName = ap.hasKey("geninfo") ? ap.getKeyValue("geninfo") : ap.getKeyValue("g");
						gen = new Genome("Genome", new File(fName), true);
					}else{
					    gen = null;
					}
				}
				
				String modelFile = Args.parseString(args, "d", null);	// read distribution file
//				if (modelFile == null){
//					defaultModel = new BindingModel(BindingModel.defaultEmpiricalDistribution);
//		        } else {
//					File pFile = new File(modelFile);
//					if(!pFile.isFile()){
//						System.err.println("\nCannot find read distribution file: "+modelFile);
//						System.exit(1);
//					}
//					defaultModel = new BindingModel(pFile);
//				}
				
				//Fixed per base read limits
				perBaseReadLimit = Args.parseInteger(args,"fixedpb",-1);
				//Do per base filtering with a Poisson Gaussian window
				poissonGaussWinPerBaseFilter = Args.parseFlags(args).contains("poissongausspb");
				//Switch off per-base read filtering?
				perBaseReadFiltering = !Args.parseFlags(args).contains("nopbfilter");
				
				//Parse command-line experiments
				String fileFormat = Args.parseString(args, "format", "SAM").toUpperCase();
				ArrayList<String> exptTags = new ArrayList<String>();
				for(String s : ap.getKeys()){
					if(!exptTags.contains(s)){
			        	if(s.contains("expt") || s.contains("ctrl")){
			        		String type = "";
			        		String condrep = "";
			        		boolean signal = true;
			        		if(s.contains("ctrl"))
			        			signal = false;
			        		if(s.startsWith("rdb")){
			        			type = "READDB";
			        			condrep = s.replaceFirst("rdbexpt", "");
			        			condrep = condrep.replaceFirst("rdbctrl", "");
			        		}else{
			        			type = fileFormat;
			        			condrep = s.replaceFirst("expt", "");
			        			condrep = condrep.replaceFirst("ctrl", "");
			        		}
			        		//Parse the condition & replicate names
			        		String cond = condrep.split("-")[0];
			        		if(cond.length()==0 || cond.equals(""))
			        			cond = signal ? "EXPERIMENT" : "DEFAULT";
			        		String rep = condrep.split("-").length>1 ? condrep.split("-")[1] : "";
			        		if(rep.length()==0 || rep.equals("")){
			        			rep = signal ? "Rep1" : "DEFAULT";
			        		}
			        		//Parse the file/rdb names
			        		Collection<String> sourceNames = Args.parseStrings(args, s);
			        		List<Pair<String,String>> sources = new ArrayList<Pair<String,String>>();
			        		for(String n : sourceNames)
			        			sources.add(new Pair<String,String>(n,type));
			        		
//			        		expts.add(new ExptDescriptor(cond, rep, signal, sources, defaultModel, perBaseReadLimit));
			        		
			        		exptTags.add(s);
			        	}
		        	}
				}
		     
				//Parse experiment design file
				// Format: (tab separated)
				// Signal/Control   SrcName   Type   Condition   Replicate	[Read Distribution]	[per-base max]
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
				            	String cond="", rep="";
//				            	BindingModel currModel=null;
					            
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
						            if(words.length>=4 && !words[3].equals("")){
						            	cond = words[3];
						            }else
						            	cond = signal ? "EXPERIMENT" :"DEFAULT";
						            
						            if(words.length>=5 &&  !words[4].equals("")){
						            	rep = words[4];
						            }else{
						            	rep = signal ? "Rep1" : "DEFAULT";
						            }
						            
						            //Read distribution in field 6
						            if(words.length>=6 && !words[5].equals("") && !words[5].equals("NULL")){
						            	String currModFile = words[5];
						            	File pFile = new File(currModFile);
										if(!pFile.isFile()){
											System.err.println("\nCannot find read distribution file: "+currModFile);
											System.exit(1);
										}
	//									currModel = new BindingModel(pFile);
						            }
	//					            if(currModel==null && signal)
		//				            	currModel = defaultModel;
						            
						            //Per-base read limit in field 7
						            float currCondPerBaseReads = perBaseReadLimit;
						            if(words.length>=7 && !words[6].equals("")){
						            	if(words[6].equals("P"))
						            		currCondPerBaseReads=0;
						            	else
						            		currCondPerBaseReads = new Float(words[6]);
						            }
						            
						            //Check if we have other entries for this experiment
						            boolean found=false;
						            for(ExptDescriptor e : expts){
						            	if(e.signal==signal && e.condition.equals(cond) && e.replicate.equals(rep)){
						            		found = true;
						            		e.sources.add(src);
						            	}
						            }
						            if(!found){
//						            	expts.add(new ExptDescriptor(cond, rep, signal, src, currModel, currCondPerBaseReads));
						            }
					            }
				            }else{
				            	System.err.println("Incorrectly formatted line in design file:\n\t"+line+"\n");
				            }
			            }
			    	}
				}
				
				
			} catch (NotFoundException e) {
				e.printStackTrace();
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
	public boolean getNonUnique(){return nonUnique;}
	public List<ExptDescriptor> getExperiments(){return expts;}
	public boolean helpWanted(){return printHelp;}
	
	//Some accessors to allow modification of options after config .
	public void setPerBaseReadFiltering(boolean pbrf){perBaseReadFiltering = pbrf;}
	
	
	
	/**
	 * returns a string describing the arguments handled by this parser. 
	 * @return String
	 */
	public String getArgsList(){
		return(new String("" +
				"Genome:" +
				"\t--species <Organism;Genome>\n" +
				"\tOR\n" +
				"\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" +
				"Experiments:\n" +
				"\t--design <file name>\n" +
				"\tOR\n" +
				"\t--expt <read file name> AND --format <SAM/BED/IDX/BOWTIE/NOVO/ELAND>\n" +
				"AND/OR" +
				"\t--rdbexpt <ReadDB experiment identifier>\n" +
				"Miscellaneous Experiment Loading Args:\n" +
				"\t--fixedpb <fixed per base limit>\n" +
				"\t--poissongausspb <filter per base using a Poisson threshold parameterized by a local Gaussian sliding window>\n" +
				""));
	}
}
