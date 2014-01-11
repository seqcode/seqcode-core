package edu.psu.compbio.seqcode.projects.akshay.bayesments.framework;

import java.io.*;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExptDescriptor;


public class Config {
	
	protected Genome gen;
	protected List<ExptDescriptor> expts =  new ArrayList<ExptDescriptor>();
	protected String[] args;
	protected boolean printHelp = false;
	protected boolean poissonGaussWinPerBaseFilter = false;
	protected float perBaseReadLimit;
	protected boolean nonUnique=false;
	protected int scalingSlidingWindow = 10000;
	protected boolean estimateScaling = true;
	protected boolean scalingByMedian =false; // default is scaling by regression
	
	public Config(String[] arguments){
		this.args = arguments;
		ArgParser ap =  new ArgParser(args);
		if(ap.hasKey("h") || args.length == 0){
			printHelp = true;
		}
		else{
			try{
				//Load genomes
				if(ap.hasKey("species") || ap.hasKey("genome") || ap.hasKey("gen")){
					Pair<Organism,Genome> orggen = Args.parseGenome(this.args);
					this.gen = orggen.cdr();
				}
				else{
					if(ap.hasKey("geneinfo") || ap.hasKey("g")){
						String infofilename = ap.hasKey("geneinfo") ? ap.getKeyValue("geneinfo") : ap.getKeyValue("g"); 
						this.gen = new Genome("Genome",new File(infofilename),true);
					}
					else{
						this.gen = null;
					}
				}
				
				//Reading per base pair read limits
				// having the --poissongausspb overrides --fixedpb <float>
				//having --poissongausspb is same as having --fixedpb 0
				this.poissonGaussWinPerBaseFilter = Args.parseFlags(this.args).contains("poissongausspb");
				this.perBaseReadLimit = Args.parseFloat(this.args,"fixedpb" , -1);
				
				//Load expts from design file
				//Format: (tab separated)
				//Signal/Control   SrcName   Type   Condition   Replicate [Chromatin/Factor default is chromatin state] [per-base max]
				if(ap.hasKey("design")){
					String dfile = ap.getKeyValue("design");
					File df = new File(dfile);
					BufferedReader reader = new BufferedReader(new FileReader(df));
				    String line;
			        while ((line = reader.readLine()) != null) {
			        	if(!line.startsWith("#")){
				            line = line.trim();
				            String[] words = line.split("\\t");
				            if(words.length >=6){
				            	String cond="", rep="", feature="";
					            boolean signal = words[0].toUpperCase().equals("SIGNAL") ? true : false;
					            Pair<String, String> src = new Pair<String, String>(words[1], words[2]);
					            if(!words[3].equals("")){
					            	cond = words[3];
					            }else
					            	cond = signal ? "EXPERIMENT" :"DEFAULT";
					            
					            if(!words[4].equals("")){
					            	rep = words[4];
					            }else{
					            	rep = signal ? "Rep1" : "DEFAULT";
					            }
					            if(!words[5].equals("")){
					            	feature = words[5].toUpperCase();
					            }else{
					            	feature = signal ? "CHROMATIN" : "DEFAULT";
					            }
					            //Per-base read limit in field 6
					            float currCondPerBaseReads = this.perBaseReadLimit;
					            if(words.length>=6 && !words[6].equals("")){
					            	if(words[5].equals("P"))
					            		currCondPerBaseReads=0;
					            	else
					            		currCondPerBaseReads = new Float(words[6]);
					            }
					            //Check if we have other entries for this experiment
					            boolean found=false;
					            for(ExptDescriptor e : this.expts){
					            	if(e.signal==signal && e.condition.equals(cond) && e.replicate.equals(rep)){
					            		found = true;
					            		e.sources.add(src);
					            	}
					            }
					            if(!found){
					            	expts.add(new ExptDescriptor(cond, rep, feature,  signal, src, currCondPerBaseReads));
					            }
				            }else{
				            	System.err.println("Error in design file. Cannot parse line: \n"+line);
				            }
			            }
			    	}
			        reader.close();
				}
				
				/****Miscellaneous arguments****/
				
				//Turn off scaling estimation
				estimateScaling = Args.parseFlags(args).contains("noscaling") ? false : true;
				//Scale by median or regression
				scalingByMedian = Args.parseFlags(args).contains("medianscale") ? true : false;
				
				
				
				
			}
			catch (NotFoundException e){
				e.printStackTrace();
			}
			catch (FileNotFoundException e){
				e.printStackTrace();
			}
			catch (IOException e){
				e.printStackTrace();
			}
			
			
		}
	}
	
	
	
	// List of all Accessors
	
	public Genome getGenome(){return this.gen;}
	public List<ExptDescriptor> getExperiments(){return expts;}
	public boolean helpWanter(){return printHelp;}
	public float getPerBaseReadLimit(){return this.perBaseReadLimit;}
	public boolean doPoissonGaussWinPerBaseFiltering(){return this.poissonGaussWinPerBaseFilter;}
	public boolean getNonUnique(){return nonUnique;}
	public int getScalingSlidingWindow(){return scalingSlidingWindow;}
	public boolean getEstimateScaling(){return estimateScaling;}
	public boolean getScalingByMedian(){return scalingByMedian;}
	
	//Some accessors to allow modification of options after config
	public void setScalingSlidingWindow(int ssw){scalingSlidingWindow = ssw;}
	
	
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
		this.gen =new Genome("Genome", chrLenMap);
		return this.gen;		
	}
	
	
	public String getArgsList(){
		return(new String("" +
				"Genome:" +
				"\t--species|genome|gen <Organism;Genome>\n" +
				"\tOR\n" +
				"\t--geninfo|g <genome info file> AND --seq <fasta seq directory>\n"+
				"\t--fixedpb <fixed per base limit>\n"+
				"\t--design <file name>\n"+
				"\t--noscaling [flag to turn off signal vs control scaling]\n" +
				"\t--medianscale [flag to use scaling by median (default = regression)]\n" +
				"\t--poissongausspb <flag to filter per base using Poisson Gaussian sliding window> (overrides --fixedpb)"));
	}
	

}
