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
	protected double perBaseLogConf;
	protected boolean nonUnique=false;
	protected int scalingSlidingWindow = 10000;
	protected boolean estimateScaling = true;
	protected boolean scalingByMedian =false; // default is scaling by regression
	protected int factorWindow;
	protected int chromatinWindow;
	protected double mappableGenome = 0.8;
	protected File locations;
	protected int numChromStates;
	protected int numBindingStates;
	protected boolean plotEM;
	protected int numItrs;        // number of EM training steps
	protected String out_name;
	protected String out_base;
	protected File out_dir=null, images_dir=null;
	protected boolean loadReads;
	
	
	
	
	public final int W_MARGIN = 80;
	public final int H_MARGIN = 60;
	public final int W = 800;
	public final int H = 800;
	
	
	
	
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
				
				this.plotEM = Args.parseFlags(args).contains("plotEM");
				
				//Reading the global win sizes if provided in the command line
				//These values will be over-ridded if they are also provided in the design file
				this.chromatinWindow = Args.parseInteger(this.args, "chromwin", 500);
				this.factorWindow = Args.parseInteger(this.args, "facwin", 50);
				
				this.numChromStates = Args.parseInteger(args, "nChrmStates", 4);
				this.numBindingStates = Args.parseInteger(args, "nFacStates", 3);
				
				this.numItrs = Args.parseInteger(args, "numItrs", 50);
				
				
				
				//Load expts from design file
				//Format: (tab separated)
				//Signal/Control   SrcName   Type   Condition   Replicate [Chromatin/Factor default is chromatin state] [per-base max] [window size]
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
					            // Window regions to consider for this feature
					            int winsize = words[4].toUpperCase().equals("CHROMATIN") ? this.chromatinWindow : this.factorWindow;
					            if(words.length>=7  && !words[7].equals("")){
					            	winsize = new Integer(words[7]);
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
					            	expts.add(new ExptDescriptor(cond, rep, feature,  signal, src, currCondPerBaseReads,winsize));
					            }
				            }else{
				            	System.err.println("Error in design file. Cannot parse line: \n"+line);
				            }
			            }
			    	}
			        reader.close();
				}
				
				//Checking if the peaks/locations file exists or not (stops the program if it does not exist)
				String peaksFileName = ap.getKeyValue("peaks");
				this.locations = new File(peaksFileName);
				if(!locations.exists()){
					System.err.println("Peaks File does not exist\n");
					System.exit(1);    
				}
				
				/****Miscellaneous arguments****/
				
				//Turn off scaling estimation
				estimateScaling = Args.parseFlags(args).contains("noscaling") ? false : true;
				//Scale by median or regression
				scalingByMedian = Args.parseFlags(args).contains("medianscale") ? true : false;
				//bacground parameters 
				perBaseLogConf = Args.parseDouble(args,"pblogconf",-7);
				
				this.loadReads = Args.parseFlags(args).contains("noreads") ? false : true;
				
				out_name = Args.parseString(args, "out", "bayesments_out");
				out_dir =  new File(out_name); //Output directory
				out_base = out_dir.getName(); //Last part of name
				
				
				
				
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
	public double getPerBaseLogConf(){return this.perBaseLogConf;}
	public double getMappableGenomeProp(){return mappableGenome;}
	public File getPeaksFile(){return this.locations;}
	public int getNumChrmStates(){return this.numChromStates;}
	public int getNumFacStates(){return this.numBindingStates;}
	public boolean doEMplot(){return this.plotEM;}
	public int getNumItrs(){return this.numItrs;}
	public String getOutName(){return out_name;}
	public String getOutBase(){return out_base;}
	public File getOutputParentDir(){return out_dir;}
	public File getOutputImagesDir(){return images_dir;}
	public boolean loadReads(){return loadReads;}
	
	
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
	
	public void makeGPSOutputDirs(boolean makeImageDirs){
		//Test if output directory already exists. If it does,  recursively delete contents
		out_dir =  new File(out_name);
		if(out_dir.exists())
			deleteDirectory(out_dir);
		out_base = out_dir.getName();
		//(re)make the output directory
		out_dir.mkdirs();
		if(makeImageDirs){
			//Make the image results output directory
			images_dir = new File(out_dir.getAbsolutePath()+File.separator+"images");
			images_dir.mkdirs();
		}
	}
	
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
				"\t--chromwin <widow size around factor binding site for a chromatin expt>\n"+
				"\t--facwin <window size around the factor binding site>\n"+
				"\t--numItrs <num of EM steps>\n"+
				"\t--nchrmStates <number of chromatin states to be learned>\n"+
				"\t--nfacStates <number of factor states to be learned>\n"+
				"\t--pblogconf <per base log confidence for the background model>\n"+
				"\t--peaks <file location for the genomic locations>\n"+
				"\t--plotEM <flag to plot the EM steps>\n"+
				"\t--out <name of the output directory>\n"+
				"\t--noreads <flag to turn off loading rads>\n"+
				"\t--poissongausspb <flag to filter per base using Poisson Gaussian sliding window> (overrides --fixedpb)"));
	}
	

}
