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
	protected File Top_Locations;
	protected int numChromStates;
	protected int numBindingStates;
	protected boolean plotEM;
	protected int numItrs;        // number of EM training steps
	protected String out_name;
	protected String out_base;
	protected File out_dir=null, images_dir=null;
	protected boolean loadReads;
	protected boolean doashin;
	protected boolean doSeqAshin;
	protected File inter_dir=null;
	protected String genomeSequencePath;
	protected boolean onlyChrom;
	protected boolean regularize;
	
	protected String MEMEpath="";
	protected String MEMEargs=" -dna -mod zoops -revcomp -nostatus ";
	protected int MEMEminw=6;
	protected int MEMEmaxw=18;
	protected int MEMEnmotifs;
	
	protected double chromatin_weightage;
	protected double sequence_weightage;
	protected boolean boundSigma;
	protected boolean Simulate_reads = false;
	protected double lambda;
	protected int buffer_Sigma;
	
	protected int N;
	
	protected int minChromStates;
	protected int maxChromStates;
	
	protected int minFacStates;
	protected int maxFacStates;

	
	
	
	
	public final int W_MARGIN = 80;
	public final int H_MARGIN = 60;
	public final int W = 800;
	public final int H = 800;
	public final double Motif_EVal_Cutoff = Double.valueOf("1.0e-10");
	public final int MOTIF_FINDING_NEGSEQ=5000; //Number of negative sequences for motif significance tests
    public final double MOTIF_FINDING_ALLOWED_REPETITIVE = 0.2; //Percentage of the examined sequence window allowed to be lowercase or N
    public final double MOTIF_MIN_ROC = 0.7; //Motif prior is used only if the ROC is greater than this .
	
	
	
	
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
				//IMPORTANT: The factorWindow is also chosen for extracting the sequence  around the peak pair locaitons
				this.chromatinWindow = Args.parseInteger(this.args, "chromwin", 500);
				this.factorWindow = Args.parseInteger(this.args, "facwin", 50);
				
				this.numChromStates = Args.parseInteger(args, "nChrmStates", 4);
				this.numBindingStates = Args.parseInteger(args, "nFacStates", 3);
				
				this.numItrs = Args.parseInteger(args, "numItrs", 50);
				
				//Load expts from design file
				//Format: (tab separated)
				//Signal/Control   SrcName   Type   Condition   Replicate [Chromatin/Factor default is chromatin state] [per-base max] [window size]
				// If the desing option is not given then it looks for the simulate flag... (If the simulate and design options are not given the program stops)
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
				else if(Args.parseFlags(args).contains("simulate")){
					this.Simulate_reads = true;
				}
				else{
					System.err.println("You have to either provide a desing file using the \"--design\" options or turn on the simulate datasets flag using \"--simulate\"");
					System.exit(1);
				}
				
				//Checking if the peaks/locations file exists or not (stops the program if it does not exist)
				if(!this.Simulate_reads){
					String peaksFileName = ap.getKeyValue("peaks");
					this.locations = new File(peaksFileName);
					if(!locations.exists()){
						System.err.println("Peaks File does not exist\n");
						System.exit(1);    
					}
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
				this.doashin = Args.parseArgs(args).contains("ashin") ? true : false;
				this.doSeqAshin = Args.parseFlags(args).contains("ashinSeq")? true : false;
				onlyChrom = Args.parseFlags(args).contains("onlyChrom") ? true : false;
				regularize = Args.parseFlags(args).contains("regularize") ? true : false;
				buffer_Sigma = Args.parseInteger(args, "bufferSigma", 0);
				this.minChromStates = Args.parseInteger(args, "minChromStates", 1);
				this.maxChromStates = Args.parseInteger(args, "maxChromStates", 10);
				this.minFacStates = Args.parseInteger(args, "minFacStates", 1);
				this.maxFacStates = Args.parseInteger(args, "maxFacStates", 5);
				
				
				this.lambda = Args.parseDouble(args, "regWeight", 10);
				
				if(!onlyChrom){ //Do we need to load sequences?
					genomeSequencePath = ap.hasKey("seq") ? ap.getKeyValue("seq") : null;
					if(genomeSequencePath==null){
						System.err.println("You have requested to include sequence features, but no genome sequence data was provided with --seq");
						System.exit(1);
					}
					
					//Read the locations to do meme search provided by the user
					String top_locs = ap.getKeyValue("meme_search_regions");
					this.Top_Locations = new File(top_locs);
					if(!this.Top_Locations.exists()){
						System.err.println("Top locations File does not exist\n");
						System.exit(1);    
					}
				}
				
				this.chromatin_weightage = Args.parseDouble(args, "chromWeight", 1.0);
				this.sequence_weightage = Args.parseDouble(args, "seqWeight", 1.0);
				
				
				
				boundSigma = Args.parseFlags(args).contains("capSigma") ? true : false;
				
				
				//MEME path
				MEMEpath = Args.parseString(args, "memepath", MEMEpath);
				//MEME args
				MEMEargs = Args.parseString(args, "memeargs", MEMEargs);
				//MEME minw
				MEMEminw = Args.parseInteger(args, "mememinw", MEMEminw);
				//MEME maxw
				MEMEmaxw = Args.parseInteger(args, "mememaxw", MEMEmaxw);
				//MEME nmotifs option
				MEMEnmotifs = Args.parseInteger(args,"memenmotifs", 3);
				MEMEargs = MEMEargs + " -nmotifs "+MEMEnmotifs + " -minw "+MEMEminw+" -maxw "+MEMEmaxw;
				
				
				
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
	public File getOutputInterDir(){return inter_dir;}
	public boolean loadReads(){return loadReads;}
	public boolean doChipAshin(){return this.doashin;}
	public boolean doSeqAshin(){return this.doSeqAshin;}
	public int getSeqWinSize(){return factorWindow;}
	public String getGenomeSeqPath(){return genomeSequencePath;}
	public boolean runOnlyChrom(){return this.onlyChrom;}
	public String getMEMEpath(){return this.MEMEpath;}
	public String getMEMEargs(){return this.MEMEargs;}
	public double getChromWeight(){return this.chromatin_weightage;}
	public double getSeqWeight(){return this.sequence_weightage;}
	public boolean capSigma(){return this.boundSigma;}
	public boolean doRegularization(){return this.regularize;}
	public double getLambda(){return this.lambda;}
	public int getBufferSigmaVal(){return this.buffer_Sigma;}
	public boolean doSimulation(){return this.Simulate_reads;}
	public int getMaxChromStates(){return this.maxChromStates;}
	public int getMinChromStates(){return this.minChromStates;}
	public int getMaxFacStates(){return this.maxFacStates;}
	public int getMinFacStates(){return this.minFacStates;}
	public double getMotifCuttof(){return this.Motif_EVal_Cutoff;}
	public File getTopLocations(){return this.Top_Locations;}
	
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
		
		inter_dir = new File(out_dir.getAbsolutePath()+File.separator+"intermediate-results");
		inter_dir.mkdirs();
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
				"\t--numItrs <num of EM steps: if onlyChrom flag is off, this option is invalid>\n"+
				"\t--numChromItrs <num of chromatin feature interations after which seq features are to be included: invalid if onlyChrom flag is on>\n"+
				"\t--numSeqItrs <num of iterations to run after the sequence features are included: invalid if onlyChrom flag is on>\n"+
				"\t--nchrmStates <number of chromatin states to be learned>\n"+
				"\t--nfacStates <number of factor states to be learned>\n"+
				"\t--pblogconf <per base log confidence for the background model>\n"+
				"\t--peaks <file location for the genomic locations>\n"+
				"\t--plotEM <flag to plot the EM steps>\n"+
				"\t--out <name of the output directory>\n"+
				"\t--noreads <flag to turn off loading rads>\n"+
				"\t--ashin <flag to turn on ashin transcormation for chip experiments>\n"+
				"\t--ashinSeq <flag to turn on ashin transformation on seq scores>\n"+
				"\t--seq<Path to seq fasta dir>\n"+
				"\t--regWeight <(lambda) weight for regularaization>\n"+
				"\t--onlyChrom<flag to turn only only chromatin tracks as features>\n"+
				"\t--regularize<flag to turn on regularization of features>\n"+
				"\t--memepath <path to the meme bin dir (default: meme is in $PATH)>\n" +
				"\t--memenmotifs <number of motifs MEME should find for each condition>\n" +
				"\t--mememinw <minw arg for MEME (default="+MEMEminw+")>\n"+
				"\t--mememaxw <maxw arg for MEME (default="+MEMEmaxw+")>\n"+
				"\t--memenmotifs <number of motifs you want to include>\n"+
				"\t--capSigma <flag to turn on an upperbound cap on sigma for each guassian>\n"+
				"\t--chromWeight <weightage for chromatin features to take responsibility for state>\n"+
				"\t--seqWeight<weightage for seq features to take responsibility for state>\n"+
				"\t--bufferSigma<int for testing purpose>\n"+
				"\t--minChromStates<min no of chrom states to test>\n"+
				"\t--maxChromStates<max no of chrom states to test>\n"+
				"\t--minFacStates<>\n"+
				"\t--maxFacStates<>\n"+
				"\t--poissongausspb <flag to filter per base using Poisson Gaussian sliding window> (overrides --fixedpb)"));
	}
	

}
