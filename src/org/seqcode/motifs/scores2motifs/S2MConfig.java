package org.seqcode.motifs.scores2motifs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.RepeatMaskedRegion;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gsebricks.verbs.location.RepeatMaskedGenerator;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;


/**
 * @author akshaykakumanu
 * @twitter ikaka89
 * @email auk262@psu.edu
 */
public class S2MConfig implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public static String version = "0.1.3";
	
	// General options
	protected GenomeConfig gcon;
	protected SequenceGenerator<Region> seqgen = null;
	protected RepeatMaskedGenerator<Region> repMask;
	protected String[] args;
	/** Number K-mers in the model */
	protected int numK;
	protected int minK=4;
	protected int maxK=5;
	
	/** Execution option 1: load a list of scored points */
	protected String peaksFile=null; 
		
	/** Minimum length to consider for motif finding */
	protected int minM=6;
	/** Maximum length for motif finding */
	protected int maxM=10;
	/** The base name of the output directory */
	protected String outbase;
	/** The output directory file */
	protected File outdir;
	/** The minimum value of model scan score to consider form motif finding */
	protected double thresold_hills = 0.1;
	/** The number of hills in each cluster to consider for MEME motif finding */
	public static final int NUM_HILLS = 750;
	/** No of iterations of clustering */
	public static final int ITRS_CLUS=10;
	/** minimum ROC required to report motif*/
	protected double motifMinROC=0.7;
	/** Flag to mask repeats */
	protected boolean screenReps=false;
	protected double repPropLimit=0.5;
	protected int repMaskWin=100; //window around points to look for repeats. 
	
	// Meme parameters
	protected String MEMEpath;
	protected String MEMEargs = " -dna -mod zoops -revcomp -nostatus ";
	protected int MEMEminw = 6;
	protected int MEMEmaxw = 11;
	protected int MEMEnmotifs = 3;
	protected int MEMEwin = 16;
	public static final double MOTIF_FINDING_ALLOWED_REPETITIVE = 0.2;
	
	public Genome getGenome(){return gcon.getGenome();}
	public SequenceGenerator<Region> getSeqGen(){return seqgen;}
	public String getPeaksFile(){return peaksFile;}
	public int getKmin(){return minK;}
	public int getKmax(){return maxK;}
	public int getNumK(){return numK;}
	public File getOutDir(){return outdir;}
	public String getOutbase(){return outbase;}
	public String getMemeArgs(){String memeargs = MEMEargs+" -nmotifs "+MEMEnmotifs + " -minw "+MEMEminw+" -maxw "+MEMEmaxw; return memeargs;}
	public String getMemePath(){return MEMEpath;}
	public int getMemeSearchWin(){return MEMEwin;}
	public int getMinM(){return minM;}
	public int getMaxM(){return maxM;}
	public double getMotifMinROC(){return motifMinROC;}
	public int getKmerBaseInd(String kmer){
		int baseInd = 0;
		for(int k=minK; k<kmer.length(); k++){
			baseInd += (int)Math.pow(4, k);
		}
		return baseInd;
	}
	public double getHillsThresh(){return thresold_hills;}
	public boolean getScreenReps(){return screenReps;}
	public double getRepPropLimit(){return repPropLimit;}
	public int getRepMaskWin(){return repMaskWin;}
	
	public RepeatMaskedGenerator<Region> getRepMask(){return repMask;}
	public String getVersion(){return version;}

	public S2MConfig(String[] arguments) throws IOException {
		// Loading general options
		args = arguments;
		ArgParser ap = new ArgParser(args);
		if(ap.hasKey("h") || ap.hasKey("help") || args.length == 0){
			System.err.println(S2MConfig.getArgsList());
			System.exit(1);
		}
		gcon = new GenomeConfig(args);
		seqgen = gcon.getSequenceGenerator();
		
		// First load all options needed for making the arff file
		minK = Args.parseInteger(args, "mink", 4);
		maxK = Args.parseInteger(args, "maxk", 5);

		// Get outdir and outbase and make them; delete dirs that exist with the same
		outbase = Args.parseString(args, "out", "seqUnwinder_out");
		outdir = new File(outbase);
		if(!ap.hasKey("fasta"))
			makeOutputDirs();

		numK = 0;
		for(int k=minK; k<=maxK; k++ ){
			numK += (int)Math.pow(4, k);
		}

		// Load peaks and annotations
		if(!ap.hasKey("points")){
			System.err.println("Please provide genomic locations with scores and try again !!");
			S2MConfig.getArgsList();
			System.exit(1);
		}else if(ap.hasKey("points")){
			// Reading peaks files and storing annotations
			peaksFile = ap.getKeyValue("points");
		}
		//support for other input types would go here. 
		
		
		if(ap.hasKey("screenrepeats")){
			screenReps = true;
			repMask = new RepeatMaskedGenerator<Region>(gcon.getGenome()); 
		}
	
		

		// Load all MEME arguments
		// Path to MEME binary
		MEMEpath = Args.parseString(args, "memepath", "");
		MEMEargs = Args.parseString(args, "memeargs", MEMEargs);
		
		MEMEminw = Args.parseInteger(args, "mememinw", 6);
		MEMEmaxw = Args.parseInteger(args, "mememaxw", 13);
		//Size of the focussed meme search win
		MEMEwin = Args.parseInteger(args, "memesearchwin", 16);
		MEMEnmotifs = Args.parseInteger(args, "memenmotifs", 3);
		//Minimum ROC to report motifs
		motifMinROC =  Args.parseDouble(args, "motifminROC", 0.7);

		// Load arguments for motif-finding analysis
		minM = Args.parseInteger(args, "minscanlen", 6);
		maxM = Args.parseInteger(args, "maxscanlen", 14);
		thresold_hills = Args.parseDouble(args, "hillsthresh", 0.1);

	}
	


	public void makeOutputDirs(){
		//Test if output directory already exists. If it does,  recursively delete contents
		if(outdir.exists())
			deleteDirectory(outdir);
		//(re)make the output directory
		outdir.mkdirs();
		
	}

	public static String getArgsList(){
		return(new String("Scores2Motifs" +
				"\n OPTIONS:\n" +
				" General:\n"+
				"\t--out <prefix>: Ouput file prefix. All output will be put into a directory with the prefix name\n" +
				"\t--threads <n>: Use n threads to train SeqUnwinder model. Default is 5 threads\n" +
				"\t--debug: Flag to run in debug mode; prints extra output\n" +
				"\t--memepath <path>: path to the meme bin dir (default: meme is in $PATH)\n" +
				" Specify the genome:\n" +
				"\t--geninfo <genome info file> This file should list the lengths of all chromosomes on separate lines using the format chrName<tab>chrLength\n" + 
				"\t\tAND\n" +  
				"\t--seq <path>: A directory containing fasta format files corresponding to every named chromosome is required\n" +
				"\t--mink <int>: Minimum length of k-mer (default = 4)\n" + 
				"\t--maxk <int>: Maximum length of k-mer (default = 5)\n" + 
				" Other SeqUnwinder options (Highly recommend using defaul options): \n"+
				"\t--minscanlen <value>: Minimum length of the window to scan K-mer models. Default=8.\n"+
				"\t--maxscanlen <value>: Maximum length of the window to scan K-mer models. Default=14.\n"+
				"\t--hillsthresh <value>: Scoring threshold to identify hills. Default=0.1.\n"+
				"\t--mememinw <value>: minw arg for MEME. Default=6.\n"+
				"\t--mememaxw <value>: maxw arg for MEME. Default=13. This value should always be less than \"maxscanlen\".\n"+
				"\t--memenmotifs <int>: Number of motifs MEME should find in each condition (default=3)\n" +
				"\t--memeargs <args> : Additional args for MEME (default:  -dna -mod zoops -revcomp -nostatus)\n"+
				"\t--memesearchwin <value>: Window around hills to search for discriminative motifs. Default=16. (Only applicable when run with \"genregs\").\n"+
				"\t--motifminROC <value>: minimum class-specific ROC required to report motif. Default=0.7." +
				""));
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

	
}
