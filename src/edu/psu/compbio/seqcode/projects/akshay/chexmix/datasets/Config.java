package edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets;

import java.io.File;

import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;

public class Config {
	
	protected int blrange; // default value set to range+-10
	protected int range; //default value set to 30
	protected double pcc_seed_cutoff; // default value set to 0.8
	protected double pcc_cutoff; // default value is set to 0.6
	protected String tagsfile;
	protected String tagsfiletype; // default value is set to IDX
	protected String genome_name;
	protected String scheme_name;
	protected String genome_path;
	protected String peak_location; 
	protected int smoothing; //default value is set to 0
	protected int no_top_bl; // default value is set to 50
	protected boolean printHelp=false;
	protected boolean useCenterApproach=true;
	protected final boolean useNonUnique = false;
	protected String out_name; //default value is out
	protected String out_base;
	protected File outDir = null;
	protected int no_of_cycles; //default value is 3
	protected final int depth_to_search_patters=6;
	protected final int factor_to_add_one_by_one=4;
	protected final int factor_to_get_new_topbls = 10;
	protected int percentage_list_to_consider;
	protected final int factor_to_refine_profiles = 5;
	protected String MEMEPath = "";
	protected String MEMEargs=" -dna -mod zoops -revcomp -nostatus ";
	protected int MEMEminw=6;
	protected int MEMEmaxw=18;
	protected int MEMEnmotifs;
	
	
	protected String[] args;
	
	public Config(String[] arguments) {
		this.args=arguments;
		ArgParser ap = new ArgParser(args);
		if(args.length ==0 || ap.hasKey("h")){
			this.printHelp = true;
		}
		else{
			range = Args.parseInteger(args, "IR", 30);
			int offset = Args.parseInteger(args, "offset", 10);
			no_top_bl = Args.parseInteger(args, "TBL", 50);
			blrange = range+offset;
			pcc_seed_cutoff = Args.parseDouble(args, "PSC", 0.8);
			pcc_cutoff= Args.parseDouble(args, "PC", 0.6);
			smoothing =  Args.parseInteger(args, "S", 0);
			tagsfiletype = Args.parseString(args, "tagstype", "IDX").toUpperCase();
			out_name = Args.parseString(args, "Out", "out");
			outDir = new File(out_name);
			out_base = outDir.getName();
			
			no_of_cycles = Args.parseInteger(args, "NI", 3);
			percentage_list_to_consider = Args.parseInteger(args, "PL", 60);
			this.MEMEPath = Args.parseString(args, "memepath", MEMEPath);
			MEMEargs = Args.parseString(args, "memeargs", MEMEargs);
			//MEME minw
			MEMEminw = Args.parseInteger(args, "mememinw", MEMEminw);
			//MEME maxw
			MEMEmaxw = Args.parseInteger(args, "mememaxw", MEMEmaxw);
			//MEME nmotifs option
			MEMEnmotifs = Args.parseInteger(args,"memenmotifs", 3);
			MEMEargs = MEMEargs + " -nmotifs "+MEMEnmotifs + " -minw "+MEMEminw+" -maxw "+MEMEmaxw;
			
			for(String s: ap.getKeys()){
				if(s.equals("tags")){
					this.tagsfile = ap.getKeyValue(s);
				}
				if(s.equals("G")){
					this.genome_name = ap.getKeyValue(s);
				}
				if(s.equals("BM")){
					this.scheme_name = ap.getKeyValue(s);
				}
				if(s.equals("seq")){
					this.genome_path = ap.getKeyValue(s);
				}
				if(s.equals("P")){
					this.peak_location = ap.getKeyValue(s);
				}
			}
			
			for(String flag: Args.parseFlags(args)){
				if(flag.equals("O")){
					this.useCenterApproach = false;
				}
			}
			
		}
		
		
	}
	
	public void makeChexmixOuputDirs(){
		if(outDir.exists()){
			deleteDirectory(outDir);
		}
		outDir.mkdirs();
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
	
	//Accessories
	public boolean helpWanted(){return printHelp;}
	public String getTagsPath(){return this.tagsfile;}
	public String getTagsFormat(){return this.tagsfiletype;}
	public int getSmoothSize(){return this.smoothing;}
	public int getIntSize(){return this.range;}
	public boolean useCenter(){return useCenterApproach;}
	public double getSeedCutoff(){return this.pcc_seed_cutoff;}
	public int getBlsize(){return this.blrange;}
	public double getPccCutoff(){ return this.pcc_cutoff;}
	public boolean useNonUnique(){return this.useNonUnique;}
	public String getPeaksFilePath(){return this.peak_location;}
	public int getNoTopBls(){return this.no_top_bl;}
	public String getSchemename(){return this.scheme_name;}
	public String getOutName(){return this.out_name;}
	public String getOutBase(){return this.out_base;}
	public File getOutParentDir(){return this.outDir;}
	public int getNoOfCycles(){return this.no_of_cycles;}
	public int getHowDeepToSearch(){return this.depth_to_search_patters;}
	public int getFactorToAddIteratively(){return this.factor_to_add_one_by_one;}
	public int getFactorToGetNewTopBls(){return this.factor_to_get_new_topbls;}
	public int getListPercentageToCosider(){return this.percentage_list_to_consider;}
	public int getFactorToRefineSeedProfiles(){return this.factor_to_refine_profiles;}
	public String getGenomeName(){return this.genome_name;}
	public String getMEMEPath(){return this.MEMEPath;}
	public String getMEMEargs(){return this.MEMEargs;}
	
	
	public String getArgsList(){
		return(new String("" +
				"\t--tags <tags file>\n" +
				"\t--tagstype <tags file type; default value set to IDX>\n" +
				"\t--seq <Path to genome fasta files>\n" +
				"\t--P <Peak locations in bed format>\n" +
				"\t--G <name of the genome; hg18, hg19..>\n" +
				"\t--BM <Build method; scheme1, scheme2..>\n" +
				"\t--TBL <Top binding locations to consider for building seed>\n" +
				"\t--PSC <PCC seed cutoff; default value 0.8>\n"+
				"\t--PC <PCC cutoff; default value 0.6>\n" +
				"\t--offset <displacement to consider, +- this value will be considered; default value set to 10>\n"+
				"\t--IR <Interval Range; *2 this values will be the interval size; default value set to 30>\n" +
				"\t--S <Smoothing value; 0 for no smoothing; default value set to 0>\n"+
				"\t--O <turn Off center approach; by default it is on>\n" +
				"\t--Out <output tag name; files will be written to in the working directory; default value is out>\n"+
				"\t--NI <No of iterations of the whole mehtod; default value is set to 3>\n"+
				"\t--memepath <path to the meme bin dir (default: meme is in $PATH)>\n"+
				"\t--mememinw <minw arg for MEME (default="+MEMEminw+")>\n"+
				"\t--mememaxw <maxw arg for MEME (default="+MEMEmaxw+")>\n"+
				"\t--memenmotifs <number of motifs MEME should find for each condition>\n" +
				"\t--PL <Percentage of the list to consider>"));
				
	}
}
