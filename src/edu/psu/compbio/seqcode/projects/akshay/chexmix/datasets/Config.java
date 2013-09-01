package edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets;

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
				"\t--O <turn Off center approach; by default it is on>\n"));
				
	}
}
