package edu.psu.compbio.seqcode.projects.shaun.rnaseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.deepseq.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.ewok.verbs.PointParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RegionParser;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.shaun.rnaseq.genemodels.GeneTUnit;

public abstract class RNASeqAnalyzer {

	protected boolean headers = true;
	protected Genome gen=null;
	protected boolean dbconnected = false;
	protected double genomeLen = 0;
	protected ArrayList<DeepSeqExpt> experiments = new ArrayList<DeepSeqExpt>();
	private ArrayList<String> conditionNames = new ArrayList<String>();
	protected String outputName="out.txt";
	protected GTFReader gtfReader;
	protected List<GeneTUnit> knownGenes = new ArrayList<GeneTUnit>();
	private ArrayList<Region> regionsOfInterest = new ArrayList<Region>();
	
	public RNASeqAnalyzer(String[] args){
		try{
			if(args.length==0){
				printError();System.exit(1);
			}
			ArgParser ap = new ArgParser(args);
			
			if (Args.parseFlags(args).contains("noheaders")) {
				headers = false;
			}
			
			//Load genome
			if(ap.hasKey("species")){
				Pair<Organism, Genome> pair = Args.parseGenome(args);
				if(pair != null){
					gen = pair.cdr();
					dbconnected=true;
					genomeLen = gen.getGenomeLength();
				}
			}else{
				//Make fake genome... chr lengths provided???
				if(ap.hasKey("geninfo") || ap.hasKey("g")){
					String fName = ap.hasKey("geninfo") ? ap.getKeyValue("geninfo") : ap.getKeyValue("g");
					gen = new Genome("Genome", new File(fName), true);
					genomeLen = gen.getGenomeLength();
				}else{
				    gen = null;
				}
			}
			
			//Load GTF file of transcripts
			if(ap.hasKey("gtf")){
				String gtfFile = ap.getKeyValue("gtf");
				gtfReader = new GTFReader(new File(gtfFile), gen);
				knownGenes = gtfReader.loadGenes();
				//Could define an initial fake genome here for the case where no genome was defined
			}
			
			//Load file of coords
			if(ap.hasKey("regions")){
				String regFile = ap.getKeyValue("regions");
				setRegionsOfInterest(loadRegionsFromFile(regFile, -1));
			}			
			
			//Output file
			if(ap.hasKey("out")){
				outputName = ap.getKeyValue("out");
			}
			
			//Load experiments
			Vector<String> exptTags=new Vector<String>();
			for(String s : ap.getKeys())
	        	if(s.contains("expt"))
	        		if(!exptTags.contains(s))
	        			exptTags.add(s);
			
			if(exptTags.size()==0){
			    System.err.println("Error: No experiments provided.\nUse the --expt option.");
			    printError();
			    System.exit(1);
			}
	        // each tag represents a condition
	        for(String tag : exptTags){
	        	String name="";
	        	if(tag.startsWith("rdb")){
	        		name = tag.replaceFirst("rdbexpt", ""); 
	        		getConditionNames().add(name);
	        	}else{
	        		name = tag.replaceFirst("expt", ""); 
	        		getConditionNames().add(name);
	        	}

	        	if(name.length()>0)
	        		System.out.println("    loading condition: "+name);
	        	
	        	List<ChipSeqLocator> rdbexpts = Args.parseChipSeq(args,"rdbexpt"+name);
	        	List<File> expts = Args.parseFileHandles(args, "expt"+name);
	        	boolean nonUnique = Args.parseFlags(args).contains("nonunique") ? true : false;
	        	String fileFormat = Args.parseString(args, "f", "TOPSAM").toUpperCase();

	        	if(expts.size()>0 && rdbexpts.size()==0){
	        		experiments.add(new DeepSeqExpt(gen, expts, nonUnique, fileFormat, -1));		        	
		        }
		        else if(rdbexpts.size()>0 && expts.size() == 0){
		        	if(gen==null){
	        			System.err.println("Error: the genome must be defined in order to use ReadDB."); 
	        			System.exit(1);
	        		}
		        	experiments.add(new DeepSeqExpt(gen, rdbexpts, "readdb", -1));
		        }
		        else{
		        	System.err.println("\n\nMust provide either an aligner output file or Gifford lab DB experiment name for the signal experiment (but not both)");
		        	printError();
		        	System.exit(1);
		        }
	        }
	        //We will need to combine all fake genomes here for the case where the genome hasn't been defined. 
		} catch (NotFoundException e) {
			e.printStackTrace();
			for(DeepSeqExpt x : experiments)
				x.closeLoaders();
		}
	}
	
	//Load a set of regions from a peak file
	protected ArrayList<Region> loadRegionsFromFile(String filename, int win){
		ArrayList<Region> regs = new ArrayList<Region>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line = "";
	        if (headers) {
	        	line = reader.readLine(); //Ignore first line
	        }
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            
	            if(words.length>=3 && win!=-1){
	                PointParser pparser = new PointParser(gen);
	            	Point p = pparser.execute(words[2]);
	            	int rstart = p.getLocation()-(win/2)<1 ? 1:p.getLocation()-(win/2);
                	int rend = p.getLocation()+(win/2)>gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+(win/2)-1;
                	Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
                	regs.add(r);
                }else if(words.length>=1){
	            	RegionParser parser = new RegionParser(gen);
	            	Region q = parser.execute(words[0]);
	            	if(win!=-1){
	            		int rstart = q.getMidpoint().getLocation()-(win/2)<1 ? 1:q.getMidpoint().getLocation()-(win/2);
	                	int rend = q.getMidpoint().getLocation()+(win/2)>gen.getChromLength(q.getChrom()) ? gen.getChromLength(q.getChrom()):q.getMidpoint().getLocation()+(win/2)-1;
	                	Region r = new Region(q.getGenome(), q.getChrom(), rstart, rend);
	                	if(r!=null){regs.add(r);}
	            	}else{
	            		if(q!=null){regs.add(q);}
	            	}
	            }
	        }reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(regs);
	}
	
	public void cleanup(){
		for(DeepSeqExpt x : experiments)
			x.closeLoaders();
	}
	
	public abstract void printError();
	public void printError(String programName){
		System.err.println(programName +
				"\nUsing with Gifford Lab ReadDB:" +
                "\n\t--rdbexpt <solexa expt>" +
                "\n\t--rdbctrl <background expt>" +
                "\nUsing with flat-files:\n" +
                "\n\t--expt <aligned reads file for expt>" +
                "\n\t--format <TOPSAM/ELAND/NOVO/BOWTIE/BED/SAM (default TOPSAM)>" +
                "\n\t--nonunique [use nonunique reads]" +
                "\nRequired:"+
                "\n\t--species <organism name;genome version>\n\tOR"+
                "\n\t--geninfo <file with chr name/length pairs>" +
                "\nOptions:" +
                "\n\t--gtf <GTF file of known transcripts>" +
                "\n\t--out <output filename>" +
				"\n");
	}

	public void setConditionNames(ArrayList<String> conditionNames) {
		this.conditionNames = conditionNames;
	}

	public ArrayList<String> getConditionNames() {
		return conditionNames;
	}

	public void setRegionsOfInterest(ArrayList<Region> regionsOfInterest) {
		this.regionsOfInterest = regionsOfInterest;
	}

	public ArrayList<Region> getRegionsOfInterest() {
		return regionsOfInterest;
	}
}
