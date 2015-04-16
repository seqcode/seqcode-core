package edu.psu.compbio.seqcode.genome;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

/**
 * GenomeConfig:
 * A config parser that loads genome objects from the command-line or config files. 
 * You can also use the Args class directly to load Genomes from the command-line. 
 * However, GenomeConfig allows convenient loading of cached sequences as well, 
 * and fits the schema of the other config parser classes.   
 * 
 * @author mahony
 *
 */
public class GenomeConfig {
	private Genome gen=null;
	private String genomeSequencePath=null; //Path to sequence data file directories
	private SequenceGenerator<Region> seqgen=null;
	private boolean sequenceAvailable=false;
	private boolean printHelp=false;
	
	private String[] args;
	public String getArgs(){
		String a="";
		for(int i=0; i<args.length; i++)
			a = a+" "+args[i];
		return a;
	}
	
	public GenomeConfig(String [] arguments){
		this.args=arguments;
		ArgParser ap = new ArgParser(args);
		seqgen = new SequenceGenerator<Region>();
		
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
				
				//Load genome
				if(ap.hasKey("species") || ap.hasKey("genome") || ap.hasKey("gen")){
					Pair<Organism, Genome> pair = Args.parseGenome(args);
					if(pair != null){
						gen = pair.cdr();
						sequenceAvailable=true;
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
				
				//Cache genome sequence
				if(ap.hasKey("seq")){
					genomeSequencePath = ap.getKeyValue("seq");
					seqgen.setGenomePath(genomeSequencePath);
					seqgen.useCache(true);
					seqgen.useLocalFiles(true);
					sequenceAvailable=true;
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
	public SequenceGenerator getSequenceGenerator(){return seqgen;}
	public String getGenomeSequencePath(){return genomeSequencePath;}
	public boolean sequenceAvailable(){return sequenceAvailable;}
	public boolean helpWanted(){return printHelp;}
	
	/**
	 * Close database connections in Genome object 
	 */
	public void close(){
		if(!gen.isClosed())
			gen.close();
	}
	
	/**
	 * Returns a string describing the arguments handled by this config parser. 
	 * @return String
	 */
	public String getArgsList(){
		return(new String("" +
				"Genome:" +
				"\t--species <Organism;Genome>\n" +
				"\tOR\n" +
				"\t--geninfo <genome info file>" +
				"Genome Sequence Caching:" +
				"\t--seq <fasta seq directory>\n" +
				""));
	}
}
