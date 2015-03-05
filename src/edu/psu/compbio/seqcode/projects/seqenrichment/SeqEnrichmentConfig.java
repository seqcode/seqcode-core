package edu.psu.compbio.seqcode.projects.seqenrichment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import edu.psu.compbio.seqcode.gse.utils.ArgParser;

public class SeqEnrichmentConfig {

	protected boolean printHelp=false;
	protected String[] args;
	public String getArgs(){
		String a="";
		for(int i=0; i<args.length; i++)
			a = a+" "+args[i];
		return a;
	}
	
	public SeqEnrichmentConfig(String[] args){
		this.args = args;
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
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	/**
	 * returns a string describing the arguments handled by this parser. 
	 * @return String
	 */
	public String getArgsList(){
		return(new String("" +
				"\t--config <config file: all options can be specified in a name<space>value text file, over-ridden by command-line args>\n"+
				""));
	}
				
}
