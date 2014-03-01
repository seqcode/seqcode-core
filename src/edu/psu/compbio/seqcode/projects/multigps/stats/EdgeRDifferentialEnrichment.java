package edu.psu.compbio.seqcode.projects.multigps.stats;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import edu.psu.compbio.seqcode.projects.multigps.framework.Config;

/**
 * EdgeRDifferentialEnrichment: Construct a script and run EdgeR from the command-line.
 * Obviously, this won't work if R and the EdgeR library are not installed. 
 * Tested only on R v2.14 and above
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class EdgeRDifferentialEnrichment extends DifferentialEnrichment{

	protected Config config;
	protected CountsDataset data;
	
	//Constructor
	public EdgeRDifferentialEnrichment(Config c){
		super();
		config = c;
	}
	
	/**
	 * execute: make & run an appropriate EdgeR script 
	 */ 
	public CountsDataset execute(CountsDataset dat) {
		this.data = dat;
		String repCountsFilename = config.getOutputParentDir()+File.separator+config.getOutBase()+".replicates.counts";
		String scriptFilename = config.getOutputParentDir()+File.separator+"call_DE_GLM"+fileIDname+".R";
		
		//Write the R script
		try {
    		FileWriter fout = new FileWriter(scriptFilename);
    		fout.write(	"#!/usr/bin/Rscript --vanilla \n"+
    					"# Rscript call_DE_GLM.R repCounts.txt 0.15 \n" +
    					"# Dataset: "+data.getCondName(data.getFocalCondition())+"\n" +
    					"library(edgeR) \n" +
    					"args <- commandArgs(TRUE) \n"+
    					"raw <- read.delim(args[1], row.names=\"Point\") \n"+
    					"y <- DGEList(raw)"+
    					"\n");

    		String factorString = "c(";
    		int[] design = data.getDesignArray();
    		for(int d=0; d<design.length; d++){
    			factorString = factorString+"\""+design[d]+"\"";
    			if(d<design.length-1)
    				factorString = factorString+",";
    		}
    		factorString = factorString+")";
    		fout.write(	"group <- factor("+factorString+")\n");

    		fout.write(	"design <- model.matrix(~group) \n" +
    					"design \n"+
    					"y <- calcNormFactors(y, method=\"TMM\") # TMM method \n"+
    					"y <- estimateGLMCommonDisp(y, design, method=\"deviance\", robust=TRUE, subset=NULL) # this only works on later versions of edgeR, which you get with R >= 2.14 \n"+
    					"print(\"Best-fit overdispersion parameter:\") \n"+
    					"print(y$common.dispersion) \n"+
    					"y$common.dispersion = as.numeric(args[2]) \n"+
    					"print(\"In tests, using overdispersion parameter:\") \n"+
    					"print(y$common.dispersion) \n"+
    					"fit <- glmFit(y, design, dispersion=y$common.dispersion) \n"+
    					"\n");

    		int ref = data.getFocalCondition();
			for(int x=0; x<data.getNumConditions(); x++){
    			if(x!=ref){
    				fout.write(	"# "+data.getCondName(x)+" vs "+data.getCondName(ref)+"\n");
    				String contrastString = "contrast=c(0,";
    				for(int c=1; c<data.getNumConditions(); c++){ 
    					if(c==ref){contrastString = contrastString+"-1";}
    					else if(c==x){contrastString = contrastString+"1";}
    					else{contrastString = contrastString+"0";}
    					if(c<data.getNumConditions()-1){contrastString = contrastString+",";}
    					else{contrastString = contrastString+")";}
    				}
    				String diffFilename = config.getOutputParentDir()+File.separator+config.getOutBase()+".overdisp"+config.getEdgeROverDisp()+
    						"."+data.getCondName(x)+"vs"+data.getCondName(ref)+".edgeR_GLM_DE.txt";
    				fout.write(	"calls <- glmLRT(fit, "+contrastString+") # later versions of edgeR (R >= 2.15) dropped the y argument in the call to glmLRT \n"+
    							"all_tags = topTags(calls, n=1000000)$table \n"+
    							"write.table(all_tags, file=\""+diffFilename+"\", sep=\"\\t\") \n"+
    							"\n");
    			}
			}
    		fout.close();

    		//Run the R script
    		String Rscriptcmd = config.getRpath()+"Rscript ";
    		System.err.println("Running: "+Rscriptcmd+" "+scriptFilename+" "+repCountsFilename+" "+config.getEdgeROverDisp());
    		Process proc = Runtime.getRuntime().exec(Rscriptcmd+" "+scriptFilename+" "+repCountsFilename+" "+config.getEdgeROverDisp());
    		// any error message? 
			StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "R_ERR", true); 
			// any output? 
			StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "R_OUT", true); 
			// kick them off 
			errorGobbler.start(); 
			outputGobbler.start(); 
			int exitVal = proc.waitFor(); 
    		System.err.println("R ExitValue: " + exitVal);
    		proc.destroy();
    		
    		//Import the EdgeR results
    		for(int x=0; x<data.getNumConditions(); x++){
    			if(x!=ref){
    				String diffFilename = config.getOutputParentDir()+File.separator+config.getOutBase()+".overdisp"+config.getEdgeROverDisp()+
    						"."+data.getCondName(x)+"vs"+data.getCondName(ref)+".edgeR_GLM_DE.txt";
    				File dfFile = new File(diffFilename);
    				if(!dfFile.isFile()){System.err.println("ERROR: Differential enrichment file not generated: "+diffFilename);}
    		        BufferedReader reader = new BufferedReader(new FileReader(dfFile));
    		        //First line should have column labels
    		        String line= reader.readLine();
    		        while ((line = reader.readLine()) != null) {
    		            line = line.trim();
    		            String[] words = line.split("\\s+");
    		            //Edit for Double correctness
    		            if(words[1].equals("Inf")){words[1]="Infinite";} if(words[1].equals("-Inf")){words[1]="-Infinite";}
    		            if(words[2].equals("Inf")){words[2]="Infinite";} if(words[2].equals("-Inf")){words[2]="-Infinite";}
    		            if(words[5].equals("Inf")){words[5]="Infinite";} if(words[5].equals("-Inf")){words[5]="-Infinite";}
    		            //Format: 
    		            //Point logFC logCPM LR PValue FDR
    		            String pointStr = words[0].replaceAll("\"", "");
    		            Integer pointIndex = data.getUnitID(pointStr);
    		            Double logFC = new Double(words[1]);
    		            if(logFC>config.LOG_FC_LIMIT){ logFC=config.LOG_FC_LIMIT; }
    		            if(logFC<-config.LOG_FC_LIMIT){ logFC=-config.LOG_FC_LIMIT; }
    		            Double logCPM = new Double(words[2]);
    		            Double FDR = new Double(words[5]);
    		            
    		            int u = data.getUnitID(pointStr);
    		            data.setDEpval(u, x, FDR);
    		            data.setCondFold(u, x, logFC);
    		            data.setCondMean(u,x,logCPM);
    		        }reader.close();
    			}
    		}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		
		return data;
	}
	
}
