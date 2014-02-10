package edu.psu.compbio.seqcode.projects.akshay.bayesments.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.stats.StreamGobbler;

/**
 * Plotter class that should be used to plot cumulative plots of count data
 * @author akshaykakumanu
 *
 */
public class CoutPlotter {
	// List of counts to be plotted
	public float[] counts;
	public Config config;
	//name tag to name the output file
	public String name_tag;
	
	/**
	 * Constructor class
	 * @param counts
	 * @param config
	 * @param name
	 */
	public CoutPlotter(float[] counts, Config config, String name) {
		this.counts = counts;
		this.config =config;
		this.name_tag = name;
	}
	
	/**
	 * Plots the data
	 */
	public void plot(){
		// Writing the data to the intermediate file
		try{
			String fdata = config.getOutputInterDir().getAbsolutePath()+File.separator
					+name_tag+"_data.tab";
			String imageName = config.getOutputImagesDir().getAbsolutePath()+File.separator
					+name_tag+"_boxplot.png";
			FileWriter fout = new FileWriter(fdata);
			for(int i=0; i<counts.length; i++){
				fout.write(Double.toString(counts[i])+"\n");
			}
			fout.close();
			
			//Writing the R script
			
			String rscript = config.getOutputInterDir().getAbsolutePath()+File.separator+
					"rscript.R";
			FileWriter rout = new FileWriter(rscript);
			rout.write("dat = read.table("+"\""+fdata+"\""+")\n"+
						"dat_sorted = sort(dat[,1])\n"+
						"maxx = max(dat_sorted)\n"+
						"breaks = seq(0,maxx,10)\n"+
						"cuts = cut(dat_sorted,breaks,right=FALSE)\n"+
						"freqs=table(cuts)\n"+
						"cumfreq0 = c(0, cumsum(freqs))\n"+
						"png("+"\""+imageName+"\""+")\n"+
						"plot(breaks,cumfreq0,xlab=\"number of tags\", ylab=\"number of events\")");
			
			rout.close();
			
			//Running the rscript
			
			Process proc = Runtime.getRuntime().exec("Rscript "+rscript);
			StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "R_ERR", true); 
			// any output? 
			StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "R_OUT", true); 
			// kick them off 
			errorGobbler.start(); 
			outputGobbler.start(); 
			int exitVal = proc.waitFor(); 
    		System.err.println("R ExitValue: " + exitVal);
    		proc.destroy();
		}
		catch(IOException e){
			e.printStackTrace();
		}catch (InterruptedException e) {
			e.printStackTrace();
		}
		
		
		
		
	}
	

}
