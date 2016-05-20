package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import javax.swing.JFrame;

import org.math.plot.FrameView;
import org.math.plot.Plot2DPanel;
import org.math.plot.plots.ColoredScatterPlot;

import com.jujutsu.tsne.FastTSne;
import com.jujutsu.tsne.MatrixOps;
import com.jujutsu.tsne.TSne;

import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;

public class TSNE_testing {

	private double[][] data=null;
	private int numCols=0;
	private int numRows=0;
	private List<String> colNames=new ArrayList<String>();
	private List<String> labels=new ArrayList<String>();
	private List<String> coords=new ArrayList<String>();
	private double perplexity=20;
	private int initialDims=10;
	private int iters=1000;
	
	public TSNE_testing(String datafile, double perplexity, int initial_dims, int iters){
		this.perplexity = perplexity;
		this.initialDims = initial_dims;
		this.iters = iters;
		loadDataFile(datafile);
	}
	
	/**
	 * Train tSNE
	 */
	public void execute(String outFilename){
		if(data!=null){
			TSne tsne = new FastTSne();
	    	System.out.println("Running " + iters + " iterations of t-SNE");
	        System.out.println("Shape is: " + data.length + " x " + data[0].length);
	        
	        System.out.println("Starting TSNE: " + new Date());
	        double [][] Y = tsne.tsne(data, 2, initialDims, perplexity, iters);
	        System.out.println("Finished TSNE: " + new Date());
	        System.out.println("Result is = " + Y.length + " x " + Y[0].length);
	        
	        
	        saveFile(new File(outFilename), MatrixOps.doubleArrayToString(Y));
	        Plot2DPanel plot = new Plot2DPanel();
	        
	        String[] labelArr = (String[]) labels.toArray();
	        ColoredScatterPlot dataPlot = new ColoredScatterPlot("data", Y, labelArr);
	        plot.plotCanvas.setNotable(true);
	        plot.plotCanvas.setNoteCoords(true);
	        plot.plotCanvas.addPlot(dataPlot);
	                
	        FrameView plotframe = new FrameView(plot);
	        plotframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	        plotframe.setVisible(true);
		}
	}
	
	private void saveFile(File file, String text) {
		saveFile(file,text,false);
	}
	
	private void saveFile(File file, String text, boolean append) {
        try{ 
        	FileWriter fw = new FileWriter(file, append);
        	BufferedWriter bw = new BufferedWriter(fw);
            bw.write(text);
            bw.close();
        } catch (IOException e) {
            System.err.println(e.toString());
        }
	}
	
	/**
	 * Load the data file
	 * @param infile
	 */
	private void loadDataFile(String infile){
		List<String> dataLines = new ArrayList<String>();
		String header;
		
		try{
			//Read the file
			File pFile = new File(infile);
			if(!pFile.isFile()){System.err.println("Invalid file name "+infile);System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        header = reader.readLine(); 
	        String line;
	        while ((line = reader.readLine()) != null) {
	            dataLines.add(line.trim());
	            numRows++;
	        }reader.close();
	        

			//Process out the labels and data
	        String[] hWords = header.split("\\t");
			numCols = hWords.length-2;
			data = new double[numRows][numCols];
			for(int c=2; c<hWords.length; c++)
				colNames.add(hWords[c]);
			int r=0;
			for(String s : dataLines){
				String[] sWords = s.split("\\t");
				labels.add(sWords[0]);
				coords.add(sWords[1]);
				for(int c=2; c<sWords.length; c++)
					data[r][c-2] = new Double(sWords[c]);
				r++;
			}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * main
	 * @param args
	 */
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
		if(args.length==0 || !ap.hasKey("data") || ap.hasKey("h")){
			System.err.println("TSNE_testing:");
			System.err.println("" +
					"\t--data <filename>\n" +
					"\t--out <filename>\n" +
					"\t--perplexity <value>\n" +
					"\t--initialdims <value>\n" +
					"\t--iters <value>\n"
					);
		}else{
			String dataName = Args.parseString(args, "data", null);
			String outName = Args.parseString(args, "out", "out");
			double perplexity = Args.parseDouble(args, "perplexity", 20);
			int dims = Args.parseInteger(args, "initialdims", 10);
			int iters = Args.parseInteger(args, "iters", 1000);
		
			TSNE_testing tester = new TSNE_testing(dataName, perplexity, dims, iters);
			tester.execute(outName);
		}
	}

}
