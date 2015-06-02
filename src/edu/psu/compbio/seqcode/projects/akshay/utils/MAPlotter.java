package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.awt.BasicStroke;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;

import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.viz.scatter.ScatterPlot;

import Jama.Matrix;

public class MAPlotter {
	protected Matrix DEPval;
	protected Matrix CondFold;
	protected Matrix CondMean;
	protected HashMap<String, Integer> unitToIndex;
	protected HashMap<Integer, String> indextoUnit;
	
	public final double LOG_FC_LIMIT = 10;
	
	// Import edgeR results
	public void loadMatices(String diffFilename) throws IOException{
		File dfFile = new File(diffFilename);
		BufferedReader reader = new BufferedReader(new FileReader(dfFile));
		//First line should have column labels
		String line= reader.readLine();
		int u = 0;
		unitToIndex = new HashMap<String, Integer>();
		indextoUnit = new HashMap<Integer, String>();
		while ((line = reader.readLine()) != null) {
			line = line.trim();
			String[] words = line.split("\\s+");
			 String pointStr = words[0].replaceAll("\"", "");
			unitToIndex.put(pointStr, u);
			indextoUnit.put(u, pointStr);
			u++;
		}reader.close();
		
		DEPval = new Matrix(u,1);
		CondFold = new Matrix(u,1);
		CondMean = new Matrix(u,1);
		
		BufferedReader reader1 = new BufferedReader(new FileReader(dfFile));
		line= reader1.readLine();
		
		u=0;
		
        while ((line = reader1.readLine()) != null) {
            line = line.trim();
            String[] words = line.split("\\s+");
            //Edit for Double correctness
            if(words[1].equals("Inf")){words[1]="Infinite";} if(words[1].equals("-Inf")){words[1]="-Infinite";}
            if(words[2].equals("Inf")){words[2]="Infinite";} if(words[2].equals("-Inf")){words[2]="-Infinite";}
            if(words[5].equals("Inf")){words[5]="Infinite";} if(words[5].equals("-Inf")){words[5]="-Infinite";}
            //Format: 
            //Point logFC logCPM LR PValue FDR
           
            
            Double logFC = new Double(words[1]);
            if(logFC>LOG_FC_LIMIT){ logFC=LOG_FC_LIMIT; }
            if(logFC<-LOG_FC_LIMIT){ logFC=-LOG_FC_LIMIT; }
            Double logCPM = new Double(words[2]);
            Double FDR = new Double(words[5]);
            
            
            DEPval.set(u, 0, FDR);
            CondFold.set(u, 0, logFC);
            CondMean.set(u,0,logCPM);
            u++;
        }reader1.close();
        
	}
	
	public void plot(String DEevents, String ConstantEvents) throws IOException{
		
		List<Pair<Double,Double>> DEhighlight = new ArrayList<Pair<Double,Double>>();
		List<Pair<Double,Double>> ConstantHighlight = new ArrayList<Pair<Double,Double>>();
		List<Pair<Double,Double>> otherMA = new ArrayList<Pair<Double,Double>>();
		
		HashMap<String, Integer> inDE = new HashMap<String, Integer>();
		HashMap<String, Integer> isConst = new HashMap<String, Integer>();
		
		
		BufferedReader DEreader = new BufferedReader(new FileReader(DEevents));
		String line;
		while ((line = DEreader.readLine()) != null) {
			line  = line.trim();
			String[] words = line.split("\\s+");
			inDE.put(words[0], 0);
		}DEreader.close();
		
		BufferedReader Constreader = new BufferedReader(new FileReader(ConstantEvents));
		while ((line = Constreader.readLine()) != null) {
			line  = line.trim();
			String[] words = line.split("\\s+");
			isConst.put(words[0], 0);
		}Constreader.close();
		
		double A_min=1;
		
		for(int d=0; d<unitToIndex.keySet().size(); d++){
			double fold=CondFold.get(d, 0), avg = CondMean.get(d,0);
			if(avg<A_min)
				avg = A_min;
			
			if(inDE.containsKey(indextoUnit.get(d)))
				DEhighlight.add(new Pair<Double,Double>(fold,avg));
			else if(isConst.containsKey(indextoUnit.get(d)))
				ConstantHighlight.add(new Pair<Double,Double>(fold,avg));
			else
				otherMA.add(new Pair<Double,Double>(fold,avg));
		}
		
		// make matrices
		Matrix maMatrixDE = new Matrix(DEhighlight.size(),2);
		Matrix maMatrixConst = new Matrix(ConstantHighlight.size(),2);
		Matrix maMatrixOther = new Matrix(otherMA.size(),2);
		int count=0;
		for(Pair<Double,Double> v : DEhighlight){
			maMatrixDE.set(count, 0, v.cdr());
			maMatrixDE.set(count, 1, v.car());
			count++;
		}
		count=0;
		for(Pair<Double,Double> v : ConstantHighlight){
			maMatrixConst.set(count, 0, v.cdr());
			maMatrixConst.set(count, 1, v.car());
			count++;
		}
		count=0;
		for(Pair<Double,Double> v : otherMA){
			maMatrixOther.set(count, 0, v.cdr());
			maMatrixOther.set(count, 1, v.car());
			count++;
		}
		
		
		
		
		ThreeClassScatterPlotMaker plotter = new ThreeClassScatterPlotMaker("Main");
		plotter.saveMAplot(maMatrixOther, maMatrixDE, maMatrixConst,  0.0, "ScatterPlot", true);
		
		
		
		
	}
	
	
	public class ThreeClassScatterPlotMaker extends ScatterPlot{
		public ThreeClassScatterPlotMaker(String title) {
			super(title);
		}
		
		public void saveMAplot(Matrix datapoints, Matrix datapoints_highlight, Matrix datapoints_constant, Double yLine, String outFilename, boolean rasterImage){
			this.setWidth(800);
			this.setHeight(800);
			this.setXLogScale(false);
			this.setYLogScale(false);
			this.addDataset("other", datapoints, new Color(75,75,75,80), 3);
			if(datapoints_highlight!=null)
				this.addDataset("DE", datapoints_highlight, new Color(0,0,255,80), 3);
			if(datapoints_constant !=null)
				this.addDataset("Constant", datapoints_constant, new Color(255,80,0,0), 3);
			this.setXAxisLabel("A");
			this.setYAxisLabel("M");
			this.setXRangeFromData();
			this.setYRange(-10.5,10.5);

			//Set the tick units according to the range
			double xUpper = daxis.getRange().getUpperBound();
			double xLower = daxis.getRange().getLowerBound();
	    	//if(daxis instanceof org.jfree.chart.axis.NumberAxis)
	    	//	((NumberAxis)daxis).setTickUnit(new NumberTickUnit(5));
	    	double yUpper = raxis.getRange().getUpperBound();
			double yLower = raxis.getRange().getLowerBound();
	    	if(raxis instanceof org.jfree.chart.axis.NumberAxis)
	    		((NumberAxis)raxis).setTickUnit(new NumberTickUnit(3));
			//Draw a line along y = yLine
	    	if(yLine!=null){
	    		XYLineAnnotation lineAnnot = new XYLineAnnotation(xLower, yLine, xUpper, yLine, new BasicStroke(1), new Color(0,0,0));
	    		this.plot.addAnnotation(lineAnnot);
	    	}
	    	
	    	try {
				this.saveImage(new File(outFilename), width, height, rasterImage);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
	}
	
	
	public static void main(String[] args) throws IOException{
		
		ArgParser ap = new ArgParser(args);
		
		String diffFilename = ap.getKeyValue("EdgeRmat");
		String DEevents = ap.getKeyValue("DEevents");
		String ConstantEvents = ap.getKeyValue("ConstantEvents");
		
		MAPlotter ma = new MAPlotter();
		
		ma.loadMatices(diffFilename);
		ma.plot(DEevents, ConstantEvents);
	}
	
	
	

}
