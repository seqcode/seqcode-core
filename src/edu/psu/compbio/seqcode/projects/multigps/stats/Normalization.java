package edu.psu.compbio.seqcode.projects.multigps.stats;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;

import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.shaun.viz.ScatterPlotConfigured;


/**
 * Normalization: parent class for normalization methods that calculate scaling factors. 
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public abstract class Normalization {

	protected double [] scalingFactors; //Overall scaling factors
	protected double [] depthScaling;   //Seq depth scaling factor
	protected double [] propScaling;	//RNA composition proportion scaling
	protected final double LOG_2 = Math.log(2.0);
	
	public Normalization(int numSamples){
		scalingFactors = new double[numSamples];
		depthScaling = new double[numSamples];
		propScaling = new double[numSamples];
		for(int c=0; c<numSamples; c++){
			scalingFactors[c]=1.0;
			depthScaling[c]=1.0;
			propScaling[c] =1.0;
		}
	}
	
	/**
	 * normalize: method that must be implemented by Normalization classes.
	 * Input a CountsDataset object, output the scaling factors array. 
	 * @param data
	 * @return
	 */
	public abstract double[] normalize(CountsDataset data);
	
	/**
	 * Return a string describing the normalization method. 
	 * @return
	 */
	public abstract String selfDescriptor();
	
	/**
	 * Get scaling factors array
	 * @return
	 */
	public double[] getScalingFactors(){return scalingFactors;}
	
	/**
	 * Get scaling factor for experiment x
	 * @param x
	 * @return
	 */
	public double getScalingFactor(int x){return scalingFactors[x];}
	
	/**
	 * Print pairwise MA values to files
	 * @param data
	 */
	public void printPairwiseMAData(CountsDataset data){
		int ref = 0;
		//Set one sample as the reference (the deepest sequenced sample in the focal condition)
		double maxTotal=0;
		for(int s=0; s<data.numSamples; s++)
			if(data.design[s] == data.focalCondition)
				if(data.totals[s]>maxTotal){
					ref=s; maxTotal=data.totals[s];
				}
		
		//Scale all samples against the reference
		for(int x=0; x<data.numSamples; x++){
			if(x!=ref){
				//Make an MAval set from the current sample
				List<MAval> MA = new ArrayList<MAval>();
				for(int d=0; d<data.numUnits; d++)
					if(data.getCount(d,x)>0 && data.getCount(d,ref)>0)
						MA.add(new MAval(data.getCount(d,x), data.getCount(d,ref), data.getTotal(x), data.getTotal(ref)));
				
				//Print to file
				Pair<String,String> refName = data.getExptName(ref);
				Pair<String,String> currName = data.getExptName(x);
				String fileName = refName.car()+"-"+refName.cdr()+"_vs_"+currName.car()+"-"+currName.cdr()+".ma.txt";
				try {
					FileWriter fw = new FileWriter(fileName);
					fw.write("Point\tA\tM\n");
					int count=0;
					for(MAval v : MA){
						fw.write(count+"\t"+v.A+"\t"+v.M+"\n"); count++;
					}
					fw.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	/**
	 * Generate MA images. May be over-ridden by subclasses.
	 * @param data
	 */
	public void savePairwiseMAPlots(CountsDataset data, boolean rasterImage){
		int ref = 0;
		//Set one sample as the reference (the deepest sequenced sample in the focal condition)
		double maxTotal=0;
		for(int s=0; s<data.numSamples; s++)
			if(data.design[s] == data.focalCondition)
				if(data.totals[s]>maxTotal){
					ref=s; maxTotal=data.totals[s];
				}
		
		//Scale all samples against the reference
		for(int x=0; x<data.numSamples; x++){
			if(x!=ref){
				//Make an MAval set from the current sample
				List<MAval> MA = new ArrayList<MAval>();
				for(int d=0; d<data.numUnits; d++)
					if(data.getCount(d,x)>0 && data.getCount(d,ref)>0)
						MA.add(new MAval(data.getCount(d,x), data.getCount(d,ref), data.getTotal(x), data.getTotal(ref)));
				
				//Make the MA matrix
				Matrix maMatrix = new Matrix(MA.size(),2);
				int count=0;
				for(MAval v : MA){
					maMatrix.set(count, 0, v.A);
					maMatrix.set(count, 1, v.M);
					count++;
				}
				
				//Image name
				Pair<String,String> refName = data.getExptName(ref);
				Pair<String,String> currName = data.getExptName(x);
				String fileName = refName.car()+"-"+refName.cdr()+"_vs_"+currName.car()+"-"+currName.cdr()+".ma";
				if(rasterImage)
					fileName = fileName+".png";
				else
					fileName = fileName+".svg";
				
				//Generate image
				ScatterPlotConfigured plotter = new ScatterPlotConfigured("MA plot");
				plotter.saveMAplot(maMatrix, null, Math.log(propScaling[x]/propScaling[ref])/LOG_2, fileName, rasterImage);
			}
		}
	}
}
