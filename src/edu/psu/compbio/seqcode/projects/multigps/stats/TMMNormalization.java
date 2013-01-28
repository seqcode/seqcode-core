package edu.psu.compbio.seqcode.projects.multigps.stats;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import Jama.Matrix;

import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.shaun.viz.ScatterPlotConfigured;

/**
 * TMMNormalization: a Normalization class that implements the "trimmed mean of M values" (TMM) normalization
 * as proposed by Robinson & Oshlack, Genome Biology, 2010
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class TMMNormalization extends Normalization{

	int ref=0;
	double Mtrim = 0.3;
	double Atrim = 0.05;
	double M_min=-Double.MAX_VALUE, M_max=Double.MAX_VALUE;
	double A_min=-Double.MAX_VALUE, A_max=Double.MAX_VALUE;
	
	
	/**
	 * Constructor: provide the number of samples, the fraction of extreme M values to trim, and the fraction of lowest A values to trim
	 * @param numSamples
	 * @param Mtrim
	 * @param Atrim
	 */
	public TMMNormalization(int numSamples) {this(numSamples, 0.3, 0.05);}
	public TMMNormalization(int numSamples, double Mtrim, double Atrim) {
		super(numSamples);
		this.Mtrim = Mtrim;
		this.Atrim = Atrim;
	}

	@Override
	public double[] normalize(CountsDataset data) {
		
		//Set one sample as the reference (the deepest sequenced sample in the focal condition)
		double maxTotal=0;
		for(int s=0; s<data.numSamples; s++)
			if(data.design[s] == data.focalCondition)
				if(data.totals[s]>maxTotal){
					ref=s; maxTotal=data.totals[s];
				}
		
		//Scale all samples against the reference
		for(int x=0; x<data.numSamples; x++){
			if(x==ref){
				scalingFactors[x]=1;
			}else{
				//Make an MAval set from the current sample
				List<MAval> MA = new ArrayList<MAval>();
				for(int d=0; d<data.numUnits; d++)
					if(data.getCount(d,x)>0 && data.getCount(d,ref)>0)
						MA.add(new MAval(data.getCount(d,x), data.getCount(d,ref), data.getTotal(x), data.getTotal(ref)));
				
				//Trim the dataset
				MA = trimMA(MA);
				
				//Estimate scaling
				double sumW=0, sumWM=0; 
				for(MAval v : MA){
					sumW+=v.w;
					sumWM+=v.w*v.M;
				}
				double TMM = Math.pow(2.0, sumWM/sumW);
				depthScaling[x]= (data.totals[x] / data.totals[ref]);
				propScaling[x] = TMM;
				scalingFactors[x] = depthScaling[x]*propScaling[x];
				
				double logTMM = sumWM/sumW;
				//System.err.println("logTMM: "+logTMM+"\tTMM: "+TMM+"\tTotalRatio:"+(data.totals[x] / data.totals[ref])+"\tScaling:"+scalingFactors[x]);				
			}
		}
		data.setScalingFactors(scalingFactors);
		
		return scalingFactors;
	}

	//Over-riding to highlight the trimmed MA value set
	@Override 
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
				
				List<MAval> highlightMA = new ArrayList<MAval>();
				List<MAval> otherMA = new ArrayList<MAval>();
				for(MAval v : MA){
					if(v.M>M_min && v.M<M_max && v.A>A_min && v.A<A_max)
						highlightMA.add(v);
					else
						otherMA.add(v);
				}
					
				//Make the MA matrices
				Matrix maMatrixHighlight = new Matrix(highlightMA.size(),2);
				Matrix maMatrixOther = new Matrix(otherMA.size(),2);
				int count=0;
				for(MAval v : highlightMA){
					maMatrixHighlight.set(count, 0, v.A);
					maMatrixHighlight.set(count, 1, v.M);
					count++;
				}
				count=0;
				for(MAval v : otherMA){
					maMatrixOther.set(count, 0, v.A);
					maMatrixOther.set(count, 1, v.M);
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
				plotter.saveMAplot(maMatrixOther, maMatrixHighlight, Math.log(propScaling[x]/propScaling[ref])/LOG_2, fileName, rasterImage);
			}
		}
	}
	
	/**
	 * trimMA: trim the MA values list according to the Mtrim and Atrim fractions
	 * @param input
	 * @return
	 */
	protected List<MAval> trimMA(List<MAval> input){
		List<MAval> trimmed = new ArrayList<MAval>();
		
		//sort on M values
		if(Mtrim>0){
			Collections.sort(input, new Comparator<MAval>(){
				public int compare(MAval a, MAval b){
					return a.compareByM(b);
				}
			});
			M_min = input.get((int)(input.size()*(Mtrim/2))).M;
			M_max = input.get((int)(input.size()*(1-(Mtrim/2)))).M;
		}
		
		//sort on A values
		if(Atrim>0){
			Collections.sort(input, new Comparator<MAval>(){
				public int compare(MAval a, MAval b){
					return a.compareByA(b);
				}
			});
			A_min = input.get((int)(input.size()*(Atrim))).A;
		}
		
		for(MAval v : input)
			if(v.M>M_min && v.M<M_max && v.A>A_min && v.A<A_max)
				trimmed.add(v);
		return trimmed;
	}
	
	@Override
	public String selfDescriptor() {
		return "Trimmed mean of M values (TMM) normalization as proposed by Robinson & Oshlack, Genome Biology, 2010 (used in edgeR).";
	}
}
