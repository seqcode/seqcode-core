package edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import edu.psu.compbio.seqcode.gse.utils.probability.NormalDistribution;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.GenomicLocations;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;

/**
 * Class that performs  MAP assignment using the parameters of the trained Bayeseian network using the EM framework
 * @author akshaykakumanu
 *
 */

public class MAPassignment {
	
	//All the parameters of the Bayesian network
	public double[] PIj;
	public double[][] MUc;
	public double[][] MUf;
	public double[][] SIGMAc;
	public double[][] SIGMAf;
	public double[][] Bjk;
	public GenomicLocations trainingdata;
	public EMtrain model;
	public Config conf;
	
	// A 2-d array that stores the map-assignment values. Rows as training example and colums as a list of size "2".
	//The first element being he chromatin assignment and the second one being the factor assignment
	public double[][] MapAssignment;
	
	/**
	 * Constructor method
	 * @param model
	 * @param conf
	 */
	public MAPassignment(EMtrain model, Config conf) {
		this.model = model;
		this.conf = conf;
		PIj = model.getPIj();
		MUc = model.getMUc();
		MUf = model.getMUf();
		SIGMAc = model.getSIGMAc();
		SIGMAf = model.getSIGMAf();
		Bjk = model.getBjk();
		trainingdata = model.getTrainingData();
		MapAssignment = new double[trainingdata.getNumTrainingExamples()][2];
	}
	
	/**
	 * Performs the MAP assignment
	 */
	public void execute(){
		int N = trainingdata.getNumTrainingExamples();
		int numChromStates = conf.getNumChrmStates();
		int numFacState = conf.getNumFacStates();
		int C = trainingdata.getNumChromatinCons();
		int F = trainingdata.getNumFacCons();
		float[][] Xc = trainingdata.getChromatinCounts();
		float[][] Xf = trainingdata.getFactorCounts();
		
		for(int i=0; i<N; i++){ // over all training examples
			double[] assignment = new double[2];
			double maxLiklehood = 0.0;
			for(int j=0; j<numChromStates; j++){
				for(int k=0; k<numFacState; k++){
					double liklehood = 0.0;
					double chromGausssianProd=0.0;
					double facGaussianProd = 0.0;
					for(int c=0; c<C; c++){
						NormalDistribution gaussian = new NormalDistribution(MUc[j][c],Math.pow(SIGMAc[j][c], 2.0));
						chromGausssianProd = (c==0 ? gaussian.calcProbability((double) Xc[i][c]): chromGausssianProd* gaussian.calcProbability((double) Xc[i][c]));
					}
					for(int f=0; f< F; f++){
						NormalDistribution gaussian = new NormalDistribution(MUf[k][f],Math.pow(SIGMAf[k][f], 2.0));
						facGaussianProd = (f == 0 ? gaussian.calcProbability((double) Xf[i][f]): facGaussianProd* gaussian.calcProbability((double) Xf[i][f]));
					}
					// If the guassian prodocuts are NaN, make them 0.0
					chromGausssianProd = ( Double.isNaN(chromGausssianProd)) ? 0.0 : chromGausssianProd;
					facGaussianProd = (Double.isNaN(facGaussianProd)) ? 0.0: facGaussianProd;
					
					liklehood = PIj[j]*chromGausssianProd*Bjk[j][k]*facGaussianProd;
					liklehood = Double.isNaN(liklehood) ? 0 : liklehood;
					if(liklehood > maxLiklehood){
						maxLiklehood = liklehood;
						assignment[0] = j; //Chromatin assignment
						assignment[1] = k; //Factor assignment
					}
				}
			}
			this.MapAssignment[i][0] = assignment[0];
			this.MapAssignment[i][1] = assignment[1];
		}
		
		// Printing the locations that were assigned to each configuration
		for(int j=0; j<numChromStates; j++){ // over all chromatin states
			for(int k=0; k<numFacState; k++){ // over all factor states
				File outfile = new File(conf.getOutputInterDir().getAbsolutePath()+File.separator+Integer.toString(j)+"chromatin_"+
			Integer.toString(k)+"_factor.tab");
				try {
					FileWriter fw = new FileWriter(outfile);
					for(int i=0; i<N; i++){
						if((int) MapAssignment[i][0] == j && (int) MapAssignment[i][1] == k){
							fw.write(trainingdata.getLocations().get(i).getLocationString()+"\n");
						}
					}
					fw.close();
				}catch(IOException e){
					e.printStackTrace();
				}
			}
		}
	}
	

}
