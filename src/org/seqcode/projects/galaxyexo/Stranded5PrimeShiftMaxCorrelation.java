package org.seqcode.projects.galaxyexo;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.Vector;

import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.gseutils.Args;
import org.seqcode.viz.metaprofile.BinningParameters;
import org.seqcode.viz.metaprofile.MetaConfig;
import org.seqcode.viz.metaprofile.MetaUtils;
import org.seqcode.viz.metaprofile.PointProfile;
import org.seqcode.viz.metaprofile.PointProfiler;
import org.seqcode.viz.metaprofile.Stranded5PrimeProfiler;

/**
 * Utility to compute a base pair shift that maximizes the correlation between stranded tags. 
 * 
 * Input:
 * 		- Genome
 * 		- Signal experiment
 * 		- Reference point
 * Output:
 * 		- Base pair shift that maximizes stranded tag correlation
 * 
 * @author naomi yamada
 */
public class Stranded5PrimeShiftMaxCorrelation {
	
	private GenomeConfig gconfig;
	private MetaConfig mconfig;
	private ExptConfig econfig;
	private ExperimentManager manager;
	private MetaUtils utils;
	private BinningParameters params; 
	private int minShift=6;
	private int maxShift=100;
	private boolean isLog=true;
	
	public Stranded5PrimeShiftMaxCorrelation(GenomeConfig g, MetaConfig m, ExptConfig e){
		gconfig = g;
		mconfig = m;
		econfig = e;
		manager = new ExperimentManager(econfig, true);
		utils = new MetaUtils(gconfig.getGenome());
		params = new BinningParameters(mconfig.winLen, mconfig.bins);
		System.out.println("Binding Parameters:\tWindow size: "+params.getWindowSize()+"\tBins: "+params.getNumBins());
	}
	
	public void setMinShift(int shift){minShift=shift;}
	public void setMaxShift(int shift){maxShift=shift;}
	public void setLog(boolean hasLog){ isLog=hasLog;}
	
	public void run() throws IOException{
		
//		File file_mat = new File(mconfig.outName+"_matrix.peaks");

//		FileWriter fw_mat = new FileWriter(file_mat.getAbsoluteFile());
//		BufferedWriter br_mat = new BufferedWriter(fw_mat);
		
		String peakFile = mconfig.peakFiles.get(0);
		Vector<Point> points = utils.loadPoints(new File(peakFile));
		
		double maxCorr=0;
		double bestShift=0;
		for (int i=minShift; i<maxShift; i++){
			double currCorr= correlationBpShift(points,i);
			if (currCorr > maxCorr){
				maxCorr = currCorr;
				bestShift =i;			
			}
		}	
		
		System.out.println("MaximumCorrelation : "+maxCorr+"\tMaximumShift : "+bestShift);
	}
	
	public double correlationBpShift(Vector<Point> points, int shift){
		
		PointProfiler posProfiler= new Stranded5PrimeProfiler(gconfig, params, manager, '+', shift, mconfig.baseLimit, mconfig.baseLimitRelPosition);
		PointProfiler negProfiler= new Stranded5PrimeProfiler(gconfig, params, manager, '-', shift, mconfig.baseLimit, mconfig.baseLimitRelPosition);

		double[][] pos_mat = null; double[][] neg_mat = null;
		double sum_pos = 0; double sum_neg = 0;
		double total_len = 0;
		for(int k=0; k<points.size(); k++){
			PointProfile posProfile = (PointProfile) posProfiler.execute(points.get(k));
			PointProfile negProfile = (PointProfile) negProfiler.execute(points.get(k));
			if(k==0){
				pos_mat = new double[points.size()][posProfile.length()];
				neg_mat = new double[points.size()][negProfile.length()];
				for(int j=0; j< posProfile.length(); j++){
					pos_mat[k][j] = posProfile.value(j);
					neg_mat[k][j] = negProfile.value(j);
					if (isLog){
						if (pos_mat[k][j]==0){ pos_mat[k][j]=1;} //pseudo counts
						if (neg_mat[k][j]==0){ neg_mat[k][j]=1;} 
						pos_mat[k][j] = Math.log(pos_mat[k][j]);
						neg_mat[k][j] = Math.log(neg_mat[k][j]);
					}
					sum_pos += pos_mat[k][j];
					sum_neg += neg_mat[k][j];
					total_len++;
				}
			}
			else{
				for(int j=0; j< posProfile.length(); j++){
					pos_mat[k][j] = posProfile.value(j);
					neg_mat[k][j] = negProfile.value(j);
					if (isLog){
						if (pos_mat[k][j]==0){ pos_mat[k][j]=1;} //pseudo counts
						if (neg_mat[k][j]==0){ neg_mat[k][j]=1;} 
						pos_mat[k][j] = Math.log(pos_mat[k][j]);
						neg_mat[k][j] = Math.log(neg_mat[k][j]);
					}
					sum_pos += pos_mat[k][j];
					sum_neg += neg_mat[k][j];
					total_len++;
				}
			}
		}
		double ave_pos = sum_pos/total_len;
		double ave_neg = sum_neg/total_len;
		double cov = 0;
		double var_pos = 0; double var_neg = 0;
		for (int k=0; k< points.size();k++){
			for (int j=0; j < pos_mat[k].length; j++){
				double posj = pos_mat[k][j] - ave_pos;
				double negj = neg_mat[k][j] - ave_neg;
				cov += posj*negj;
				var_pos += posj*posj;
				var_neg += negj*negj;
			}
		}
		return cov/(Math.sqrt(var_pos)*Math.sqrt(var_neg));	
	}
	
	public static void main(String[] args) throws IOException {
		GenomeConfig gconfig = new GenomeConfig(args);
		ExptConfig econfig = new ExptConfig(gconfig.getGenome(), args);
		MetaConfig mconfig = new MetaConfig(args);
		
		System.out.println("winwow size set? "+mconfig.winLen);
		
		Stranded5PrimeShiftMaxCorrelation maxcorr = new Stranded5PrimeShiftMaxCorrelation(gconfig, mconfig, econfig);
		boolean isLog= Args.parseFlags(args).contains("notlog") ? false: true;
		maxcorr.setLog(isLog);
		
		maxcorr.run();
	}

}
