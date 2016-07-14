package org.seqcode.genome.sequence.seqfunctions;

import java.util.ArrayList;
import java.util.List;

import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.ScoredPoint;
import org.seqcode.genome.location.ScoredStrandedPoint;
import org.seqcode.genome.sequence.ScoredSequence;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.genome.sequence.SequenceUtils;

/**
 * Create a composite scoring function of a particular type of SeqFunction
 * 
 * @author mahony
 *
 * @param <F> SeqFunction
 */
public class CompositeSeqFunction<F extends SeqFunction> {

	F function;
	int width;
	double [][] means;
	double [][] variances;
	List<ScoredSequence> scoredSeqs;
	int zeroOffset=0;
	
	public CompositeSeqFunction(F fn, int width, int zeroOffset){
		function = fn;
		this.width = width;
		means = new double[function.scoreDimension()][width];
		variances = new double[function.scoreDimension()][width];
		for(int i=0; i<function.scoreDimension(); i++)
			for(int j=0; j<width; j++){
				means[i][j]=0;
				variances[i][j]=0;
			}
		this.scoredSeqs = new ArrayList<ScoredSequence>();
		this.zeroOffset = zeroOffset;
	}
	
	//Accessors
	public int getWidth(){return width;}
	public double[][] getMeans(){return means;}
	public double[][] getVariances(){return variances;}
	public int getNumSeqs(){return scoredSeqs.size();}
	public F getFunction(){return function;}
	public int getZeroOffset(){return zeroOffset;}
	
	/**
	 * Add sequences (possibly weighted) to the composite list
	 * @param points
	 */
	public void addSequencesFromPoints(List<Point> points, SequenceGenerator seqgen){
		for(Point p : points){
			double weight = 1;
			char strand = '+';
			if(p instanceof ScoredStrandedPoint){
				strand = p.getStrand();
				weight = ((ScoredPoint)p).getScore();
			}
			
			Region r = new Region(p.getGenome(), p.getChrom(), p.getLocation()-(width/2)+1, p.getLocation()+(width/2));
			String seq = seqgen.execute(r);
			//Handle minus strand regions
	        if(strand=='-')
	        	seq = SequenceUtils.reverseComplement(seq);
	        
	        //Ignore regions that are too short (presumably because they are at the edge of the chromosome)
	        if(seq.length()==width){
	        	scoredSeqs.add(new ScoredSequence(seq, weight));
	        }
		}
	}
	
	/**
	 * Add equally weighted sequences to the composite list
	 * @param seqs
	 */
	public void addSequences(List<String> seqs){
		for(String s : seqs){
			if(s.length()==width){
				scoredSeqs.add(new ScoredSequence(s, 1.0));
			}
		}
	}
	
	/**
	 * Add weighted sequences to the composite list
	 * @param seqs
	 */
	public void addScoredSequences(List<ScoredSequence> sseqs){
		scoredSeqs.addAll(sseqs);
	}
	
	/**
	 * Populate the composite profile with regions of a fixed width around points
	 * @param points
	 */
	public void populate(){
		try{
			double normTotal=0;

			//Weighted mean
			for(ScoredSequence s : scoredSeqs){
				if(s.getSeq().length()==width){
		        	normTotal+=s.getScore();
		        	
		        	double[][] tmpScores = function.score(s.getSeq());
					for(int i=0; i<function.scoreDimension(); i++)
		    			for(int j=0; j<width; j++)
		    				means[i][j]+=(tmpScores[i][j]*s.getScore());
				}
			}
			for(int i=0; i<function.scoreDimension(); i++)
				for(int j=0; j<width; j++)
					means[i][j] = means[i][j]/normTotal;
			
			//Weighted variance
			for(ScoredSequence s : scoredSeqs){
				if(s.getSeq().length()==width){
		        	double[][] tmpScores = function.score(s.getSeq());
					for(int i=0; i<function.scoreDimension(); i++)
		    			for(int j=0; j<width; j++)
		    				variances[i][j]+=s.getScore()*(tmpScores[i][j]-means[i][j])*(tmpScores[i][j]-means[i][j]);
				}
			}
			for(int i=0; i<function.scoreDimension(); i++)
				for(int j=0; j<width; j++)
					variances[i][j] = variances[i][j]/normTotal;
			
			
		} catch (SeqFunctionException e) {
			e.printStackTrace();
		}
	}
	
	
	public void printFunctionMeans(){
		String[] labels = function.dimensionLabels();
		for(int i=0; i<function.scoreDimension(); i++){
			System.out.print(labels[i]);
			for(int j=0; j<width; j++)
				System.out.print("\t"+means[i][j]);
			System.out.println("");
		}
	}
	
	public void printFunctionVariances(){
		String[] labels = function.dimensionLabels();
		for(int i=0; i<function.scoreDimension(); i++){
			System.out.print(labels[i]);
			for(int j=0; j<width; j++)
				System.out.print("\t"+variances[i][j]);
			System.out.println("");
		}
	}

}
