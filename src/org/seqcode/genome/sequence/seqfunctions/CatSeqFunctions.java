package org.seqcode.genome.sequence.seqfunctions;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.ScoredPoint;
import org.seqcode.genome.location.ScoredStrandedPoint;
import org.seqcode.genome.sequence.ScoredSequence;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.genome.sequence.SequenceUtils;

/**
 * Concatenates the results of a set of SeqFunctions for each sequence (fixed-width) in a collection
 *  
 * @author mahony
 *
 */
public class CatSeqFunctions {

	List<SeqFunction> seqFuncs;
	int seqWidth; //Fixed width of all sequences that will be input
	int catLength; //Length of the concatenated vector created for each input sequence
	String [] positionLabels; //labels associated with each entry in the concatenated vector
	int positionLabelOffset; //Offset used for labels only (e.g define some position in each sequence to be the "0")
	List<ScoredSequence> scoredSeqs = new ArrayList<ScoredSequence>();
	Map<ScoredSequence, Double[]> seqResults = new HashMap<ScoredSequence, Double[]>();

	public CatSeqFunctions(List<SeqFunction> fns, int seqWidth){this(fns, seqWidth, 0);}
	public CatSeqFunctions(List<SeqFunction> fns, int seqWidth, int positionLabelOffset){
		this.seqWidth = seqWidth;
		seqFuncs = fns;
		
		//Count the length of the full concatenated vector
		catLength = 0;
		for(SeqFunction f : seqFuncs){
			catLength+=f.scoreDimension()*seqWidth;
		}
		
		//Initialize the positionLabels
		positionLabels = new String[catLength];
		int i=0;
		for(SeqFunction f : seqFuncs){
			String posType = f.isBetweenNucleotides() ? "step" : "pos";
			for(int p=0; p<seqWidth; p++){
				int position = p-positionLabelOffset;
				for(int d=0; d<f.scoreDimension(); d++){
					
					positionLabels[i] = f.dimensionLabels()[d]+"_at_"+posType+"_"+position;
					
					i++;
				}
			}
		}
	}
	
	/**
	 * Print labels to a file
	 * @param outFile
	 */
	public void printLabels(Path outFile){
		List<String> labels = new ArrayList<String>();
		for(int i=0; i<catLength; i++){
			labels.add((i+1)+":"+positionLabels[i]);
		}
		//Output file
		try {
			Files.write(outFile, labels, Charset.forName("UTF-8"));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Print data to a file in libSVM format
	 * @param outFile
	 */
	public void printLibSVMData(Path outFile){
		List<String> data = new ArrayList<String>();
		for(ScoredSequence s : seqResults.keySet()){
			Double[] scores = seqResults.get(s);
			String currLine = new Double(s.getScore()).toString(); 
			for(int i=0; i<catLength; i++){
				if(scores[i]!=0)
					currLine=currLine+"\t"+(i+1)+":"+String.format("%.3f", scores[i]);
			}
			data.add(currLine);
		}
		//Output file
		try {
			Files.write(outFile, data, Charset.forName("UTF-8"));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Print sequences to a file
	 * @param outFile
	 */
	public void printSequences(Path outFile){
		List<String> seqs = new ArrayList<String>();
		for(ScoredSequence s : scoredSeqs)
			seqs.add(s.getScore()+"\t"+s.getSeq());
		
		//Output file
		try {
			Files.write(outFile, seqs, Charset.forName("UTF-8"));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
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
			
			Region r = new Region(p.getGenome(), p.getChrom(), p.getLocation()-(seqWidth/2)+1, p.getLocation()+(seqWidth/2));
			String seq = seqgen.execute(r);
			//Handle minus strand regions
	        if(strand=='-')
	        	seq = SequenceUtils.reverseComplement(seq);
	        
	        //Ignore regions that are too short (presumably because they are at the edge of the chromosome)
	        if(seq.length()==seqWidth){
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
			if(s.length()==seqWidth){
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
	 * Populate the concatenated SeqFunctions set
	 * @param points
	 */
	public void execute(){
		try {		
			for(ScoredSequence s : scoredSeqs){
				if(s.getSeq().length()==seqWidth){
					Double [] catScores = new Double[catLength];
					int i=0;
					for(SeqFunction function : seqFuncs){
						double min = function.getMinScore();
						double span = (function.getMaxScore()-function.getMinScore());
						double[][] tmpScores = function.score(s.getSeq());
						for(int p=0; p<seqWidth; p++){
							for(int d=0; d<function.scoreDimension(); d++){
			    				catScores[i]=(tmpScores[d][p]-min)/span;
			    				
			    				i++;
			    			}
						}
					}
					seqResults.put(s, catScores);
				}
			}
		} catch (SeqFunctionException e) {
			e.printStackTrace();
		}
	}
	
}
