package edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder.motifs;

import java.io.IOException;
import java.text.ParseException;
import java.util.HashMap;
import java.util.HashSet;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder.framework.SeqUnwinderConfig;


public class ScoreMotif {
	
	protected SeqUnwinderConfig seqConfig;
	protected HashMap<String,double[]> scores = new HashMap<String,double[]>();
	
	public ScoreMotif(SeqUnwinderConfig seqcon) {
		seqConfig = seqcon;
	}
	
	public void execute(){
		for(WeightMatrix wm : seqConfig.getDiscrimMotifs()){
			scores.put(wm.getName(), new double[seqConfig.getDiscrimMotifs().size()]);
		}
		
		for(WeightMatrix mot : seqConfig.getDiscrimMotifs()){
			HashSet<String> motKmers = WeightMatrix.getConsensusKmers(mot, seqConfig.getKmin(), seqConfig.getKmax());
			int j=0;
			for(String modName : seqConfig.getMNames()){
				double weight=0.0;
				for(String s : motKmers){
					String revS = SequenceUtils.reverseComplement(s);
					int indS = RegionFileUtilities.seq2int(s);
					int indRevS = RegionFileUtilities.seq2int(revS);
					int KmerInd  = indS<indRevS ? indS : indRevS;
					int ind = seqConfig.getKmerBaseInd(s) + KmerInd;
					weight = weight + seqConfig.getKmerWeights().get(modName)[ind];
				}
				scores.get(mot.getName())[j] = weight;
				j++;
			}
		}
		
		seqConfig.setDiscrimMotifScores(scores);
		
	}
	
	
	
	
	public static void main(String[] args) throws IOException, ParseException {
		
	}

}
