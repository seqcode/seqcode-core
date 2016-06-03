package org.seqcode.projects.sequtils;

import java.util.ArrayList;
import java.util.List;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import org.seqcode.gse.utils.sequence.SequenceUtils;


/**
 * Filter a list of StrandedBaseCounts by genomic base at a defined position
 * @author mahony
 *
 */
public class StrandedBaseCountFilterByBase {
	private SequenceGenerator seqgen=null;
	char base; //Base required
	int baseRelPosition; //Position base is reqd to appear at, wrt the five prime position of the SBC. 
	
	public StrandedBaseCountFilterByBase(GenomeConfig genConfig, char base, int baseRelPosition){
		seqgen = genConfig.getSequenceGenerator();
		this.base = base;
		this.baseRelPosition = baseRelPosition;
		
	}
	
	/**
	 * This method filters an input list of StrandedBaseCounts by the base rules defined in the constructor.
	 * Assumes that the input list of reads is defined by a single Region on the genome. 
	 * @param currentReg
	 * @param reads
	 * @return
	 */
	public List<StrandedBaseCount> execute(Region currentReg, List<StrandedBaseCount> reads){
		List<StrandedBaseCount> output = new ArrayList<StrandedBaseCount>();
		//Get sequence if required
		char [] seq=null;
		if(seqgen!=null)
			seq = seqgen.execute(currentReg).toCharArray();
		
		for(StrandedBaseCount sbc : reads){
			if(seq!=null && getBaseAtPosition(sbc, baseRelPosition, currentReg, seq)==base){
				output.add(sbc);
			}
		}
		
		return output;
	}
	
	private char getBaseAtPosition(StrandedBaseCount a, int position, Region queryReg, char[] seq){
		char b = '.';
		int wantedPos = a.getStrand()=='+' ? 
				a.getCoordinate()+position : 
				a.getCoordinate()-position;
		if(wantedPos>=queryReg.getStart() && wantedPos<queryReg.getEnd()){
			b = a.getStrand()=='+' ? seq[wantedPos-queryReg.getStart()] : SequenceUtils.complementChar(seq[wantedPos-queryReg.getStart()]);
		}
		return b;
	}
}
