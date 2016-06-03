package org.seqcode.projects.shaun.teqseq.core;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;

public class TUnitJunction extends TUnit{

	public TUnitJunction(Genome g, String c, int s, int e, int i, char str) {
		super(g, c, s, e, i, str);
	}
	public TUnitJunction(Region r, int i, char str){
		super(r,i, str);
	}

	/**
	 * Overrides the parent method, because start and end refer to the end & start of two blocks in TUnitJunction
	 */
	public boolean overlaps(AlignHit hit){
		boolean oLap = false;
		if(coords.getChrom().equals(hit.getChrom())){
			AlignBlock[] blocks = hit.getAlignmentBlocks();
			for(int i=0; i<blocks.length-1; i++){
				if(blocks[i].getReferenceEnd() == coords.getStart() && blocks[i+1].getReferenceStart() == coords.getEnd()){
					oLap=true;
					break;
				}
			}
		}
		return oLap;
	}
	
	public boolean overlaps(TUnitJunction j){
		boolean oLap = false;
		if(coords.getChrom().equals(j.getCoords().getChrom()) && coords.getStart()==j.getCoords().getStart() && coords.getEnd()==j.getCoords().getEnd()){
			oLap=true;
		}
		return oLap;
	}
}
