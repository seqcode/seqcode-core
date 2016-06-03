package org.seqcode.projects.shaun;


public class ConsensusSequenceScoreProfile {
	private ConsensusSequence consensus;
	private int[] forward, reverse;
	
	public ConsensusSequenceScoreProfile(ConsensusSequence c, int[] f, int[] r) {
		consensus = c;
		if(f.length != r.length) { throw new IllegalArgumentException(); }
		forward = (int[])f.clone();
		reverse = (int[])r.clone();
	}
	
	public int length() { return forward.length; }
	public ConsensusSequence getConsensus() { return consensus; }
	public int[] getForwardScores() { return forward; }
	public int[] getReverseScores() { return reverse; }
	
	/** Strand of lowest mismatch over the whole sequence*/
	public char getLowestMismatchStrand() { 
		int min = Integer.MAX_VALUE; 
		int minScore = consensus.getMaxMismatch();
		char str = '+';
		for(int i = 0; i < forward.length; i++) { 
			int ms = getLowestMismatch(i);
			if(min == Integer.MAX_VALUE || ms < minScore) { 
				minScore= ms; 
				min = i;
				str = getLowestMismatchStrand(i);
			}
		}
		return str;
	}
	
	/** Strand of lowest mismatch (from 2 strands) at this position */
	public char getLowestMismatchStrand(int i) { 
		if(forward[i] <= reverse[i]) { 
			return '+';
		} else { 
			return '-';
		}
	}
	
	/** Lowest mismatch (from 2 strands) at this position */
	public int getLowestMismatch(int i) { 
		return Math.min(forward[i], reverse[i]); 
	}
	/** Lowest mismatch over the whole sequence*/
	public int getLowestMismatchBoundedWin(int l, int r){
		int min = consensus.getMaxMismatch();
		if(l>0 && r>=0 && l<forward.length && r<forward.length &&  l<r){
			for(int i=l; i<=r; i++){
				if(getLowestMismatch(i)<min)
					min = getLowestMismatch(i);
			}
			return min;
		}else{
			return -1;
		}
	}
	/** Lowest mismatch over the whole sequence*/
	public int getLowestMismatch(){
		return(getLowestMismatch(getLowestMismatchIndex()));
	}
	/** Index of lowest mismatch level. If multiple equivalent, leftmost is returned**/
	public int getLowestMismatchIndex() { 
		int min = Integer.MAX_VALUE; 
		int minScore = consensus.getMaxMismatch();
		for(int i = 0; i < forward.length; i++) { 
			int ms = getLowestMismatch(i);
			if(min == Integer.MAX_VALUE || ms<minScore) { 
				minScore = ms;
				min = i;
			}
		}
		return min;
	}
}
