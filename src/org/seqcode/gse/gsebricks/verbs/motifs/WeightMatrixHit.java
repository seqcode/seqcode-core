/**
 * Created 2/15/08
 */
package org.seqcode.gse.gsebricks.verbs.motifs;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.ScoredStrandedRegion;
import org.seqcode.gse.datasets.motifs.WeightMatrix;

/**
 * @author tdanford
 *
 */
public class WeightMatrixHit extends ScoredStrandedRegion {
	
	private WeightMatrix matrix;

	public WeightMatrixHit(Genome g, String c, int start, int end, double score, char strand, WeightMatrix wm) {
		super(g, c, start, end, score, strand);
		matrix = wm;
	}
	
	public WeightMatrix getMatrix() { return matrix; }
	
	public String toString() { 
		String s = super.toString();
		return String.format("%s (%s)", s, matrix.toString());
	}
	
	public int hashCode() { 
		int code = super.hashCode();
		code += matrix.hashCode(); code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof WeightMatrixHit)) { 
			return false;
		}
		WeightMatrixHit h = (WeightMatrixHit)o;
		if(!matrix.equals(h.matrix)) { return false; }
		return super.equals(h);
	}
}
