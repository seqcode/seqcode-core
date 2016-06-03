package org.seqcode.projects.akshay.regulatorydomains;

import java.util.List;

import org.seqcode.genome.location.Point;
import org.seqcode.gse.datasets.motifs.WeightMatrix;

/**
 * 
 * @author akshaykakumanu
 *
 */
public class BSRegulatoryRegion extends RegulatoryRegion{
	/**
	 * 
	 * @param p
	 * @param pStrength
	 * @param w
	 * @param motifs
	 * @param s
	 */
	public BSRegulatoryRegion(Point p, double pStrength, int w,
			List<WeightMatrix> motifs,  List<Double> motifMarkovThresholds, String s) {
		super(p, pStrength, w, motifs, motifMarkovThresholds, s);
		// TODO Auto-generated constructor stub
	}
	
	/**
	 * 
	 * @param p
	 * @param pStrength
	 * @param pDynamics
	 * @param w
	 * @param motifs
	 * @param s
	 */
	public BSRegulatoryRegion(Point p, double pStrength, double pDynamics, int w,
			List<WeightMatrix> motifs,  List<Double> motifMarkovThresholds, String s) {
		super(p, pStrength, pStrength,w, motifs, motifMarkovThresholds, s);
		// TODO Auto-generated constructor stub
	}

	@Override
	public int compareTo(RegulatoryRegion o) {
		if(bindingStrength<o.bindingStrength){return(-1);}
  		else if(bindingStrength>o.bindingStrength){return(1);}
  		else{return(0);}
	}

}
