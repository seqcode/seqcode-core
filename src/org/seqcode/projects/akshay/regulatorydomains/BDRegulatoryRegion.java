package org.seqcode.projects.akshay.regulatorydomains;

import java.util.List;

import org.seqcode.genome.location.Point;
import org.seqcode.gse.datasets.motifs.WeightMatrix;

/**
 * 
 * @author akshaykakumanu
 *
 */
public class BDRegulatoryRegion extends RegulatoryRegion{

	public BDRegulatoryRegion(Point p, double pStrength, double pDynamics,
			int w, List<WeightMatrix> motifs, List<Double> motifMarkovThresholds,String seq) {
		super(p, pStrength, pDynamics, w, motifs,motifMarkovThresholds, seq);
		// TODO Auto-generated constructor stub
	}

	@Override
	public int compareTo(RegulatoryRegion o) {
		if(bindingDynamicsIndex<o.bindingDynamicsIndex){return(-1);}
  		else if(bindingDynamicsIndex>o.bindingDynamicsIndex){return(1);}
  		else{return(0);}
	}

}
