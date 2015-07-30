package edu.psu.compbio.seqcode.projects.akshay.regulatorydomains;

import java.util.List;

import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
/**
 * 
 * @author akshaykakumanu
 *
 */
public class BDRegulatoryRegion extends RegulatoryRegion{

	public BDRegulatoryRegion(Point p, double pStrength, double pDynamics,
			int w, List<WeightMatrix> motifs, String seq) {
		super(p, pStrength, pDynamics, w, motifs, seq);
		// TODO Auto-generated constructor stub
	}

	@Override
	public int compareTo(RegulatoryRegion o) {
		if(bindingDynamicsIndex<o.bindingDynamicsIndex){return(-1);}
  		else if(bindingDynamicsIndex>o.bindingDynamicsIndex){return(1);}
  		else{return(0);}
	}

}
