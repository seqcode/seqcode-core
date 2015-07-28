package edu.psu.compbio.seqcode.projects.akshay.regulatorydomains;

import java.util.List;

import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;

public class BSRegulatoryRegion extends RegulatoryRegion{

	public BSRegulatoryRegion(Point p, double pStrength, int w,
			List<WeightMatrix> motifs, String s) {
		super(p, pStrength, w, motifs, s);
		// TODO Auto-generated constructor stub
	}

	@Override
	public int compareTo(RegulatoryRegion o) {
		if(bindingStrength<o.bindingStrength){return(-1);}
  		else if(bindingStrength>o.bindingStrength){return(1);}
  		else{return(0);}
	}

}
