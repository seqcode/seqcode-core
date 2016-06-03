package org.seqcode.projects.shaun;


import java.util.Collection;
import java.util.ArrayList;

import org.seqcode.genome.location.Gene;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;

public class ChipSeqPeak implements Comparable<ChipSeqPeak>{
	public Region coords;
	public Point peak;
	public double ipHits;
	public double backHits;
	public double Z;
	public double overrep;
	public char strand; //In here to support CLIP 
	public Gene nearestGene=null;
	public int distToGene=Integer.MAX_VALUE;
    public Collection<Region> annotations;
	
	public ChipSeqPeak(){
		coords=null;
		peak=null;
		ipHits=0;
		backHits=0;
		Z=0;
		strand='+';
	}
	
    public void addAnnotation(Region r) {
        if (annotations == null) {
            annotations = new ArrayList<Region>();
        }
        annotations.add(r);
    }

	public int compareTo(ChipSeqPeak p) {
		if(Z>p.Z){return(-1);}
		else if(Z<p.Z){return(1);}
		else{return(0);}
	}
}