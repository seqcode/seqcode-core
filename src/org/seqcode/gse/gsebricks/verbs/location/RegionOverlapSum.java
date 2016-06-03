/*
 * Author: tdanford
 * Date: Apr 2, 2009
 */
package org.seqcode.gse.gsebricks.verbs.location;

import java.util.Collection; 
import java.util.ArrayList;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.utils.Interval;
import org.seqcode.gse.utils.OverlapSum;


public class RegionOverlapSum extends OverlapSum {
	
	private Region totalRegion;

	public RegionOverlapSum(Region r) { 
		totalRegion = r;
	}
	
	public RegionOverlapSum(Genome g, String chr) { 
		totalRegion = new Region(g, chr, 1, g.getChromLength(chr));
	}
	
	public Collection<Region> collectRegions(int threshold) { 
		Collection<Interval> intvs = super.collect(threshold);
		
		ArrayList<Region> rs = new ArrayList<Region>();
		for(Interval intv : intvs) { 
			rs.add(new Region(totalRegion.getGenome(), totalRegion.getChrom(),
					intv.start, intv.end));
		}
		return rs;
	}
	
	public Region getTotalRegion() { return totalRegion; }
	
    public void addRegion(Region r) {
        if(!totalRegion.getGenome().equals(r.getGenome())) { 
        	throw new IllegalArgumentException(r.getGenome().toString()); 
        }
        if(!totalRegion.getChrom().equals(r.getChrom())) { 
        	throw new IllegalArgumentException(r.getChrom()); 
        }
        
        int start = Math.max(totalRegion.getStart(), r.getStart());
        int end = Math.min(totalRegion.getEnd(), r.getEnd());
        addInterval(start, end);
    }
}
