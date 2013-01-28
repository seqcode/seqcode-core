/**
 * 
 */
package edu.psu.compbio.seqcode.gse.conservation;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.binding.*;
import edu.psu.compbio.seqcode.gse.utils.Closeable;

/**
 * @author tdanford
 *
 */
public class ProbedRegionFilter<X extends Probe> 
	implements Filter<Region,Region>, Closeable { 

	private RegionProber<X> prober;
	
	public ProbedRegionFilter(RegionProber<X> pr) { 
		prober = pr;
	}
	
	public boolean isClosed() { 
		return prober == null;
	}
	
	public void close() { 
		if(prober instanceof Closeable) { 
			Closeable c = (Closeable)prober;
			if(!c.isClosed()) { c.close(); }
		}
		prober = null;
	}
	
	public Region execute(Region a) {
		/*
		String chrom = a.getChrom();
		if(chrom.startsWith("Un_") || chrom.endsWith("random")) { 
			return null;
		}
		*/
		
		Iterator<X> probes = prober.execute(a);
		boolean hasProbes = probes.hasNext();
		if(probes instanceof Closeable) { 
			Closeable c = (Closeable)probes;
			if(!c.isClosed()) { 
				c.close();
			}
		}
		
		if(hasProbes) { 
			return a;
		} else { 
			return null;
		}
	}
	
}
