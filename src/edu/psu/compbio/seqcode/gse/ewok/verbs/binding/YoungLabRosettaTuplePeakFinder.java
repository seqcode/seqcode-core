/*
 * Created on Mar 20, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.binding;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.binding.ProbedBindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.MSPProbe;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.utils.*;

/**
 * @author tdanford
 */
public class YoungLabRosettaTuplePeakFinder implements Filter<Vector<MSPProbe>,BindingEvent> {
    
    private Genome genome;
    private String key;
    private double p3Cutoff, pCutoff, nCutoff, twoCutoff;
    
    public YoungLabRosettaTuplePeakFinder(Genome g, String k, double p3, double p1, double n, double two) { 
        genome = g;
        key = k;
        
        p3Cutoff = p3;
        pCutoff = p1;
        nCutoff = n;
        twoCutoff = two;        
    }

    public YoungLabRosettaTuplePeakFinder(Genome g, String k) {
        genome = g;
        key = k;
        
        p3Cutoff = 0.001;
        pCutoff = 0.001;
        nCutoff = 0.1;
        twoCutoff = 0.005;
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Distiller#execute(null)
     */
    public BindingEvent execute(Vector<MSPProbe> p) {
        return distillCurrent(p);
    }
    
    private boolean currentMeetsPredicate(Vector<MSPProbe> triple) {
        if(triple.size() != 3) { throw new IllegalArgumentException(); }
        if(triple.get(0) == null) { return false; }

		if(!triple.get(0).containsKey(key)) { 
			//return false; 
			throw new IllegalStateException();
		}
        
        for(int i = 1; i < 3; i++) { 
			if(!triple.get(i).containsKey(key)) { 
				throw new IllegalStateException();
				//return false; 
			}

            if(!triple.get(i).getChrom().equals(triple.get(0).getChrom())) { 
                //return false; 
				throw new IllegalStateException();
            }

			if(triple.get(i).getLocation() < triple.get(i-1).getLocation()) { 
				throw new IllegalStateException();
			}
        }

        MSPProbe left = triple.get(0), middle = triple.get(1), right = triple.get(2);
        
        boolean firstCond = middle.getPval3(key) <= p3Cutoff;
        if(!firstCond) { return false; }
        
        boolean sc1 = middle.getPval(key) <= pCutoff;
        boolean sc2 = left.getPval(key) <= nCutoff || right.getPval(key) <= nCutoff;
        boolean secondCond = sc1 && sc2;
        if(secondCond) { return true; }
        
        int count = 0;
        for(int i = 0; i < 3; i++) { 
            if(triple.get(i).getPval(key) <= twoCutoff) { 
                count += 1; 
            }
        }
        boolean thirdCond = count >= 2;
        return thirdCond;
    }
    
    private BindingEvent distillCurrent(Vector<MSPProbe> triple) { 
        if(currentMeetsPredicate(triple)) { 
            int ns = triple.get(0).getLocation(), ne = triple.get(2).getLocation();
            String c = triple.get(1).getChrom();
            LinkedList<Probe> probes = new LinkedList<Probe>();
            probes.addLast(triple.get(1));
            return new ProbedBindingEvent(genome, c, ns, ne, 
                    triple.get(1).getRatio(key), triple.get(1).getPval3(key), key, probes);
        }
        
        return null;        
    }

}
