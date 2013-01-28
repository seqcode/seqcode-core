/*
 * Created on Mar 10, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.binding;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;

/**
 * @author tdanford
 */
public class BestBindingFilter implements Filter<Region,BindingEvent> {
    
    private PeakCaller caller;
    private BindingEventAnalyzer analyzer;
    
    public BestBindingFilter(PeakCaller pc, BindingEventAnalyzer a) { 
        caller = pc;
        analyzer = a;
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Filter#execute(java.lang.Object)
     */
    public BindingEvent execute(Region a) {
        Iterator<BindingEvent> evtItr = caller.execute(a);
		BindingEvent ret = null;
        double maxScore = -Double.MAX_VALUE;
        
        while(evtItr.hasNext()) { 
            BindingEvent e = evtItr.next();
            double score = analyzer.execute(e);
            if(ret == null || score > maxScore) { 
                maxScore = score;
                ret = e;
            }
        }

		return ret;
    }
}
