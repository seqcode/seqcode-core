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
public class BindingFilter implements Filter<Region,Region> {
    
    private PeakCaller caller;
    
    public BindingFilter(PeakCaller pc) { 
        caller = pc;
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Filter#execute(java.lang.Object)
     */
    public Region execute(Region a) {
        Iterator<BindingEvent> evtItr = caller.execute(a);
		Region ret = null;
        if(evtItr.hasNext()) { ret = a; }

		/*
		System.out.println("\n" + a.toString() + " --> {");
		while(evtItr.hasNext()) { 
			System.out.println(evtItr.next());
		}
		System.out.println("}");
		*/

		return ret;
    }
}
