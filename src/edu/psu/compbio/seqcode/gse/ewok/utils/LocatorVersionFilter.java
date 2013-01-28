/*
 * Created on Apr 3, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.utils;

import edu.psu.compbio.seqcode.gse.datasets.locators.ExptLocator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Filter;

/**
 * @author tdanford
 */
public class LocatorVersionFilter implements Filter<ExptLocator,ExptLocator> {
    
    private String version;

    public LocatorVersionFilter(String v) {
        version = v;
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Filter#execute(null)
     */
    public ExptLocator execute(ExptLocator a) {
        if(a.getNameVersion().version.equals(version)) { 
            return a;
        } else { 
            return null;
        }
    }

}
