/**
 * 
 */
package edu.psu.compbio.seqcode.gse.seqview.model;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;

/**
 * @author tdanford
 */
public class BindingEventModel extends RegionExpanderModel<BindingEvent> {
    public BindingEventModel(Expander<Region,BindingEvent> exp) { 
        super(exp);
    }
}
