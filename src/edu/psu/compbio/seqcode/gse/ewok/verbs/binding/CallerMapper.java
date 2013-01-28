/*
 * Created on Nov 27, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.binding;

import java.util.Map;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingExtent;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.locators.ExptLocator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper;

public interface CallerMapper extends Mapper<ExptLocator,Expander<Region,BindingExtent>> {
    
    public Map<String,String> getLastParams();
}