/*
 * Created on Feb 2, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.psu.compbio.seqcode.gse.datasets.locators;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.*;
import edu.psu.compbio.seqcode.gse.datasets.general.CellLine;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptCondition;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptTarget;

public class ExptLocatorFilterOptions {
    
    public static final int AGILENT_TYPE = 0;
    public static final int MSP_TYPE = 1;
    public static final int BAYES_TYPE = 2;

    private Map<String,Object> optsMap;
    
    public ExptLocatorFilterOptions() { 
        optsMap = new HashMap<String,Object>();
    }
    
    public Set<String> getOptionKeys() { return new HashSet<String>(optsMap.keySet()); }
    public boolean hasOptionValue(String k) { return optsMap.containsKey(k); }
    public Object getOptionValue(String k) { return optsMap.get(k); }
    public void setOptionValue(String k, String v) { optsMap.put(k, v); }
    
    public void clearOptionValue(String k) { optsMap.remove(k); }
    
    public void setCells(CellLine c) { optsMap.put("cells", c); }
    public void setFactor(ExptTarget f) { optsMap.put("factor", f); }
    public void setCondition(ExptCondition c) { optsMap.put("condition", c); }
    
    public void setChipChipType() { optsMap.put("expt_type", AGILENT_TYPE); }
    public void setMSPType() { optsMap.put("expt_type", MSP_TYPE); }
    public void setBayesType() { optsMap.put("expt_type", BAYES_TYPE); }
    
    public CellLine getCells() { return optsMap.containsKey("cells") ? (CellLine)optsMap.get("cells") : null; }
    public ExptTarget getFactor() { return optsMap.containsKey("factor") ? (ExptTarget)optsMap.get("factor") : null; }
    public ExptCondition getCondition() { return optsMap.containsKey("condition") ? (ExptCondition)optsMap.get("condition") : null; }
    
    public int getExptType() { return (Integer)optsMap.get("expt_type"); }
    
    public void clearCells() { optsMap.remove("cells"); }
    public void clearCondition() { optsMap.remove("condition"); }
    public void clearFactor() { optsMap.remove("factor"); }
}
