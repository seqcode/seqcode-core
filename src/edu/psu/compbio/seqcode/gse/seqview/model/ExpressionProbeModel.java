/**
 * 
 */
package edu.psu.compbio.seqcode.gse.seqview.model;

import edu.psu.compbio.seqcode.gse.datasets.expression.*;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;

/**
 * @author tdanford
 */
public class ExpressionProbeModel extends RegionExpanderModel<LocatedExprMeasurement> {
    public ExpressionProbeModel(Expander<Region,LocatedExprMeasurement> exp) { 
        super(exp);
    }
}
