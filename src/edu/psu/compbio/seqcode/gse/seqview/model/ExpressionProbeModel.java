/**
 * 
 */
package edu.psu.compbio.seqcode.gse.seqview.model;

import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.datasets.expression.*;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Expander;

/**
 * @author tdanford
 */
public class ExpressionProbeModel extends RegionExpanderModel<LocatedExprMeasurement> {
    public ExpressionProbeModel(Expander<Region,LocatedExprMeasurement> exp) { 
        super(exp);
    }
}
