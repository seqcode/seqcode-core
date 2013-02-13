/**
 * 
 */
package edu.psu.compbio.seqcode.gse.seqview.model;

import edu.psu.compbio.seqcode.gse.datasets.chippet.ChipPetDatum;
import edu.psu.compbio.seqcode.gse.datasets.expression.*;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper;
import edu.psu.compbio.seqcode.gse.utils.Interval;

/**
 * @author tdanford
 */
public class ScoreTrackModel extends RegionExpanderModel<ChipPetDatum> {
    public ScoreTrackModel(Expander<Region,ChipPetDatum> exp) { 
        super(exp);
    }
}
