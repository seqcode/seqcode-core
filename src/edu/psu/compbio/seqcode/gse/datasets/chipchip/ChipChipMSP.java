package edu.psu.compbio.seqcode.gse.datasets.chipchip;

import java.util.*;

import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.io.parsing.*;

public interface ChipChipMSP extends GenericExperiment {
    public float getRatio(int i);
    public float getX(int i);
    public float getPval(int i);
    public float getPval3(int i);
    public float getMedianOfRatios(int i);


}
