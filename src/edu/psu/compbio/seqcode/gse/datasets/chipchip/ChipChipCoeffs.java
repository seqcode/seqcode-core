package edu.psu.compbio.seqcode.gse.datasets.chipchip;

import java.io.*;
import java.util.*;

import edu.psu.compbio.seqcode.gse.utils.io.parsing.*;

/* represents a set of coefficients from the fragment shear size distribution . */

public interface ChipChipCoeffs {
    // get the relative intensity and distance d, which
    // can be positive or negative (allowing for asymmetric
    // distributions)
    public double getCoeff(int d);
    public int getMaxDist();
}
