package edu.psu.compbio.seqcode.gse.datasets.chipchip;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.*;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseException;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseFactory;

public class SQLDataWCEVals extends SQLData {

    public SQLDataWCEVals (String exptname, String exptversion, int genomeid, Set<String> replicates) throws NotFoundException {
        super(exptname,exptversion,genomeid,replicates);
    }

    public double getMax() {
        return 65535;
    }
    public double getMax(String chrom, int start, int stop) throws NotFoundException {
        return 65535;
    }

    public double getValue(int i, int j) {
        return getWCE(i,j);
    }

}
