/*
 * Created on Mar 21, 2007
 */
package edu.psu.compbio.seqcode.gse.datasets.chipchip;

import edu.psu.compbio.seqcode.gse.datasets.general.CellLine;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptCondition;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptTarget;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

/**
 * @author tdanford
 */
public class ChipChipMetadataQuery {

    private Genome genome;
    private CellLine cells;
    private ExptCondition cond;
    private ExptTarget factor;

    public ChipChipMetadataQuery() {  
        genome = null;
        cells = null;
        cond = null;
        factor = null;
    }

    public Genome getGenome() { return genome; }
    public void setGenome(Genome g) { genome = g; }

    public CellLine getCells() { return cells; }
    public void setCells(CellLine c) { cells = c; }

    public ExptCondition getCondition() { return cond; }
    public void setCondition(ExptCondition c) { cond = c; }

    public boolean equals(Object o) { 
        if(!(o instanceof ChipChipMetadataQuery)) { return false; }
        return true;
    }

    public int hashCode() { 
        int code = 17;
        code += genome != null ? genome.hashCode() : 0; code *= 37;
        code += cells != null ? cells.hashCode() : 0; code *= 37;
        return code;
    }
}
