package edu.psu.compbio.seqcode.gse.seqview.components;

import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.seqview.WarpOptions;

public interface PainterContainer {
    public void addPaintersFromOpts(WarpOptions opts);
    public Genome getGenome();
}
