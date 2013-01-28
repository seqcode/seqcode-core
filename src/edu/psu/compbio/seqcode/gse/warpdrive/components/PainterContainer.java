package edu.psu.compbio.seqcode.gse.warpdrive.components;

import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.warpdrive.WarpOptions;

public interface PainterContainer {
    public void addPaintersFromOpts(WarpOptions opts);
    public Genome getGenome();
}
