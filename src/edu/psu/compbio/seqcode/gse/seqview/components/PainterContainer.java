package edu.psu.compbio.seqcode.gse.seqview.components;

import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.seqview.SeqViewOptions;

public interface PainterContainer {
    public void addPaintersFromOpts(SeqViewOptions opts);
    public Genome getGenome();
}
