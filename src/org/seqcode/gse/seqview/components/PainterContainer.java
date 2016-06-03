package org.seqcode.gse.seqview.components;

import org.seqcode.genome.Genome;
import org.seqcode.gse.seqview.SeqViewOptions;

public interface PainterContainer {
    public void addPaintersFromOpts(SeqViewOptions opts);
    public Genome getGenome();
}
