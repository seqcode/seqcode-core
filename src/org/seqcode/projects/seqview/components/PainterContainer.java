package org.seqcode.projects.seqview.components;

import org.seqcode.genome.Genome;
import org.seqcode.projects.seqview.SeqViewOptions;

public interface PainterContainer {
	public void addPaintersFromOpts(SeqViewOptions opts);

	public Genome getGenome();
}
