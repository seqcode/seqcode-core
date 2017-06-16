package org.seqcode.deepseq;

import org.seqcode.genome.Genome;

public class ExtReadHit extends ReadHit {

	public ExtReadHit(Genome g, ReadHit r, int startShift, int fivePrimeExt, int threePrimeExt) {
		super(r.getChrom(),
				r.getStrand() == '+' ? Math.max(1, r.getStart() + startShift - fivePrimeExt)
						: Math.max(1, r.getStart() - startShift - threePrimeExt),
				r.getStrand() == '+' ? Math.min(g.getChromLength(r.getChrom()), r.getEnd() + startShift + threePrimeExt)
						: Math.min(g.getChromLength(r.getChrom()), r.getEnd() - startShift + fivePrimeExt),
				r.getStrand(), r.getWeight());
	}

}
