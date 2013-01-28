package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.util.Iterator;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Gene;
import edu.psu.compbio.seqcode.gse.ewok.nouns.GeneDomainData;

public class DomainAnalysis<DOMAIN extends Region> 
	implements Mapper<Gene,GeneDomainData> {

	private int window;
	private Expander<Region,DOMAIN> domainCaller;
	
	public DomainAnalysis(Expander<Region,DOMAIN> dc, int w) { 
		domainCaller = dc;
		window = w;
	}

	public GeneDomainData execute(Gene g) {
		GeneDomainData gdd = new GeneDomainData(g, window);
		Iterator<DOMAIN> itr = domainCaller.execute(gdd.getWindow());
		while(itr.hasNext()) { 
			gdd.addDomain(itr.next());
		}
		return gdd;
	}
	
	
}
