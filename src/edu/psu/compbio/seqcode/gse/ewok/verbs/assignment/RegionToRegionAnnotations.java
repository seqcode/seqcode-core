package edu.psu.compbio.seqcode.gse.ewok.verbs.assignment;

import java.util.*;

import javax.swing.JProgressBar;

import edu.psu.compbio.seqcode.gse.datasets.binding.*;
import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Gene;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.binding.*;

/* Doesn't provide annotations itself, but the Items are the results of
   applying Expander<Region,Region> evts to the whole genome
*/

public class RegionToRegionAnnotations extends CachedAnnotations<Region,Region> {
    
    private RefGeneGenerator generator;
    private JProgressBar progressBar;
    
    public RegionToRegionAnnotations(Genome g, Expander<Region,Region> evts) {
        super();
        progressBar = null;
        init(evts,g);
    }

    public RegionToRegionAnnotations(Genome g, Expander<Region,Region> evts, JProgressBar pbar) { 
        super();
        progressBar = pbar;
        init(evts,g);
    }
    public void init(Expander<Region,Region> evts, Genome g) {

        ChromRegionIterator chroms = new ChromRegionIterator(g);
        Iterator<Region> regionChroms = 
            new MapperIterator<NamedRegion,Region>(new CastingMapper<NamedRegion,Region>(), chroms);
        Iterator<Region> evtItr = new ExpanderIterator<Region,Region>(evts, regionChroms);
        
        LinkedList<Region> totalEvents = new LinkedList<Region>();
        while(evtItr.hasNext()) { 
            totalEvents.addLast(evtItr.next());
        }
        
        addItems(totalEvents);
        
        if(progressBar != null) { 
            progressBar.setMaximum(getNumItems());
        }
        
        generator = new RefGeneGenerator<Region>(g);        
        addAnnotations("genes", generator);
    }
    
    public void markProgress() { 
        if(progressBar != null) {
            progressBar.setValue(progressBar.getValue() + 1);
        }
    }
    
    public Expander<Region,Gene> getGeneGenerator() { return generator; }
}
