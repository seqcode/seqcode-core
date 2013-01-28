package edu.psu.compbio.seqcode.gse.ewok.verbs.probers;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.*;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipBayes;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.SQLBayes;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.SQLData;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.locators.BayesLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.*;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;
import edu.psu.compbio.seqcode.gse.ewok.verbs.binding.RegionProber;
import edu.psu.compbio.seqcode.gse.utils.*;

public class BayesImmediateProbeGenerator implements RegionProber<Probe> {
    
    private ChipChipBayes data;
    private Genome genome;
    private double postThresh, strThresh;

    public BayesImmediateProbeGenerator(Genome g, BayesLocator loc) {
        data = loc.createObject();
        genome = g;
        postThresh = -1.0;
        strThresh = -1.0;
    }

    public BayesImmediateProbeGenerator(Genome g, BayesLocator loc, double pt, double st) {
        data = loc.createObject();
        genome = g;
        postThresh = pt;
        strThresh = st;
    }

    public BayesImmediateProbeGenerator(Genome g, ChipChipBayes d) {
        data = d;
        genome = g;
    }
    
    public Iterator<edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe> execute(Region r) {
        try {
            if(strThresh > 0.0 && postThresh > 0.0) { 
                data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            } else { 
                data.window(r.getChrom(), r.getStart(), r.getEnd(), postThresh, strThresh);
            }
            
            LinkedList<edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe> probes = 
                new LinkedList<edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe>();
            
            for(int i = 0; i < data.getCount(); i++) { 
                probes.addLast(createProbe(data, i, genome, r.getChrom()));
            }
            return probes.iterator();

        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
    }
    
    private static edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe createProbe(ChipChipBayes data, int index, Genome genome, String chrom) { 
        edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe p = 
            new edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe(genome, chrom, data.getPos(index));
        
        double[] total = new double[2];
        
        total[0] = data.getStrength(index);
        total[1] = data.getPosterior(index);
        
        String n = data.getName();
        p.addValue(n, total);
        
        return p;
    }
}
