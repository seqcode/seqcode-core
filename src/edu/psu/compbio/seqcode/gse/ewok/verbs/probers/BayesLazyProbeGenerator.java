package edu.psu.compbio.seqcode.gse.ewok.verbs.probers;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.*;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipBayes;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.SQLData;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.locators.BayesLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.*;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;
import edu.psu.compbio.seqcode.gse.ewok.verbs.binding.RegionProber;
import edu.psu.compbio.seqcode.gse.utils.*;

public class BayesLazyProbeGenerator implements RegionProber<Probe> {
    
    private BayesLocator locator;
    private Genome genome;
    private double strThresh, postThresh;

    public BayesLazyProbeGenerator(Genome g, BayesLocator loc) {
        locator = loc;
        genome = g;
        strThresh = -1.0;
        postThresh = -1.0;
    }

    public BayesLazyProbeGenerator(Genome g, BayesLocator loc, double pt, double st) {
        locator = loc;
        genome = g;
        strThresh = pt;
        postThresh = st;
    }

    public Iterator<edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe> execute(Region r) {
        try {
            ChipChipBayes data = locator.createObject();
            if(strThresh < 0.0 || postThresh < 0.0) { 
                data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            } else { 
                data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd(), postThresh, strThresh);
            }
            return new LazyIterator(data, genome, r.getChrom());
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
    
    private static class LazyIterator implements Iterator<edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe> {
        
        private ChipChipBayes data;
        private Genome genome;
        private String chrom;
        private int index;
        
        private LazyIterator(ChipChipBayes d, Genome g, String chrom) { 
            data = d;
            index = 0;
            genome = g;
            this.chrom = chrom;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#hasNext()
         */
        public boolean hasNext() {
            if(data != null && index >= data.getCount()) { 
                ((SQLData)data).close();
                data = null;                
            }
            
            return data != null;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#next()
         */
        public edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe next() {
            if(!hasNext()) { throw new IllegalStateException(); }
            
            edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe p = createProbe(data, index, genome, chrom);
            
            //System.out.println(index + ": " + p);
            
            index += 1;
            if(index >= data.getCount()) { 
                ((SQLData)data).close();
                data = null;
            }
            
            return p;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#remove()
         */
        public void remove() {
            throw new UnsupportedOperationException();
        } 
        
    }

}
