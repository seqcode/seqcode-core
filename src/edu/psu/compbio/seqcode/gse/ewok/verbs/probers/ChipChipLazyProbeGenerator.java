package edu.psu.compbio.seqcode.gse.ewok.verbs.probers;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.*;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipData;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.SQLData;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.locators.ChipChipLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.*;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;
import edu.psu.compbio.seqcode.gse.utils.*;

public class ChipChipLazyProbeGenerator implements Expander<Region,Probe> {
    
    private ChipChipLocator locator;
    private Genome genome;

    public ChipChipLazyProbeGenerator(Genome g, ChipChipLocator loc) {
        locator = loc;
        genome = g;
    }

    public Iterator<edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe> execute(Region r) {
        try {
            ChipChipData data = locator.createObject();
            data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            return new LazyIterator(data, genome, r.getChrom());
        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
    }
    
    private edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe createProbe(ChipChipData data, int index, Genome genome, String chrom) { 
        edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe p = 
            new edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe(genome, chrom, data.getPos(index));
        
        double total = 0.0;
        for(int i = 0; i < data.getReplicates(index); i++) { 
            total += data.getRatio(index, i);
        }
        if(data.getReplicates(index) > 0) { 
            total /= (double)data.getReplicates(index);
        }
        
        //String n = data.getName();
        String n = locator.name + "," + locator.version;
        p.addValue(n, total);
        
        return p;
    }
    
    private class LazyIterator implements Iterator<Probe>, Closeable {
        
        private ChipChipData data;
        private Genome genome;
        private String chrom;
        private int index;
        
        private LazyIterator(ChipChipData d, Genome g, String chrom) { 
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
            if(isClosed()) { throw new IllegalStateException(); }
            if(!hasNext()) { throw new IllegalStateException(); }
            
            edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe p = createProbe(data, index, genome, chrom);
            
            //System.out.println(index + ": " + p);
            
            index += 1;
            if(index >= data.getCount()) { close(); }
            
            return p;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#remove()
         */
        public void remove() {
            throw new UnsupportedOperationException();
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.utils.Closeable#close()
         */
        public void close() {
            ((SQLData)data).close();
            data = null;
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.utils.Closeable#isClosed()
         */
        public boolean isClosed() {
            return data == null;
        } 
        
    }
}
