package edu.psu.compbio.seqcode.gse.ewok.verbs.probers;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.*;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipMSP;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.MSPProbe;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.SQLData;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.SQLMSP;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.locators.*;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.*;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;
import edu.psu.compbio.seqcode.gse.utils.*;

public class MSPImmediateProbeGenerator<X extends Region> implements Expander<X,MSPProbe>, Closeable {
    
    private MSPLocator locator;
    private ChipChipMSP data;
    private Genome genome;

    public MSPImmediateProbeGenerator(Genome g, MSPLocator loc) {
        locator = loc;
        data = loc.createObject();
        genome = g;
        
        System.out.println("++++ ChipChipMSP Opened.");
    }

    public Iterator<MSPProbe> execute(X r) {
        if(isClosed()) { throw new IllegalStateException(); }
        try {
            data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            
            LinkedList<MSPProbe> probes = 
                new LinkedList<MSPProbe>();
            for(int i = 0; i < data.getCount(); i++) { 
                probes.addLast(createProbe(data, i, genome, r.getChrom()));
            }
            return probes.iterator();

        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
    }
    
    private MSPProbe createProbe(ChipChipMSP data, int index, 
            Genome genome, String chrom) {
        
        String key = locator.name + "," + locator.version;
        MSPProbe p = 
            new MSPProbe(genome, chrom, data.getPos(index), key, data, index);
        
        return p;
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.utils.Closeable#close()
     */
    public void close() {
        ((SQLMSP)data).close();
        data = null;
        System.out.println("---- ChipChipMSP Closed.");
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.utils.Closeable#isClosed()
     */
    public boolean isClosed() {
        return data == null;
    }
}
