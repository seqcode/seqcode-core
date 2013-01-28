package edu.psu.compbio.seqcode.gse.tools.fragdist;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.database.*;

public class DumpFragDist {

    public static void main(String args[]) throws Exception {
        String fd = Args.parseString(args,"fragdist",(String)null);
        if (fd == null) {
            throw new IllegalArgumentException("Must supply --fragdist 'name;version' on the command line");
        }
        String pieces[] = fd.split(";");
        if (pieces.length != 2) {
            throw new IllegalArgumentException("Must supply --fragdist 'name;version' on the command line");            
        }
        ChipChipMetadataLoader loader = new ChipChipMetadataLoader();
        FragDist dist = loader.loadFragDist(pieces[0],pieces[1]);
        for (int i = 0; i < dist.getLength(); i++) {
            System.out.println(i + "\t" + dist.getValue(i));
        }
    }
}