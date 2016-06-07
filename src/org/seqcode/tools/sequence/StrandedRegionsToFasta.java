package org.seqcode.tools.sequence;

import java.io.*;

import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gsebricks.verbs.location.StrandedRegionParser;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.Pair;


/**
 *  reads a list of stranded regions (one per line) on STDIN.
    Produces on stdout a fasta file containing those regions.
    The genome is specified on the command line as
    --species "Mus musculus;mm8"
*/

public class StrandedRegionsToFasta {

    public static void main(String args[]) {
        try {
            Pair<Species,Genome> pair = Args.parseGenome(args);
            Species organism = pair.car();
            Genome genome = pair.cdr();
            StrandedRegionParser parser = new StrandedRegionParser(genome);
            SequenceGenerator seqgen = new SequenceGenerator();
            boolean cache = Args.parseFlags(args).contains("cache");
            seqgen.useCache(cache);
            String line;
            BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                StrandedRegion r = parser.execute(line);
                String seq = seqgen.execute(r);
                if(r.getStrand()=='-'){
                	seq = SequenceUtils.reverseComplement(seq);
                }
                System.out.println(">"+r.toString()+"\n"+seq);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}
