package edu.psu.compbio.seqcode.gse.tools.sequence;

import java.io.*;

import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.StrandedRegionParser;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;

/* reads a list of stranded regions (one per line) on STDIN.
   Produces on stdout a fasta file containing those regions.
   The genome is specified on the command line as
   --species "Mm;mm8"
*/

public class StrandedRegionsToFasta {

    public static void main(String args[]) {
        try {
            Pair<Organism,Genome> pair = Args.parseGenome(args);
            Organism organism = pair.car();
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
