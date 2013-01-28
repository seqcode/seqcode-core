package edu.psu.compbio.seqcode.gse.tools.core;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;


/**
 * Dumps gene annotations in GFF format
 *
 */

public class DumpGeneGFF {

    public static void main(String args[]) throws Exception {
        Genome genome = Args.parseGenome(args).cdr();
        RefGeneGenerator genegen = Args.parseGenes(args).get(0);
        genegen.retrieveExons(true);
        genegen.setWantAlias(true);
        Map<String,Collection<Gene>> genes = new HashMap<String,Collection<Gene>>();
        Iterator<Gene> all = genegen.getAll();
        while (all.hasNext()) {
            Gene g = all.next();
            String n = g.getName();
            if (n.equals(g.getID())) {
                Collection<String> notid = g.getNonIDNames();
                if (notid.size() > 0) {
                    Iterator<String> iter = notid.iterator();
                    n = iter.next();
                }
            }

            if (!genes.containsKey(n)) {
                genes.put(n, new ArrayList<Gene>());
            }
            genes.get(n).add(g);
        }

        for (String id : genes.keySet()) {
            String chrom = null;
            int minpos = Integer.MAX_VALUE, maxpos = 0;
            char strand = '+';
            boolean mixedchroms = false;
            for (Gene g : genes.get(id)) {
                if (chrom == null) {
                    chrom = g.getChrom();
                    strand = g.getStrand();
                } else {
                    if (!g.getChrom().equals(chrom)) {
                        mixedchroms = true;
                    }

                }

                minpos = Math.min(minpos, g.getStart());
                maxpos = Math.max(maxpos, g.getEnd());
            }
            if (!chrom.matches("^.*")) {
                chrom = "chr" + chrom;
            }

            System.out.println(String.format("%s\tprotein_coding\tgene\t%d\t%d\t.\t%s\t.\tID=%s",
                                             chrom,minpos,maxpos,Character.toString(strand),id));

            for (Gene g : genes.get(id)) {
                System.out.println(String.format("%s\tprotein_coding\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s",
                                                 g.getChrom(),g.getStart(),g.getEnd(),Character.toString(g.getStrand()),g.getID(),id));
                if (g instanceof ExonicGene) {
                    ExonicGene exonic = (ExonicGene)g;
                    Iterator<Region> iter = exonic.getExons();
                    int count = 1;
                    while (iter.hasNext()) {
                        Region e = iter.next();
                        System.out.println(String.format("%s\tprotein_coding\texon\t%d\t%d\t.\t%s\t.\tID=%s.%d;Parent=%s",
                                                         chrom,e.getStart(),e.getEnd(),Character.toString(g.getStrand()),g.getID(),count,g.getID()));                       
                        count++;
                    }
                }          
            }
        }
    }
}