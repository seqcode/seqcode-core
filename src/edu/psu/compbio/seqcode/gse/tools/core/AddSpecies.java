package edu.psu.compbio.seqcode.gse.tools.core;

import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

/**
 * AddSpecies --species "mus musculus"
 */

public class AddSpecies {

    public static void main(String args[]) throws Exception {
        String species = Args.parseString(args,"species",null);
        if (species == null) {
            System.err.println("Must supply --species");
            System.exit(1);
        }
        try {
            Organism o = new Organism(species);
            System.err.println("Species " + species + " already exists");
            System.exit(2);
        } catch (NotFoundException e) {
            Organism.insertOrganism(species);
        }
    }

}