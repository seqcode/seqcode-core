package org.seqcode.genome.tools;

import org.seqcode.genome.Species;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.NotFoundException;

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
            Species o = new Species(species);
            System.err.println("Species " + species + " already exists");
            System.exit(2);
        } catch (NotFoundException e) {
            Species.insertSpecies(species);
        }
    }

}