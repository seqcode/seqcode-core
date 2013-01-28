package edu.psu.compbio.seqcode.gse.tools.proteinprotein;

import java.util.*;
import java.sql.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.proteinprotein.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.*;

public class StringDB {

    private static int species;
    private static edu.psu.compbio.seqcode.gse.datasets.proteinprotein.StringDB db;

    public static void main(String[] args) throws Exception {
        String speciesname = args[0];
        String cmd = args[1];
        db = new edu.psu.compbio.seqcode.gse.datasets.proteinprotein.StringDB();
        species = db.getSpeciesID(speciesname);
        if (cmd.equals("geneid")) {
            getGeneID(args[2]);
        } else if (cmd.equals("genelinks")) {
            getGeneLinks(args[2]);
        } else if (cmd.equals("geneactions")) {
            getGeneActions(args[2]);
        } else if (cmd.equals("geneslike")) {
            getGenesLike(args[2]);
        } else {
            System.err.println("Unknown command " + cmd);
            System.exit(1);
        }


    }
    private static void getGeneID(String name) throws SQLException, NotFoundException {
        String id = db.getGeneID(species, name);
        System.out.println(name + " is " + id);
    }
    private static void getGenesLike(String name)  throws SQLException, NotFoundException {
        List<String> names = db.getGenesLike(species, name);
        for (String n : names) {
            System.out.println(n);
        }

    }
    private static void getGeneLinks(String name) throws SQLException, NotFoundException {
        String id = db.getGeneID(species, name);
        List<Link> links = db.getGeneLinks(species, id);
        for (Link l : links) {            
            List<String> other = db.getGeneName(species,l.geneB,"Ensembl_EntrezGene");
            System.out.println(String.format("%s\t%s\t%f",
                                             name, other.toString(), l.score));
                                             
        }
    }
    private static void getGeneActions(String name) throws SQLException, NotFoundException {
        String id = db.getGeneID(species, name);
        List<Action> actions = db.getGeneActions(species,id);
        for (Action a : actions) {            
            List<String> other = db.getGeneName(species,a.geneB,"Ensembl_EntrezGene");
            System.out.println(String.format("%s\t%s\t%f\t%s\t%s",
                                             name, other.toString(), a.score,a.mode,a.action));
        }
    }
}
