package edu.psu.compbio.seqcode.projects.shaun;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;

public class StrandedSequenceGenerator {

	private Genome gen;
	
	public StrandedSequenceGenerator(Genome g){
		gen = g;
	}
	
	public List<String> generateSequences(List<StrandedRegion> regs){
		return RegionFileUtilities.getSequencesForStrandedRegions(regs,null);
	}
	
	public List<String> generateFASTA(List<StrandedRegion> regs){
		List<String> seqs = RegionFileUtilities.getSequencesForStrandedRegions(regs,null);
		List<String> fasta = new ArrayList<String>();
		for(int i=0; i<regs.size(); i++){
			fasta.add(new String(">seq_"+regs.get(i).toString()+"\n"+seqs.get(i)));
		}
		return fasta;
	}
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") ||(!ap.hasKey("peaks"))) { 
            System.err.println("Usage:\n " +
                               "StrandedSequenceGenerator \n" +
                               " Required: \n" +
                               "  --species <species;version> " +
                               "  --peaks <file containing stranded coordinates> \n" +
                               "  --win <window of sequence around peaks> \n"+
                               "  --out output filename\n" +
                               "");
            return;
        }
        try {
        	Pair<Species, Genome> pair = Args.parseGenome(args);
        	Species currorg = pair.car();
        	Genome currgen = pair.cdr();
        	StrandedSequenceGenerator seqGen = new StrandedSequenceGenerator(currgen);
        	
	        String peaksFile = ap.getKeyValue("peaks");
	        String outName = ap.hasKey("out") ? ap.getKeyValue("out") : "out.seq";
	    	int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():-1;
	    	
	    	List<StrandedRegion> regions = RegionFileUtilities.loadStrandedRegionsFromMotifFile(currgen, peaksFile, win);
	    	
	    	List<String> fasta = seqGen.generateFASTA(regions);
	    	
	    	FileWriter fw = new FileWriter(outName);
	    	for(String s : fasta)
	    		fw.write(s+"\n");
			fw.close();
			System.err.println("FASTA file written to "+outName);
        } catch (IOException e) {
        	e.printStackTrace();
        } catch (NotFoundException e) {
			e.printStackTrace();
		}
	}
}
