package edu.psu.compbio.seqcode.projects.shaun;

import java.io.IOException;
import java.sql.SQLException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqLoader;
import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromosomeGenerator;
import edu.psu.compbio.seqcode.gse.projects.readdb.Client;
import edu.psu.compbio.seqcode.gse.projects.readdb.ClientException;
import edu.psu.compbio.seqcode.gse.projects.readdb.PairedHit;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.RealValuedHistogram;

public class PairInsertDistribution {

	Genome gen;
	private Client client;
    private Set<ChipSeqAlignment> alignments;
    private Set<String> ids;
    
	public static void main(String[] args) throws SQLException, NotFoundException {
	    Pair<Organism,Genome> pair = Args.parseGenome(args);
        List<ChipSeqLocator> expts = Args.parseChipSeq(args,"expt");
        if (expts.size() == 0) {
            System.err.println("Usage:\n " +
                               "PairInsertDistribution " +
                               "--species <organism name;genome version> "+
                               "--expt <solexa expt> " );
            return;
        }
        
        PairInsertDistribution pid = new PairInsertDistribution(pair.cdr(), expts);
        pid.execute();
	}
	
	public PairInsertDistribution(Genome gen, List<ChipSeqLocator> expts){
		try {	
			this.gen=gen;
			this.client = new Client();
			this.alignments = new HashSet<ChipSeqAlignment>();
			ChipSeqLoader loader = new ChipSeqLoader(); 
			for(ChipSeqLocator csl : expts){
				this.alignments.addAll(csl.loadAlignments(loader, gen));
			}
			
			ids = new HashSet<String>();
	        for (ChipSeqAlignment a : alignments) {
	            ids.add(Integer.toString(a.getDBID()));
	        }
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClientException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void execute(){
		try{
			RealValuedHistogram histo = new RealValuedHistogram(0, 2000, 200);
			Iterator<Region> testRegions = new ChromosomeGenerator().execute(gen);
			while (testRegions.hasNext()) {
				Region currReg = testRegions.next();
				if(!currReg.getChrom().endsWith("_random")){
					System.err.println(currReg.getLocationString());
					for (String alignid : ids) {
						List<PairedHit> r = client.getPairedHits(alignid,
							currReg.getGenome().getChromID(currReg.getChrom()),
						    true,
						    currReg.getStart(),
						    currReg.getEnd(),
						    null,
						    null);
						
						for (PairedHit h : r) {
							if (h.leftPos > h.rightPos)
	                            h.flipSides();
	                        if(h.leftStrand && !h.rightStrand){
								int insert = h.rightPos - h.leftPos;
								histo.addValue(insert);
	                        }
						}
					}
				}
			}
			System.out.println("Bin\tCount");
			histo.printContents();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClientException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
