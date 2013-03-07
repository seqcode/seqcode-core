package edu.psu.compbio.seqcode.projects.shaun;

import java.util.*;
import java.io.File;
import java.sql.SQLException;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExptHandler;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.deepseq.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromosomeGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
/**
 Simple class to count the numbers of reads in an experiment/experiments
 */
public class ChipSeqCountReads {
	
	private Genome gen;
	private double genomeLen=0;
	private ArrayList<SeqExptHandler> IPhandles;
	private double iphittot;
		
	/* command-line driver */
	public static void main(String[] args) throws SQLException, NotFoundException {
		boolean  metaPeak=false;
        Pair<Organism,Genome> pair = Args.parseGenome(args);
        List<SeqLocator> expts = Args.parseChipSeq(args,"expt");
        if (expts.size() == 0) {
            System.err.println("Usage:\n " +
                               "ChipSeqCountReads " +
                               "--species <organism name;genome version> "+
                               "--expt <solexa expt> " );
            return;
        }
        
        //Initialize the peak finder
        ChipSeqCountReads counter = new ChipSeqCountReads(pair.cdr(), expts);
                		
	}
	
	/* Constructor requires a species, genome version and lists of ip and background solexa experiments */
	public ChipSeqCountReads(Genome gen, List<SeqLocator> ips) throws SQLException, NotFoundException {
        this.gen = gen;
      //  System.out.println("Initializing the ChIP-Seq Peak Finder");
        iphittot=0; 
        IPhandles = new ArrayList<SeqExptHandler>();
       // try {
			genomeLen = gen.getGenomeLength(); 
			
			DeepSeqExpt signal = new DeepSeqExpt(gen, ips, "db", 32);
            //IPhandles.add(curr);
			iphittot = signal.getHitCount();
			
		/*	Iterator<Region> testRegions=new ChromosomeGenerator().execute(gen);
			while (testRegions.hasNext()) {
				Region r = testRegions.next();
				Region x = new Region(gen, r.getChrom(),r.getStart(), r.getEnd()+1000000);
				for(StrandedRegion s : signal.loadHits(x)){
					if(s.getFivePrime()>=gen.getChromLength(s.getChrom())){
						System.out.println(s.getLocationString());
					}
				}
			}
		*/
			System.out.print(String.format("%.0f reads loaded\n", iphittot));
        //}
	}
}
