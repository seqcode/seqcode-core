package org.seqcode.projects.shaun;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.datasets.seqdata.SeqAlignment;
import org.seqcode.gse.datasets.seqdata.SeqDataLoader;
import org.seqcode.gse.datasets.seqdata.SeqLocator;
import org.seqcode.gse.projects.readdb.Client;
import org.seqcode.gse.projects.readdb.ClientException;
import org.seqcode.gse.projects.readdb.PairedHit;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.NotFoundException;
import org.seqcode.gse.utils.Pair;
import org.seqcode.gse.utils.RealValuedHistogram;
import org.seqcode.projects.shaun.teqseq.core.GenomeLoader;
import org.seqcode.projects.shaun.teqseq.geneloaders.ARegion;
import org.seqcode.projects.shaun.teqseq.geneloaders.AnnotationFilter;
import org.seqcode.projects.shaun.teqseq.geneloaders.AnnotationLoader;
import org.seqcode.projects.shaun.teqseq.geneloaders.GTFAnnotationLoader;


public class PairedReadInserts {

	Genome gen;
	private GenomeLoader genLoader;
	private Client client;
    private Set<SeqAlignment> alignments;
    private Set<String> ids;
    private AnnotationLoader annotLoader;
    private AnnotationFilter annotFilter= new AnnotationFilter();
    private int minUTR=5000;
    private int maxFrag=1000;
    
	public static void main(String[] args) throws SQLException, NotFoundException {
	    Pair<Species,Genome> pair = Args.parseGenome(args);
        List<SeqLocator> expts = Args.parseSeqExpt(args,"expt");
        String gtfFile = Args.parseString(args, "gtf", null);
        int minUTR = Args.parseInteger(args, "minutr", 5000);
        int maxFrag = Args.parseInteger(args, "maxfrag", 1000);
        if (expts.size() == 0 || pair==null || gtfFile==null) {
            System.err.println("Usage:\n " +
                               "PairedReadInserts " +
                               "--species <organism name;genome version> \n" +
                               "--gtf <annotation file>\n" +
                               "--minutr <min 3'UTR length>\n"+
                               "--maxfrag <max frag length length>\n"+
                               "--expt <solexa expt> " );
            return;
        }
        
        PairedReadInserts pid = new PairedReadInserts(pair.cdr(), expts, gtfFile, minUTR, maxFrag);
        pid.execute();
	}
	
	public PairedReadInserts(Genome gen, List<SeqLocator> expts, String gtfFile, int utr_min, int max_frag){
		try {	
			this.gen=gen;
			genLoader = new GenomeLoader(gen);
			this.client = new Client();
			this.alignments = new HashSet<SeqAlignment>();
			SeqDataLoader loader = new SeqDataLoader(); 
			for(SeqLocator csl : expts){
				this.alignments.addAll(loader.loadAlignments(csl, gen));
			}
			this.minUTR = utr_min;
			this.maxFrag = max_frag;
			annotLoader = new GTFAnnotationLoader(new File(gtfFile), genLoader);
			
			ids = new HashSet<String>();
	        for (SeqAlignment a : alignments) {
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
			RealValuedHistogram histo = new RealValuedHistogram(0, maxFrag, maxFrag/10);
			ArrayList<Integer> dists = new ArrayList<Integer>();
			double sum=0, count=0;
			double mean, stdDev;
			
			System.err.println("Parsing annotation");
			Collection<ARegion> testRegions = annotFilter.getSingleIsoform3UTRs(annotLoader.loadGenes(), minUTR);
			System.err.println("Analyzing paired read distribution over "+testRegions.size()+" single isoform 3'UTRs.");

			for(ARegion a : testRegions){
				Region currReg = a.getCoords();
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
							if(insert<maxFrag){
								histo.addValue(insert);
								dists.add(insert);
								sum += insert;
								count++;
							}
                        }
					}
				}
			}
			
			//Histogram of inserts
			System.out.println("Bin\tCount");
			histo.printContents();
			
			//Mean & std dev
			mean = sum/count;
			double v = 0;
			for(Integer i : dists)
				v += (i-mean)*(i-mean);
			stdDev = Math.sqrt(v/(double)dists.size());
			
			System.out.println("\nMean:\t"+mean+"\nStdDev:\t"+stdDev);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClientException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
