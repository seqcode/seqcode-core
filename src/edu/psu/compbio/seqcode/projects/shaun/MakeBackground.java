package edu.psu.compbio.seqcode.projects.shaun;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.motifs.BackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.CountsBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.FrequencyBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.species.Gene;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ExpanderIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RefGeneGenerator;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.io.DatasetsGeneralIO;
import edu.psu.compbio.seqcode.gse.utils.io.motifs.BackgroundModelIO;
import edu.psu.compbio.seqcode.gse.utils.io.parsing.FASTAStream;

public class MakeBackground {

	/**
	 * @param args
	 * @throws NotFoundException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws NotFoundException, IOException {
		BackgroundModel backMod;
		CountsBackgroundModel cbm = null;
		
		//hardcoded args for a one-off run
		//args = new String[] {"--species", "Mus musculus;mm5", "--genome", "mm5","--model", "count", "--k", "3", "--out", "testcount3.back", "--regions", "freqTestRegions.txt"};
		
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("genome")||!ap.hasKey("k")) { 
            System.err.println("Usage:\n" +
                               "MakeBackground " +
                               "--model [count|freq|markov] " +
                               "--species <organism name> " +
                               "--genome <genome version> "+
                               "--k <model length> "+
                               "--out <output file name> "+
                               "--seq <optional FASTA file> "+
                               "--regions <optional region coords file> "+
                               "--promoters "+
                               "--geneannot <annotation name: default refGene> "+
                               "--tssup <length of upstream sequence> "+
                               "--tssdown <length of downstream sequence> "+
                               "--stranded "+
                               "--outseq <filename for background sequences>");
            return;
        }
        String species = ap.getKeyValue("species");
        String genome = ap.getKeyValue("genome");
        int modLen = new Integer(ap.getKeyValue("k")).intValue();
        String modType = ap.hasKey("model") ? ap.getKeyValue("model") : "markov";
        String outFile = ap.hasKey("out") ? ap.getKeyValue("out") : "out.back";
        boolean usingRFile = ap.hasKey("regions");
        boolean usingSeqs = ap.hasKey("seq");
        String inSeqFile = ap.hasKey("seq") ? ap.getKeyValue("seq") : null;
        boolean printSeqs = ap.hasKey("outseq");
        String seqFile = ap.hasKey("outseq") ? ap.getKeyValue("outseq") : null;
        boolean stranded = ap.hasKey("stranded");
        boolean usingGenes = ap.hasKey("promoters");
        int tssup = 10000; int tssdown = 5000;
        String annotation = ap.hasKey("geneannot") ? ap.getKeyValue("geneannot") : "refGene";
        if(usingGenes){
        	if(ap.hasKey("tssup")){tssup = new Integer(ap.getKeyValue("tssup")).intValue();}
        	if(ap.hasKey("tssdown")){tssdown = new Integer(ap.getKeyValue("tssdown")).intValue();}
        }
        
        
        if(!modType.equals("count") && !modType.equals("freq") && !modType.equals("markov")){
        	System.err.println(modType+" is an unknown model type");
        	System.exit(1);
        }
               
        if(usingRFile){
        	String rFile = ap.getKeyValue("regions");
        	Genome gen = Organism.findGenome(genome);
        	cbm = CountsBackgroundModel.modelFromRegionList(gen, DatasetsGeneralIO.readRegionsFromFile(gen, rFile));
        }else if(usingGenes){
        	Genome g;
			try {
				g = Organism.findGenome(genome);
				RefGeneGenerator<NamedRegion> geneGen = new RefGeneGenerator<NamedRegion>(g, annotation);
				cbm = CountsBackgroundModel.modelFromRegionList(g, MakeBackground.getTSSRegions(geneGen, tssup, tssdown));
			} catch (NotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}        	
        }else if(usingSeqs){
        	try {
				FASTAStream faFile = new FASTAStream(new File(inSeqFile));
				cbm = CountsBackgroundModel.modelFromFASTAStream(faFile);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        }else{
          cbm = CountsBackgroundModel.modelFromWholeGenome(Organism.findGenome(genome));
        }System.out.println("Counted");
        
        if (!stranded) {
          cbm.degenerateStrands();
        }
        if (modType.equals("freq")) {
          backMod = new FrequencyBackgroundModel(cbm);
        }
        else if (modType.equals("markov")) {
          backMod = new MarkovBackgroundModel(cbm);
        }
        else {
          backMod = cbm;
        }
        		
        
        BackgroundModelIO.printProbsToFile(backMod, outFile);
        System.out.println("Printed to: "+outFile);
	}


  public static List<Region> getTSSRegions(RefGeneGenerator<NamedRegion> geneGen, int up, int down)
      throws NotFoundException {
    ArrayList<Region> regList = new ArrayList<Region>();
    Genome gen = geneGen.getGenome();
    ChromRegionIterator chroms = new ChromRegionIterator(gen);
    Iterator<Gene> genes = new ExpanderIterator<NamedRegion, Gene>(geneGen, chroms);
    while (genes.hasNext()) {
      Gene gene = genes.next();

      int tss = 0, start = 0, stop = 0;
      char strand = gene.getStrand();
      if (strand == '+') {
        tss = gene.getStart();
        start = tss - up;
        stop = tss + down;
      }
      else {
        tss = gene.getEnd();
        start = tss - down;
        stop = tss + up;
      }
      if (start < 1) {
        start = 1;
      }
      if (stop > gen.getChromLength(gene.getChrom())) {
        stop = gen.getChromLength(gene.getChrom());
      }

      regList.add(new Region(gen, gene.getChrom(), start, stop));
    }
    return regList;

  }

}

