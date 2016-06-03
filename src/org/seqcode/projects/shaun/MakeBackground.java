package org.seqcode.projects.shaun;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Gene;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.datasets.motifs.BackgroundModel;
import org.seqcode.gse.datasets.motifs.CountsBackgroundModel;
import org.seqcode.gse.datasets.motifs.FrequencyBackgroundModel;
import org.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import org.seqcode.gse.gsebricks.verbs.ExpanderIterator;
import org.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gse.gsebricks.verbs.location.RefGeneGenerator;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.NotFoundException;
import org.seqcode.gse.utils.io.BackgroundModelIO;
import org.seqcode.gse.utils.io.DatasetsGeneralIO;
import org.seqcode.gse.utils.io.parsing.FASTAStream;


public class MakeBackground {

	/**
	 * @param args
	 * @throws NotFoundException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws NotFoundException, IOException {
		BackgroundModel backMod;
		CountsBackgroundModel cbm = null;
				
		ArgParser ap = new ArgParser(args);
		GenomeConfig gConfig = new GenomeConfig(args);
        if(!ap.hasKey("species") ||!ap.hasKey("k")) { 
            System.err.println("Usage:\n" +
                               "MakeBackground " +
                               "--model [count|freq|markov] " +
                               "--species <organism;version> " +
                               "--seq <genome sequence directory> "+
                               "--k <model length> "+
                               "--out <output file name> "+
                               "--fasta <optional FASTA file> "+
                               "--regions <optional region coords file> "+
                               "--promoters "+
                               "--geneannot <annotation name: default refGene> "+
                               "--tssup <length of upstream sequence> "+
                               "--tssdown <length of downstream sequence> "+
                               "--stranded "+
                               "--outseq <filename for background sequences>");
            return;
        }
        Genome gen = gConfig.getGenome();
        int modLen = new Integer(ap.getKeyValue("k")).intValue();
        String modType = ap.hasKey("model") ? ap.getKeyValue("model") : "markov";
        String outFile = ap.hasKey("out") ? ap.getKeyValue("out") : "out.back";
        boolean usingRFile = ap.hasKey("regions");
        boolean usingSeqs = ap.hasKey("fasta");
        String inSeqFile = ap.hasKey("fasta") ? ap.getKeyValue("fasta") : null;
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
        	cbm = CountsBackgroundModel.modelFromRegionList(gen, DatasetsGeneralIO.readRegionsFromFile(gen, rFile), modLen);
        }else if(usingGenes){
        	RefGeneGenerator<NamedRegion> geneGen = new RefGeneGenerator<NamedRegion>(gen, annotation);
        	cbm = CountsBackgroundModel.modelFromRegionList(gen, MakeBackground.getTSSRegions(geneGen, tssup, tssdown), modLen);
        }else if(usingSeqs){
        	try {
				FASTAStream faFile = new FASTAStream(new File(inSeqFile));
				cbm = CountsBackgroundModel.modelFromFASTAStream(faFile, modLen);
			} catch (IOException e) {
				e.printStackTrace();
			}
        }else{
          cbm = CountsBackgroundModel.modelFromWholeGenome(gen, modLen);
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

