package edu.psu.compbio.seqcode.projects.shaun.teqseq.exptprops;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.GenomeLoader;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.GenomeSegmenter;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.ReadLoader;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.SAMReadLoader;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.AGene;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.GTFAnnotationLoader;

public class ExperimentPropertyAnalyzer {

	private GenomeLoader genLoader=null;
	private ArrayList<ReadLoader> experiments;
	Collection<AGene> genes;
	private int seqBiasK = 7;
	private Map<String, List<Region>> subChrRegions;
	
	public ExperimentPropertyAnalyzer(GenomeLoader gl, ArrayList<ReadLoader> expts, Collection<AGene> genes, Map<String, List<Region>> subChrRegions){
		genLoader = gl;
		experiments = expts;
		this.genes = genes;
		this.subChrRegions = subChrRegions;
	}
	
	/**
	 * Estimate the following experiment properties:
	 *   1. Sequence read start biases
	 *TODO   2. Histogram of read starts over long, single-isoform 5' and 3' UTRs
	 *TODO   3. Known single-isoform gene read counts
	 *TODO   4. Known single-isoform exon read counts
	 *TODO   5. Contiguous locus boundaries
	 */
	public void execute(){
		
		//First estimate sequencing biases
		HashMap<String, SeqBiasModel> biasModels = new HashMap<String, SeqBiasModel>();
		for(ReadLoader currLoader : experiments){
			
			String currName = currLoader.getSourceName();
			currName = currName.replaceAll(".bam", "");
			currName = currName.replaceAll(".sam", "");
			String outName = currName+"."+seqBiasK+"mer.seqbias";
			
			SeqBiasModel model = new ObsExpSeqBiasModel(genLoader, seqBiasK, subChrRegions);
			biasModels.put(currName, model);
			model.execute(currLoader, genes);
			model.save(outName);
			model.print();
		}
	}
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("gtf")) { 
            System.err.println("Usage:\n " +
                               "  --species <species;genome>\n" +
                               "  --gtf <GTF file>\n" +
                               "  --sam <SAM file of reads>"+
                               "");
            return;
        }
        try {
        	GenomeLoader gLoad=null;
        	if(ap.hasKey("species")){
				Pair<Species, Genome> pair = Args.parseGenome(args);
				Genome currgen = pair.cdr();
				gLoad = new GenomeLoader(currgen);
        	}else{
        		String faFile = ap.getKeyValue("fa");
        		gLoad = new GenomeLoader(new File(faFile), true);
        	}
			String gtfFile = ap.getKeyValue("gtf");		
			List<File> samFiles = Args.parseFileHandles(args, "sam");
			
			//Load genes
			GTFAnnotationLoader reader = new GTFAnnotationLoader(new File(gtfFile), gLoad);
			List<AGene> geneSet = reader.loadGenes();
			
			GenomeSegmenter segmenter = new GenomeSegmenter(gLoad, 5000000, 100000);
			Map<String, List<Region>> subChrRegions = segmenter.segmentWithGenes(geneSet);

			//Load experiments
			ArrayList<ReadLoader> expts = new ArrayList<ReadLoader>();
			for(File f : samFiles)
				expts.add(new SAMReadLoader(f, "", gLoad.getGenome(), gLoad.getNameTranslator()));
			
			ExperimentPropertyAnalyzer analyzer = new ExperimentPropertyAnalyzer(gLoad, expts, geneSet, subChrRegions);
			
			
						
        } catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
