package org.seqcode.projects.shaun.teqseq.core;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.NotFoundException;
import org.seqcode.gse.utils.Pair;
import org.seqcode.projects.shaun.teqseq.geneloaders.AGene;
import org.seqcode.projects.shaun.teqseq.geneloaders.GTFAnnotationLoader;


public class ReadLoaderTester {

	//Basics
	protected ExptCollection experiments;
	protected GenomeLoader gLoad;
	protected Collection<AGene> genes;

	public ReadLoaderTester(GenomeLoader gLoad, ExptCollection expts, Collection<AGene> geneSet) {
		experiments = expts;
		this.gLoad = gLoad;
		this.genes = geneSet;
	}
	
	public void execute(){
		/*
		//Test region of interest loading
		while(experiments.hasNextRegions()){
			Pair<Region, HashMap<String, List<AlignHit>>> p = experiments.getAllNextRegionHit(false);
			if(p!=null && p.car()==null){
				Region currRegion = p.car();
				HashMap<String, List<AlignHit>> allHits = p.cdr();
				
				System.out.println("Loaded region: "+currRegion.getLocationString());
				for(String e : allHits.keySet()){
					System.out.println("\t"+e+": "+allHits.get(e).size());
				}
			}
		}
				
		//Test reset 
		experiments.resetHitExtractors();
		*/
		//Test sub-chromosomal region loading
		while(experiments.hasNextRegions()){ //Iterate through regions
			Pair<Region, HashMap<String, List<AlignHit>>> p = experiments.getAllNextSubChrHit(false);
			if(p!=null && p.car()!=null){
				Region currRegion = p.car();
				HashMap<String, List<AlignHit>> allHits = p.cdr();
				
				System.out.println("Loaded sub-chr: "+currRegion.getLocationString());
				for(String e : allHits.keySet()){
					System.out.println("\t"+e+": "+allHits.get(e).size());
				}
			}
		}
	}
	
	/**
	 * Main method for program execution
	 * @param args Command-line arguments
	 */
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
		if((!ap.hasKey("fa") &&!ap.hasKey("fai"))  || !ap.hasKey("gtf")) { 
            System.err.println("Usage:\n" +
                               "  --fa <FASTA file> OR --fai <FASTA index>\n" +
                               "  --gtf <GTF file>\n"+
                               "  --exptNAME <SAM file of reads>\n" +
                               "  --exptname <experiment collection name>\n" +
                               "  --minsubsize <sub-chr size>\n" +
                               "  --minsubspace <sub-chr empty spacing>\n" +
                               "");
            return;
        }
        try {
        	//Process command-line options
        	GenomeLoader gLoad=null;
        	if(ap.hasKey("species")){
				Pair<Species, Genome> pair = Args.parseGenome(args);
				Genome currgen = pair.cdr();
				gLoad = new GenomeLoader(currgen);
        	}else if(ap.hasKey("fa")){
        		String faFile = ap.getKeyValue("fa");
        		gLoad = new GenomeLoader(new File(faFile), true);
        	}else if(ap.hasKey("fai")){
        		String faiFile = ap.getKeyValue("fai");
        		gLoad = new GenomeLoader(new File(faiFile), false);
        	}
			Args.parseFileHandles(args, "sam");
			String gtfFile = ap.getKeyValue("gtf");
			String exptName = Args.parseString(args, "exptname", "RNAseqExpt");
			Integer minSubChrSize =Args.parseInteger(args, "minsubsize", 5000000);
			Integer minSubChrSpace = Args.parseInteger(args, "minsubspace", 500000);
			
			//Load the experiments
			HashMap<String,List<File>> samFiles = new HashMap<String,List<File>>();
			ArrayList<String> conditionNames = new ArrayList<String>();
			Vector<String> exptTags=new Vector<String>();
			for(String s : args)
	        	if(s.contains("--expt"))
	        		if(!exptTags.contains(s)){
	        			exptTags.add(s);
	        			String name = s.replaceFirst("--expt", ""); 
		        		conditionNames.add(name);
	        		}
			if(exptTags.size()==0){
			    System.err.println("Error: No experiments provided.\nUse the --expt option.");
			    System.exit(1);
			}
			for(String name : conditionNames){
				if(!samFiles.containsKey(name)){
        			samFiles.put(name, new ArrayList<File>());
        		}
        		samFiles.get(name).addAll(Args.parseFileHandles(args, "expt"+name));
			}	        	
			
			//Load genes
			GTFAnnotationLoader reader = new GTFAnnotationLoader(new File(gtfFile), gLoad);
			List<AGene> geneSet = reader.loadGenes();
			GenomeSegmenter segmenter = new GenomeSegmenter(gLoad, minSubChrSize, minSubChrSpace);
			Map<String, List<Region>> regionsOfInterest = segmenter.segmentWithGenes(geneSet);
			
			//Initialize experiment collection
			int numReps=0;
			ExptCollection expts = new ExptCollection(exptName, regionsOfInterest);
			HashMap<String, ExptCondition> conds = new HashMap<String, ExptCondition>();
			ArrayList<ReadLoader> loaders = new ArrayList<ReadLoader>();
			for(String name : conditionNames){
				ExptCondition c = new ExptCondition(name);
				for(File sf : samFiles.get(name)){
					ReadLoader rl = new SAMReadLoader(sf, name, gLoad.getGenome(), gLoad.getNameTranslator());
					ExptReplicate rep = new ExptReplicate(sf.getName(), numReps, gLoad, rl, regionsOfInterest, geneSet, gLoad.seqIsAvailable());
					loaders.add(rl);
					c.addReplicate(rep);
					numReps++;
				}
				conds.put(name, c);
				expts.addCondition(c);
			}
			
			ReadLoaderTester tester = new ReadLoaderTester(gLoad, expts, geneSet);
			tester.execute();
			
	    } catch (NotFoundException e) {
			e.printStackTrace();
		}
	}
}
