package edu.psu.compbio.seqcode.projects.shaun;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExptHandler;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

/* A constrained version of the ChipSeqPeakFinder with some options enforced
 * 
 */
public class CLIPSeqPeakFinder extends ChipSeqPeakFinder{
	
	
	/* command-line driver */
	public static void main(String[] args) throws SQLException, NotFoundException {
		boolean  metaPeak=false;
		Pair<Organism,Genome> pair = Args.parseGenome(args);
		if(pair==null){printError();return;}
		List<SeqLocator> expts = Args.parseSeqExpt(args,"expt");
        List<SeqLocator> back = Args.parseSeqExpt(args,"back");
        if (expts.size() == 0) {
            printError();
            return;
        }
        double rLen = Args.parseDouble(args,"readlen",SeqExptHandler.defaultReadLength);
        double rExt = Args.parseDouble(args,"readextend",0);
        double rShift=0;
        
        //Initialize the peak finder
        CLIPSeqPeakFinder finder = new CLIPSeqPeakFinder(pair.cdr(), expts, back, rLen, rExt, rShift);
        
        //Load options for peak calling 
        finder.setPeakFileName(Args.parseString(args,"outpeak",finder.outPeakName));
        finder.setSeqFileName(Args.parseString(args,"outseq",finder.outSeqName));
        finder.setSeqwin(Args.parseInteger(args,"seqwin",finder.seqwin));
        finder.setBinWidth(Args.parseDouble(args,"binwidth",finder.winWidth));
        finder.setBinStep(Args.parseDouble(args,"binstep",finder.winStep));
        finder.setMinZ(Args.parseDouble(args,"minz",finder.sigThres));
        finder.setPoissonThres(Args.parseDouble(args,"poisson",finder.poissThres));
        finder.setTowerFilter(Args.parseFlags(args).contains("filtertowers"));
        finder.setNeedleFilter(Args.parseFlags(args).contains("filterneedles"));
        finder.setGeneOverlap(Args.parseFlags(args).contains("geneOverlap"));
        finder.setMaxGeneDistance(Args.parseInteger(args,"maxgenedist",finder.maxGeneDistance));
        finder.peakCenter = Args.parseFlags(args).contains("peakcenter");
        Collection<String> genes = Args.parseStrings(args,"genes");
        if (genes.size() > 0) {
            finder.annots = new String[genes.size()];
            int i = 0;
            Iterator<String> iter = genes.iterator();
            while (iter.hasNext()) {
                finder.annots[i++] = iter.next();
            }
        } 
        finder.namedRegions = Args.parseStrings(args,"namedregion");
        finder.namedStrandedRegions = Args.parseStrings(args,"namedstrandedregion");
        finder.namedTypedRegions = Args.parseStrings(args,"namedtypedregion");
        finder.repeatMaskerAnnots = Args.parseFlags(args).contains("repeatmasker");
        finder.setGeneScanning(Args.parseFlags(args).contains("scangenesonly"));
        if(Args.parseArgs(args).contains("dynpoisson")){finder.addPoissonModel(Args.parseInteger(args,"dynpoisson",5000));}
        else{finder.addPoissonModel(5000);}
        
        //Run the peak finder
        ArrayList<ChipSeqPeak> enriched = finder.execute();
        System.out.println("Printing");
		finder.printEnrichedPeaks(enriched);
		finder.printPeakSequences(enriched);
		
		System.out.println("Finished!");		
	}
	public CLIPSeqPeakFinder(Genome gen, List<SeqLocator> ips,
			List<SeqLocator> backs, double len, double ext, double shift) {
		super();
		winWidth=75;
		winStep = 10;
		sigThres = 2.33; //p<0.01?
		poissThres=-9; //10th power for Poission threshold
		geneOverlap=true; // if true, then only return genes that overlap the bound region; don't do the closest gene thing
		addStrandedness=true;
	    maxGeneDistance=50000;
		towerfiltering=false; needlefiltering=false;
		annots = new String[]{
			"sgdGene"
		};
		genePeaksOnly=false;
		
		
		System.out.println("Initializing the Peak Finder");
		this.gen = gen;
        readLength = len;
        readExtension = ext;
        readShift = shift;
        if(readShift > 0)
        	shiftTags=true;
        iphittot=0; backhittot=0;
                
        //Loading experiments
        loadExperiments(ips, backs);
        
        //Initialize Poisson Models
        addPoissonModel(-1);
	}
	
	public static void printError(){
		System.err.println("Usage:\n " +
                "CLIPSeqPeakFinder \n" +
                "--species <organism name;genome version> "+
                "--expt <solexa expt> " +
                "--back <background expt> \n"+
                "--outpeak <output file> "+
                "--outseq <output file> "+
                "--seqwin <sequence length> \n"+
                "--binwidth <width of bins> "+
                "--binstep <offset of bins> "+
                "--minz <minimum z score> "+
                "--poisson <Poisson threshold (10 to the power of)> \n" +
                "--dynpoisson <dynamic Poisson threshold (Xbp) default 5Kbp> \n" +
                "--readlen <length> " +
                "--readextend <extension> \n"+
                "--filtertowers [filter towers] " +
                "--filterneedles [filter needles] " +
                "--printpeakdistrib [prints meta-peak] \n" +
                "--genes <name of gene annotation to examine> " +
                "--scangenesonly <only look for peaks within genes (CLIP-seq specific)> " +
                "--maxgenedist <distance from gene> " +
                "--geneOverlap [only look at genes that overlap a peak] \n" +
                "--shifttags <distance to shift tags (overrides extension)> \n");
	}
}
