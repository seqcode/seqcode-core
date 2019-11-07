package org.seqcode.projects.galaxyexo;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

/**
 * Utility to access peak enrichment at a set of genomic regions.  Statistical significance
 * of peak enrichment over the set of sites is accessed using Poisson model.  The null 
 * distribution comes from randomly placed peaks throughout a genome.  
 * 
 * Input:
 * 		- Genome
 * 		- Peak locations in gff/points/bed
 *		- A set of genomic regions to test peak enrichment
 * Output:
 * 		- A text file indicating the number of overlap and Poisson p-value.
 * 
 * @author naomi yamada
 */

public class PointEnrichmentTester {
	protected GenomeConfig gconfig;	
	protected String outbase;
	protected List<Point> gff;
	protected List<Region> regions; // region that you want to test significance in overlap
	protected int ext; //distance to expand so that I don't double count points
	protected int numItr = 1000;
	protected Poisson poisson;
	protected boolean printRandOverlap = false; // flag to print number of random overlap
	protected String nulldistrib = "";
	
	public PointEnrichmentTester(String base, GenomeConfig gcon,List<Point> g, List<Region> r, int numTest){
		outbase=base;
		gconfig=gcon;
		gff=g;
		regions=r;	
		poisson = new Poisson(1, new DRand());
		numItr = numTest;
	}
	
	// set pseudo counts
	public void setExpansion(int e){ext=e;}
	public void printRandOverlap(){printRandOverlap=true;}
	
	public void execute() throws FileNotFoundException{
		
		int totalRegionSize = 0;
		for (Region reg : regions){totalRegionSize+=reg.getWidth();}	// total size of regions
		// expand regions for 20bp so that I don't double counts the peaks nearby
		List<Region> expandedGff = new ArrayList<Region>();
		for (Point point : gff){
			expandedGff.add(point.expand(ext));
		}
		List <Region> mergedGff = Region.mergeRegions(expandedGff);	
		int totalOverlap = 0; // number of total overlap between two regions
		for (Region gff : mergedGff){
			for (Region reg : regions){
				if (gff.overlaps(reg))
					totalOverlap ++;
			}
		}
		
		File outFile = new File(outbase+File.separator+"point_enrichment.txt");
		outFile.getParentFile().mkdirs();
		PrintWriter writer = new PrintWriter(outFile);		
		writer.println("total number of non-overlapping peaks : "+mergedGff.size());
		writer.println("number of overlap between peaks and regions with size "+totalRegionSize+" : "+totalOverlap);	
		
		double maxPval = 0;
		// Determined p-val based on Poisson distributions for numItr times and return the max p-val
		for (int i=0 ; i < numItr ; i++){
			double pValuePoisson =1;
			// produce random hits through genome
			List<Region> randomRegions = randomRegionPick(gconfig.getGenome(), null, mergedGff.size(),1);
			int numRandOverlaps=0;
			for (Region randRegion : randomRegions){
				for (Region reg : regions){
					if (randRegion.overlaps(reg)){
						numRandOverlaps++;
						break;
					}}}
			
			if (i ==1){writer.println("number of overlap with random regions for iteration one : "+numRandOverlaps);}
			nulldistrib+=numRandOverlaps+",";
			if (numRandOverlaps >totalOverlap){
				pValuePoisson=1;
			}else{
				poisson.setMean(numRandOverlaps);
				int cA = (int)Math.ceil(totalOverlap);
				pValuePoisson = 1 - poisson.cdf(cA) + poisson.pdf(cA);           	
			}
			if (pValuePoisson >maxPval) {maxPval = pValuePoisson;}		
		}
		writer.println("Poisson p-val : "+maxPval);
		
		if (printRandOverlap)
			writer.println(nulldistrib);
		
		writer.close();
	}
	
	public static List<Region> randomRegionPick(Genome gen, List<Region> blackList, int numSamples, int sampleSize){
		List<Region> regs = new ArrayList<Region>();
		Random rand = new Random();
		int validSamples=0;
		
		//First see how big the genome is:
		int numChroms=0;
		long genomeSize=0;
		long [] chromoSize = new long[gen.getChromList().size()];
		String [] chromoNames = new String[gen.getChromList().size()];
		Iterator<NamedRegion> chroms = new ChromRegionIterator(gen);
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			genomeSize += (double)currentChrom.getWidth();
			chromoSize[numChroms]=currentChrom.getWidth();
			chromoNames[numChroms]=currentChrom.getChrom();
			numChroms++;				
		}

		//Now, iteratively generate random positions and check if they are valid and not overlapping repeats. 
		while(validSamples<numSamples){
			Region potential;				
			long randPos = (long)(1+(rand.nextDouble()*genomeSize));
			//find the chr
			boolean found=false;
			long total=0;
			for(int c=0; c<numChroms && !found; c++){
				if(randPos<total+chromoSize[c]){
					found=true;
					if(randPos+sampleSize<total+chromoSize[c]){
						potential = new Region(gen, chromoNames[c], (int)(randPos-total), (int)(randPos+sampleSize-total));
						
						//is this region in the blacklist? 
						boolean valid=true;
						if(blackList!=null){
							for(Region r : blackList){
								if(potential.overlaps(r)){valid=false;}
							}
						}
						if(valid){
							validSamples++;
							regs.add(potential);
						}
					}
				}total+=chromoSize[c];
			}
		}
		return(regs);
	}
	
	public static void main(String[] args) throws FileNotFoundException{
		
		ArgParser ap = new ArgParser(args);		
		GenomeConfig gconf = new GenomeConfig(args);
		Genome genome = gconf.getGenome();
		
		if (!ap.hasKey("gff") && !ap.hasKey("points") && !ap.hasKey("bed")){
            System.err.println("Usage:\n " +
                    "PointEnrichmentTester\n " +
                    "--geninfo <genome info file> \n " +
                    "--gff <gff containing site coordinates> OR --points <points containing site coordinates> OR --bed <bed containing site coordinates>\n " +
                    "--region <region of the genome for enrichment test> OR --regbed <region in bed format for enrichment test>\n " +
                    "\nOPTIONS:\n " +
                    "--out <output directory (default = working directory)> \n " +
                    "--pseudo <pseudocounts to suppress telomere enrichment (default=0) > \n " +
                    "--ext <window size to merge gff points to prevent event double counts (default=20) > \n " +
                    "--print <flag to print number of random overlaps with region> \n " +
                    "--numItr <number of test conducted> \n " +
                    "");
			System.exit(0);
		}
		
		// get peak positions from file
		List<Point> points = new ArrayList<Point>();
		if (ap.hasKey("gff")){
			points = RegionFileUtilities.loadPointsFromGFFFile(ap.getKeyValue("gff"),genome);
		}else if (ap.hasKey("points")){
			points = RegionFileUtilities.loadPointsFromFile(ap.getKeyValue("points"), genome);
		}else if (ap.hasKey("bed")){
			List<Region> peakRegs = RegionFileUtilities.loadRegionsFromBEDFile(genome, ap.getKeyValue("bed"),-1);
			for (Region p : peakRegs){points.add(p.getMidpoint());} // convert bed to points
		}else{
			System.err.println("peak files have zero hits.");
			System.exit(0);
		}
		
		// get regions from file
		List<Region> reg= null;
		if (ap.hasKey("region")){
			reg = RegionFileUtilities.loadRegionsFromFile(ap.getKeyValue("region"),gconf.getGenome(),-1);
		}else if (ap.hasKey("regbed")){
			reg = RegionFileUtilities.loadRegionsFromBEDFile(genome, ap.getKeyValue("regbed"),-1);
		}else{
			System.err.println("region files have zero hits.");
			System.exit(0);
		}
		
		int expand = Args.parseInteger(args,"ext", 20);
		int itr = Args.parseInteger(args,"numItr", 1000);
		// Get outdir and outbase and make them;
		String outbase = Args.parseString(args, "out", System.getProperty("user.dir"));
			
		PointEnrichmentTester tester = new PointEnrichmentTester(outbase,gconf,points,reg,itr);
		
		tester.setExpansion(expand);
		if (ap.hasKey("print")){tester.printRandOverlap();}		
		tester.execute();
	}	
}
