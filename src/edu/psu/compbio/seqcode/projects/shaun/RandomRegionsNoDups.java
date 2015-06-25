package edu.psu.compbio.seqcode.projects.shaun;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.RepeatMaskedRegion;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RepeatMaskedGenerator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class RandomRegionsNoDups {

	private int numSamples = 1000;
	private int validSamples=0;
	private int sampleSize=200;
	private Genome gen;
	private RepeatMaskedGenerator repMask;
	private double genomeSize=0;
	private long [] chromoSize;
	private String [] chromoNames;
	private int numChroms=0;
	private ArrayList<Region> regList = new ArrayList<Region>();
	private Random rand = new Random();
	private SequenceGenerator seqgen = new SequenceGenerator();
	private double repPropLimit=0.5;
	private boolean screenRepeats=false;
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
		
		if(!ap.hasKey("species") || !ap.hasKey("genome")){
			System.out.println("RandomRegionsNoReps Usage:\n" +
					"--species <species name> " +
					"--genome <genome version>\n" +
					"--len <length of each sequence>\n" +
					"--num <number of sequences>\n" +
					"--seqout <output file name>\n" +
					"--regout <output file name>\n" +
					"--peakout <output file name>\n" +
					"--screenrepeats\n");
		}
		String species = ap.getKeyValue("species");
	    String genome = ap.getKeyValue("genome");
	    String seqFile = ap.hasKey("seqout") ? ap.getKeyValue("seqout") : null;
	    String regFile = ap.hasKey("regout") ? ap.getKeyValue("regout") : null;
	    String peakFile = ap.hasKey("peakout") ? ap.getKeyValue("peakout") : null;
	    boolean sr = Args.parseFlags(args).contains("screenRepeats");
        
        try{
        	Species org = Species.getSpecies(species);
        	Genome g = new Genome(org, genome);
        	
        	RandomRegionsNoDups rrnd = new RandomRegionsNoDups(g);
        	rrnd.setScreenRepeats(sr);
        	rrnd.setLen(ap.hasKey("len") ? new Integer(ap.getKeyValue("len")).intValue() : 200);
        	rrnd.setNum(ap.hasKey("num") ? new Integer(ap.getKeyValue("num")).intValue() : 10);
        	rrnd.execute();
        	
        	if(seqFile!=null)
        		rrnd.printSeqsToFile(seqFile);
        	if(regFile!=null)
        		rrnd.printRegionsToFile(regFile);
        	if(peakFile!=null)
        		rrnd.printPeaksToFile(peakFile);
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public RandomRegionsNoDups(Genome g){
		gen = g;
		repMask = new RepeatMaskedGenerator(gen);
	}
	//Mutators
	public void setNum(int n){numSamples=n;}
	public void setLen(int l){sampleSize=l;}
	public void setScreenRepeats(boolean s){screenRepeats=s;}
	
	public List<Region> execute(){
		//First see how big the genome is:
		chromoSize = new long[gen.getChromList().size()];
		chromoNames = new String[gen.getChromList().size()];
		Iterator<NamedRegion> chroms = new ChromRegionIterator(gen);
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			genomeSize += (double)currentChrom.getWidth();
			chromoSize[numChroms]=currentChrom.getWidth();
			chromoNames[numChroms]=currentChrom.getChrom();
			//System.out.println(chromoNames[numChroms]+"\t"+chromoSize[numChroms]);
			numChroms++;				
		}//System.out.println(genomeSize);
		
		//Now, iteratively generate random positions and check if they are valid 
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
						potential = new Region(gen, chromoNames[c], (int)(randPos-total), (int)(randPos+sampleSize-total-1));
						boolean regionOK = true;
						
						//screen repeats
						if(screenRepeats){
							//is this overlapping a repeat?
							double repLen=0;
							Iterator<RepeatMaskedRegion> repItr = repMask.execute(potential);
							while(repItr.hasNext()){
								RepeatMaskedRegion currRep = repItr.next();
								if(currRep.overlaps(potential)){
									repLen +=(double)currRep.getWidth();
								}
							}if(repLen/(double)potential.getWidth() >repPropLimit)
								regionOK=false;
							
							//Is the sequence free from N's?
							String potSeq=seqgen.execute(potential);
							if(potSeq.indexOf('N')>=0){regionOK=false;}
						}
						//Screen dupicates
						for(Region r : regList){
							if(potential.overlaps(r))
								regionOK=false;
						}
						
						if(regionOK){
							validSamples++;
							regList.add(potential);
							System.out.println(potential.getChrom()+":"+potential.getStart()+"-"+potential.getEnd());
						}
					}
				}total+=chromoSize[c];
			}				
		}
		return(regList);
	}
	//Print the list to a file
	public void printRegionsToFile(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			for(Region currRe : regList){
				fout.write(currRe.getChrom()+":"+currRe.getStart()+"-"+currRe.getEnd()+"\n");
			}
			fout.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	//Print the list as pseudo-peaks to a file
	public void printPeaksToFile(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			for(Region currRe : regList){
				int mid = (currRe.getStart()+currRe.getEnd())/2;
				int off = mid - currRe.getStart();
				fout.write(currRe.getChrom()+":"+currRe.getStart()+"-"+currRe.getEnd()+"\t"+currRe.getWidth()+"\t"+currRe.getChrom()+":"+mid+"\t"+off+"\n");
			}
			fout.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	//Print the sequences to a file
	public void printSeqsToFile(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			for(Region currRe : regList){
				String seq=seqgen.execute(currRe);
				String name = String.format(">%s:%d-%d\n", currRe.getChrom(),currRe.getStart(),currRe.getEnd());
				fout.write(name+seq+"\n");
			}
			fout.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
