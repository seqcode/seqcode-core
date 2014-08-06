package edu.psu.compbio.seqcode.projects.shaun.teqseq.core;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.NoSuchElementException;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.IllegalSymbolException;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;

/**
 * GenomeLoader: Load a Genome from one of:
 * <ul>
 * <li>A FASTA file</li>
 * <li>The Gifford Lab DB</li>
 * <li>A genome index (name/length pairs)</li>
 * </ul>
 * Sequence is available from FASTA or GiffordDB. However, be aware that the GiffordDB chromosome names do not contain "chr" 
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class GenomeLoader {

	protected Genome gen=null;
	protected boolean seqAvailable=false;
	protected boolean giffordDB=false;
	protected boolean fastaFiles=false; 
	protected Map<String,Integer> chrLengths = new HashMap<String,Integer>();
	protected SequenceDB seqDB; 			 //Used by FASTA loading
	protected BufferedInputStream faStream;  //Used by FASTA loading
	protected SequenceGenerator seqgen=null; //Used by Gifford DB loading
	protected Map<String, String> chromNameTranslator = new HashMap<String,String>();
	
	/**
	 * Constructor: Load a Genome from a FASTA file or a genome index file. 
	 * Uses the BioJava methods for FASTA loading. 
	 * @param inFile File, either a FASTA file of a tab-delim pairing of chr names and lengths
	 * @param isFasta boolean, true if the file is a FASTA
	 */
	public GenomeLoader(File inFile, boolean isFasta){
		if(isFasta){
			//Load a FASTA file containing the genome
			try {
				System.err.println("Loading genome sequence");
				faStream = new BufferedInputStream(new FileInputStream(inFile));
				Alphabet alpha = AlphabetManager.alphabetForName("DNA");
				seqDB = SeqIOTools.readFasta(faStream, alpha);
				seqAvailable=true;
				fastaFiles=true;
				SequenceIterator seqIter = seqDB.sequenceIterator();
				while(seqIter.hasNext()){
					Sequence curr = seqIter.nextSequence();
					chromNameTranslator.put(curr.getName().replaceFirst("^chr", ""), curr.getName());
					chrLengths.put(curr.getName().replaceFirst("^chr", ""), curr.length());
				}
				gen = new Genome("myGenome", chrLengths);
			}catch (BioException ex) {
				//not in fasta format or wrong alphabet
				ex.printStackTrace();
			}catch (NoSuchElementException ex) {
				//no fasta sequences in the file
				ex.printStackTrace();
			}catch (FileNotFoundException ex) {
				//problem reading file
				ex.printStackTrace();
			}
		}else{
			//Assume that this is a tab delimited file of chr names & their lengths
			seqAvailable=false;
			BufferedReader reader;
			try {
				reader = new BufferedReader(new FileReader(inFile));
				String line;
				while ((line = reader.readLine()) != null) {
					line = line.trim();
					String[] words = line.split("\\s+");
					if(words.length>=2){
						String chr = words[0].replaceFirst("^chr", "");
						chromNameTranslator.put(chr, words[0]);
						Integer len = Integer.parseInt(words[1]);
						chrLengths.put(chr, len);
					}	
				}		
				gen = new Genome("myGenome", chrLengths);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (NumberFormatException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Constructor: Load a Genome from an initialized Genome that connects to the Gifford DB. 
	 * @param g Genome
	 */
	public GenomeLoader(Genome g){
		if(g.getDBID()!=-1){
			gen = g; 
			seqgen = new SequenceGenerator();
			chrLengths = gen.getChromLengthMap();
			seqAvailable=true;
			giffordDB=true;
			for(String chr : gen.getChromList())
				chromNameTranslator.put(chr, chr);
		}else{
			//TODO: Throw an exception
			System.err.println("Genome not connected to Gifford Lab DB");
		}
	}
	
	/**
	 * Constructor: Load a Genome from a Collection of ReadLoaders.
	 * When the ReadLoader is a SAMReadLoader, this will have the effect of loading the chromosome names and lengths from the sequence dictionary or will estimate from the read locations.
	 * @param readLoaders A Collection of ReadLoaders
	 */
	public GenomeLoader(Collection<ReadLoader> readLoaders){
		seqAvailable=false;
		Map<String,Integer> tmpChrLengths = new HashMap<String,Integer>();
		for(ReadLoader rl : readLoaders){
			Genome currEstGen = rl.estimateGenome();
			Map<String, Integer> currMap = currEstGen.getChromLengthMap();
			for(String chr : currMap.keySet()){
				if(tmpChrLengths.containsKey(chr)){
					if(tmpChrLengths.get(chr)<currMap.get(chr)){
						tmpChrLengths.put(chr, currMap.get(chr));
					}
				}else{
					tmpChrLengths.put(chr, currMap.get(chr));
				}
			}
		}
		for(String chr : tmpChrLengths.keySet()){
			chrLengths.put(chr.replaceFirst("^chr", ""), tmpChrLengths.get(chr));
			chromNameTranslator.put(chr.replaceFirst("^chr", ""), chr);
		}
		
		for(ReadLoader rl : readLoaders)
			rl.setNameTranslator(chromNameTranslator);
	}
	
	
	//Accessors
	public Genome getGenome(){return gen;}
	public Map<String,Integer> getChrMap(){return chrLengths;}
	public boolean isGiffordDB(){return(giffordDB);}
	public boolean isFastaLoaded(){return(fastaFiles);}
	public boolean seqIsAvailable(){return(seqAvailable);}
	public Map<String,String> getNameTranslator(){return chromNameTranslator;}
	
	/**
	 * Get a sequence corresponding to a region
	 * @param r Region
	 * @return String
	 */
	public String getSequence(Region r){
		String seq ="";
		if(seqAvailable){
			if(fastaFiles){
				try {
					//Hack here to keep chromosome names consistent
					Sequence s = seqDB.getSequence("chr"+r.getChrom());
					seq = s.subStr(r.getStart(), r.getEnd());
				} catch (IllegalIDException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (BioException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}else if(giffordDB){
				seq = seqgen.execute(r);
			}
		}else{
			//TODO: Throw an exception
			System.err.println("No sequence available!");
		}
		return seq;
	}
	
	/**
	 * Get a sequence corresponding to a region
	 * @param c String (chromosome)
	 * @param start int
	 * @param end int
	 * @return String
	 */
	public String getSequence(String c, int start, int end){
		String seq ="";
		if(seqAvailable){
			if(fastaFiles){
				try {
					//Hack here to keep chromosome names consistent
					Sequence s = seqDB.getSequence("chr"+c);
					seq = s.subStr(start, end);
				} catch (IllegalIDException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (BioException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}else if(giffordDB){
				seq = seqgen.execute(new Region(gen, c, start, end));
			}
		}else{
			//TODO: Throw an exception
			System.err.println("No sequence available!");
		}
		return seq;
	}
	
	/**
	 * Get a sequence corresponding to an entire chromosome
	 * @param chrom String corresponding to chromosome name
	 * @return Sequence
	 */
	public Sequence getSequence(String chrom){
		Sequence seq=null;
		if(seqAvailable){
			if(fastaFiles){
				try {
					//Hack here to keep chromosome names consistent
					seq = seqDB.getSequence("chr"+chrom);
					
				} catch (IllegalIDException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (BioException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}else if(giffordDB){
				try {
					seq = DNATools.createDNASequence(seqgen.execute(new Region(gen, chrom, 1, gen.getChromLength(chrom))), "dna");
				} catch (IllegalSymbolException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}else{
			//TODO: Throw an exception
			System.err.println("No sequence available!");
		}
		return seq;
	}
	
	/**
	 * Main method to test functions
	 * @param args Command line args
	 */
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("fa")) { 
            System.err.println("Usage:\n " +
                               "  --fa <file>\n" +
                               "");
            return;
        }
        
		String fastaFile = Args.parseString(args, "fa", null);
		GenomeLoader loader = new GenomeLoader(new File(fastaFile), true);
		
		//Test sequence query
		String testSeq = loader.getSequence("chrY", 1000, 2000);
		System.out.println(testSeq);
		//Test chromosome lengths
		Map<String, Integer> cMap = loader.getChrMap();
		for(String c : cMap.keySet()){
			System.out.println(c+"\t"+cMap.get(c));
		}
	}
}
