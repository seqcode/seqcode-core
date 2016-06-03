package org.seqcode.projects.shaun;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;

import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.NotFoundException;


/**
 * @author Shaun Mahony
 * 
 * This class aims to count the frequency of each k-mer occurring in a genome, and 
 * write the corresponding frequencies for each base in the genome to the DB (i.e. 
 * it writes how often the k-mer beginning at each base occurs in the entire genome).
 * Unique k-mers receive the value 1. Any observed k-mers containing 'N' receive the 
 * value 0.
 * 
 * This program currently uses a hash map to store the counts of each k-mer. This 
 * obviously doesn't scale very well in terms of required memory for larger values 
 * of "k". To allow the program to run on a typical computer, only a subset of all 
 * possible k-mers are stored/counted in each iteration. This choice makes the program 
 * inefficient; the genome is scanned multiple times in order to count all k-mers. 
 * However, at least the hash table now fits into memory. 
 * 
 * The other problem with the current version of this program is that the hash map uses 
 * an integer representation of each k-mer in the lookup. Since I am using 64-bit integers
 * (i.e. long), the maximum value of k allowed is 32. For longer k, there are two choices:
 * 		1. Use more bits for the integers, if such a representation exists. This option
 * 			will increase the memory footprint of the hash map. If you run into memory 
 * 			trouble, try decreasing the value of "maxKmersPerRun".
 * 		2. Re-write the code to follow Alex's original suggestion. This involved building 
 * 			linked lists of the *locations* of l-mers, where l<=k/2, and processing these 
 * 			linked lists later to reconstruct how often each k-mer occurs in the genome. 
 * 			While this is a nice idea, I began writing this class by trying it out and 
 * 			it turned out to have much bigger memory requirements than a simple hash map
 * 			(probably because you have to store two ints, chr and start, for every entry
 * 			on each linked list). 
 * 
 * Typical CountKmerOccurrence usage (on Pol2 which has ~8G RAM):
 * 	java -Xmx6144m org.seqcode.gse.shaun.CountKmerOccurrence --k 26 --species "Mus musculus" --genome "mm8"  
 */

public class CountKmerOccurrence {
	private long maxKmersPerRun=15000000;//50000000;
	private static int B=4;
	private static int bitsInType = 64;
	private static int k=26;
	private static Species org;
	private static Genome gen;
	private static String genomeName;
	private ArrayList <String> chrNameLookup; 
	private boolean strandDegeneracy=true;
	private int chrChunkSize=50000000;
	private long genomeSize=0;
	private long unique=0;
	private int numRuns=1;
	private static SequenceGenerator seqgen = new SequenceGenerator();
	private long [] kPow;
	private Map <Long, Integer> hash;
	private char [] revseq;
		
	 
	public static void main(String[] args) {
		 ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("genome") || !ap.hasKey("k")) { 
            System.err.println("Usage:\n" +
                               "CountKmerOccurrence " +
                               "--k <length of k-mers> " +
                               "--species <organism name> " +
                               "--genome <genome version> "+
            				   "--forwardonly ");
            return;
        }
        
        String species = ap.getKeyValue("species");
        String gname = ap.getKeyValue("genome");
        int userK = Integer.valueOf(ap.getKeyValue("k"));
        boolean userBothStrands = ap.hasKey("forwardonly") ? false : true; 
        
		//CountKmerOccurrence counter = new CountKmerOccurrence(26, "Monkeypox", "MPX-ZAI", true);
		//CountKmerOccurrence counter = new CountKmerOccurrence(26, "Mus musculus", "mm8", true);
		//CountKmerOccurrence counter = new CountKmerOccurrence(26, "Saccharomyces cerevisiae", "sacCer1", true);
        CountKmerOccurrence counter = new CountKmerOccurrence(userK, species, gname, userBothStrands);
		counter.run();
	}
	
	@SuppressWarnings("unchecked")	
	public CountKmerOccurrence(int k, String organism, String genome, boolean bothStrands){
		this.k=k;
		revseq = new char[k];
		if(k>(bitsInType)){
			System.err.println("Sorry: "+bitsInType+" bits are used for indexing k-mers in this program, so only k up to "+(bitsInType/2)+" is allowed.");
			System.exit(1);
		}else{
			genomeName = genome;
			strandDegeneracy = bothStrands;
			kPow = new long [k+1];
			for(int i=0; i<=k; i++){
				kPow[i]=(long)Math.pow(B, i);
			}
			
			try {
				org = Species.getSpecies(organism);
				gen = new Genome(org, genome);						
			} catch (NotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	
	public void run(){
		boolean firstRun=true;
				
		//Decide how many times to run this thing
		//Get the length of the genome 
		Iterator<NamedRegion> chroms = new ChromRegionIterator(gen);
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			genomeSize += currentChrom.getWidth();
		}
		System.out.println("Length of genome: "+genomeSize+" bp");
		numRuns = (int)((genomeSize/maxKmersPerRun)+1);
		long intervalSize = kPow[k]/numRuns;
		
		for(int r=0; r<numRuns; r++){
			hash = new HashMap <Long, Integer> (10000000);
			long intervalStart = r*intervalSize;
			long intervalStop = ((r*intervalSize)+(intervalSize-1))<kPow[k] ? (r*intervalSize)+(intervalSize-1) : kPow[k]-1;
			System.out.println("Run "+(r+1)+" of "+numRuns+" ("+intToSeq(intervalStart, k)+" to "+intToSeq(intervalStop, k)+")");
		
			//Fill the hash map
			chroms = new ChromRegionIterator(gen);
			while (chroms.hasNext()) {
				NamedRegion currentChrom = chroms.next();
				int regionStart = currentChrom.getStart();
				int regionEnd = currentChrom.getEnd();
				System.out.println("Processing chromosome: "+currentChrom.getChrom()+"\t("+currentChrom.getWidth()+" bp)");
				
				for(int x=regionStart; x<regionEnd; x+=chrChunkSize){
					int rstart = (x-(k-1))<regionStart ? regionStart : (x-(k-1));
					int rstop = (rstart+chrChunkSize)>=regionEnd ? regionEnd-1 : (rstart+chrChunkSize);
					Region currR = new Region(gen, currentChrom.getChrom(), rstart, rstop);
					
					String tmpSeq = seqgen.execute(currR);
					String regionSeq = tmpSeq.toUpperCase();
					
					for(int i=0; i<(regionSeq.length()-k+1); i++){
						String currK = regionSeq.substring(i, i+k);
						long id = seqToInt(currK, k);
						
						if(inInterval(id, intervalStart, intervalStop)){
							if(id>=0){
								if(hash.containsKey(id)){
									int count = 1+ hash.get(id);
									hash.put(id, count);
								}else{
									hash.put(id, 1);
								}						
							}
						}
					}
				}
			}
			
			//Pass 2: Scan through genome again and print landscape for each chromosome to file
			System.out.println("Adding counts to database");
			chroms = new ChromRegionIterator(gen);
			while (chroms.hasNext()) {
				NamedRegion currentChrom = chroms.next();
				int regionStart = currentChrom.getStart();
				int regionEnd = currentChrom.getEnd();
				
				for(int x=regionStart; x<regionEnd; x+=chrChunkSize){
					int rstart = (x-(k-1))<regionStart ? regionStart : (x-(k-1));
					int rstop = (rstart+chrChunkSize)>=regionEnd ? regionEnd-1 : (rstart+chrChunkSize);
					Region currR = new Region(gen, currentChrom.getChrom(), rstart, rstop);
					
					String tmpSeq = seqgen.execute(currR);
					String regionSeq = tmpSeq.toUpperCase();
					
					for(int i=0; i<(regionSeq.length()-k+1); i++){
						String currK = regionSeq.substring(i, i+k);
						long id = seqToInt(currK, k);
													
						if(id>=0){
							if(inInterval(id, intervalStart, intervalStop)){
								int cst=i+rstart;
								if(hash.containsKey(id)){
									Integer count = hash.get(id);
									if(count==1){
										unique++;
									}
									//Replace this with database entry: chr:currentChrom start:cst uniqueness:count
									//pw.println(cst+"\t"+currK +"\t"+count.toString());
								}else{
									//Replace this with database entry: chr:currentChrom start:cst uniqueness:0
									//pw.println(cst+"\t"+currK +"\t"+"0");
								}
							}
						}else if(firstRun){
							//Replace this with database entry: chr:currentChrom start:cst uniqueness:count
							//pw.println(cst+"\t"+currK +"\t"+"0");
						}
					}
				}			
			}
			firstRun=false;
		}
		
		//Uniqueness summary
		System.out.println("\n\n"+unique+" bases of "+genomeSize+" total are unique.");
	}
	
	private boolean inInterval(long x, long intStart, long intEnd){
		return((x>=intStart && x<=intEnd) ? true : false);
	}
	
	private long seqToInt(String seq, int currk){
		long bin=0;
		int len = seq.length();
		if(len != currk){System.err.println("Error: k's don't match: "+len);}
		
		for(int i=0; i<len; i++){
			long currInt = base2int(seq.charAt(i));
			if(currInt==-1){
				return(-1);
			}
			bin+=currInt*kPow[len-i-1];
		}
		
		if(strandDegeneracy){
			long revbin=0;
		    for (int i = (currk - 1); i >= 0; i--)
		      revseq[(currk - 1)-i] = complement(seq.charAt(i));
		    
			for(int i=0; i<len; i++){
				long currInt = base2int(revseq[i]);
				if(currInt==-1){
					return(-1);
				}
				revbin+=currInt*kPow[len-i-1];
			}if(revbin<bin){bin=revbin;}
		}
		
		return(bin);
	}
	
	private String intToSeq(long id, int currk){
	    long x = id;
	    StringBuffer currWord = new StringBuffer(currk);
		if(x>=kPow[currk]){
			return("");
		}else{
			for(int i=currk-1; i>=0; i--){
				long m = (long)x/kPow[i];
				currWord.append(int2base((int)m));
				x = x-(m*kPow[i]);
			}
			return(currWord.toString());
		}
	}
	
	private char complement(char c){
		char x='N';
		switch(c){
			case 'A':
				x='T'; break;
			case 'C':
				x='G'; break;
			case 'G':
				x='C'; break;
			case 'T':
				x='A'; break;
			default:
				x='N';
		}
		return(x);
	}
	
	private int base2int(char c){
		int l=-1;
		switch(c){
			case 'A':
				l=0; break;
			case 'C':
				l=1; break;
			case 'G':
				l=2; break;
			case 'T':
				l=3; break;
			default:
				l=-1;
		}
		return(l);
	}
	
	private char int2base(int x){
		char c;
		switch(x){
			case 0:
				c='A'; break;
			case 1:
				c='C'; break;
			case 2:
				c='G'; break;
			case 3:
				c='T'; break;
			default:
				c='N';
		}
		return(c);
	}
}
