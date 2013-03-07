package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.StringTokenizer;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class SeqExptProbLandscapeKmers {

	private static Organism org;
	private static Genome gen;
	private static SeqDataHandler IPhandle;
	private static SeqDataHandler backhandle;
	private static int DEFT_K = 6;
	private static int K;
	private static BufferedReader regionReader;
	private static ArrayList <Region> loadedRegList; 
	private static ArrayList <Region> regList;
	private static int windowSize= 200;
	private static int windowOff= 100;
	private static int [] kPow;
	private static int B=4;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("genome")||!ap.hasKey("expt") || !ap.hasKey("back")|| !ap.hasKey("regions")) { 
            System.err.println("Usage:\n" +
                               "SeqExptProbLandscapeKmers " +
                               "--species <organism name> " +
                               "--genome <genome version> "+
                               "--expt <solexa expt> " +
                               "--back <background expt> "+
                               "--regions <coordinate file> " +
                               "--k <length of K-mers> +" +
                               "--ecdf <expt CDF> " +
            				   "--bcdf <back CDF> ");
              return;
        }
        String species = ap.getKeyValue("species");
        String genome = ap.getKeyValue("genome");
        String exptName = ap.getKeyValue("expt");
        String backName = ap.getKeyValue("back");
        String regionFile = ap.getKeyValue("regions");
        int k = ap.hasKey("k") ? Integer.valueOf(ap.getKeyValue("k")) : DEFT_K;
        String ecdfFile = ap.hasKey("ecdf") ? ap.getKeyValue("ecdf") : "NOFILE";
        String bcdfFile = ap.hasKey("bcdf") ? ap.getKeyValue("bcdf") : "NOFILE";
        
        //Load the ip & back channels and read in the list of regions of interest
        SeqExptProbLandscapeKmers seplk = new SeqExptProbLandscapeKmers(species, genome, exptName, backName, k);
		seplk.loadFile(regionFile);
		if(!ecdfFile.equals("NOFILE")){IPhandle.loadCDFFile(ecdfFile);}
		if(!bcdfFile.equals("NOFILE")){backhandle.loadCDFFile(bcdfFile);}
		
		//Break up the loaded regions into standard-sized chunks
		regList = standardizeRegions(seplk.getLoadedRegions());
		
		//Initialize the sequence generator 
		SequenceGenerator seqgen = new SequenceGenerator();
		
		//Get the likelihood ratios for the standard-sized regions.
		printVecHeader();
		SeqDataLikelihoodRatio selr = new SeqDataLikelihoodRatio(IPhandle, backhandle);
		Iterator<Region> rit = regList.iterator();
		while(rit.hasNext()){
			Region currR = rit.next();
			
			double [] currLLR = selr.llratioRegion(currR);
			String currSeq = seqgen.execute(currR).toUpperCase();
			double [] kmervec = seq2KVec(currSeq, K);
			
			//Print the result here
			printVectors(currR, currLLR, kmervec);
			//printVectors(currR, kmervec);
		}
	}
	
	
	public SeqExptProbLandscapeKmers(String o, String g, String ex, String bk, int k){
		try {
			org = Organism.getOrganism(o);
			gen = org.getGenome(g);
			
			IPhandle = new SeqDataHandler(org, gen, ex);
			backhandle = new SeqDataHandler(org, gen, bk);
			
			loadedRegList = new ArrayList <Region>();
			
			K = k;
			initKPow(K);
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}			
	}
	
	public static void printVecHeader(){
		System.out.print("Chromosome\tStart\tEnd\tAverageLLR");
		for(int v=0; v<kPow[K]; v++){
			System.out.print("\t"+intToSeq(v, K));
		}System.out.print("\n");
	}
	public static void printVectors(Region currR, double [] llr, double [] kvec){
		System.out.print(currR.getChrom()+"\t"+currR.getStart()+"\t"+currR.getEnd()+"\t");
		
		double avLL = 0;
		for(int l=0; l<=llr.length; l++){
			avLL+=llr[l];
		}avLL = avLL/(double)llr.length;
		
		//change sign for now...
		avLL = -1*avLL;
		
		System.out.print("\t"+avLL);
		
		for(int v=0; v<kvec.length; v++){
			System.out.print("\t"+kvec[v]);
		}System.out.print("\n");
	}
	public static void printVectors(Region currR, double [] kvec){
		System.out.print(currR.getChrom()+"\t"+currR.getStart()+"\t"+currR.getEnd()+"\t");
		
		for(int v=0; v<kvec.length; v++){
			System.out.print("\t"+kvec[v]);
		}System.out.print("\n");
	}
	
	public void loadFile(String filename){
		try{
			File rFile = new File(filename);
		
			if(rFile.isFile()){
				regionReader = new BufferedReader(new FileReader(rFile));
				
				String line;
				while((line= regionReader.readLine())!=null){
					String [] tokens = line.split("[-:\\s*\\t\\r\\n\\f]");
					if(tokens.length==3){
						Region r = new Region(gen, tokens[0], new Integer(tokens[1]).intValue(), new Integer(tokens[2]).intValue());
						loadedRegList.add(r);
					}
				}
			}else{
				System.err.println("Cannot read region file!");
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public static ArrayList<Region> standardizeRegions(ArrayList<Region> rlist){
		ArrayList<Region> newList = new ArrayList <Region>();
		
		Iterator<Region> itr = rlist.iterator();
		while(itr.hasNext()){
			Region currR = itr.next();
			int stR = currR.getStart()-(windowSize/2);
			int enR = currR.getEnd()+(windowSize/2);
			for(int r=stR; r<enR; r+=windowOff){
				Region sub = new Region(currR.getGenome(), currR.getChrom(), r, r+windowSize);
				newList.add(sub);
			}
		}
		return(newList);
	}
	public ArrayList<Region> getLoadedRegions(){
		return loadedRegList;
	}
	public static double [] seq2KVec(String seq, int k){
		double [] vec = new double[kPow[k]];
		for(int i=0; i<kPow[k]; i++){
			vec[k]=0;
		}
		
		for(int i=0; i<seq.length()-k; i++){
			String mer = seq.substring(i, i+k);
			int s =seqToInt(mer, k, true);
			if(s!=-1){
				vec[s]++;
			}
		}		
		return(vec);
	}
	
	private void initKPow(int k){
		kPow = new int [k+1];
		for(int i=0; i<=k; i++){
			kPow[i]=(int)Math.pow(B, i);
		}
	}
	private static int seqToInt(String seq, int currk, boolean strandDegeneracy){
		int bin=0;
		int len = seq.length();
		if(len != currk){System.err.println("Error: k's don't match: "+len);}

		for(int i=0; i<len; i++){
			int currInt = base2int(seq.charAt(i));
			if(currInt==-1){
				return(-1);
			}
			bin+=currInt*kPow[len-i-1];
		}		
		if(strandDegeneracy){
			int revbin=0;
			char [] revseq = new char[currk];
		    for (int i = (currk - 1); i >= 0; i--)
		      revseq[(currk-1)-i] = complement(seq.charAt(i));
		    
			for(int i=0; i<len; i++){
				int currInt = base2int(revseq[i]);
				if(currInt==-1){
					return(-1);
				}
				revbin+=currInt*kPow[len-i-1];
			}if(revbin<bin){bin=revbin;}
		}
		return(bin);
	}
	
	private static String intToSeq(int id, int currk){
	    int x = id;
	    StringBuffer currWord = new StringBuffer(currk);
		if(x>=kPow[currk]){
			return("");
		}else{
			for(int i=currk-1; i>=0; i--){
				int m = (int)x/kPow[i];
				currWord.append(int2base((int)m));
				x = x-(m*kPow[i]);
			}
			return(currWord.toString());
		}
	}
	
	private static char complement(char c){
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
	
	private static int base2int(char c){
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
	
	private static char int2base(int x){
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
