package org.seqcode.genome.sequence;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.gseutils.Args;

public class WildcardKmerUtils {
	
	/** A hash map the holds all the k-mer mapping functions (should be loaded) */
	public static Map<String,List<String>> wildcardMap = new HashMap<String,List<String>>();

	/** length of the k-mers */
	public static int k;
	
	@SuppressWarnings("unchecked")
	public WildcardKmerUtils(int kmerLen) throws IOException {
		k = kmerLen;
		if(k ==8 ){
			InputStream ins = this.getClass().getResourceAsStream("wildcard_8mer_2mismatch_map_hg19_hashmap.txt");
			BufferedReader br = new BufferedReader(new InputStreamReader(ins));
			String line = null;
			while((line=br.readLine()) != null){
				String[] pieces = line.split(",");
				List<String> tmpAdd = new ArrayList<String>();
                for(int s=1; s<pieces.length; s++){
                        tmpAdd.add(pieces[s]);
                }
                wildcardMap.put(pieces[0], tmpAdd);
			}
		}
	}
	
	
	
	/**
	 * Convert a base to an int value
	 * 
	 * @param base
	 * @return
	 */
	public static int base2int(char base) {
		int intVal = -1;
		switch (base) {
		case 'A':
			intVal = 0;
			break;
		case 'C':
			intVal = 1;
			break;
		case 'G':
			intVal = 2;
			break;
		case 'T':
			intVal = 3;
			break;
		case 'N':
			intVal = 4;
		default:
			throw new IllegalArgumentException("Invalid character: " + base);
		}
		return intVal;
	}


	/**
	 * Return a base for the specified integer
	 * 
	 * @param x
	 * @return
	 */
	public static char int2base(int x) {
		char base;
		switch (x) {
		case 0:
			base = 'A';
			break;
		case 1:
			base = 'C';
			break;
		case 2:
			base = 'G';
			break;
		case 3:
			base = 'T';
			break;
		case 4:
			base = 'N';
		default:
			throw new IllegalArgumentException("Invalid int: " + x);
		}
		return (base);
	}
	
	public static List<String> map(String kmer){
		return wildcardMap.get(kmer);
	}
	
	
	public static void main(String[] args) throws IOException, ClassNotFoundException{
		WildcardKmerUtils wku = new WildcardKmerUtils(8);
		System.out.println(wku.wildcardMap.get("AAAAAAAA"));
		
	}
	
	
}
