package org.seqcode.genome.sequence;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class WildcardKmerUtils {

	/**
	 * A hash map the holds all the k-mer mapping functions (should be loaded)
	 */
	public static Map<String, Set<String>> wildcardMap = new HashMap<String, Set<String>>();

	/** length of the k-mers */
	public static int k;

	@SuppressWarnings("unchecked")
	public WildcardKmerUtils(int kmerLen) throws IOException {
		k = kmerLen;
		if (k == 8) {
			InputStream ins = this.getClass().getResourceAsStream("wildcard_8mer_2mismatch_hg19.txt");
			BufferedReader br = new BufferedReader(new InputStreamReader(ins));
			String line = null;
			while ((line = br.readLine()) != null) {
				String[] pieces = line.split(",");
				Set<String> tmpAdd = new HashSet<String>();
				for (int s = 1; s < pieces.length; s++) {
					tmpAdd.add(pieces[s]);
					if (pieces[s].contains("N")) { // is this a wild-card kmer
						if (wildcardMap.containsKey(pieces[s])) {
							wildcardMap.get(pieces[s]).add(pieces[0]);
						} else {
							wildcardMap.put(pieces[s], new HashSet<String>());
							wildcardMap.get(pieces[s]).add(pieces[0]);
						}
						// now also the rev complement
						String rev = SequenceUtils.reverseComplement(pieces[s]);
						if (wildcardMap.containsKey(rev)) {
							wildcardMap.get(rev).add(pieces[0]);
						} else {
							wildcardMap.put(rev, new HashSet<String>());
							wildcardMap.get(rev).add(pieces[0]);
						}
					}
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
	public int base2int(char base) {
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
			break;
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
	public char int2base(int x) {
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
			break;
		default:
			throw new IllegalArgumentException("Invalid int: " + x);
		}
		return (base);
	}

	public int seq2int(String seq) {
		int intVal = 0;
		int len = seq.length();

		for (int i = 0; i < len; i++) {
			long currInt = base2int(seq.charAt(i));
			if (currInt == -1) {
				return -1;
			}
			intVal = intVal * 5;
			intVal += currInt;
		}
		return intVal;
	}

	public String int2seq(long x, int kmerLen) {
		if (x >= (int) Math.pow(5, kmerLen)) {
			throw new IllegalArgumentException("Invalid int value, " + x + ", for kmerLen " + kmerLen);
		}
		StringBuffer seq = new StringBuffer(kmerLen);
		for (int i = 0; i < kmerLen; i++) {
			int baseVal = (int) (x % 5);
			seq.append(int2base(baseVal));
			x = (long) Math.floor(x / 5.0);
		}
		return seq.reverse().toString();
	}

	public Set<String> map(String kmer) {
		return wildcardMap.get(kmer);
	}

	public static void main(String[] args) throws IOException, ClassNotFoundException {
		WildcardKmerUtils wku = new WildcardKmerUtils(8);
		System.out.println(wku.wildcardMap.get("AAAAAAAA"));

	}

}
