package org.seqcode.motifs;

import java.util.*;
import java.io.*;

import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gsebricks.verbs.location.PointParser;
import org.seqcode.gsebricks.verbs.location.RegionParser;
import org.seqcode.gseutils.*;


/** Use:
 *   java org.seqcode.motifs.CountKmers --species "$SC;sacCer3" --mink 4 --maxk 6 [--outputcounts] [--includerc] [--topn 100]
 *
 *   outputcounts: output counts instead of frequencies
 *   includerc: include counts from reverse complement strand too
 *   topn: only output the top N kmers instead of all of them.
 *   table: output frequencies per region
 *
 */

public class CountKmers {

    private int mink, maxk;
    private Genome genome;
    private GenomeConfig gconfig;
    private List<Region> regions;
    private int win;
    private List<int[]> counts;
    private SequenceGenerator seqgen;
    private boolean outputCounts, includeRC;
    private int topN;
    private boolean fullTable=false;

    /* use this if you're going to feed in sequences, but not Regions */
    public void init(int mink, int maxk) {
        this.mink = mink;
        this.maxk = maxk;
        counts = new ArrayList<int[]>();
        for (int i = 0; i <= maxk; i++) {
            int[] l = new int[(4 << ((i-1) * 2))];
            for (int j = 0; j < l.length; j++) {
                l[j] = 0;
            }
            counts.add(l);
        }
    }

    public void parseArgs(String args[]) throws NotFoundException {
        gconfig = new GenomeConfig(args);
        Genome genome = gconfig.getGenome();
        win = Args.parseInteger(args,"win",-1);
        String regFile = Args.parseString(args, "regions", null);
        regions = loadRegionsFromFile(regFile, genome, win);
        outputCounts = Args.parseFlags(args).contains("outputcounts");
        includeRC = Args.parseFlags(args).contains("includerc");
        fullTable = Args.parseFlags(args).contains("table");
        topN = Args.parseInteger(args,"topn",-1);
        seqgen = gconfig.getSequenceGenerator();

        for(Region r : regions){
        	System.out.println(r.getLocationString()+"\t"+r.getWidth());
        }
        init(Args.parseInteger(args,"mink",1),
             Args.parseInteger(args,"maxk",4));
    }
    
	
	/**
	 * Loads a set of regions from the third or first column of a file
	 * (Suitable for GPS & StatisticalPeakFinder files
	 * @param filename String
	 * @param win integer width of region to impose (-1 leaves region width alone)
	 * @return
	 */
	public static List<Region> loadRegionsFromFile(String filename, Genome gen, int win){
		List<Region> regs = new ArrayList<Region>();

		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid file name: "+filename);System.exit(1);}
			BufferedReader reader = new BufferedReader(new FileReader(pFile));
			String line;
			while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            
	            if(words.length>0 && !words[0].contains("#") && !words[0].equals("Region") && !words[0].equals("Position")){
	                if(words.length>=3 && words[2].contains(":")){
		                PointParser pparser = new PointParser(gen);
		            	Point p = pparser.execute(words[2]);
		                if(win==-1 && words[0].contains(":") && words[0].contains("-")){
		                	RegionParser rparser = new RegionParser(gen);
			            	Region q = rparser.execute(words[0]);
			            	regs.add(q);
		                }else{
		                	regs.add(p.expand(win/2));
		                }
	                }else if(words.length>=1 && words[0].contains(":")){
	                	String[] coords = words[0].split(":");
		            	if(coords[1].contains("-")){
		                	RegionParser rparser = new RegionParser(gen);
			            	Region q = rparser.execute(words[0]);
			            	if(win==-1){
			                	if(q!=null){regs.add(q);}
			                }else
			                	regs.add(q.getMidpoint().expand(win/2));
		            	}else{
		            		PointParser pparser = new PointParser(gen);
			            	Point p = pparser.execute(words[0]);
			            	regs.add(p.expand(win/2));
		            	}
		            }
                }
	        }reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return(regs);
	}

    public void addToCounts(Region r) {
        String s = seqgen.execute(r);
        addToCounts(s);
        if (includeRC) {
            addToCounts(SequenceUtils.reverseComplement(s));
        }
    }
    public void addToCounts(String s) {
        char[] chars = s.toCharArray();
        for (int i = 0; i < chars.length; i++) {
            if (chars[i] == 'A' || chars[i] == 'a') {
                chars[i] = 0;
            } else if (chars[i] == 'C' || chars[i] == 'c') {
                chars[i] = 1;
            } else if (chars[i] == 'G' || chars[i] == 'g') {
                chars[i] = 2;
            } else {
                chars[i] = 3;
            }
        }
        for (int k = mink; k <= maxk; k++) {
            int[] l = counts.get(k);
            for (int i = 0; i < chars.length - k; i++) {
                int index = 0;
                for (int j = 0; j < k; j++) {
                    char c = chars[i + j];
                    if (c > 3) {
                        i += k;
                        break;
                    }
                    index = (index << 2) + c;
                }
                l[index]++;
            }
        }
    }

    public String indexToString(int index, int k) {
        char[] out = new char[k];
        int pos = k - 1;
        while (pos >= 0) {
            out[pos--] = (char)(index & 3);
            index >>= 2;
        }
        for (int i = 0; i < out.length; i++) {
            if (out[i] == 0) {
                out[i] = 'A';
            } else if (out[i] == 1) {
                out[i] = 'C';
            } else if (out[i] == 2) {
                out[i] = 'G';
            } else {
                out[i] = 'T';
            }
        }
        return new String(out);
    }

    public Set<String> getKeySet(int k) {
        HashSet<String> s = new HashSet<String>();
        int[] l = counts.get(k);
        for (int i = 0; i < l.length; i++) {
            if (l[i] > 0) {
                s.add(indexToString(i,k));
            }
        }
        return s;       
    }

    public Map<String,Integer> getCounts(int k) {
        int[] l = counts.get(k);
        Map<String,Integer> output = new HashMap<String,Integer>();
        for (int i = 0; i < l.length; i++) {
            if (l[i] > 0) {
                output.put(indexToString(i,k),
                           l[i]);
            }
        }
        return output;
    }

    public int getCount(String key, int k) throws NumberFormatException {
        char[] chars = key.toCharArray();
        for (int i = 0; i < chars.length; i++) {
            if (chars[i] == 'A' || chars[i] == 'a') {
                chars[i] = 0;
            } else if (chars[i] == 'C' || chars[i] == 'c') {
                chars[i] = 1;
            } else if (chars[i] == 'G' || chars[i] == 'g') {
                chars[i] = 2;
            } else {
                chars[i] = 3;
            }
        }
        int index = 0;
        for (int j = 0; j < k; j++) {
            char c = chars[j];
            if (c > 3) {
                throw new NumberFormatException("Invalid character in " + key);
            }
            index = (index << 2) + c;
        }
        return counts.get(k)[index];
    }

    public int getMinCount(int k) {
        if (topN < 0) {
            return 0;
        }
        ArrayList<Integer> list = new ArrayList<Integer>();
        int[] array = counts.get(k);
        for (int i = 0; i < array.length; i++) {
            list.add(array[i]);
        }
        Collections.sort(list);       
        if (topN == 0) {
            return list.get(list.size() - 1) + 1;
        }

        return list.get(list.size() - topN);
    }

    public void print(PrintWriter pw) {
        System.err.println(String.format("Printing from %d to %d", mink, maxk));
        for (int k = mink; k <= maxk; k++) {
            int[] l = counts.get(k);
            int minCount = getMinCount(k);
            System.err.println(String.format("Length at %d is %d",k,l.length));
            long total = 0;
            for (int i = 0; i < l.length; i++) {
                total += l[i];
            }
            for (int i = 0; i < l.length; i++) {
                String key = indexToString(i,k);
                if (l[i] < minCount) {
                    continue;
                }
                if (outputCounts) {
                    pw.println(key + "\t" + l[i]);
                } else {
                    pw.println(key + "\t" + (((double)l[i]) / ((double)total)));
                }
            }
        }
    }

    public static void main(String args[]) throws Exception {
        CountKmers counter = new CountKmers();
        counter.parseArgs(args);
        for (Region r : counter.regions) {
            counter.addToCounts(r);
        }
        PrintWriter pw = new PrintWriter(System.out);
        counter.print(pw);
        pw.close();
    }


}