package org.seqcode.data.io.parsing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Gene;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.Expander;

//TODO: Add GTF loader
public class TranscriptFileExpander <X extends Region> 
	implements Expander<X,Gene>{
	
	protected File f;
	protected Genome gen;
	protected List<Gene> allGenes = new ArrayList<Gene>();

	public TranscriptFileExpander(Genome g, String fileName){
		f = new File(fileName);
		gen=g;
		try{
			BufferedReader reader = new BufferedReader(new FileReader(f));
			String line;
			while ((line = reader.readLine()) != null) {
		    	line = line.trim();
		        String[] words = line.split("\\s+");
		        
		        if(words.length==5){
		        	String chr = words[0];
		        	String[] tmp = chr.split("\\.");
        			chr=tmp[0].replaceFirst("^chr", "");
		        	String name = words[1];
		        	Integer start = new Integer(words[2]);
		        	Integer end = new Integer(words[3]);
		        	char strand = words[4].charAt(0);
		        	allGenes.add(new Gene(gen, chr, start, end, name, name, strand, "file"));
		        }
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public Iterator<Gene> execute(X a) {
		List<Gene> currGenes = new ArrayList<Gene>();
		for(Gene x : allGenes){
			if(x.overlaps(a)){
				currGenes.add(x);
			}
		}return currGenes.iterator();
	}

}
