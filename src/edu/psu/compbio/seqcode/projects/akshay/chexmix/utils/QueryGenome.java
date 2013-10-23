package edu.psu.compbio.seqcode.projects.akshay.chexmix.utils;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.parsing.FASTAStream;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.*;

public abstract class QueryGenome {
	public String genomepath;
	public static Map<String,String> cache = new HashMap<String,String>();
	public abstract void fillGenomePath();
	public abstract void fillGenomePath(String path);
	public QueryGenome() {}
	
	private void cache(BindingLocation bl)  throws IOException{
		String chromid = bl.getChr();
		// if the chromosome is already cached exit method
		synchronized(cache) {
			if (cache.containsKey(chromid)) {
				return;
			}
		}
		String chromseq = null;
		File f = new File( this.genomepath + "/" + chromid + ".fa");
		if (f.exists()) {
			FASTAStream stream = new FASTAStream(f);
			while (stream.hasNext()) {
				Pair<String,String> pair = stream.next();
				if (pair.car().equals(chromid)) {
					chromseq = pair.cdr();
					break;
				}
			}
			stream.close();
		}
		else{
			System.out.print("Genome is not found at "+this.genomepath+". \n");
			System.exit(-1);
		}
		if (chromseq == null) {
			return;
		}
		synchronized(cache) {
			if (!cache.containsKey(chromid)) {
				cache.put(chromid, chromseq);
			}
		}
	}
	
	 public String execute(BindingLocation bl) throws IOException{
		 String ret = null;
		 String Chrname = bl.getChr();
		 cache(bl);
		 String chromString = null;
		 synchronized(cache) {
			 if (!cache.containsKey(Chrname)) {
				 return null;
			 }
			 chromString = cache.get(Chrname);	
		 }
		 ret = chromString.substring(bl.getCoords().get(0)-1, bl.getCoords().get(1));
		 return ret;
	 }
}
