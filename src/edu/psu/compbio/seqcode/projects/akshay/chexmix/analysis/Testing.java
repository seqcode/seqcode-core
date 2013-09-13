package edu.psu.compbio.seqcode.projects.akshay.chexmix.analysis;


import java.util.List;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.BindingLocation;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.Config;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.Node;

public class Testing {
	
	public static void main(String[] args){
		Config c = new Config(args);
		Chexmix driver = new Chexmix(c);
		if(driver.c.helpWanted()){
			System.err.println("Chexmix:");
			System.err.println(driver.c.getArgsList());
		}
		else{
			LoadTags loadTags = new LoadTags();
			loadTags.loadHits(c, c.useNonUnique());
			BindingLocation bl =  new BindingLocation(151128370,"chr1",c);
			bl.filltags(loadTags);
			//List<Integer> temp = bl.getConcatenatedTags(bl.vecpos.midpoint, bl.vecpos.range, bl.vecpos.orientation);
			Node no = new Node(bl);
			for(int i=0; i<no.composite.length;i++){
				System.out.println(i+"\t"+no.composite[i]);
			}
		}
	}

}
