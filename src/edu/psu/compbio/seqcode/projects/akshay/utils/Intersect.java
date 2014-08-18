package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;




/**
 * Similar to bedtools intersect with several added functions and access to readdb experiments
 * @author akshaykakumanu
 *
 */
public class Intersect {

	List<Point> a;
	List<Point> b;
	
	SeqLocator expta;
	SeqLocator exptb;
	
	Genome gen;
	
	int min_match_distance;
	
	
	public Intersect() {
		
	}
	/**
	 * Prints the matching overlapping win
	 */
	public void intersect(){
		
	}
	
}
