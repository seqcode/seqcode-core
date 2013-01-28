package edu.psu.compbio.seqcode.projects.shaun;

import java.util.ArrayList;

import edu.psu.compbio.seqcode.gse.utils.Pair;

public class Score2Sp {
	
	private ArrayList<Pair<Double,Double>> scores;
	
	public Score2Sp (ArrayList<Pair<Double,Double>> s){
		scores = s;
	}
	
	public Double getSp(Double score){
		Double sp=0.0;
		int i=0;
		Pair<Double,Double> p = scores.get(i);
		while(i<scores.size() && p.car()<score){
			p = scores.get(i);
			i++;
		}
		sp=p.cdr();
		return(sp);
	}
	public Double getScore(Double sp){
		Double score=0.0;
		int i=0;
		Pair<Double,Double> p = scores.get(i);
		while(i<scores.size() && p.cdr()<sp){
			p = scores.get(i);
			i++;
		}
		score=p.car();
		return(score);
	}
}
