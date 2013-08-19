package edu.psu.compbio.seqcode.projects.akshay.chexmix.utils;

import java.util.*;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.*;

public abstract class CompareVecs {
	
	public List<Integer> first;
	public List<Integer> second;
	
	public CompareVecs(List<Integer> first, List<Integer> second) {
		this.first = first;
		this.second = second;
	}
	
	public abstract double doComparision();
	

}
