package org.seqcode.projects.akshay.chexmix.utils;

import java.util.ArrayList;
import java.util.List;

import org.seqcode.projects.akshay.chexmix.datasets.Vec;


public class Pearson extends CompareVecs{

	public Pearson(List<Integer> first, List<Integer> second) {
		super(first, second);
		// TODO Auto-generated constructor stub
	}

	@Override
	public double doComparision() {
		double pcc=0.0;
		if(this.first.size() != this.second.size()){
			System.err.println("The lengths of the vectors must be equal");
			System.exit(1);
		}
		else{
			double meanVec1=0.0;
			double meanVec2=0.0;
			for(int i=0; i<this.first.size(); i++){
				meanVec1 = meanVec1+this.first.get(i);
				meanVec2 = meanVec2+this.second.get(i);
			}
			meanVec1 = meanVec1/this.first.size();
			meanVec2 = meanVec2/this.second.size();
			double num=0.0;
			double den1=0.0;
			double den2=0.0;
			for(int i=0; i<this.first.size(); i++){
				num = num + (this.first.get(i)-meanVec1)*(this.second.get(i)-meanVec2);
				den1 = den1+(this.first.get(i)-meanVec1)*(this.first.get(i)-meanVec1);
				den2 = den2+(this.second.get(i)-meanVec2)*(this.second.get(i)-meanVec2);
			}
			double den = Math.sqrt(den1*den2);
			pcc = num/den;
				
		}
	return pcc;
	}
	
	

}
