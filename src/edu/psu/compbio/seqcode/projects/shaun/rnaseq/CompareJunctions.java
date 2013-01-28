package edu.psu.compbio.seqcode.projects.shaun.rnaseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import cern.jet.stat.Probability;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;

public class CompareJunctions {
	
	private String fileA, fileB;
	private ArrayList<Junction> junctionsA, junctionsB;
	private Double totalCoverageA=0.0, totalCoverageB=0.0;
	private Double totalHitsA=0.0, totalHitsB=0.0;
	private int collapseDist = 5;
	private ArrayList<MatchedJunctions> matched=new ArrayList<MatchedJunctions>();
	private Double LOG2 = Math.log(2);
	
	public CompareJunctions(File j1, File j2, Integer c){
		fileA = j1.getName(); fileB = j2.getName();
		junctionsA = loadJunctions(j1);
		junctionsB = loadJunctions(j2);
		collapseDist = c;
		
		//Sum coverage
		totalCoverageA = sumCoverage(junctionsA);
		totalCoverageB = sumCoverage(junctionsB);
		totalHitsA = totalCoverageA;
		totalHitsB = totalCoverageB;
		
	}
	
	//Run the methods
	public void execute(){
		System.err.println("A:\t"+junctionsA.size()+"\t"+totalCoverageA+"\t"+totalHitsA);
		System.err.println("B:\t"+junctionsB.size()+"\t"+totalCoverageB+"\t"+totalHitsB);

		ArrayList<Integer> processedA = new ArrayList<Integer>();
		ArrayList<Integer> processedB = new ArrayList<Integer>();
				
		//Collapse junctions
		junctionsA = collapseJunctions(junctionsA);
		junctionsB = collapseJunctions(junctionsB);
		
		//Compare coverage levels 
		
		//First overlapping
		for(int x = 0; x<junctionsA.size(); x++){
			Junction a = junctionsA.get(x);
			boolean found=false;
			for(int y = 0; y<junctionsB.size() && !found; y++){
				Junction b = junctionsB.get(y);
				if(a.isEquivalentTo(b, collapseDist)){
					matched.add(new MatchedJunctions(a, b));
					processedA.add(x);
					processedB.add(y);
					found=true;
				}
			}
		}
		//Now set-specific
		for(int x = 0; x<junctionsA.size(); x++){
			Junction a = junctionsA.get(x);
			if(!processedA.contains(x)){
				matched.add(new MatchedJunctions(a, null));
			}
		}
		for(int y = 0; y<junctionsB.size(); y++){
			Junction b = junctionsB.get(y);
			if(!processedB.contains(y)){
				matched.add(new MatchedJunctions(null, b));
			}
		}
		
		//Sort matched here
		Collections.sort(matched);
		matched = benjaminiHochbergCorrection(matched);
		Collections.sort(matched);
		
		//Print
		System.out.println("# A = "+fileA+", B = "+fileB);
		System.out.println("#JuncA\tJuncB\tCoverageA\tCoverageB\tlog2(Ratio)\tP");
		for(MatchedJunctions m : matched){
			System.out.println(m.toString());
		}
	}
	
	//Load junctions from a BED file
	private ArrayList<Junction> loadJunctions(File jFile){
		ArrayList<Junction> juncs = new ArrayList<Junction>();
		BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(jFile));
			String line = reader.readLine(); //Ignore first line
		    while ((line = reader.readLine()) != null) {
		    	line = line.trim();
		    	String[] words = line.split("\\s+");
		    	String[] offsets = words[10].split(","); 
		    	String chr = words[0];
		    	int l = (new Integer(words[1]))+(new Integer(offsets[0]));
		    	int r = (new Integer(words[2]))-(new Integer(offsets[1]));
		    	double reads = new Double(words[4]);
		    	Junction j = new Junction(chr, l, r, reads); 
		    	juncs.add(j);
		    }
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(juncs);
	}
	
	//Collapse junctions
	private ArrayList<Junction> collapseJunctions(ArrayList<Junction> juncs){
		System.err.println(juncs.size()+" junctions loaded");
		ArrayList<Junction> collapsed = new ArrayList<Junction>();
		while(!juncs.isEmpty()){
			Junction candidate = juncs.get(0);
			juncs.remove(candidate); 
			ArrayList<Integer> removeList = new ArrayList<Integer>();
			for(int x=0; x<juncs.size() && candidate.right > juncs.get(x).left && candidate.chr.equals(juncs.get(x).chr); x++){
				if(candidate.isEquivalentTo(juncs.get(x), collapseDist)){
					candidate.addCollapsed(juncs.get(x));
					removeList.add(x);
				}
			}
			for(int y = removeList.size()-1; y>=0; y--){
				juncs.remove(removeList.get(y).intValue());
			}collapsed.add(candidate);
		}System.err.println(collapsed.size()+" junctions after collapse");
		return collapsed;		
	}
	
	//Sum coverage
	private Double sumCoverage(ArrayList<Junction> juncs){
		Double sumCov=0.0;
		for(Junction x : juncs){
			sumCov+=x.coverage;
		}
		return(sumCov);
	}
	//Set the total hits (if you want to use a full total hit count for statistical testing)
	private void setTotalHitsA(Double tha){totalHitsA=tha;}
	private void setTotalHitsB(Double thb){totalHitsB=thb;}

	//Accessors
	private Double getTotalCoverageA(){return totalCoverageA;}
	private Double getTotalCoverageB(){return totalCoverageB;}
	
	//Multiple hypothesis testing correction -- assumes peaks ordered according to p-value
	protected ArrayList<MatchedJunctions> benjaminiHochbergCorrection(ArrayList<MatchedJunctions> mj){
		double total = mj.size();
		ArrayList<MatchedJunctions> res = new ArrayList<MatchedJunctions>();
		double rank =1;
		for(MatchedJunctions m : mj){
			m.significance = m.significance*(total/rank);
			if(m.significance>1)
				m.significance=1.0;
			res.add(m);
			rank++;
		}return(res);
	}
	
	public static void main(String[] args) {
		if(args.length==0){
			System.out.println("Arguments:\n" +
					"\t--fileA: junctions BED file A\n" +
					"\t--fileB: junctions BED file B\n" +
					"\t--collapse: collapse window distance\n" +
					"\t--atotal: total hits in expt A\n" +
					"\t--btotal: total hits in expt B");
		}else{
			String fileA = Args.parseString(args,"fileA","xxx");
			String fileB = Args.parseString(args,"fileB","xxx");
			File fA = new File(fileA);
			if(!fA.isFile()){System.err.println("Invalid file name: "+fileA);System.exit(1);}
			File fB = new File(fileB);
			if(!fB.isFile()){System.err.println("Invalid file name: "+fileB);System.exit(1);}
			Integer cDist = Args.parseInteger(args,"collapse",5);
			
			CompareJunctions compare = new CompareJunctions(fA, fB, cDist);
			if(Args.parseArgs(args).contains("atotal"))
				compare.setTotalHitsA(Args.parseDouble(args, "atotal", compare.getTotalCoverageA()));
			if(Args.parseArgs(args).contains("btotal"))
				compare.setTotalHitsB(Args.parseDouble(args, "btotal", compare.getTotalCoverageB()));
			compare.execute();
		}
	}
	
	///CLASS: Junction ///
	private class Junction{
		public String chr;
		public int left;
		public int right;
		public double coverage;
		ArrayList<Integer> collapsedLefts = new ArrayList<Integer>();
		ArrayList<Integer> collapsedRights = new ArrayList<Integer>();
		
		public Junction(String c, int l, int r, double cov){
			chr=c; left=l; right=r; coverage=cov;
			collapsedLefts.add(l);
			collapsedRights.add(r);
		}
		
		public boolean isEquivalentTo(Junction j, int window){
			boolean leftEqual=false; 
			boolean rightEqual=false;
			if(chr.equals(j.chr)){
				for(Integer i : collapsedLefts){
					if(Math.abs(j.left-i)<=window)
						leftEqual=true;
				}for(Integer i : collapsedRights){
					if(Math.abs(j.right-i)<=window)
						rightEqual=true;
				}
				if(leftEqual && rightEqual)
					return true;
			}
			return false;
		}
		public String toString(){
			if(chr.equals("NOTFOUND")){
				return(String.format("%s", chr));
			}
			return(String.format("%s:%d-%d", chr, left, right));
		}
		
		public void addCollapsed(Junction j){
			collapsedLefts.add(j.left);
			collapsedRights.add(j.right);
			if(j.coverage>coverage){
				coverage=j.coverage;
				left = j.left;
				right=j.right;
			}
		}
	}

	///CLASS: MatchedJunction ///
	private class MatchedJunctions  implements Comparable<MatchedJunctions>{
		public Junction juncA;
		public Junction juncB;
		public Double logRatio=1.0;
		public Double significance = 0.0;
		
		public MatchedJunctions(Junction A, Junction B){
			juncA=A; juncB=B;
			if(A==null)
				juncA = new Junction("NOTFOUND", -1, -1, 1.0);
			if(B==null)
				juncB = new Junction("NOTFOUND", -1, -1, 1.0);
			
			logRatio = Math.log((juncA.coverage/totalCoverageA)/(juncB.coverage/totalCoverageB))/LOG2;
			significance = calcSignificance();			
		}
		
		private Double calcSignificance(){   //Binomial test
			Double s =0.0;
			s = binomialSampleEquality(juncA.coverage, juncB.coverage, totalHitsA, totalHitsB);
			return(s);
		}
		public String toString(){
			return(String.format("%s\t%s\t%.1f\t%.1f\t%.4f\t%e", juncA.toString(),juncB.toString(),juncA.coverage,juncB.coverage,logRatio,significance));
		}
		
		// Binomial test for differences between two population proportions 
		protected double binomialSampleEquality(double X1, double X2, double n1, double n2){
			double P1 = X1/n1;
			double P2 = X2/n2;
			double P = (X1+X2)/(n1+n2);
			double Z = (P1-P2)/(Math.sqrt(P*(1-P)*((1/n1)+(1/n2))));
			if(!Double.isNaN(Z))
				if(X1/n1 > X2/n2)
					return(1-Probability.normal(Z));
				else
					return(Probability.normal(Z));
			else
				return(-1);
		}
				
		//Comparable default method
		public int compareTo(MatchedJunctions m) {
			if(significance<m.significance){return(-1);}
			else if(significance>m.significance){return(1);}
			else{return(0);}
		}
	}	
}
