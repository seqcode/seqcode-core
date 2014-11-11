package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.stats.StatUtil;

public class StudentT {
	public double[] vec1;
	public double[] vec2;
	public double pval;
	public double Tstat;
	
	public StudentT(double[] vec1, double[] vec2) {
		this.vec1=vec1;
		this.vec2=vec2;
		
		double mean1 = StatUtil.mean(vec1);
		double mean2 = StatUtil.mean(vec2);
		
		double std1 = StatUtil.std(vec1);
		double std2 = StatUtil.std(vec2);
		
		double Tstat_num = Math.pow(mean1, 2.0) - Math.pow(mean2, 2.0);
		double Tstat_den = Math.sqrt((Math.pow(std1, 2.0)/vec1.length)+(Math.pow(std2, 2.0)/vec2.length));
		Tstat = Tstat_num/Tstat_den;
		
		double v_num =  Math.pow((Math.pow(std1, 2.0)/vec1.length)+(Math.pow(std2, 2.0)/vec2.length),2);
		double v_den = (Math.pow((Math.pow(std1, 2.0)/vec1.length),2.0)/(vec1.length-1))+(Math.pow((Math.pow(std2, 2.0)/vec2.length),2.0)/(vec2.length-1));
		
		double v = v_num/v_den;
		pval = StatUtil.studentTPvalue(Tstat, v);
	}
	
	public double getPval(){return this.pval;}
	public double getTstat(){return this.Tstat;}
	
	
	public static void main(String[] args){
		try{
			ArgParser ap = new ArgParser(args);
			String xFilename = ap.getKeyValue("setA");
			String yFilename = ap.getKeyValue("setB");
			double[] x;
			double[] y;
		
			File xfile = new File(xFilename);
			File yfile = new File(yFilename);
			if(!xfile.exists() || !yfile.exists()){
				throw new FileNotFoundException("Can't find file " + xfile.getName() + " or " + yfile.getName());
			}else{
				BufferedReader xreader = new BufferedReader(new FileReader(xfile));
				String line;
				List<Double> tmpX = new ArrayList<Double>();
				while((line = xreader.readLine()) != null){
					String[] pieces  = line.split("\t");
					if(pieces.length == 1){
						tmpX.add(Double.parseDouble(pieces[0]));
					}else{
						tmpX.add(Double.parseDouble(pieces[1]));
					}
				}
				xreader.close();
				BufferedReader yreader = new BufferedReader(new FileReader(yfile));
				line="";
				List<Double> tmpY = new ArrayList<Double>();
				while((line = yreader.readLine()) != null){
					String[] pieces  = line.split("\t");
					if(pieces.length == 1){
						tmpY.add(Double.parseDouble(pieces[0]));
					}else{
						tmpY.add(Double.parseDouble(pieces[1]));
					}
				}
				yreader.close();
				x = new double[tmpX.size()];
				y = new double[tmpY.size()];
				for(int i=0; i<tmpX.size(); i++){
					x[i] = tmpX.get(i);
				}
				for(int i=0; i<tmpY.size(); i++){
					y[i] = tmpY.get(i);
				}
				
				StudentT tester = new StudentT(x,y);
				System.out.println(tester.getPval());
				System.out.println(tester.getTstat());
				
			}
			
			
		}catch (FileNotFoundException e){
			e.printStackTrace();
		}catch (IOException e){
			e.printStackTrace();
		}
		
		
	}
	
	
	

}
