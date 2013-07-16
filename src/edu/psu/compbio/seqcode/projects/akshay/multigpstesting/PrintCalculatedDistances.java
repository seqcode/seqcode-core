package edu.psu.compbio.seqcode.projects.akshay.multigpstesting;

import java.io.*;
import java.util.*;

import edu.psu.compbio.seqcode.projects.akshay.multigpstesting.datatypes.*;

public class PrintCalculatedDistances {
	public List<Expts> conditions = new ArrayList<Expts>();
	public Map<Range,Integer> unionBlacklist = new HashMap<Range,Integer>();
	
	private void generateExpFromFile(String filename) throws IOException{
		String exptname;
		String[] pieceslvl1 = filename.split("/");
		String[] pieceslvl2= pieceslvl1[pieceslvl1.length-1].split("_");
		exptname = pieceslvl2[2]+"_"+pieceslvl2[3]+"_"+pieceslvl2[4];
		BufferedReader br =null;
		String currentline;
		br = new BufferedReader(new FileReader(filename));
		currentline = br.readLine();
		Expts currentExpt = new Expts();
		while(currentline != null){
			if(!currentline.startsWith("#")){
				Point currentPoint = new Point(currentline.split("\t")[0]);
				currentExpt.mapPointstoRnages(currentPoint, exptname);
			}
		}
		conditions.add(currentExpt);
	br.close();
	}
	
	public void fillUnionBlacklist(){
		for(Expts iterExpts : conditions){
			iterExpts.fillBlacklist();
			for(Range iterRange: iterExpts.getBlacklist()){
				if(!unionBlacklist.containsKey(iterRange)){
					unionBlacklist.put(iterRange, 1);
				}
			}
		}
	}
	
	public void fillAllisolatedPoints(){
		for(Expts iterExpts : conditions){
			iterExpts.fillIsolatedPoints(unionBlacklist);
		}
	}
	
	
	
	
	
	public void ScanConditions(String designfile) throws IOException{
		BufferedReader brnew = null;
		brnew = new BufferedReader(new FileReader(designfile));
		String CurrentLine = brnew.readLine();
		while(CurrentLine != null){
			generateExpFromFile(CurrentLine);
		}
		brnew.close();
	}
	
	
	

	
	
}
