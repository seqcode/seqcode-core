package edu.psu.compbio.seqcode.projects.akshay.multigpstesting;

import java.io.*;
import java.util.*;

import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.projects.akshay.multigpstesting.datatypes.*;

public class PrintDiffEnrichedFraction {
	
	public static void main(String[] args) throws IOException{
		Expts allevents = new Expts();
		List<Expts> combinations = new ArrayList<Expts>();
		Map<Point,String> diffsite = new HashMap<Point, String>();
		
		Set<String> commandlineArgs = Args.parseArgs(args);
		if(commandlineArgs.contains("help") || commandlineArgs.size() == 0){
			System.out.println("PrintDiffEnrichedFraction usage:\n" +
					"\t--design <design file>\n");
		}
		else{
			Collection<String> designfiles =  Args.parseStrings(args, "design");
			String designfilename = (String) designfiles.toArray()[0];
			BufferedReader br  = null;
			br = new BufferedReader(new FileReader(designfilename));
			String currentline = br.readLine();
			allevents.mapPointstoRanges(currentline);
			currentline = br.readLine();
			while(currentline != null){
				Expts currentexpt = new Expts();
				currentexpt.mapPointstoRanges(currentline);
				combinations.add(currentexpt);
				currentline = br.readLine();
			}
			br.close();
			for (Point currentpt : allevents.getListofPointMap().keySet()){
				boolean diff = false;
				for(Expts expt:combinations){
					diff = (diff ? true : expt.getListofPointMap().containsKey(currentpt));
				}
				if(diff && !diffsite.containsKey(currentpt)){
					diffsite.put(currentpt, "");
				}
				
			}
			System.out.println(diffsite.size());
			System.out.println(allevents.getListofPointMap().size());
		}
		
	}

}
