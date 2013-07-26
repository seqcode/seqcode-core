package edu.psu.compbio.seqcode.projects.akshay.multigpstesting;

import java.io.*;
import java.util.*;

import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.projects.akshay.multigpstesting.datatypes.*;

public class PrintCalculatedDistances {
	public List<Expts> conditions = new ArrayList<Expts>();
	public Map<Range,Integer> unionBlacklist = new HashMap<Range,Integer>();
	
	
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
	
	public void printDistances(String outfilename)throws IOException{
		File file = new File(outfilename);
		if(!file.exists()){
			file.createNewFile();
		}
		FileWriter fw = new FileWriter(file.getAbsoluteFile());
		BufferedWriter brwrite = new BufferedWriter(fw);
		for(int i=0; i<conditions.size()-1; i++){
			for(int j=i+1; j<conditions.size();j++){
				for(Point iterPoint : conditions.get(i).getIsolatedPointMap().keySet()){
					brwrite.write(Integer.toString(conditions.get(j).getNearestDistance(iterPoint))+"\n");
				}
			}
		}
		brwrite.close();
	}
	
	public void ScanConditions(String designfile) throws IOException{
		BufferedReader brnew = null;
		brnew = new BufferedReader(new FileReader(designfile));
		String CurrentLine = brnew.readLine();
		while(CurrentLine != null){
			Expts currentExpt = new Expts();
			currentExpt.mapPointstoRanges(CurrentLine);
			conditions.add(currentExpt);
			CurrentLine = brnew.readLine();
		}
		brnew.close();
	}
	
	public static void main(String[] args) throws IOException{
		Set<String> commandlineArgs = Args.parseArgs(args);
		if(commandlineArgs.contains("help") || commandlineArgs.size() == 0){
			System.out.println("PrintCalculatedDistances usage:\n" +
					"\t--design <design file>\n" +
					"\t--out <output file name>");
		}
		else{
			PrintCalculatedDistances driver = new PrintCalculatedDistances();
			Collection<String> designfiles =  Args.parseStrings(args, "design");
			Collection<String> outputfilenames = Args.parseStrings(args, "out");
			String designfilename = (String) designfiles.toArray()[0];
			String outfilename = (String) outputfilenames.toArray()[0];
			driver.ScanConditions(designfilename);
			driver.fillUnionBlacklist();
			driver.fillAllisolatedPoints();
			driver.printDistances(outfilename);
		}
	}
}
