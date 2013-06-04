package edu.psu.compbio.seqcode.projects.akshay.Querying;

import java.util.*;
import java.io.*;


public class Scanner {
	
	public Map<String,String[]> queries = new HashMap<String,String[]>();
	
	
	
	public Scanner(String filename) throws IOException {
		BufferedReader br = null;
		try{
			String CurrentLine;
			br = new BufferedReader(new FileReader(filename));
			CurrentLine = br.readLine();
			while(CurrentLine != null){
				if (!queries.containsKey(CurrentLine)){
					String[] pieces = CurrentLine.split("\t");
					queries.put(CurrentLine, pieces);
				}
				CurrentLine = br.readLine();
			}
		}
		catch (IOException e){
			e.printStackTrace();
		}
		catch (NullPointerException e){
			e.printStackTrace();
		}
		finally{
			if (br != null){
				br.close();
			}
			
		}
		}
	
	
	
	public Map<String,String[]> getMap(){
		return queries;
	}
	
	
	
	public static void main(String [ ] args) throws IOException {
		Scanner sc = new Scanner("/Users/akshaykakumanu/Desktop/temp");
		Map<String,String[]> queries = sc.getMap();
		/*System.out.println(queries);*/
		System.out.println(queries.get("kkkk	sasda")[1]);
		
		
		
	}
		
	}
	
