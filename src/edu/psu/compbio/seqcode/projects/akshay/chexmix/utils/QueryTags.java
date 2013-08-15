package edu.psu.compbio.seqcode.projects.akshay.chexmix.utils;

import java.io.*;
import java.util.*;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.*;

public abstract class QueryTags {
	public int midpoint;
	public String chr;
	public int range;
	public String tagsbedfilepath;
	public Map<Integer, Integer> tags = new TreeMap<Integer, Integer>();
	
	public QueryTags(int midpoint, int range, String chr) {
		this.midpoint = midpoint;
		this.range = range;
		this.chr = chr;
	}
	
	
	public void prepareQureybed(String orientation) throws IOException{
		String currdir = System.getProperty("user.dir");
		FileOutputStream fop = null;
		File file;
		file = new File(currdir+"/temp/"+"tempQueryInterval.bed");
		if(file.exists()){
			file.delete();
		}
		file.createNewFile();

		fop = new FileOutputStream(file);
		String content="";
		for(int i=this.midpoint-this.range/2; i<this.midpoint+this.range/2;i++){
			if(content != ""){
				content = content+"\n";
			}
			content = content + this.chr+"\t"+Integer.toString(i)+"\t"+Integer.toString(i+1)+"\t"+"*"+"\t"+"*"+"\t"+orientation;
		}
		
		byte[] contentinbytes = content.getBytes();
		fop.write(contentinbytes);
		fop.flush();
		fop.close();
		
	}
	
	public abstract void fillTagsBedPath(String tagspath);

	public Vec getTags(String orientation){
		Vec ret = null;
		String currdir = System.getProperty("user.dir");
		List<String> command = new ArrayList<String>();
		command.add("bedtools");
		command.add("intersect");
		command.add("-c");
		command.add("-s");
		command.add("-a");
		command.add(currdir+"/temp/"+"tempQueryInterval.bed");
		command.add("=b");
		command.add(this.tagsbedfilepath);
		command.add("temp.fa");
		ProcessBuilder pb = new ProcessBuilder(command);
		try{
			Process shell = pb.start();
			shell.waitFor();
	        BufferedReader br = new BufferedReader(new InputStreamReader(shell.getInputStream()));
	        String line = null, previous = null;
	        while ((line = br.readLine()) != null){
	            if (!line.equals(previous)) {
	                previous = line;
	                String[] pieces = line.split("\t");
	                this.tags.put(Integer.parseInt(pieces[1]),Integer.parseInt(pieces[6]));
	             }
	        }
	        
	        ret = new Vec(this.range, this.chr, orientation, false, 1, this.tags);
	        
		}
		catch (IOException e) {
			System.out.println("Error occured while executing Linux command. Error Description: "
			+ e.getMessage());
		} catch (InterruptedException e) {
			System.out.println("Error occured while executing Linux command. Error Description: "
					+ e.getMessage());
		}
		return ret;
	}
}
		

