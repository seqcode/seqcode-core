package edu.psu.compbio.seqcode.projects.akshay.chexmix.utils;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.*;

public abstract class QueryGenome {
	public String genomepath;
	public String chr;
	public int midpoint;
	public int range;
	public abstract void fillGenomePath();
	public abstract void fillGenomePath(String path);
	public QueryGenome(String chr, int midpoint, int range) {
		this.chr = chr;
		this.midpoint = midpoint;
		this.range = range;
	}
	public Seq getSeq(String orientation) throws IOException{
		String currdir = System.getProperty("user.dir");
		FileOutputStream fop = null;
		File file;
		file = new File(currdir+"/temp/"+"tempSeqQuery.bed");
		file.createNewFile();
		fop = new FileOutputStream(file);
		String content = this.chr+"\t"+Integer.toString(this.midpoint-this.range)+"\t"+Integer.toString(this.midpoint+this.range)+"\t"+"*"+"\t"+"*"+"\t"+orientation;
		byte[] contentinbytes = content.getBytes();
		fop.write(contentinbytes);
		fop.flush();
		fop.close();
		List<String> command = new ArrayList<String>();
		command.add("bedtools");
		command.add("getfasta");
		command.add("-s");
		command.add("-fi");
		command.add(this.genomepath+"/"+this.chr+".fa");
		command.add("-bed");
		command.add(currdir+"/temp/tempSeqQuery.bed");
		command.add("-fo");
		command.add("temp.fa");
		ProcessBuilder pb = new ProcessBuilder(command);
		try{
			Process shell = pb.start();
			shell.waitFor();
		}
		catch (IOException e) {
			System.out.println("Error occured while executing Linux command. Error Description: "
			+ e.getMessage());
		} catch (InterruptedException e) {
			System.out.println("Error occured while executing Linux command. Error Description: "
					+ e.getMessage());
		}
		BufferedReader br = null;
		br = new BufferedReader(new FileReader(currdir+"/temp/"+"temp.fa"));
		String currentline = br.readLine();
		String tempseq = null;
		while(!currentline.startsWith(">") && currentline != null){
			tempseq = currentline;
			currentline = br.readLine();
		}
		br.close();
		Seq ret = new Seq(this.midpoint, this.range, this.chr, orientation, tempseq);
		return ret;
	}
}
