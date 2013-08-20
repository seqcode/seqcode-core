package edu.psu.compbio.seqcode.projects.akshay.chexmix.utils;

import java.io.*;

public class QueryTagsIdx extends QueryTags{

	public QueryTagsIdx(int midpoint, int range, String chr) {
		super(midpoint, range, chr);
	}

	@Override
	public void fillTagsBedPath(String tagspath) {
		String currdir = System.getProperty("user.dir");
		File temptagsbed = new File(currdir+"/temp/temptagsbed");
		if(temptagsbed.exists()){
			this.tagsbedfilepath = currdir+"/temp/temptagsbed";
		}
		else{
			try{
				BufferedReader br = null;
				br = new BufferedReader(new FileReader(tagspath));
				String currentline = br.readLine();
				FileWriter fstream = new FileWriter(currdir+"/temp/temptagsbed", false);
		        BufferedWriter out = new BufferedWriter(fstream);
		        System.out.println(currentline);
				while(!currentline.startsWith("chrom") && currentline != null){
					
					String[] pieces = currentline.split("\t");
					for(int i=0; i< Integer.parseInt(pieces[2]); i++){
						String content = pieces[0]+"\t"+Integer.toString(Integer.parseInt(pieces[0])+1)+"\t"+"*"+"\t"+"*"+"+";
						out.write(content+"\n");
					}
					for(int i=0; i< Integer.parseInt(pieces[3]); i++){
						String content =  pieces[0]+"\t"+Integer.toString(Integer.parseInt(pieces[0])+1)+"\t"+"*"+"\t"+"*"+"-";
						out.write(content+"\n");
					}
					currentline = br.readLine();
				}
				br.close();
				out.close();
				this.tagsbedfilepath = currdir+"/temp/temptagsbed";	
				
			}
			catch (IOException e){
				System.out.println("No such IDX file. Error Description: "
						+ e.getMessage());
			}
		}
	}
			
}
		


