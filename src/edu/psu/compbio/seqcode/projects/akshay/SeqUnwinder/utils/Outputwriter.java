package edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder.framework.SeqUnwinderConfig;
import edu.psu.compbio.seqcode.projects.multigps.stats.StreamGobbler;
import edu.psu.compbio.seqcode.projects.multigps.utilities.Utils;

public class Outputwriter {
	SeqUnwinderConfig seqConfig;
	
	public Outputwriter(SeqUnwinderConfig seqcon) {
		seqConfig = seqcon;
	}
	
	public void writeDiscrimMotifs() throws IOException{
		//First write a transfac file with all the identified discrim motifs
		File motifsTransfac = new File(seqConfig.getOutDir().getAbsoluteFile()+File.separator+"Discrim_motifs.transfac");
		FileWriter fw = new FileWriter(motifsTransfac);
		BufferedWriter bw = new BufferedWriter(fw);

		for(int m =0; m<seqConfig.getDiscrimMotifs().size(); m++){
			String out = WeightMatrix.printTransfacMatrix(seqConfig.getDiscrimMotifs().get(m),seqConfig.getDiscrimMotifs().get(m).getName());
			bw.write(out);
		}
		bw.close();

		// Next, print the ROC values
		File motifsROC = new File(seqConfig.getOutDir().getAbsoluteFile()+File.separator+"Discrim_motifs_roc.tab");
		fw = new FileWriter(motifsROC);
		bw = new BufferedWriter(fw);
		for(String s : seqConfig.getDiscrimMotifRocs().keySet()){
			bw.write(s+"\t"+Double.toString(seqConfig.getDiscrimMotifRocs().get(s))+"\n");
		}
		bw.close();

		// Now print the logos
		File motifLogos = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"motif_logos");
		motifLogos.mkdirs();
		// Finally, draw the motif logos
		for(WeightMatrix fm : seqConfig.getDiscrimMotifs()){
			File motifFileName = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"motif_logos"+File.separator+fm.getName()+".png");
			Utils.printMotifLogo(fm, motifFileName, 150, fm.getName(), true);
			File motifFileNameRC = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"motif_logos"+File.separator+fm.getName()+"_rc.png");
			Utils.printMotifLogo(WeightMatrix.reverseComplement(fm), motifFileNameRC, 150, fm.getName(), true);
		}
	}

	public void writeClssifierOutput() throws IOException{
		File classifierOutFile = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"classifier.out");
		FileWriter fw = new FileWriter(classifierOutFile);
		BufferedWriter bw= new BufferedWriter(fw);
		bw.write(seqConfig.getClassifierOutput());
		bw.close();
	}
	
	public void writeDiscrimScore() throws IOException{
		File discrimScoreFile = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"Discrim_motifs.scores");
		FileWriter fw = new FileWriter(discrimScoreFile);
		BufferedWriter bw = new BufferedWriter(fw);
		StringBuilder sb = new StringBuilder();
		sb.append("MotifName");sb.append("\t");
		//First header
		for(String s : seqConfig.getMNames()){
			sb.append(s.replaceAll("#", ""));sb.append("\t");
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append("\n");
		for(String s : seqConfig.getDiscrimMotsScore().keySet()){
			sb.append(s.replaceAll("#", ""));sb.append("\t");
			for(double scr : seqConfig.getDiscrimMotsScore().get(s)){
				sb.append(scr);sb.append("\t");
			}
			sb.deleteCharAt(sb.length()-1);
			sb.append("\n");
		}
		bw.write(sb.toString());
		bw.close();
		
		// Now simpler version
		File discrimScoreFileSimple = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"Discrim_motifs_simple.scores");
		fw = new FileWriter(discrimScoreFileSimple);
		bw = new BufferedWriter(fw);
		sb = new StringBuilder();
		//First header
		for(String s:seqConfig.getMNames()){
			String[] pieces = s.split("#");
			if(pieces.length < 2){
				sb.append(s.replaceAll("#", ""));sb.append("\t");
			}
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append("\n");
		for(String s : seqConfig.getDiscrimMotsScore().keySet()){
			String[] piecesMotName = s.split("#");
			if(piecesMotName.length>1)
				continue;
			sb.append(s.replaceAll("#", ""));sb.append("\t");
			int i=0;
			for(double scr : seqConfig.getDiscrimMotsScore().get(s)){
				String[] piecesModName = seqConfig.getMNames().get(i).split("#");
				i++;
				if(piecesModName.length>1)
					continue;
				sb.append(scr);sb.append("\t");
			}
			sb.deleteCharAt(sb.length()-1);
			sb.append("\n");
		}
		bw.write(sb.toString());
		bw.close();	
	}
	
	public void makeDiscrimHeatmaps() throws IOException, InterruptedException{
		//First make the Rscript to plot heatmaps
		File rScript = new File(seqConfig.getOutDir().getAbsolutePath()+File.separator+"plotHeatmap.R");
		FileWriter fw = new FileWriter(rScript);
		BufferedWriter bw = new BufferedWriter(fw);
		bw.write("library(gplots) \n");
		bw.write("args <- commandArgs(TRUE) \n");
		bw.write("data <- as.matrix(read.table(args[1], header=TRUE, sep = \"\t\", row.names = 1,as.is=TRUE)) \n");
		bw.write("colors = c(seq(0,0.5,length=6),seq(0.51,3,length=6)) \n");
		bw.write("my_palette <- colorRampPalette(c(\"lightblue\",\"red\"))(n = 11) \n");
		bw.write("png(file=args[2],width=600, height=800, bg = \"transparent\") \n");
		bw.write("heatmap.2(data,Rowv=FALSE,Colv=FALSE,breaks=colors,col=my_palette,symm=F,symkey=F,symbreaks=T,dendrogram=\"none\",trace=\"none\") \n");
		bw.write("dev.off()");
		bw.close();
		
		// Now run the R script
		String Rscriptcmd = seqConfig.getRpath()+"Rscript ";
		Process proc = Runtime.getRuntime().exec(Rscriptcmd+" "+seqConfig.getOutDir().getAbsolutePath()+"/plotHeatmap.R"+" Discrim_motifs.scores"+" Discrim_motifs_heatmap.png");
		// any error message? 
		StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "R_ERR", true); 
		// any output? 
		StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "R_OUT", true); 
		// kick them off 
		errorGobbler.start(); 
		outputGobbler.start(); 
		int exitVal = proc.waitFor(); 
		System.err.println("R ExitValue: " + exitVal);
		proc.destroy();
		
		proc = Runtime.getRuntime().exec(Rscriptcmd+" "+seqConfig.getOutDir().getAbsolutePath()+"/plotHeatmap.R"+" Discrim_motifs_simple.scores"+" Discrim_motifs_simple_heatmap.png");
		// any error message? 
		errorGobbler = new StreamGobbler(proc.getErrorStream(), "R_ERR", true); 
		// any output? 
		outputGobbler = new StreamGobbler(proc.getInputStream(), "R_OUT", true); 
		// kick them off 
		errorGobbler.start(); 
		outputGobbler.start(); 
		exitVal = proc.waitFor(); 
		System.err.println("R ExitValue: " + exitVal);
		proc.destroy();
	}
	
	public void writeHTMLfile() throws IOException{
		// Write the HTML file
		File htmlOut = new File(seqConfig.getOutDir().getAbsoluteFile()+File.separator+"SeqUnwinder_results.html");
		FileWriter fout = new FileWriter(htmlOut);
    	fout.write("<html>\n" +
    			"\t<head><title>SeqUnwinder results ("+seqConfig.getOutbase()+")</title></head>\n" +
    			"\t<style type='text/css'>/* <![CDATA[ */ table, th{border-color: #600;border-style: solid;} td{border-color: #600;border-style: solid;} table{border-width: 0 0 1px 1px; border-spacing: 0;border-collapse: collapse;} th{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} td{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} /* ]]> */</style>\n" +
    			"\t<script language='javascript' type='text/javascript'><!--\nfunction motifpopitup(url) {	newwindow=window.open(url,'name','height=75');	if (window.focus) {newwindow.focus()}	return false;}// --></script>\n" +
    			"\t<script language='javascript' type='text/javascript'><!--\nfunction fullpopitup(url) {	newwindow=window.open(url,'name');	if (window.focus) {newwindow.focus()}	return false;}// --></script>\n" +
    			"\t<body>\n" +
    			"\t<h1>SeqUnwinder results ("+seqConfig.getOutbase()+")</h1>\n" +
    			"");
    	DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
    	Date date = new Date();
    	fout.write("\t<p>SeqUnwinder version "+SeqUnwinderConfig.version+" run completed on: "+dateFormat.format(date));
    	fout.write(" with arguments:\n "+seqConfig.getArgs()+"\n</p>\n");

    	// Write the error on training dataset
    	fout.write("\t<h2>Error on training data</h2>\n");
    	fout.write("\t<table>\n");
    	fout.write("\t\t<tr>" +
    			"\t\t<th>Class</th>\n" +
    			"\t\t<th>TP Rare</th>\n" +
    			"\t\t<th>FP Rare</th>\n" +
    			"\t\t<th>Precision</th>\n" +
    			"\t\t<th>Recall</th>\n" +
    			"\t\t<th>F-Measure</th>\n" +
    			"\t\t<th>MCC</th>\n" +
    			"\t\t<th>ROC Area</th>\n" +
    			"\t\t<th>PRC Area</th>\n");
    	fout.write("\t\t</tr>\n");
    	for(String s : seqConfig.getTrainSetStats().keySet()){
    		fout.write("\t\t<tr>" +
    				"\t\t<td>"+s+"</th>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[0])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[1])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[2])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[3])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[4])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[5])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[6])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTrainSetStats().get(s)[7])+"</td>\n");
    		fout.write("\t\t</tr>\n");
    	}
    	fout.write("\t</table>\n");

    	// Cross-validation error
    	fout.write("\t<h2>Stratified cross-validation</h2>\n");
    	fout.write("\t<table>\n");
    	fout.write("\t\t<tr>" +
    			"\t\t<th>Class</th>\n" +
    			"\t\t<th>TP Rare</th>\n" +
    			"\t\t<th>FP Rare</th>\n" +
    			"\t\t<th>Precision</th>\n" +
    			"\t\t<th>Recall</th>\n" +
    			"\t\t<th>F-Measure</th>\n" +
    			"\t\t<th>MCC</th>\n" +
    			"\t\t<th>ROC Area</th>\n" +
    			"\t\t<th>PRC Area</th>\n");
    	fout.write("\t\t</tr>\n");
    	for(String s : seqConfig.getTestSetStats().keySet()){
    		fout.write("\t\t<tr>" +
    				"\t\t<td>"+s+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[0])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[1])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[2])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[3])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[4])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[5])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[6])+"</td>\n" +
    				"\t\t<td>"+Double.toString(seqConfig.getTestSetStats().get(s)[7])+"</td>\n");
    		fout.write("\t\t</tr>\n");
    	}
    	fout.write("\t</table>\n");
    	
    	// Print the de novo identified motifs
    	fout.write("\t<h2>De novo identified motifs</h2>\n");
    	fout.write("\t<table>\n");
    	fout.write("\t\t<tr>" +
    			"\t\t<th>Class</th>\n" +
    			"\t\t<th>Motif</th>\n");
    	fout.write("\t\t</tr>\n");

    	for(String s: seqConfig.getMNames()){
    		fout.write("\t\t<tr>" +
    				"\t\t<td>"+s+"</td>\n");
    		List<String> selectedMotifNames = new ArrayList<String>();
    		for(String motName : seqConfig.getDiscrimMotifRocs().keySet()){
    			if(motName.contains(s+"_c"))
    				selectedMotifNames.add(motName);
    		}
    		if(selectedMotifNames.size()>0){
    			fout.write("\t\t<td>");
    			for(String selecMotName : selectedMotifNames){
    				fout.write("\t\t<img src='motif_logos/"+selecMotName+".png'><a href='#' onclick='return motifpopitup(\"motif_logos/"+selecMotName+"_rc.png\")'>rc</a>");
    				fout.write("<br>\n");
    			}
    			fout.write("</td>\n");
    		}else{
    			fout.write("\t\t<td>No motif found</td>\n");
    		}
    		fout.write("\t\t</tr>\n");
    	}
    	fout.write("\t</table>\n");
    	
    	// Print motif score heatmaps for both complete and simple versions
    	fout.write("\t<h2>Label specific score of identified de novo motifs</h2>\n");
    	fout.write("\t<table>\n");
    	fout.write("\t\t<tr>\n");
		fout.write("\t\t<td><a href='#' onclick='return fullpopitup(\"Discrim_motifs_heatmap.png\")'><img src='Discrim_motifs_heatmap.png' height='200'></a></td>\n");
		fout.write("\t\t<td><a href='#' onclick='return fullpopitup(\"Discrim_motifs_simple_heatmap.png\")'><img src='Discrim_motifs_simple_heatmap.png' height='200'></a></td>\n");
		fout.write("\t\t</tr>\n");
		fout.write("\t</table>\n");
		
		fout.write("\t<p><a href='Discrim_motifs.transfac'> Discriminative motifs.</a> Try inputting these motifs into <a href='http://www.benoslab.pitt.edu/stamp/'>STAMP</a> for validation.</p>\n");
		
		fout.write("\t</body>\n</html>\n");
    	fout.close();
	}

}
