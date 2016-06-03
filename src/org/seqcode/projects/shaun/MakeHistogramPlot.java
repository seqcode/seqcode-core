package org.seqcode.projects.shaun;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.Random;

import javax.swing.JFrame;

import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.NotFoundException;
import org.seqcode.gse.utils.RealValuedHistogram;
import org.seqcode.gse.viz.paintable.PaintableFrame;
import org.seqcode.projects.shaun.viz.HistogramPaintable;


public class MakeHistogramPlot {
	private HistogramPaintable painter;
	private PaintableFrame plotter;
	private int screenSizeX=1000, screenSizeY=900;
	private String inputFile;
	private Species org;
	private Genome gen;
	private final boolean display = true;
	private final boolean useBackground = true;

	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("in")) { 
            System.err.println("Usage:\n " +
                               "MakeHistogramPlot " +
                               "--in <file name> ");
            return;
        }
        String dfile = ap.getKeyValue("in");
		MakeHistogramPlot maker = new MakeHistogramPlot(dfile);		
	}
	
	public MakeHistogramPlot(String df){
		try {
			org = Species.getSpecies("Mus musculus");
			gen = new Genome(org, "mm8");
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		inputFile = df;
		RealValuedHistogram histo, back;
		
		//Background load
		if(useBackground){
			ChipSeqMetaPeakMaker bm = new ChipSeqMetaPeakMaker(org, gen, 5000, 500);		
			bm.loadPeaks(inputFile); 
			back = bm.execute(ppg_back);
		}
		
		//Experiment load
		for(String[] experiment : ppg_exp){ 
			ChipSeqMetaPeakMaker hm = new ChipSeqMetaPeakMaker(org, gen, 5000, 500);		
			hm.loadPeaks(inputFile);
			String [][] ex = new String[1][];
			ex[0]=experiment;
			histo = hm.execute(ex);
			
			if(useBackground){
				histo.subtractHistogram(back);
			}
			
			/*histo = new RealValuedHistogram(-500, 500, 50);
			Random rgen = new Random(); 
			for(int i=0; i<100000; i++){
				histo.addValue(rgen.nextGaussian()*100, 0.0001);
			}//histo.printContents();
			*/
			
			String title = new String(experiment[0]+"_"+experiment[1]);
			String saveFile = new String(title+".jpg");
			String saveFile2 = new String(title+".svg");
			painter = new HistogramPaintable(histo, title, Color.blue);
			//painter.zeroTheY(true);
			try {
				painter.saveImage(new File(saveFile2), 800, 800, false);
				painter.saveImage(new File(saveFile), 800, 800, true);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			if(display){
				plotter = new PaintableFrame("Histogram", painter);
				plotter.setSize(screenSizeX, screenSizeY);
				plotter.setVisible(true);
				plotter.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			}
		}
	}
	
	private static String [][] ppg_back, yl_back, mik_back, ppg_exp, yl_exp, mik_exp, sigova, rarChipSeq, backChipSeq, chipSeqExperiments;	
	static {
		ppg_back = new String[][]{
				{"PPG_Solexa_WCE_2+1", "ELAND_unique"},
				{"PPG_Solexa_WCE_ES+2d", "ELAND_unique"}
		};
		ppg_exp = new String[][]{
				{"PPG_Solexa_RAR_8hr", "ELAND_unique"},
				{"PPG_Solexa_RAR_2+1", "ELAND_unique"},
				{"PPG_Solexa_RAR_ES+2d", "ELAND_unique"},
				{"PPG_Solexa_WCE_2+1", "ELAND_unique"},
				{"PPG_Solexa_WCE_ES+2d", "ELAND_unique"}
		};
		yl_back = new String [][]{
				{"YoungLab_Solexa_WCE", "ELAND_unique"}
		};
		yl_exp = new String [][]{
				{"YoungLab_Solexa_H3K36me3", "ELAND_unique"},
				{"YoungLab_Solexa_H3K4me3", "ELAND_unique"},
				{"YoungLab_Solexa_H3K79me2", "ELAND_unique"},
				{"YoungLab_Solexa_Nanog", "ELAND_unique"},
				{"YoungLab_Solexa_Oct4", "ELAND_unique"},
				{"YoungLab_Solexa_Sox2", "ELAND_unique"},
				{"YoungLab_Solexa_TCF3", "ELAND_unique"},
				{"YoungLab_Solexa_WCE", "ELAND_unique"}
		};
		mik_back=new String[][]{
				{"Mikkelsen07_Solexa_WCE_ES", "Unique"},
				{"Mikkelsen07_Solexa_WCE_MEF", "Unique"}
		};
		mik_exp=new String[][]{
				{"Mikkelsen07_Solexa_H3K20me3_ES", "Unique"},
				{"Mikkelsen07_Solexa_H3K27me3_ES", "Unique"},
				{"Mikkelsen07_Solexa_H3K27me3_MEF", "Unique"},
				{"Mikkelsen07_Solexa_H3K27me3_NP", "Unique"},
				{"Mikkelsen07_Solexa_H3K36me3_ES", "Unique"},
				{"Mikkelsen07_Solexa_H3K36me3_ESHyb", "Unique"},
				{"Mikkelsen07_Solexa_H3K36me3_MEF", "Unique"},
				{"Mikkelsen07_Solexa_H3K4me3_ES", "Unique"},
				{"Mikkelsen07_Solexa_H3K4me3_ESHyb", "Unique"},
				{"Mikkelsen07_Solexa_H3K4me3_MEF", "Unique"},
				{"Mikkelsen07_Solexa_H3K4me3_NP", "Unique"},
				{"Mikkelsen07_Solexa_H3K9me3_ES", "Unique"},
				{"Mikkelsen07_Solexa_H3K9me3_ESHyb", "Unique"},
				{"Mikkelsen07_Solexa_H3K9me3_MEF", "Unique"},
				{"Mikkelsen07_Solexa_H3K9me3_NP", "Unique"},
				{"Mikkelsen07_Solexa_panH3_ES", "Unique"},
				{"Mikkelsen07_Solexa_RNAPolII_ES", "Unique"},
				{"Mikkelsen07_Solexa_WCE_ES", "Unique"},
				{"Mikkelsen07_Solexa_WCE_MEF", "Unique"},
				{"Mikkelsen08_Solexa_H3K27me3_BIV1", "Unique"},
				{"Mikkelsen08_Solexa_H3K27me3_iPS", "Unique"},
				{"Mikkelsen08_Solexa_H3K27me3_MCV6", "Unique"},
				{"Mikkelsen08_Solexa_H3K27me3_MCV8", "Unique"},
				{"Mikkelsen08_Solexa_H3K27me3_MEF", "Unique"},
				{"Mikkelsen08_Solexa_H3K4me3_BIV1", "Unique"},
				{"Mikkelsen08_Solexa_H3K4me3_iPS", "Unique"},
				{"Mikkelsen08_Solexa_H3K4me3_MCV6", "Unique"},
				{"Mikkelsen08_Solexa_H3K4me3_MCV8", "Unique"},
				{"Mikkelsen08_Solexa_H3K4me3_MEF", "Unique"}
		};
		sigova=new String[][]{
				{"Sigova_Solexa_MES_Transcriptome_no_rRNA", "ELAND_iterative"},
				{"Sigova_Solexa_MES_Transcriptome_polyA_RNAase", "ELAND_iterative"},
				{"Sigova_Solexa_MES_Transcriptome_totalRNA_noRT", "ELAND_iterative"}
		};
		chipSeqExperiments = new String[][]{
				{"PPG_Solexa_RAR_2+1", "ELAND_unique"},
				{"PPG_Solexa_RAR_8hr", "ELAND_unique"},
				{"PPG_Solexa_RAR_ES+2d", "ELAND_unique"},
				{"PPG_Solexa_WCE_2+1", "ELAND_unique"},
				{"PPG_Solexa_WCE_ES+2d", "ELAND_unique"},
				{"YoungLab_Solexa_H3K36me3", "ELAND_unique"},
				{"YoungLab_Solexa_H3K4me3", "ELAND_unique"},
				{"YoungLab_Solexa_H3K79me2", "ELAND_unique"},
				{"YoungLab_Solexa_Nanog", "ELAND_unique"},
				{"YoungLab_Solexa_Oct4", "ELAND_unique"},
				{"YoungLab_Solexa_Sox2", "ELAND_unique"},
				{"YoungLab_Solexa_TCF3", "ELAND_unique"},
				{"YoungLab_Solexa_WCE", "ELAND_unique"},
				{"Mikkelsen07_Solexa_H3K20me3_ES", "Unique"},
				{"Mikkelsen07_Solexa_H3K27me3_ES", "Unique"},
				{"Mikkelsen07_Solexa_H3K27me3_MEF", "Unique"},
				{"Mikkelsen07_Solexa_H3K27me3_NP", "Unique"},
				{"Mikkelsen07_Solexa_H3K36me3_ES", "Unique"},
				{"Mikkelsen07_Solexa_H3K36me3_ESHyb", "Unique"},
				{"Mikkelsen07_Solexa_H3K36me3_MEF", "Unique"},
				{"Mikkelsen07_Solexa_H3K4me3_ES", "Unique"},
				{"Mikkelsen07_Solexa_H3K4me3_ESHyb", "Unique"},
				{"Mikkelsen07_Solexa_H3K4me3_MEF", "Unique"},
				{"Mikkelsen07_Solexa_H3K4me3_NP", "Unique"},
				{"Mikkelsen07_Solexa_H3K9me3_ES", "Unique"},
				{"Mikkelsen07_Solexa_H3K9me3_ESHyb", "Unique"},
				{"Mikkelsen07_Solexa_H3K9me3_MEF", "Unique"},
				{"Mikkelsen07_Solexa_H3K9me3_NP", "Unique"},
				{"Mikkelsen07_Solexa_panH3_ES", "Unique"},
				{"Mikkelsen07_Solexa_RNAPolII_ES", "Unique"},
				{"Mikkelsen07_Solexa_WCE_ES", "Unique"},
				{"Mikkelsen07_Solexa_WCE_MEF", "Unique"},
				{"Mikkelsen08_Solexa_H3K27me3_BIV1", "Unique"},
				{"Mikkelsen08_Solexa_H3K27me3_iPS", "Unique"},
				{"Mikkelsen08_Solexa_H3K27me3_MCV6", "Unique"},
				{"Mikkelsen08_Solexa_H3K27me3_MCV8", "Unique"},
				{"Mikkelsen08_Solexa_H3K27me3_MEF", "Unique"},
				{"Mikkelsen08_Solexa_H3K4me3_BIV1", "Unique"},
				{"Mikkelsen08_Solexa_H3K4me3_iPS", "Unique"},
				{"Mikkelsen08_Solexa_H3K4me3_MCV6", "Unique"},
				{"Mikkelsen08_Solexa_H3K4me3_MCV8", "Unique"},
				{"Mikkelsen08_Solexa_H3K4me3_MEF", "Unique"},
				{"Sigova_Solexa_MES_Transcriptome_no_rRNA", "ELAND_iterative"},
				{"Sigova_Solexa_MES_Transcriptome_polyA_RNAase", "ELAND_iterative"},
				{"Sigova_Solexa_MES_Transcriptome_totalRNA_noRT", "ELAND_iterative"}
		};
		rarChipSeq = new String[][]{
	    	//	{"PPG_Solexa_RAR_ES+2d", "ELAND_unique"},
        		{"PPG_Solexa_RAR_8hr", "ELAND_unique"}        		
        };
        backChipSeq = new String[][]{
        		{"PPG_Solexa_WCE_2+1", "ELAND_unique"},
        		{"PPG_Solexa_WCE_ES+2d", "ELAND_unique"}
        };
	}
}
