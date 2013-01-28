package edu.psu.compbio.seqcode.projects.shaun;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipData;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipDataset;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.ExptNameVersion;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.SQLData;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class ChipIDReporter {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String [][][] experiments = new String[][][] { h3k79Versions, h3k4Versions, h3k27Versions, suz12Versions, ring1bVersions};
		try {
			Organism org = new Organism("Mus musculus");
			Genome gen = org.getGenome("mm8");
			
			for(int i=0; i<experiments.length; i++){
				for(int t=0; t<experiments[i].length; t++){
					String [] exp = experiments[i][t];
					ChipChipDataset ccset = new ChipChipDataset(gen);
			        ExptNameVersion expt1=null;
			        if(exp.length==2)
			        	expt1= new ExptNameVersion(exp[0], exp[1]);
			        else if(exp.length==3)
			        	expt1= new ExptNameVersion(exp[0], exp[1], exp[2]);
			        
			        SQLData data1 = (SQLData)ccset.getData(expt1);
			        data1.window("1",100, 10000000);
			        int ID = data1.getExptID(1, 0);
			        System.out.println(exp[0]+"\t"+exp[1]+"\t"+exp[2]+"\t"+ID);
			        
				}
			}
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

	}
	public static String[][] h3k4Versions, h3k27Versions, h3k79Versions, ring1bVersions, suz12Versions;
	static { 
        h3k4Versions = new String[][] { 
                { "Mm H3K4me3:HBG3:ES Stage vs H3:HBG3:ES Stage", "median linefit, quantile norm", "1" },
                //{ "Mm H3K4me3:HBG3:ES Stage vs H3:HBG3:ES Stage", "median linefit, quantile norm", "2" },
                { "Mm H3K4me3:HBG3:ES Stage vs H3:HBG3:ES Stage", "median linefit, quantile norm", "3" },
                { "Mm H3K4me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA", "median linefit, quantile norm", "1" },
                { "Mm H3K4me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA", "median linefit, quantile norm", "2" },
                { "Mm H3K4me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm", "1" },
                { "Mm H3K4me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm", "2 (2/1/08)" },
                { "Mm H3K4me3:HBG3:2+1 day vs H3:HBG3:2+1 day", "median linefit, quantile norm", "1"},
                //{ "Mm H3K4me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm", "3" },
                { "Mm H3K4me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm", "chip1" },
                { "Mm H3K4me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm", "chip2" },
                { "Mm H3K4me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm", "1" },
                { "Mm H3K4me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm", "2" },
                //{ "Mm H3K4me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm", "3" }
        };

        h3k27Versions = new String[][] { 
                { "Mm H3K27me3:HBG3:ES Stage vs H3:HBG3:ES Stage", "median linefit","1" },
                { "Mm H3K27me3:HBG3:ES Stage vs H3:HBG3:ES Stage", "median linefit","2" },
                { "Mm H3K27me3:HBG3:ES Stage vs H3:HBG3:ES Stage", "median linefit","3" },
                { "Mm H3K27me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA", "median linefit", "1" },
                { "Mm H3K27me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA", "median linefit", "2" },
                { "Mm H3K27me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA", "median linefit", "1" },
                { "Mm H3K27me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA", "median linefit", "2 (2/1/08)" },
                { "Mm H3K27me3:HBG3:2+1 day vs H3:HBG3:2+1 day", "median linefit", "1"},
                { "Mm H3K27me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit","chip1" },
                { "Mm H3K27me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit","chip2" },
                //{ "Mm H3K27me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm","3" },
                { "Mm H3K27me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm", "1" },
                //{ "Mm H3K27me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm", "2" }
        }; 
        
        h3k79Versions = new String[][] { 
                { "Mm H3K79me2:HBG3:mES vs H3:HBG3:mES", "median linefit, quantile norm 6tp", "1" },
                { "Mm H3K79me2:HBG3:mES vs H3:HBG3:mES", "median linefit, quantile norm 6tp","2 (10/26/07, Hox Array)" },
                { "Mm H3K79me2:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage", "median linefit, quantile norm 6tp", "1 (5/20/08)" },
                { "Mm H3K79me2:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm 6tp",  "1 (5/20/08)"},
                { "Mm H3K79me2:HBG3:2+1 day vs WCE:HBG3:2+1 day", "median linefit, quantile norm 6tp", "1 (5/20/08)" },
                { "Mm H3K79me2:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm 6tp", "1" },
                { "Mm H3K79me2:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm 6tp", "2 (10/26/07, Hox Array)" },
                { "Mm H3K79me2:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm 6tp", "1" },
                { "Mm H3K79me2:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm 6tp", "2 (10/26/07, Hox Array)" }
        };
        ring1bVersions = new String[][] {
        		{"Mm Ring1b:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage", "median linefit, quantile norm", "1 (7/24/08)"},	
        		{"Mm Ring1b:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA","median linefit, quantile norm", "1 (7/24/08)"},	
        		{"Mm Ring1b:HBG3:2+1 day vs WCE:HBG3:2+1 day", "median linefit, quantile norm", "1 (7/24/08)"}        		
        };
        suz12Versions = new String[][] {
        		{"Mm Suz12:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage","median linefit", "1 (5/20/08)"},	
        		{"Mm Suz12:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA","median linefit", "1 (5/20/08)"},	
        		{"Mm Suz12:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA","median linefit", "2 (7/24/08)"},
        		{"Mm Suz12:HBG3:2+1 day vs WCE:HBG3:2+1 day","median linefit", "1 (5/20/08)"},
        		{"Mm Suz12:HBG3:2+1 day vs WCE:HBG3:2+1 day","median linefit", "2 (7/24/08)"}
        };
	}
}
