package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipData;
import edu.psu.compbio.seqcode.gse.datasets.locators.ChipChipLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Gene;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RefGeneGenerator;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class PrintMetaGenes {
	private static double [][] metagene;
	private static double [][] metacount;
	private static String [][] experiment;
	private static Organism org;
	private static Genome gen;
	private static int windowSize=10000;
	private static int smoothWin=500;
	private static int smoothOff=250;
	private static int numBins;
	private static int numTotalGenes=0;
	private static int numFoundGenes=0;
	
	public PrintMetaGenes(String [][] expts, String species, String genome){
		numBins = (windowSize/smoothOff)-1;
		metagene = new double[expts.length][numBins];
		metacount = new double[expts.length][numBins];
		for(int i=0; i<expts.length; i++){
			for(int j=0; j<numBins; j++){metagene[i][j]=0;metacount[i][j]=0;}
		}
		
		experiment = expts;
		try {
			org = Organism.getOrganism(species);
			gen = org.getGenome(genome);
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void calcMetagenes(String [] genelist){
		numTotalGenes = genelist.length;
		RefGeneGenerator refgen = new RefGeneGenerator(gen);
		for(int x=0; x<experiment.length; x++){
			ChipChipLocator loc =	new ChipChipLocator(gen, experiment[x][0], experiment[x][1]);
			ChipChipData data = loc.createObject();
			
			System.out.println(experiment[x][2]);
			for(int g=0; g<genelist.length; g++){
				Iterator<Gene> iter = refgen.byName(genelist[g]);
				if(iter.hasNext()){
					if(x==0){numFoundGenes++;}
					Gene head = iter.next();
					System.out.println(head.toString());
					String c = head.getChrom();
					int st, ed;
					if(head.getStrand()=='+'){
						st = head.getStart()-(windowSize/2);
						ed = head.getStart()+(windowSize/2);
					}else{
						st = head.getEnd()-(windowSize/2);
						ed = head.getEnd()+(windowSize/2);
					}
					
					try {
						data.window(c, st, ed);
						
						double[] array = new double[data.getCount()];
						int[] positions = new int[data.getCount()];
						
						for(int i = 0; i < data.getCount(); i++) {
							positions[i] = data.getPos(i);
							array[i] = getMeanRatio(data, i);	
						}
						
						//Data loaded, now smooth over overlapping windows of sequence
						int l=0; 
						for(int i = st; i <= ed-smoothWin; i+=smoothOff) {
							double count =0;
							double smoothVal=0;
							for(int j=0; j<array.length; j++){
								if(positions[j]>=i &&  positions[j]<=i+smoothWin){
									smoothVal+=array[j];
									count+=1;
								}							
							}
							if(count>=1){
								if(head.getStrand()=='+'){
									metagene[x][l]+=smoothVal/count;
									metacount[x][l]++;
								}else{
									metagene[x][numBins-l-1]+=smoothVal/count;
									metacount[x][numBins-l-1]++;
								}
							}
							l++;
						}
					} catch (NotFoundException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			for(int j=0; j<numBins; j++){if(metacount[x][j]>0){metagene[x][j]=metagene[x][j]/metacount[x][j];}}
		}
	}
	
	public void printMetagenes(){
		System.out.println(numFoundGenes+" of "+numTotalGenes+" listed genes included in metagene.");
		System.out.print("Experiment");
		int pos = (-1*(windowSize/2))+(smoothWin/2);
		for(int j=0; j<numBins; j++){
			System.out.print("\t"+pos);
			pos+=smoothOff;
		}System.out.print("\n");
		for(int i=0; i<metagene.length; i++){
			System.out.print(experiment[i][2]);
			for(int j=0; j<numBins; j++){
				System.out.print("\t"+metagene[i][j]);
			}System.out.print("\n");
		}
	}
	
	public double getMeanRatio(ChipChipData d, int i) { 
		double sum = 0.0;
		double count = 0;
		for(int j = 0; j < d.getReplicates(i); j++) { 
			if(!Double.isNaN(d.getRatio(i, j))){
				sum += d.getRatio(i, j);
				count += 1;
			}
		}
		return count > 0 ? sum / count : 0.0;
	}
	
	public static String[] loadGeneFile(String filename){
		ArrayList<String> geneL = new ArrayList<String>();
		try {
			File rFile = new File(filename);
			if(rFile.isFile()){
				BufferedReader reader = new BufferedReader(new FileReader(rFile));
				String line;
				while((line= reader.readLine())!=null){
					String [] tokens = line.split("[\\s*\\t\\r\\n\\f]");
					geneL.add(tokens[0]);
				}
				reader.close();
			}else{
				System.err.println("Gene file not valid");
				System.exit(0);
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		String[] geneS = new String[geneL.size()];
		int g=0;
		Iterator<String> it = geneL.iterator();
		while(it.hasNext()){
			String n = it.next();
			geneS[g]=n;
			g++;
		}
		return(geneS);
	}
	
	private static String [][] h3k4Versions = new String[][] { 
            { "Mm H3K4me3:HBG3:ES Stage vs H3:HBG3:ES Stage", "median linefit, quantile norm", "ES" },
            { "Mm H3K4me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA", "median linefit, quantile norm", "ES+2" },
            { "Mm H3K4me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm", "RA+8hrs" },
            { "Mm H3K4me3:HBG3:2+1 day vs H3:HBG3:2+1 day", "median linefit, quantile norm", "RA+1d" },
            { "Mm H3K4me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm", "Olig2" },
            { "Mm H3K4me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm", "Hb9" }
    };

	private static String [][] h3k27Versions = new String[][] { 
            { "Mm H3K27me3:HBG3:ES Stage vs H3:HBG3:ES Stage", "median linefit, quantile norm", "ES" },
            { "Mm H3K27me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA", "median linefit, quantile norm", "ES+2" },
            { "Mm H3K27me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm", "RA+8hrs" },
            { "Mm H3K27me3:HBG3:2+1 day vs H3:HBG3:2+1 day", "median linefit, quantile norm", "RA+1d" },
            { "Mm H3K27me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm", "Olig2" },
            { "Mm H3K27me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm", "Hb9" }
    }; 
    
	private static String [][] h3k79Versions = new String[][] { 
            { "Mm H3K79me2:HBG3:mES vs H3:HBG3:mES", "median linefit, quantile norm", "ES" },
            { "Mm H3K79me2:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm", "Olig2" },
            { "Mm H3K79me2:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm", "Hb9" }
    };
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("genome")||!ap.hasKey("expt")||!ap.hasKey("genes")) { 
            System.err.println("Usage:\n" +
                               "PrintMetaGenes" +
                               "--species <organism name> " +
                               "--genome <genome version> "+
                               "--expt <expt name> "+
                               "--genes <file name> ");
            return;
        }
        String species = ap.getKeyValue("species");
        String genome = ap.getKeyValue("genome");
        String exptName = ap.getKeyValue("expt");
        String geneFile = ap.getKeyValue("genes");
        
        String [][] currExpt=h3k4Versions;
        if(exptName.equals("K4")){
        	currExpt = h3k4Versions;
        }else if(exptName.equals("K27")){
        	currExpt = h3k27Versions;
        }else if(exptName.equals("K79")){
        	currExpt = h3k79Versions;
        }
        
        String[] genelist;// = new String[]{"rarb", "olig2", "hoxa1"};
        genelist = loadGeneFile(geneFile);
        
        PrintMetaGenes pmg = new PrintMetaGenes(currExpt, species, genome);
        
        pmg.calcMetagenes(genelist);
        pmg.printMetagenes();
	}
}

