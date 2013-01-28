package edu.psu.compbio.seqcode.projects.shaun.rnaseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.shaun.rnaseq.genemodels.GeneTUnit;

/*
 * Reads a differential expression table, appends the binding site that is closest to any of a gene's TSSs 
 */
public class AddBindingToExprTable {

	protected Genome gen;
	protected GTFReader gtfReader;
	protected HashMap<String, GeneTUnit> knownGenes = new HashMap<String, GeneTUnit>();
	protected HashMap<String, String> exprTable = new HashMap<String, String>();
	protected ArrayList<PeakItem> peaks = new ArrayList<PeakItem>(); 
	
	public AddBindingToExprTable(Genome gen, String gtfFile, String exprFile, int exprIDCol, String bindingFile, String bindingType){
		this.gen=gen;
		try {
				
			//Load genes
			gtfReader = new GTFReader(new File(gtfFile), gen);
			for(GeneTUnit g : gtfReader.loadGenes()){
				knownGenes.put(g.getID(), g);
			}
			System.err.println("Loaded: "+knownGenes.size()+" known genes");
			
			//Load expression table
			BufferedReader ereader = new BufferedReader(new FileReader(new File(exprFile)));
			String line;
			while ((line = ereader.readLine()) != null) {
	        	line = line.trim();
	            String[] words = line.split("\\t");
	            String[] tmp = words[exprIDCol].split(":");
	            String ID = tmp[0];
	            exprTable.put(ID, line);
			}
			System.err.println("Loaded: "+exprTable.size()+" expression entries");
			
			//Load peaks
			BufferedReader preader = new BufferedReader(new FileReader(new File(bindingFile)));
			line = preader.readLine(); //Skip first line
			while ((line = preader.readLine()) != null) {
	        	line = line.trim();
	            String[] words = line.split("\\t");
	            String ID = words[exprIDCol];
	            Point currCoord=null;
	            double currIP=0.0; double currCtrl=0.0;
	            double currP = 0.0;
	            if(bindingType.equals("STAT")){
	            	String scoord = words[2];
	            	String[] pieces = scoord.split(":");
	            	if(pieces.length==2){
	            		currCoord = new Point(gen, pieces[0], new Integer(pieces[1]).intValue());
	            	}
	            	currIP = new Double(words[7]).doubleValue();
	            	currCtrl = new Double(words[8]).doubleValue();
	            	currP = new Double(words[6]).doubleValue();
	            }else if(bindingType.equals("GPS")){
	            	String scoord = words[0];
	            	String[] pieces = scoord.split(":");
	            	if(pieces.length==2){
	            		currCoord = new Point(gen, pieces[0], new Integer(pieces[1]).intValue());
	            	}
	            	currIP = new Double(words[1]).doubleValue();
	            	currCtrl = new Double(words[2]).doubleValue();
	            	currP = new Double(words[5]).doubleValue();
	            }
	            if(currCoord!=null){
	            	peaks.add(new PeakItem(currCoord, currIP, currCtrl, currP));
	            }
			}
			System.err.println("Loaded: "+peaks.size()+" peaks");
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		
	}
	
	public void execute(){
		for(String id : exprTable.keySet()){
			
			if(knownGenes.containsKey(id)){
				GeneTUnit currGene = knownGenes.get(id);
				Collection<Point> gTSSs = currGene.getTSSs();
				int minDist = 100000000;
				PeakItem closestPeak=null;
				Point closestTSS=null;
				for(PeakItem p : peaks){
					if(p.coord.getChrom().equals(currGene.getCoords().getChrom())){
						for(Point tss : gTSSs)
							if(tss.distance(p.coord)<minDist){
								minDist = tss.distance(p.coord);
								closestPeak=p;
								closestTSS = tss;
							}
					}
				}
				if(closestPeak==null || closestTSS==null)
					System.out.println(exprTable.get(id)+"\tNOTFOUND");
				else
					System.out.println(exprTable.get(id)+"\t"+closestTSS.getLocationString()+"\t"+minDist+"\t"+closestPeak.toString());
			}else{
				System.out.println(exprTable.get(id)+"\tGENENOTKNOWN");
			}
		}
	}
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("gtf") || !ap.hasKey("exprtable") ||(!ap.hasKey("statpeak")&&!ap.hasKey("gps"))) { 
            System.err.println("Usage:\n " +
                               "  --species <species;genome>\n" +
                               "  --gtf <GTF file>\n"+
                               "  --exprtable <expression table>" +
                               "  --exprid <expression ID column>\n"+
                               "  --statpeak <stat peak calls>\n"+
                               "  --gps <GPS peak calls>\n" +
                               "  Option:\n" +
                               "  --exprid <ID column in expression table>\n"+
                               "");
            return;
        }
        try {
			Pair<Organism, Genome> pair = Args.parseGenome(args);
			Genome currgen = pair.cdr();
			String gtfFile = ap.getKeyValue("gtf");			
			String exprFile = ap.getKeyValue("exprtable");
			String peakFile=""; String peakType="STAT";
			if(ap.hasKey("statpeak")){
				peakFile = ap.getKeyValue("statpeak");
				peakType="STAT";
			}else if(ap.hasKey("gps")){
				peakFile = ap.getKeyValue("gps");
				peakType="GPS";
			}
			int exprIDCol = 1;
			if(ap.hasKey("exprid")){
				exprIDCol = new Integer(ap.getKeyValue("exprid")).intValue();
			}
			AddBindingToExprTable adder = new AddBindingToExprTable(currgen, gtfFile, exprFile, exprIDCol, peakFile, peakType);
			adder.execute();
        } catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public class PeakItem{
		public Point coord;
		public double ipHits;
		public double ctrlHits;
		public double pval;
		
		public PeakItem(Point p, double ip, double ctrl, double pv){
			coord=p; ipHits=ip; ctrlHits=ctrl; pval=pv;
		}
		
		public String toString(){
			return(String.format("%s\t%.1f\t%.1f\t%e", coord.getLocationString(), ipHits, ctrlHits, pval));
		}
	}
}
