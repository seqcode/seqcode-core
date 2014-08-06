package edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Organism;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.GenomeLoader;

public class GTFAnnotationLoader extends AnnotationLoader{

	private File inFile;
	public GTFAnnotationLoader(File f, GenomeLoader gl){
		super(gl);
		inFile=f;
	}
	
	public HashMap<String, List<AGene>> loadGeneHashMap(Collection<AGene> allGenes){
		HashMap<String, List<AGene>> geneMap = new HashMap<String, List<AGene>>();
		
		//List<AGene> allGenes = loadGenes();
		//////
		for (AGene gene : allGenes){
			String chr = gene.getCoords().getChrom();
			if (geneMap.containsKey(chr)){
				geneMap.get(chr).add(gene);
			}
			else{
				List<AGene> newList = new ArrayList<AGene>();
				newList.add(gene);
				geneMap.put(chr, newList);
			}
		}
		
		return(geneMap);		
	}
	
	public List<AGene> loadGenes(){
		BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(inFile));
			String line;
			String activeGeneID="UNDEF";
			String activeTransID="UNDEF";
			AIsoform currTrans = null;
			AGene currGene = null;
			while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	        	if(line.charAt(0)!='#'){
		            String[] words = line.split("\\t");
		            String chr = words[0];
		            String[] tmp = chr.split("\\.");
		           	chr=tmp[0].replaceFirst("^chr", "");
		            chr=chr.replaceFirst(">", "");
		            
		            if(chr.equals("MT")){chr="M";}
	            	if(gen.containsChromName(chr)){//Ignore chromosomes not defined in genome
	            		
		            	String gType = words[1];
		            	String unitType = words[2];
		            	ARegionType aType = ARegion.translateToARegionType(unitType);
		            	int start = new Integer(words[3]).intValue();
		            	int stop = new Integer(words[4]).intValue();
		            	char strand = words[6].charAt(0);
		            	String[] tags = words[8].split(";");
		            	String gName="", gID="", tName="", tID="", eNum="", pID="";
		            	for(String x : tags){ 
		            		String[] pair = x.split("\\s+");
		            		int pl = pair.length;
		            		if(pl>0)
		            			pair[pl-1]=pair[pl-1].replaceAll("\"", "");
		            		if(pair[pl-2].startsWith("gene_id"))
		            			gID = pair[pl-1];
		            		else if(pair[pl-2].startsWith("transcript_id"))
		            			tID =pair[pl-1];
		            		else if(pair[pl-2].startsWith("gene_name"))
		            			gName = pair[pl-1];
		            		else if(pair[pl-2].startsWith("transcript_name"))
		            			tName = pair[pl-1];
		            		else if(pair[pl-2].startsWith("exon_number"))
		            			eNum = pair[pl-1];
		            		else if(pair[pl-2].startsWith("protein_id"))
		            			pID = pair[pl-1];
		            	}
		            	String exonID =tID+"_"+unitType+"_"+eNum;
		            	
		            	Region currR = new Region(gen, chr, start, stop);
		            	ARegion aReg = new ARegion(currR, strand, exonID, exonID, gType, aType);  
		            	if(aType != null){
			            	//Load the entry
			            	if(activeGeneID.equals(gID)){
			            		if(activeTransID.equals(tID)){
			            			
			            		}else{//New transcript, same gene
			            			if(currTrans!=null)
				            			currGene.addIsoform(currTrans);
			            			currTrans = new AIsoform(currR, strand, tName, tID, gType, ARegionType.isoform);
			            		}
			            	}else{//New gene, new transcript
			            		if(currTrans!=null && currGene!=null)
		            				currGene.addIsoform(currTrans);
			            		if(currGene!=null)
			            			genes.add(currGene);
			            		currGene = new AGene(currR, strand, gName, gID, gType, ARegionType.gene);
			            		currTrans = new AIsoform(currR, strand, tName, tID, gType, ARegionType.isoform);
			            	}
			            	currTrans.addPart(aReg);
			            	activeGeneID = gID;
			            	activeTransID = tID;
		            	}
	            	}
	        	}
			}
			if(currTrans!=null && currGene!=null)
				currGene.addIsoform(currTrans);
			if(currGene!=null){
				currGene.sortIsoforms();
    			genes.add(currGene);
			}
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.err.println(genes.size()+" genes loaded");
		Collections.sort(genes);
		return(genes);
	}
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("gtf")) { 
            System.err.println("Usage:\n " +
                               "  --species <species;genome>\n" +
                               " OR\n" +
                               "  --fai <genome index>\n" +
                               " OR\n" +
                               "  --seq <FASTA file>\n\n" +
                               "  --gtf <GTF file>\n"+
                               "");
            return;
        }
        try {
        	GenomeLoader gl=null;
        	if(ap.hasKey("species")){
        		Pair<Organism, Genome> pair = Args.parseGenome(args);
        		Genome currgen = pair.cdr();
        		gl = new GenomeLoader(currgen);
        	}else if(ap.hasKey("fai")){
        		gl = new GenomeLoader(new File(ap.getKeyValue("fai")), false);
        	}else if(ap.hasKey("seq")){
        		gl = new GenomeLoader(new File(ap.getKeyValue("seq")), true);
        	}
			
			String gtfFile = ap.getKeyValue("gtf");			
			
			//Load genes
			GTFAnnotationLoader reader = new GTFAnnotationLoader(new File(gtfFile), gl);
			Collection<AGene> geneSet = reader.loadGenes();
						
			//Test loaded genes
			/*for(AGene g : geneSet){
				System.out.println(g.toString());
				for(AIsoform i : g.getIsoforms())
					System.out.println("\t"+i.toString());
			}*/
			
			reader.geneLengths();
			reader.exonExonDistances();
        } catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	

}
