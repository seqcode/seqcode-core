package edu.psu.compbio.seqcode.projects.shaun.rnaseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.shaun.rnaseq.genemodels.GeneTUnit;
import edu.psu.compbio.seqcode.projects.shaun.rnaseq.genemodels.SplicedTUnit;
import edu.psu.compbio.seqcode.projects.shaun.rnaseq.genemodels.TUnit;


public class GTFReader {
	private File inFile;
	private Genome gen;
	public GTFReader(File f, Genome g){
		inFile=f;
		gen=g;
	}
	
	public List<GeneTUnit> loadGenes(){
		HashMap<String, GeneTUnit> genes = new HashMap<String, GeneTUnit>();
		BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(inFile));
			String line;
			while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	        	if(line.charAt(0)!='#'){
		            String[] words = line.split("\\t");
		            String chr = words[0];
		            String[] tmp = chr.split("\\.");
	            	chr=tmp[0].replaceFirst("chr", "");
	            	chr=chr.replaceFirst("^>", "");
	            	if(chr.equals("MT")){chr="M";}
	            	if(gen.containsChromName(chr)){//Ignore chromosomes not defined in genome
	            		
		            	String gType = words[1];
		            	String unitType = words[2];
		            	int start = new Integer(words[3]).intValue();
		            	int stop = new Integer(words[4]).intValue();
		            	char strand = words[6].charAt(0);
		            	String[] tags = words[8].split(";");
		            	String gName="", gID="", tName="", tID="";
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
		            	}
		            	
		            	Region currR = new Region(gen, chr, start, stop);
		            	TUnit tu = new TUnit(currR, strand);  
		            	if(unitType.equals("exon") || unitType.equals("CDS")){
		            		if(!genes.containsKey(gID))
		            			genes.put(gID, new GeneTUnit(gName, gID, gType, currR, strand));
		            		if(!genes.get(gID).hasIsoform(tID))
		            			genes.get(gID).addIsoform(new SplicedTUnit(tName, tID, gType, currR, strand));
		            		
		            		if(unitType.equals("exon"))
		            			genes.get(gID).addExon(tID, tu);
	            			else if(unitType.equals("CDS"))
	            				genes.get(gID).addCDS(tID, tu);
		            	}
	            	}
	        	}
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
		List<GeneTUnit> geneList = new ArrayList<GeneTUnit>();
		geneList.addAll(genes.values());
		return geneList;
	}
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("gtf")) { 
            System.err.println("Usage:\n" +
                               "  --species <species;genome>\n" +
                               "  --gtf <GTF file>\n" +
                               "  --tss [print all TSS coordinates]\n" +
                               "  --singleisotss [print single isoform gene TSSs]\n"+
                               "  --juncs [print all juncs]\n"+
                               "  --junctions [print all junctions with gene names]\n" +
                               "  --test [test the loaded annotation]\n"+
                               "");
            return;
        }
        try {
			Pair<Species, Genome> pair = Args.parseGenome(args);
			Genome currgen = pair.cdr();
			String gtfFile = ap.getKeyValue("gtf");			
			
			//Load genes
			GTFReader reader = new GTFReader(new File(gtfFile), currgen);
			List<GeneTUnit> genes = reader.loadGenes();
			
			
			if(Args.parseFlags(args).contains("tss")){
				System.out.println("#TSS-coord\tGeneName\tGeneID\tGeneType\tIsoformName\tIsoformID");
				for(GeneTUnit g : genes){
					for(SplicedTUnit i : g.getIsoforms()){
					    if(i.getTSS()!=null)
						System.out.println(i.getTSS().getLocationString()+":"+i.getStrand()+"\t"+g.getName()+"\t"+g.getID()+"\t"+g.getType()+"\t"+i.getName()+"\t"+i.getID());
					}
				}
			}
			if(Args.parseFlags(args).contains("singleisotss")){
				System.out.println("#TSS-coord\tGeneName\tGeneID\tGeneType\tIsoformName\tIsoformID");
				for(GeneTUnit g : genes){
					if(g.getIsoforms().size()==1){
						for(SplicedTUnit i : g.getIsoforms()){
						    if(i.getTSS()!=null)
							System.out.println(i.getTSS().getLocationString()+":"+i.getStrand()+"\t"+g.getName()+"\t"+g.getID()+"\t"+g.getType()+"\t"+i.getName()+"\t"+i.getID());
						}
					}
				}
			}
			if(Args.parseFlags(args).contains("juncs")){
				for(GeneTUnit g : genes){
					for(SplicedTUnit i : g.getIsoforms()){
						Iterator<TUnit> inIter = i.getIntronIterator();
						while(inIter.hasNext()){
							TUnit intron = inIter.next();
							System.out.println(intron.getCoords().getChrom()+"\t"+intron.getCoords().getStart()+"\t"+intron.getCoords().getEnd()+"\t"+intron.getStrand());
						}
					}
				}
			}
			if(Args.parseFlags(args).contains("junctions")){
				for(GeneTUnit g : genes){
					for(SplicedTUnit i : g.getIsoforms()){
						Iterator<TUnit> inIter = i.getIntronIterator();
						while(inIter.hasNext()){
							TUnit intron = inIter.next();
							System.out.println(intron.getCoords()+":"+intron.getStrand()+"\t"+g.getName()+"\t"+g.getID()+"\t"+i.getID());
						}
					}
				}
			}
			
			if(Args.parseFlags(args).contains("test")){
				//Test loaded genes
				for(GeneTUnit g : genes){
					System.out.println(g.getName()+"\t"+g.getID()+"\t"+g.getType()+"\t"+g.getCoords().getLocationString()+":"+g.getStrand()+"\t"+g.getNumExons()+"\t"+g.getExonLength()+"\t"+g.getNumIsoforms());
					for(SplicedTUnit i : g.getIsoforms()){
						//System.out.println("\t"+i.getName()+"\t"+i.getID()+"\t"+i.getCoords().getLocationString()+":"+i.getStrand()+"\t"+i.getNumExons()+"\t"+i.getExonLength());
					}
				}
			}
        } catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
