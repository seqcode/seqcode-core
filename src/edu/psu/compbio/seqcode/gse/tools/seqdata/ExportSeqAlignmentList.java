package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.sql.SQLException;
import java.util.HashMap;

import edu.psu.compbio.seqcode.gse.datasets.core.MetadataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExpt;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;

/**
 * Export existing SeqExpt & SeqAlignment descriptions to a file. 
 * This exports a file like deepseq.list for backing up the seqdata mysql database.
 * Since a couple of deepseq.list fields are not stored in the database, we unfortunately 
 * must provide the old deepseq.list file as input. 
 *  
 * @author mahony
 *
 * Usage: ExportSeqAlignmentList --list "deepseq.list" --out "outname"
 * 
 * The assumed input & output files are in the deepseq.list format, with the following fields:
 * 
 *0) ReadDBID
 *1) ExptType
 *2) Lab
 *3) ExptCondition
 *4) ExptTarget
 *5) CellLine
 *6) Replicate
 *7) Aligner
 *8) Genome
 *9) Permissions
 *10) PubSource
 *11) PublicDBID
 *12) CollabExptID
 *13) CollabAlignID
 *14) ReadType
 *15) AlignType
 *16) ReadLength
 *17) TotalReads
 *18) AlignedHits
 *19) UniquelyAlignedHits
 *20) DBLoadedHits
 *21) DBLoadedWeight
 *22) DBLoadedPairs
 *23) DBLoadedPairWeight
 *24) ReadsFile
 *25) AlignDir
 *26) AlignFile
 *27) IDXFile
 *28) AlignParamFile
 *29) ExptNote
 *30) LoadDate
 *31) ExptName
 * 
 */
public class ExportSeqAlignmentList {
	
	private String outFile, filename;
	
	public static void main(String args[]) {
		String outFile="out";
		String filename = Args.parseString(args,"list",null);
		ArgParser ap = new ArgParser(args);
		if(ap.hasKey("h") || filename==null){
			System.err.println("ExportSeqAlignmentList:\n" +
					"\t--list <deepseq.list format file>\n"+
					"\t--out <output filename>\n"
					);
		}else{
			if(ap.hasKey("out"))
				outFile = ap.getKeyValue("out");
			
			
			ExportSeqAlignmentList esal = new ExportSeqAlignmentList(filename, outFile);
			esal.execute();
		}
	}
		
	public ExportSeqAlignmentList(String infile, String outfile){
		this.outFile = outfile;
		this.filename = infile;
	}
	
	public void execute(){
	    try{
	    	SeqDataLoader loader = new SeqDataLoader();
	        MetadataLoader core = new MetadataLoader();
	        HashMap<Integer, DeepSeqEntry> oldList = new HashMap<Integer, DeepSeqEntry>();
	    
        	FileWriter fw = new FileWriter(outFile);
        	
        	//Load the old deepseq.list file
            BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
            String line = null;
    		while ((line = reader.readLine()) != null) {
    			String[] fields = line.split("\t");
    			if(!fields[0].equals("ReadDBID") && !fields[0].startsWith("#")){//skip the first line in the deepseq.list file
    				
    				DeepSeqEntry entry = new DeepSeqEntry();
    				//Variables
    				entry.dbid = new Integer(fields[0]);
    				entry.etypestring = fields[1];
    				entry.labstring = fields[2];
    				entry.conditionstring = fields[3];
    				entry.targetstring = fields[4];
    				entry.cellsstring = fields[5];
    				entry.replicate = fields[6];
    				entry.aligner = fields[7];
    				entry.genomestring = fields[8];
    				entry.permissions = fields[9];
    				entry.publicsource = fields[10];
    				entry.publicdbid = fields[11];
    				entry.collabexptid = fields[12];
    				entry.collabalignid = fields[13];
    				entry.rtypestring = fields[14];
    				entry.atypestring = fields[15];
    				entry.readlength = new Integer(fields[16].split(" ")[0]);
    				entry.numreadsStr = new String(fields[17]);
    				entry.alignedhitsStr = new String(fields[18].split(" ")[0]);
    				entry.uniqalignedhitsStr = new String(fields[19].split(" ")[0]);
    				entry.numhits = new Integer(fields[20].split(" ")[0]);
    				entry.totalweight = new Float(fields[21].split(" ")[0]);
    				entry.numpairs = new Integer(fields[22].split(" ")[0]);
    				entry.totalpairweight = new Float(fields[23].split(" ")[0]);
    				entry.fqfile = fields[24];
    				entry.aligndir = fields[25];
    				entry.alignfile = fields[26];
    				entry.idxfile = fields[27];
    				entry.paramsfname = fields[28];
    				entry.exptnote = fields[29];
    				entry.loaddate = fields[30];
    				entry.alignname= fields[31];
    				
    				oldList.put(entry.dbid, entry);
    			}
    		}
    		System.out.println(oldList.size()+" entries loaded from old list");
    		
	        
	        //Now load what is actually in the database, filling in missing entries with the old file
    		int count=0;
	        fw.write("#ReadDBID\tExptType\tLab\tExptCondition\tExptTarget\tCellLine\tReplicate\tAligner\tGenome\tPermissions\tPubSource\tPublicDBID\tCollabExptID\tCollabAlignID\tReadType\tAlignType\tReadLength\tTotalReads\tAlignedHits\tUniquelyAlignedHits\tDBLoadedHits\tDBLoadedWeight\tDBLoadedPairs\tDBLoadedPairWeight\tReadsFile\tAlignDir\tAlignFile\tIDXFile\tAlignParamFile\tExptNotes\tLoadDate\tExptName\n");
    		for(SeqExpt expt : loader.loadAllExperiments()){
	        	for(SeqAlignment align : loader.loadAllAlignments(expt)){
	        		String permissions ="";
	        		for(String p : align.getPermissions())
	        			permissions = permissions+p+";";
	        		if(!oldList.containsKey(align.getDBID()))
	        			System.out.println("NOTFOUND\t"+align.getDBID()+"\t"+expt.getName()+";"+expt.getReplicate()+";"+align.getName()+";"+align.getGenome().getVersion());
	        		else{
	        			fw.write(align.getDBID()+"\t"+
	        				expt.getExptType().getName()+"\t"+
	        				expt.getLab().getName()+"\t"+
	        				expt.getExptCondition().getName()+"\t"+
	        				expt.getExptTarget().getName()+"\t"+
	        				expt.getCellLine().getName()+"\t"+
	        				expt.getReplicate()+"\t"+
	        				align.getName()+"\t"+
	        				align.getGenome().getVersion()+"\t"+
	        				permissions+"\t"+
	        				expt.getPublicSource()+"\t"+
	        				expt.getPublicDBID()+"\t"+
	        				expt.getCollabID()+"\t"+
	        				align.getCollabAlignID()+"\t"+
	        				expt.getReadType().getName()+"\t"+
	        				align.getAlignType().getName()+"\t"+
	        				expt.getReadLength()+"\t"+
	        				expt.getNumRead()+"\t"+
	        				oldList.get(align.getDBID()).alignedhitsStr+"\t"+
	        				oldList.get(align.getDBID()).uniqalignedhitsStr+"\t"+
	        				align.getNumHits()+"\t"+
	        				align.getTotalWeight()+"\t"+
	        				align.getNumPairs()+"\t"+
	        				align.getTotalPairWeight()+"\t"+
	        				(expt.getFQFile().equals("") ? oldList.get(align.getDBID()).fqfile : expt.getFQFile())+"\t"+
	        				align.getAlignDir()+"\t"+
	        				align.getAlignFile()+"\t"+
	        				align.getIDXFile()+"\t"+
	        				oldList.get(align.getDBID()).paramsfname+"\t"+
	        				expt.getExptNote()+"\t"+
	        				oldList.get(align.getDBID()).loaddate+"\t"+
	        				expt.getName()+";"+expt.getReplicate()+";"+align.getName()+"\n"
	        				);
	        			count++;
	        		}
	        	}
	        }
    		System.out.println(count+" entries exported to new list");
	        fw.close();
	        core.close();
	        loader.close();
        } catch (IOException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}
	
	public class DeepSeqEntry{
		public Integer dbid;
		public String etypestring;
		public String labstring;
		public String conditionstring;
		public String targetstring;
		public String cellsstring;
		public String replicate;
		public String aligner;
		public String genomestring;
		public String permissions;
		public String publicsource;
		public String publicdbid;
		public String collabexptid;
		public String collabalignid;
		public String rtypestring;
		public String atypestring;
		public Integer readlength;
		public String numreadsStr;
		public String alignedhitsStr;
		public String uniqalignedhitsStr;
		public Integer numhits;
		public Float totalweight;
		public Integer numpairs;
		public Float totalpairweight;
		public String fqfile;
		public String aligndir;
		public String alignfile;
		public String idxfile;
		public String paramsfname;
		public String exptnote;
		public String loaddate;
		public String alignname;
	}
}
