package org.seqcode.projects.akshay.multigpstesting;

import java.sql.*;
import java.util.*;

import org.seqcode.gse.utils.database.DatabaseConnectionManager;
import org.seqcode.gse.utils.database.DatabaseException;




public class Pullencodetier1 {
	private static final String seqdata = "seqdata";
	private static final String core = "core";
	public static final String celllineIds = "^36$|^401$|^34$|^135$|^77$";
	private Map<String,String> FacNames = new HashMap<String,String>();
	public Map<String,List<String>> Exptrecords = new HashMap<String,List<String>>();
	
	/*
	 * the constructor class established connections to the core and seqdata databases
	 */
	
	public Pullencodetier1() throws SQLException {
		
	}
	
	/*
	 * looks at the entire expttarget SQL table for factor aliases for an input list of factors
	 * return a hash map with keys as readdb id (expttarget) in string format and value as the fator name 
	 */
	
	public void getFactorAliases(Map<String,String[]> input) throws SQLException{
		Connection cxnCore = null;
		try{
			cxnCore = DatabaseConnectionManager.getConnection(core);
			for (String facname: input.keySet()){
				String querytarget = "select * from expttarget where name regexp "+'"'+"^"+facname+'"';
				PreparedStatement pstarget = cxnCore.prepareStatement(querytarget);
				ResultSet rstarget = pstarget.executeQuery();
				while(rstarget.next()){
					FacNames.put(rstarget.getString("id"),facname);
				}
			}
		} finally {
			if(cxnCore!=null) try {cxnCore.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+core, ex); }
		}
	}
	
	
	
	public void getExperiments() throws SQLException{
		Connection cxnSeqdata = null;
		try{
			cxnSeqdata = DatabaseConnectionManager.getConnection(seqdata);
			String queryexpt = "select * from seqexpt where cellline regexp "+'"'+celllineIds+'"';
			PreparedStatement psexpt = cxnSeqdata.prepareStatement(queryexpt);
			ResultSet rsexpt = psexpt.executeQuery();
			while(rsexpt.next()){
				if(FacNames.containsKey(rsexpt.getString("expttarget"))){
					String design = rsexpt.getString("name")+";"+rsexpt.getString("replicate")+";bowtie_unique";
					if(Exptrecords.containsKey(FacNames.get(rsexpt.getString("expttarget")))){
						Exptrecords.get(FacNames.get(rsexpt.getString("expttarget"))).add(design);
					}
					else{
						List<String> temp = new LinkedList<String>(Arrays.asList(design));
						Exptrecords.put(FacNames.get(rsexpt.getString("expttarget")),temp);
					}
				}
			}
		} finally {
			if(cxnSeqdata!=null) try {cxnSeqdata.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+seqdata, ex); }
		}
		
	}
	
	public void printExptMap(){
		for(String fac:Exptrecords.keySet()){
			System.out.println(fac);
			for(String record: Exptrecords.get(fac)){
				String[] superpiece = record.split(";");
				String[] pieces = superpiece[0].split(" ");
				
				System.out.println("Signal"+"\t"+record+"\t"+"READDB"+"\t"+fac+"_"+pieces[pieces.length-1]+"_"+superpiece[1]
						+"\t"+superpiece[1]);
			}
			System.out.println("Control"+"\t"+"ENCh-Broad-Bernstein GM12878 Input-std-ChipSeq GM12878;1;bowtie_unique"+"\t"
			+"READDB"+"\t"+"Cntrl_GM12878"+"\t"+"1");
			System.out.println("Control"+"\t"+"ENCh-Broad-Bernstein H1-hESC Input-std-ChipSeq H1-hESC;1;bowtie_unique"+"\t"
					+"READDB"+"\t"+"Cntrl_H1"+"\t"+"1");
			System.out.println("Control"+"\t"+"ENCh-Broad-Bernstein K562 Input-std-ChipSeq K562;1;bowtie_unique"+"\t"
					+"READDB"+"\t"+"Cntrl_k562"+"\t"+"1");
			
		}
	}


}
