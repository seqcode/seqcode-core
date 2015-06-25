package edu.psu.compbio.seqcode.projects.akshay.query;

import java.io.*;
import java.util.*;
import java.sql.*;
import edu.psu.compbio.seqcode.gse.utils.*;

public class Printexpttable {
	
	public static void main(String[] args) throws IOException, SQLException, NotFoundException{
		ArgParser ap = new ArgParser(args);
		String com_file = ap.getKeyValue("commands");
		Scanner sc = new Scanner(com_file);
		Map<String, String[]> map = sc.getMap();
		for(String[] command : map.values()){
			Pullexpttable gexp = new Pullexpttable(command);
			gexp.executeSQLExptCommand();
			gexp.printTable();
		}
	}
}
