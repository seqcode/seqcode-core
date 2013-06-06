package edu.psu.compbio.seqcode.projects.akshay.Querying;

import java.io.*;
import java.util.*;
import java.sql.*;
import edu.psu.compbio.seqcode.gse.utils.*;

public class printexpttable {
	public static void main(String[] args) throws IOException, SQLException, NotFoundException{
		Scanner sc = new Scanner("/Users/akshaykakumanu/Desktop/temp");
		Map<String, String[]> map = sc.getMap();
		for(String[] command : map.values()){
			getAlign galign = new getAlign(command);
			galign.executeSQLExptCommand();
			galign.executeSQLAlignCommand();
			galign.printTable();
			
		}
		
				
	}

}
