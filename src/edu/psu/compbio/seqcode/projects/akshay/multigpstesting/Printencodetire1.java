package edu.psu.compbio.seqcode.projects.akshay.multigpstesting;

import java.io.*;
import java.sql.*;
import edu.psu.compbio.seqcode.projects.akshay.query.*;

public class Printencodetire1 {

	public static void main(String[] args) throws SQLException, IOException{
		Pullencodetire1 getlist = new Pullencodetire1();
		Scanner sc = new Scanner("/Users/akshaykakumanu/Desktop/tt");
		getlist.getFactorAliases(sc.getMap());
		getlist.getExperiments();
		getlist.printExptMap();
	}
}
