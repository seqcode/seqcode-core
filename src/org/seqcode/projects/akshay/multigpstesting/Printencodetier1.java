package org.seqcode.projects.akshay.multigpstesting;

import java.io.*;
import java.sql.*;

import org.seqcode.projects.akshay.query.*;

public class Printencodetier1 {

	public static void main(String[] args) throws SQLException, IOException{
		Pullencodetier1 getlist = new Pullencodetier1();
		Scanner sc = new Scanner("/Users/akshaykakumanu/PSU/multigpstesting/encode_tf_list");
		getlist.getFactorAliases(sc.getMap());
		getlist.getExperiments();
		getlist.printExptMap();
	}
}
