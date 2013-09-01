package edu.psu.compbio.seqcode.projects.akshay.chexmix.analysis;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.BindingLocation;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.Config;

public class Chexmix {
	public Config c;
	public Chexmix(Config conf) {
		this.c = conf;
	}

	public static void main(String[] args) throws IOException{
		Config c = new Config(args);
		Chexmix driver = new Chexmix(c);
		if(driver.c.helpWanted()){
			System.err.println("Chexmix:");
			System.err.println(driver.c.getArgsList());
		}
		else{
			System.currentTimeMillis();
			System.out.println("\n============================ Loading Tags/Reads ============================");
			LoadTags tagsloader = new LoadTags();
			tagsloader.loadHits(driver.c, driver.c.useNonUnique());
			System.currentTimeMillis();
			System.out.println("\n============================ Filling Binding Locations ============================");
			BufferedReader brpeaks = new BufferedReader(new FileReader(driver.c.getPeaksFilePath()));
			ArrayList<BindingLocation> allbls = new ArrayList<BindingLocation>();
			String currentline = brpeaks.readLine();
			while(currentline != null){
				String[] pieces = currentline.split("\t");
				int tempMidpoint = Integer.parseInt(pieces[3])+((Integer.parseInt(pieces[4])-Integer.parseInt(pieces[3]))/2);
				String tempChr = pieces[0];
				BindingLocation temploc = new BindingLocation(tempMidpoint, tempChr, driver.c);
				//debug line
				System.out.println(temploc.getName());
				temploc.filltags(tagsloader);
				allbls.add(temploc);
				currentline = brpeaks.readLine();
			}
			brpeaks.close();
			System.currentTimeMillis();
			System.out.println("\n============================ Building Seed Profile ============================");
			ArrayList<BindingLocation> topbls = new ArrayList<BindingLocation>();
			for(int i=0; i<driver.c.getNoTopBls(); i++){
				topbls.add(allbls.get(i));
			}
			BuildSeed seedbuilder = new BuildSeed(topbls, driver.c);
			int[] profile=null ;
			if(driver.c.getSchemename().equals("scheme1")){
				profile= seedbuilder.executeScheme1(driver.c);
			}
			if(driver.c.getSchemename().equals("scheme2")){
				profile=seedbuilder.executeScheme2(driver.c);
			}
			System.currentTimeMillis();
			System.out.println("No of Binding Locations selected to build seed are: "+seedbuilder.getNoInSeed());
			System.out.println("Composite of seed:");
			System.out.println(profile);
			System.out.println("\n============================ Scanning the entire list of binding locations ============================");
			LocationsScanner scanner = new LocationsScanner(allbls, driver.c, profile);
			System.out.println(scanner.getListofPCCvalues());
			System.currentTimeMillis();
		}
		
		
	}

}
