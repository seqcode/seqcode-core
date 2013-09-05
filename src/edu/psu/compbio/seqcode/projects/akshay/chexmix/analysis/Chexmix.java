package edu.psu.compbio.seqcode.projects.akshay.chexmix.analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import cern.colt.Arrays;

import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.BindingLocation;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.Config;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.utils.ChexmixSandbox;

public class Chexmix {
	public Config c;
	public Chexmix(Config conf) {
		this.c = conf;
	}

	public static void main(String[] args) throws IOException{
		Config c = new Config(args);
		ChexmixSandbox sandbox = new ChexmixSandbox();
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
			List<BindingLocation> allbls = new ArrayList<BindingLocation>();
			String currentline = brpeaks.readLine();
			while(currentline != null){
				String[] pieces = currentline.split("\t");
				int tempMidpoint = Integer.parseInt(pieces[3])+((Integer.parseInt(pieces[4])-Integer.parseInt(pieces[3]))/2);
				String tempChr = pieces[0];
				BindingLocation temploc = new BindingLocation(tempMidpoint, tempChr, driver.c);
				temploc.filltags(tagsloader);
				allbls.add(temploc);
				currentline = brpeaks.readLine();
			}
			brpeaks.close();
			System.currentTimeMillis();
			
			int i=1;
			List<BindingLocation> totalbls = allbls;
			
			//debug lines
			System.out.println("Printing the composite of the entire list of bls");
			List<BindingLocation> temp = new ArrayList<BindingLocation>();
			for(int k=0; k<5000; k++){
				temp.add(allbls.get(k));
			}
			int[] allblscomposite = sandbox.getCompositeFromBlLisr(temp);
			for(int k=0; k<allblscomposite.length; k++){
				System.out.println(k+"\t"+allblscomposite[k]);
			}
			
			
			while(i<=driver.c.getNoOfCycles() && totalbls.size()>driver.c.getNoOfCycles()){
				List<BindingLocation> selectedbls = new ArrayList<BindingLocation>();
				for(int j=0; j<driver.c.getNoTopBls(); j++){
					selectedbls.add(totalbls.get(j));
				}
				System.out.println("\n============================ Building Seed Profile - "+i+" ============================");
				BuildSeed seedbuilder = new BuildSeed(selectedbls, driver.c);
				int[] profile=null ;
				if(driver.c.getSchemename().equals("scheme1")){
					profile= seedbuilder.executeScheme1(driver.c);
				}

				if(driver.c.getSchemename().equals("scheme2")){
					profile=seedbuilder.executeScheme2(driver.c);
				}
				System.currentTimeMillis();
				System.out.println("No of Binding Locations selected to build seed "+i+" are: "+seedbuilder.getNoInSeed());
				System.out.println("Composite of seed "+i+":");
				System.out.println(Arrays.toString(profile));
				File file = new File(driver.c.getOutTagname()+"_seed_profile_"+i+".tab");
				if(!file.exists()){
					file.createNewFile();
				}
				FileWriter fw = new FileWriter(file.getAbsoluteFile());
				BufferedWriter br_profile = new BufferedWriter(fw);
				for(int j=0; j<profile.length; j++){
					br_profile.write(Integer.toString(j+1)+"\t"+Integer.toString(profile[j])+"\n");
				}
				br_profile.close();
				
				System.out.println("\n============================ Scanning the entire list of binding locations - "+i+" ============================");
				LocationsScanner scanner = new LocationsScanner(totalbls, driver.c, profile);
				File file_pcc =  new File(driver.c.getOutTagname()+"_list_pcc_"+i+".tab");
				if(!file_pcc.exists()){
					file_pcc.createNewFile();
				}
				FileWriter fw_pcc = new FileWriter(file_pcc.getAbsoluteFile());
				BufferedWriter br_pcc = new BufferedWriter(fw_pcc);
				for(int j=0; j<scanner.getEntireListOfPccValues().length; j++){
					br_pcc.write(Double.toString(scanner.getEntireListOfPccValues()[j])+"\n");
				}
				br_pcc.close();
				System.currentTimeMillis();
				
				totalbls = scanner.getListOfBlsThatDoNotPassCuttOff();
				i++;
				
			}
		}
	}
}
