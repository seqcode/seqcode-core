package edu.psu.compbio.seqcode.projects.akshay.chexmix.analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
			Map<BindingLocation, String> motif_orientation = new HashMap<BindingLocation, String>();
			while(currentline != null){
				String[] pieces = currentline.split("\t");
				int tempMidpoint = Integer.parseInt(pieces[3])+((Integer.parseInt(pieces[4])-Integer.parseInt(pieces[3]))/2);
				String tempChr = pieces[0];
				BindingLocation temploc = new BindingLocation(tempMidpoint, tempChr, driver.c);
				temploc.filltags(tagsloader);
				motif_orientation.put(temploc, pieces[6]);
				allbls.add(temploc);
				currentline = brpeaks.readLine();
			}
			brpeaks.close();
			System.currentTimeMillis();
			
			int i=1;
			List<BindingLocation> totalbls= allbls;
			
			System.out.println("Total no of Binding Locations in the input peak file are: "+allbls.size() );
			
			File file_entire_composite = new File(driver.c.getOutTagname()+"_entire_locations_composite.tab");
			if(!file_entire_composite.exists()){
				file_entire_composite.createNewFile();
			}
			FileWriter fw_entire_composite = new FileWriter(file_entire_composite.getAbsoluteFile());
			BufferedWriter br_entire_composite = new BufferedWriter(fw_entire_composite);
			int[] allblscomposite = ChexmixSandbox.getCompositeFromBlLisr(allbls,motif_orientation);
			for(int k=0; k<allblscomposite.length; k++){
                br_entire_composite.write(k+"\t"+allblscomposite[k]+"\n");
			}
			br_entire_composite.close();
			
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
				if(driver.c.getSchemename().equals("scheme3")){
					profile = seedbuilder.executeScheme3(driver.c);
				}
				System.currentTimeMillis();
				// debug methods start
				if(i==1){
					//for(BindingLocation bl: seedbuilder.seed){
					//	System.out.println(bl.getName());
					//}
					//int[] tempout = ChexmixSandbox.getCompositeFromBlLisr(seedbuilder.seed, motif_orientation);
					//for (int l=0; l<tempout.length; l++){
					//	System.out.println(Integer.toString(l+1)+"\t"+tempout[l]);
					//}
					List<BindingLocation> onlycenter = new ArrayList<BindingLocation>();
					onlycenter.add(seedbuilder.center);
					int[] tempout = ChexmixSandbox.getCompositeFromBlLisr(onlycenter, motif_orientation);
					for (int l=0; l<tempout.length; l++){
						System.out.println(Integer.toString(l+1)+"\t"+tempout[l]);
					}
					System.out.println(seedbuilder.center.getName());
				}
				//debug method ends
				
				System.out.println("No of Binding Locations selected to build seed "+i+" are: "+seedbuilder.getNoInSeed());
				System.out.println("Composite of seed "+i+":");
				System.out.println(Arrays.toString(profile));
				File file = new File(driver.c.getOutTagname()+"_seed_profile_composite_"+i+".tab");
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
				System.out.println("No of locations that match the seed profile "+i+" are:"+scanner.getListOfBlsThatPassCuttoff().size());
				File file_pcc =  new File(driver.c.getOutTagname()+"_complete_list_pcc_"+i+".tab");
				if(!file_pcc.exists()){
					file_pcc.createNewFile();
				}
				FileWriter fw_pcc = new FileWriter(file_pcc.getAbsoluteFile());
				BufferedWriter br_pcc = new BufferedWriter(fw_pcc);
				for(int j=0; j<scanner.getEntireListOfPccValues().length; j++){
					br_pcc.write(Double.toString(scanner.getEntireListOfPccValues()[j])+"\n");
				}
				br_pcc.close();
				File file_tagsPass = new File(driver.c.getOutTagname()+"_scan_pass_profile_composite_"+i+".tab");
				if(!file_tagsPass.exists()){
					file_tagsPass.createNewFile();
				}
				FileWriter fw_tagsPass = new FileWriter(file_tagsPass.getAbsoluteFile());
				BufferedWriter br_tagsPass = new BufferedWriter(fw_tagsPass);
				int[][] temp = scanner.getTagsThatPassCuttoff(c);
				int[] composite_passed = new int[temp[0].length];
				for(int j=0; j<temp.length; j++){
					for(int k=0; k<temp[j].length; k++){
						if(j==0){
							composite_passed[k] = temp[j][k];
						}
						else{
							composite_passed[k] = composite_passed[k] + temp[j][k];
						}
						
					}
				}
				for(int j=0; j<composite_passed.length; j++){
					br_tagsPass.write(j+"\t"+composite_passed[j]+"\n");
				}
				br_tagsPass.close();
				System.currentTimeMillis();
				
				totalbls = scanner.getListOfBlsThatDoNotPassCuttOff();
				
				int[] tempcomposite = ChexmixSandbox.getCompositeFromBlLisr(totalbls, motif_orientation);
				File file_remaining_all_composite = new File(driver.c.getOutTagname()+"_remaining_all_composite_"+i+".tab");
				if(!file_remaining_all_composite.exists()){
					file_remaining_all_composite.createNewFile();
				}
				FileWriter fw_remaining_all_composite = new FileWriter(file_remaining_all_composite.getAbsoluteFile());
				BufferedWriter br_remaining_all_composite = new BufferedWriter(fw_remaining_all_composite);
				for(int j=0; j< tempcomposite.length; j++){
					br_remaining_all_composite.write(j+"\t"+tempcomposite[j]+"\n");
				}
				br_remaining_all_composite.close();
				i++;
				
			}
		}
	}
}
