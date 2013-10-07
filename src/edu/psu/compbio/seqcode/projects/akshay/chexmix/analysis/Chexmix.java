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
import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.CustomReturn;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.datasets.Membership;
import edu.psu.compbio.seqcode.projects.akshay.chexmix.utils.ChexmixSandbox;

public class Chexmix {
	public Config c;
	public List<int[]> profiles_in_the_dataset = new ArrayList<int[]>();
	public List<Membership> cluster_assignment = new ArrayList<Membership>();
	public double[][] cluster_assessment;
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
//========================================================================================LOADING TAGS================================================================================================
			System.out.println("\n============================ Loading Tags/Reads ============================");
			LoadTags tagsloader = new LoadTags();
			tagsloader.loadHits(driver.c, driver.c.useNonUnique());
//========================================================================================FILLING LOCATIONS=================================================================================================
			System.out.println("\n============================ Filling Binding Locations ============================");
			BufferedReader brpeaks = new BufferedReader(new FileReader(driver.c.getPeaksFilePath()));
			List<BindingLocation> allbls = new ArrayList<BindingLocation>();
			String currentline = brpeaks.readLine();
			Map<BindingLocation, String> motif_orientation = new HashMap<BindingLocation, String>();
			Map<BindingLocation, Double> location_coverage = new HashMap<BindingLocation, Double>();
			while(currentline != null){
				String[] pieces = currentline.split("\t");
				int tempMidpoint = Integer.parseInt(pieces[3])+((Integer.parseInt(pieces[4])-Integer.parseInt(pieces[3]))/2);
				String tempChr = pieces[0];
				BindingLocation temploc = new BindingLocation(tempMidpoint, tempChr, driver.c);
				temploc.filltags(tagsloader);
				motif_orientation.put(temploc, pieces[6]);
				location_coverage.put(temploc, Double.parseDouble(pieces[5]));
				allbls.add(temploc);
				currentline = brpeaks.readLine();
			}
			brpeaks.close();
			System.currentTimeMillis();
			
			int size_to_consider =  (int) (driver.c.getListPercentageToCosider()*0.01*allbls.size());
			List<BindingLocation> totalbls= new ArrayList<BindingLocation>();
			for(int k=0; k< size_to_consider; k++){
				totalbls.add(allbls.get(k));
			}
//===========================================================================================BUILDING SEEDS===============================================================================================
			int i=1;
			
			
			System.out.println("Total no of Binding Locations in the input peak file are: "+allbls.size() );
			System.out.println("Total no of Binding Locations that are being considered: "+totalbls.size() );
			
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
			
			while(i<=driver.c.getNoOfCycles() && totalbls.size()>driver.c.getHowDeepToSearch()*driver.c.getNoTopBls()){
				int counter = 0;
				System.out.println("\n============================ Building Seed Profile - "+i+" ============================");
				List<BindingLocation> selectedbls = new ArrayList<BindingLocation>();
				for(counter=0; counter<driver.c.getNoTopBls(); counter++){
					selectedbls.add(totalbls.get(counter));
				}
				
				BuildSeed seedbuilder = new BuildSeed(selectedbls, driver.c);
				int[] profile=null;
				int Final_number_in_profile=0;
				CustomReturn profile_cr = null;
				
				
				if(driver.c.getSchemename().equals("scheme4")){
					profile_cr = seedbuilder.executeScheme4(driver.c);
				}
				
				
				while(counter < driver.c.getHowDeepToSearch()*driver.c.getNoTopBls() && counter < totalbls.size()){
					//debug line
					System.out.println(counter);
					System.out.println(profile_cr.no_in_seed);
					//end
					if(profile_cr.no_in_seed >= driver.c.getNoTopBls()/4){
						profile = profile_cr.profile;
						Final_number_in_profile = profile_cr.no_in_seed;
						break;
					}
					else if(profile_cr.no_in_seed < driver.c.getNoTopBls()/4 && profile_cr.no_in_seed > driver.c.getNoTopBls()/10){
						profile = profile_cr.profile;
						int no_bl_already_in= profile_cr.no_in_seed;
						while(no_bl_already_in <= driver.c.getNoTopBls()/4 && counter < driver.c.getHowDeepToSearch()*driver.c.getNoTopBls() ){
							BindingLocation to_scan = totalbls.get(counter);
							CustomReturn temp_cr = to_scan.scanConcVecWithBl(profile, driver.c.getIntSize());
							//debug line
							System.out.println(temp_cr.pcc);
							//end
							if(temp_cr.pcc > driver.c.getSeedCutoff()){
								List<Integer> addtoprofile = to_scan.getConcatenatedTags(temp_cr.maxvec.midpoint, temp_cr.maxvec.range, temp_cr.maxvec.orientation);
								for(int k=0; k< addtoprofile.size(); k++){
									profile[k] = profile[k]+addtoprofile.get(k);
								}
								no_bl_already_in++;
							}
							counter++;
						}
						Final_number_in_profile = no_bl_already_in;
						//debug
						System.out.println(Final_number_in_profile);
						//end 
						break;
					}
					else if(profile_cr.no_in_seed <= driver.c.getNoTopBls()/10){
						int temp_count=0;
						selectedbls=new ArrayList<BindingLocation>();
						while(temp_count<driver.c.getNoTopBls()){
							selectedbls.add(totalbls.get(counter));
							counter++;
							temp_count++;
						}
						seedbuilder = new BuildSeed(selectedbls, driver.c);
						if(driver.c.getSchemename().equals("scheme4")){
							profile_cr = seedbuilder.executeScheme4(driver.c);
						}
					
					}
				}
				if(Final_number_in_profile < driver.c.getNoTopBls()/driver.c.getFactorToAddIteratively()){
					int no_of_porfiles_in_dataset = i-1;
					System.out.println("There are "+no_of_porfiles_in_dataset+" profiles in this dataset");
					break;
				}
//=============================================================================== FINISHED BUILDING SEEDS ========================================================================			
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
				
				driver.profiles_in_the_dataset.add(profile);
//==================================================================================== SCANNING THE SEED PROFILE THROUGH THE ENTIRE DATASET ==============================================				
				System.out.println("\n============================ Scanning the entire list of binding locations - "+i+" ============================");
				LocationsScanner scanner = new LocationsScanner(totalbls, driver.c, profile);
				System.out.println("No of locations that match the seed profile "+i+" are:"+scanner.getListOfBlsThatPassCuttoff().size());
				
				List<Membership> temp_add_to_cluster_assignment = scanner.getMembershipsForThoseThatPassCuttoff(i, driver.c);
				for(Membership temp :  temp_add_to_cluster_assignment){
					driver.cluster_assignment.add(temp);
				}
				
				File file_pcc =  new File(driver.c.getOutTagname()+"_complete_list_pcc_"+i+".tab");
				if(!file_pcc.exists()){
					file_pcc.createNewFile();
				}
				FileWriter fw_pcc = new FileWriter(file_pcc.getAbsoluteFile());
				BufferedWriter br_pcc = new BufferedWriter(fw_pcc);
				Map<BindingLocation, Double> tempmap = scanner.getmapofAllblscan();
				for(BindingLocation tempbl :tempmap.keySet()){
					br_pcc.write(location_coverage.get(tempbl)+"\t"+tempmap.get(tempbl)+"\n");
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
				
			} //end of while
			
			// extending the seed profiles
			System.out.println("Expanding the seed profiles for each cluster group");
			for(int l=0; l< driver.profiles_in_the_dataset.size(); l++){
				int add_till = 20;
				int k=0;
				while(add_till < driver.c.getFactorToRefineSeedProfiles()*driver.c.getNoTopBls()){
					if(driver.cluster_assignment.get(k).membership == l+1){
						List<Integer> add_to_profile = driver.cluster_assignment.get(k).bl.getConcatenatedTags(driver.cluster_assignment.get(k).cr.maxvec.midpoint,
								driver.cluster_assignment.get(k).cr.maxvec.range, driver.cluster_assignment.get(k).cr.maxvec.orientation);
						for(int m =0; m< add_to_profile.size(); m++){
							driver.profiles_in_the_dataset.get(l)[m] = driver.profiles_in_the_dataset.get(l)[m] + add_to_profile.get(m);
						}
						add_till++;
						//debug line 
						System.out.println(add_till);
					}
					k++;
				}
			}
			
			// re assigning cluster membership
			
			System.out.println("Reassigning cluster membership");
			List<Membership> temp = ChexmixSandbox.perfromReAssignmentOfMembership(driver.cluster_assignment, driver.profiles_in_the_dataset);
			driver.cluster_assignment = temp;
			temp = null;
			
			// do cluster assessment
			
			System.out.println("Cluster assessment");
			driver.cluster_assessment = ChexmixSandbox.getClusterAssesmsent(driver.cluster_assignment, driver.profiles_in_the_dataset);
			for(int m = 0; m< driver.cluster_assessment.length; m++){
				String out ="";
				for(int n=0; n< driver.cluster_assessment[m].length; n++){
					out = out + "\t" + Double.toString(driver.cluster_assessment[m][n]);
				}
				System.out.println(out);
			}
			
			
		} // end of else
	} // end of main
} // end of class
