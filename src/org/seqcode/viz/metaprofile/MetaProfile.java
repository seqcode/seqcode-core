package org.seqcode.viz.metaprofile;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * MetaProfile: a collection of profiles.
 * @author: tdanford
 * Date: Aug 12, 2008
 */
public class MetaProfile implements Profile, ProfileListener{

	protected String name; 
	protected BinningParameters params;
	protected double[] values;
	protected Double normalization;
	protected Vector<Profile> profiles;
	protected double max, min;
	protected boolean stranded=false;
	private LinkedList<ProfileListener> listeners;

		

	public MetaProfile(String n, BinningParameters bps) { 
		name = n;
		params = bps;
		values = new double[params.getNumBins()];
		max = min = 0.0;
		profiles = new Vector<Profile>();
		normalization = null;
		listeners = new LinkedList<ProfileListener>();
	}
	public void saveToFile(String fileName){
		if(profiles.size()>0){
			try {
				FileWriter fout = new FileWriter(fileName);
				int start = (-1*(params.getWindowSize()/2))+params.getBinSize()/2;
				int step = params.getWindowSize()/params.getNumBins();
				
				fout.write(name+"\n");
				int k= start;
				for(int i=0; i<values.length; i++){
					fout.write(k+"\t"+values[i]+"\n");
					k+=step;
				}			
				fout.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}else{
			System.err.println("Empty MetaProfile: nothing to write to file");
		}
	}
	public void saveProfilesToFile(String fileName){
		if(profiles.size()>0){
			try {
				FileWriter fout = new FileWriter(fileName);
				for(Profile p : profiles) { 
					for(int i = 0; i < values.length; i++) {
						fout.write(String.format("%.2f", p.value(i)));
						if(i<values.length-1)
							fout.write("\t");
					}
					fout.write("\n");
				}			
				fout.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}else{
			System.err.println("Empty MetaProfile: nothing to write to file");
		}
	}
	
	public synchronized void normalize() { 
		if(profiles.size() > 0) { 
			normalize((double)profiles.size());
		}
	}
	public synchronized void normalize(double dn) {
		max = min = 0.0;
		if(dn <= 0.0) { 
			throw new IllegalArgumentException(String.format("Can't normalize with factor %f", dn));
		}
		double norm = normalization == null ? 1.0 / dn : normalization / dn;
		for(int i = 0; i < values.length; i++) { 
			values[i] *= norm;
			max = Math.max(max, values[i]);
			min = Math.min(min, values[i]);
		}normalization = dn;
	}
	
	protected void recalculate() { 
		max = min = 0.0;
		for(int i = 0; i < values.length; i++) { 
			values[i] = 0.0;
		}
		
		for(Profile p : profiles) { 
			for(int i = 0; i < values.length; i++) { 
				values[i] += p.value(i);
			}
		}
		
		if(isNormalized()) { 
			for(int i = 0; i < values.length; i++) { 
				values[i] /= normalization;
				max = Math.max(max, values[i]);
				min = Math.min(min, values[i]);
			}
		}
	}
	
	public synchronized void clear() {
		normalization = null;
		profiles.clear();
		min = max = 0.0;
		for(int i = 0; i < values.length; i++ ){ 
			values[i] = 0.0;
		}
	}

	public int size() { return profiles.size(); }
	public int length() { return values.length; }
	public double value(int i) { return values[i]; }
	public double max() { return max; }
	public double min() { return min; }
	public void setStranded(boolean s){stranded = s;}
	public boolean isStranded(){return stranded;}

	public String getName() { return name; }
	public BinningParameters getBinningParameters() { return params; }
	public Profile profile(int i) { return profiles.get(i); }
	public boolean isNormalized() { return normalization != null; }
	
	public synchronized void addProfile(Profile p) {
		if(p.isStranded()){
			stranded=true;
		}
		if(p.length() != params.getNumBins()) { 
			throw new IllegalArgumentException(String.format("Profile length %d doesn't" +
					" match bin-length %d", p.length(), params.getNumBins()));
		}
		
		if(isNormalized()) { 
			throw new IllegalArgumentException("Can't add profile to a normalized MetaProfile");
		}
		
		if(profiles.contains(p)) { 
			/*throw new IllegalArgumentException(String.format(
					"Can't add same profile %s to MetaProfile", 
					p.getName()));*/
		}else{
		
			profiles.add(p);
			for(int i = 0; i< values.length ;i++) { 
				values[i] += p.value(i);
				max = Math.max(max, values[i]);
				min = Math.min(min, values[i]);
			}			
			dispatchChange(new ProfileEvent(this, p));
		}
	}
	
	public String toString() { return name; }
	
	public int hashCode() { return name.hashCode(); }
	
	public boolean equals(Object o) { 
		if(!(o instanceof MetaProfile)) { return false; }
		MetaProfile mp = (MetaProfile)o;
		if(!mp.name.equals(name)) { return false; }
		return true;
	}

	public int getNumProfiles() {
		return profiles.size();
	}
	
	public void addProfileListener(ProfileListener pl) { 
		listeners.add(pl);
	}
	
	public void removeProfileListener(ProfileListener pl) {  
		listeners.remove(pl);
	}
	
	protected void dispatchChange(ProfileEvent e) { 
		for(ProfileListener pl : listeners) { 
			pl.profileChanged(e);
		}
	}
	
	public synchronized void profileChanged(ProfileEvent p) {
		recalculate();
	}
}
