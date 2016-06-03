package org.seqcode.projects.shaun.teqseq.core;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.datasets.seqdata.SeqAlignment;
import org.seqcode.gse.datasets.seqdata.SeqDataLoader;
import org.seqcode.gse.datasets.seqdata.SeqLocator;
import org.seqcode.gse.projects.readdb.*;
import org.seqcode.gse.utils.NotFoundException;


/**
 * ReadDBHitLoader: load read alignment five primes from a collection of ReadDB alignments (ChipSeqLocators)
 * @author mahony
 *
 */
public class ReadDBSource{
	
	private Genome gen=null;
	private Client client=null; //ReadDB client
	private String name;
	private String condition="";
	//Single-end
	private List<String> singleExptNames =new ArrayList<String>();
	private List<SeqAlignment> singleAligns = new ArrayList<SeqAlignment>();
	private Collection<String> singleAlignIDs= new ArrayList<String>();
	//Junctions
	private List<String> junctExptNames =new ArrayList<String>();
	private List<SeqAlignment> junctAligns = new ArrayList<SeqAlignment>();
	private Collection<String> junctAlignIDs= new ArrayList<String>();
	//Paired-end
	private List<String> pairedExptNames =new ArrayList<String>();
	private List<SeqAlignment> pairedAligns = new ArrayList<SeqAlignment>();
	private Collection<String> pairedAlignIDs= new ArrayList<String>();
	
	/**
	 * Constructor: initialize the experiments
	 * @param g Genome
	 * @param locs List<ChipSeqLocator>
	 */
	public ReadDBSource(String name, Genome g, List<SeqLocator> singleLocs, List<SeqLocator> junctLocs, List<SeqLocator> pairedLocs){
		gen=g;
		this.name = name;
		try {
			client = new Client();
			
			SeqDataLoader loader = new SeqDataLoader(false, true);
			
			//Initialize Single-end ChipSeqLoaders
            for(SeqLocator locator : singleLocs){
				String exptName = locator.getExptName();
				singleExptNames.add(exptName);
				for (String replicate : locator.getReplicates()) {//Given alignment name, given replicate names
        			singleAligns.add(loader.loadAlignment(loader.loadExperiment(locator.getExptName(),
        					replicate), 
        					locator.getAlignName(),
        					g));
        		}
			}
	        for(SeqAlignment alignment : singleAligns) {
	            singleAlignIDs.add(Integer.toString(alignment.getDBID()));
	        }
	        
	        //Initialize Junction ChipSeqLoaders
            for(SeqLocator locator : junctLocs){
				String exptName = locator.getExptName();
				junctExptNames.add(exptName);
				for (String replicate : locator.getReplicates()) {//Given alignment name, given replicate names
        			junctAligns.add(loader.loadAlignment(loader.loadExperiment(locator.getExptName(),
        					replicate), 
        					locator.getAlignName(),
        					g));
        		}
			}
	        for(SeqAlignment alignment : junctAligns) {
	            junctAlignIDs.add(Integer.toString(alignment.getDBID()));
	        }
	        
	        //Initialize Paired-end ChipSeqLoaders
            for(SeqLocator locator : pairedLocs){
				String exptName = locator.getExptName();
				pairedExptNames.add(exptName);
				for (String replicate : locator.getReplicates()) {//Given alignment name, given replicate names
        			pairedAligns.add(loader.loadAlignment(loader.loadExperiment(locator.getExptName(),
        					replicate), 
        					locator.getAlignName(),
        					g));
        		}
			}
	        for(SeqAlignment alignment : pairedAligns) {
	            pairedAlignIDs.add(Integer.toString(alignment.getDBID()));
	        }
		} catch (IOException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (NotFoundException e) {
			e.printStackTrace();
		} catch (ClientException e) {
			e.printStackTrace();
		}
	}
	
	public String getName(){return name;}
	public String getConditionName(){return condition;}
	
	public void setConditionName(String c){condition=c;}
	/**
	 * Single hits overlapping a region
	 * @param r
	 * @return
	 */
	public List<SingleHit> getSingleHits(Region r){
		List<SingleHit> hits = new ArrayList<SingleHit>();
		for(SeqAlignment a : singleAligns){
				try {
					hits.addAll(client.getSingleHits(Integer.toString(a.getDBID()),
					       r.getGenome().getChromID(r.getChrom()),false, r.getStart(), r.getEnd(), null, null));
				} catch (IOException e) {
					e.printStackTrace();
				} catch (ClientException e) {
					
				}
		}
		Collections.sort(hits);
		return hits;
	}
	
	/**
	 * Junction hits overlapping a region
	 * @param r
	 * @return
	 */
	public List<PairedHit> getJunctionHits(Region r){
		List<PairedHit> hits = new ArrayList<PairedHit>();
		for(SeqAlignment a : junctAligns){
				try {
					hits.addAll(client.getPairedHits(Integer.toString(a.getDBID()),
					       r.getGenome().getChromID(r.getChrom()), true, r.getStart(), r.getEnd(), null, null));
				} catch (IOException e) {
					e.printStackTrace();
				} catch (ClientException e) {
				
				}
		}
		Collections.sort(hits);
		return hits;
	}
	
	/**
	 * Paired hits overlapping a region
	 * @param r
	 * @return
	 */
	public List<PairedHit> getPairedHits(Region r){
		List<PairedHit> hits = new ArrayList<PairedHit>();
		for(SeqAlignment a : pairedAligns){
				try {
					hits.addAll(client.getPairedHits(Integer.toString(a.getDBID()),
					       r.getGenome().getChromID(r.getChrom()), true, r.getStart(), r.getEnd(), null, null));
				} catch (IOException e) {
					e.printStackTrace();
				} catch (ClientException e) {
				
				}
		}
		Collections.sort(hits);
		return hits;
	}
	
	/**
	 * Return the total weight for the experiment, defined by the single reads
	 * @return
	 */
	public double getTotalSingleWeight(){
		double w =0;
		try {
			for(SeqAlignment alignment : singleAligns) { 
				double currWeight =(double)client.getWeight(Integer.toString(alignment.getDBID()),false, false, null, null);
				w +=currWeight;
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClientException e) {
			//Do nothing here; ClientException could be thrown because a chromosome doesn't contain any hits
			return(w);
		}
		return w;
	}
	
	/**
	 * Close the client
	 */
	public void close(){
		if(client==null)
			client.close();
	}
}
