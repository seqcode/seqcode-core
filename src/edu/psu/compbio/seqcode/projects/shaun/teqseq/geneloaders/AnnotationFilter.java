package edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.List;

import edu.psu.compbio.seqcode.genome.location.Region;

/**
 * AnnotationFilter: a set of tools to filter specific subsets from a collection of known genes
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class AnnotationFilter {

	/**
	 * Find genes overlapping a region
	 * @param currReg
	 * @param genes
	 * @return
	 */
	public static List<AGene> getOverlappingGenes(Region currReg, Collection<AGene> genes){
		List<AGene> res = new ArrayList<AGene>();
		for(AGene g : genes){
			if(g.getCoords().overlaps(currReg))
				res.add(g);
		}
		return(res);
	}
	
	/**
	 * Get all genes as regions 
	 * @param genes A set of known genes  
	 * @return A HashMap of ARegion collections indexed by chromosome name
	 */
	public static Map<String, List<Region>> getAllGeneRegions(Collection<AGene> genes){
		Map<String, List<Region>> accepted = new HashMap<String, List<Region>>();
		
		int geneCount=0;
		for(AGene currGene : genes){
			if(!accepted.containsKey(currGene.getCoords().getChrom()))
				accepted.put(currGene.getCoords().getChrom(), new ArrayList<Region>());
			accepted.get(currGene.getCoords().getChrom()).add(currGene.getCoords());
			geneCount++;
		}
		int regionCount=0;
		for(String c: accepted.keySet()){
			Collections.sort(accepted.get(c));
			regionCount+=accepted.get(c).size();
		}
		//System.err.println(regionCount+" genes");
		return(accepted);
	}
	
	/**
	 * Get all exons from all genes 
	 * @param genes A set of known genes  
	 * @return A HashMap of ARegion collections indexed by chromosome name
	 */
	public static Map<String, List<Region>> getAllExons(Collection<AGene> genes){
		Map<String, List<Region>> accepted = new HashMap<String, List<Region>>();
		
		int geneCount=0;
		for(AGene currGene : genes){
			ArrayList<AIsoform> currIsos = currGene.getIsoforms();
			for(AIsoform iso : currIsos){
				for(ARegion exon : iso.getExons()){
					if(!accepted.containsKey(exon.getCoords().getChrom()))
						accepted.put(exon.getCoords().getChrom(), new ArrayList<Region>());
					accepted.get(exon.getCoords().getChrom()).add(exon.getCoords());					
				}
			}
			geneCount++;
		}
		int regionCount=0;
		for(String c: accepted.keySet()){
			Collections.sort(accepted.get(c));
			regionCount+=accepted.get(c).size();
		}
		System.err.println(regionCount+" exons from all genes");
		return(accepted);
	}
	
	/**
	 * Get non-UTR exons from non-overlapping single-isoform genes
	 * @param genes A set of known genes  
	 * @return A HashMap of ARegion collections indexed by chromosome name
	 */
	public static Map<String, List<Region>> getNonOverlappingMiddleExons(Collection<AGene> genes){
		Map<String, List<Region>> accepted = new HashMap<String, List<Region>>();
		
		AGene lastGene=null, currGene=null;
		int geneCount=0;
		for(AGene nextGene : genes){
			//Special case for the first gene
			if(geneCount==1 && (!currGene.getCoords().overlaps(nextGene.getCoords()))){
				ArrayList<AIsoform> currIsos = currGene.getIsoforms();
				if(currIsos.size()==1){
					AIsoform iso = currIsos.get(0);
					for(ARegion exon : iso.getExons()){
						if(!exon.equals(iso.getRegion3Prime()) && !exon.equals(iso.getRegion5Prime())){
							if(!accepted.containsKey(exon.getCoords().getChrom()))
								accepted.put(exon.getCoords().getChrom(), new ArrayList<Region>());
							accepted.get(exon.getCoords().getChrom()).add(exon.getCoords());
						}
					}
				}
			}
			//Typical case
			if(geneCount>=2 && (!currGene.getCoords().overlaps(lastGene.getCoords()) && !currGene.getCoords().overlaps(nextGene.getCoords()))){
				ArrayList<AIsoform> currIsos = currGene.getIsoforms();
				if(currIsos.size()==1){
					AIsoform iso = currIsos.get(0);
					for(ARegion exon : iso.getExons()){
						if(!exon.equals(iso.getRegion3Prime()) && !exon.equals(iso.getRegion5Prime())){
							if(!accepted.containsKey(exon.getCoords().getChrom()))
								accepted.put(exon.getCoords().getChrom(), new ArrayList<Region>());
							accepted.get(exon.getCoords().getChrom()).add(exon.getCoords());
						}
					}
				}
			}
		
			//Special case for the last gene
			if(geneCount>2 && geneCount==genes.size()-1 && (!currGene.getCoords().overlaps(nextGene.getCoords()))){
				ArrayList<AIsoform> currIsos = nextGene.getIsoforms();
				if(currIsos.size()==1){
					AIsoform iso = currIsos.get(0);
					for(ARegion exon : iso.getExons()){
						if(!exon.equals(iso.getRegion3Prime()) && !exon.equals(iso.getRegion5Prime())){
							if(!accepted.containsKey(exon.getCoords().getChrom()))
								accepted.put(exon.getCoords().getChrom(), new ArrayList<Region>());
							accepted.get(exon.getCoords().getChrom()).add(exon.getCoords());
						}
					}
				}
			}
			lastGene=currGene;
			currGene = nextGene;
			geneCount++;
		}
		int regionCount=0;
		for(String c: accepted.keySet()){
			Collections.sort(accepted.get(c));
			regionCount+=accepted.get(c).size();
		}
		System.err.println(regionCount+" middle exons from non-overlapping genes");
		return(accepted);
	}
	
	/**
	 * Get 3' UTRs from a set of single-isoform genes
	 * @param genes	A set of known genes.
	 * @param minLength	Minimum length of the accepted exons.   
	 * @return A collection of ARegions
	 */
	public static Collection<ARegion> getSingleIsoform3UTRs(Collection<AGene> genes, int minLength){
		ArrayList<ARegion> accepted = new ArrayList<ARegion>();
		
		AGene lastGene=null, currGene=null;
		int geneCount=0;
		for(AGene nextGene : genes){
			if(geneCount>=2){
				ArrayList<AIsoform> currIsos = currGene.getIsoforms();
				if(currIsos.size()==1){
					AIsoform iso = currIsos.get(0);
					if(iso.getRegion3Prime().getCoords().getWidth()>=minLength)
						accepted.add(iso.getRegion3Prime());
				}
			}
			lastGene=currGene;
			currGene = nextGene;
			geneCount++;
		}
		return(accepted);
	}
	
}
