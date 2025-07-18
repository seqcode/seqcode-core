package org.seqcode.genome.sequence;

import java.io.IOException;
import java.io.File;
import java.util.*;
import java.sql.*;

import org.seqcode.data.connections.DatabaseConnectionManager;
import org.seqcode.data.connections.DatabaseException;
import org.seqcode.data.connections.UnknownRoleException;
import org.seqcode.data.io.parsing.FASTAStream;
import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gsebricks.types.*;
import org.seqcode.gsebricks.verbs.Mapper;
import org.seqcode.gseutils.*;


/** 
 * <code>SequenceGenerator</code> maps a Region to the genomic
 * sequence included in that Region.
 * 
 * 1-based genome
 */
public class SequenceGenerator<X extends Region> implements Mapper<X,String>, SelfDescribingVerb {

	private Genome genome;
    private static Map<Integer,String> cache;
    private boolean useCache = false;
    private boolean useLocalFiles = false;
    private String genomePath = null;
    private boolean genomePathIsFullGenomeFile=false;
    private int maxQuery = -1;

    private static Map<String, String[]> regionCache;
    private static Map<String, int[]> regionStarts;
    private static boolean regionIsCached = false;
    public boolean isRegionCached(){return regionIsCached;}
    public boolean usingLocalFiles(){return useLocalFiles;}
    
    /**
     * SequenceGenerator
     * maxQuery is here to avoid calling the sequence generator if the Region width is greater 
     * than a given size. This is really only useful for SeqView painters.   
     * @param g
     */
    public SequenceGenerator (Genome g) {this(g, -1);}
    public SequenceGenerator (Genome g, int maxQuery) {
    	genome = g;
    	this.maxQuery=maxQuery;
    }
    public SequenceGenerator() {}
    public void useCache(boolean b) {
        if (b && cache == null) {
            cache = new HashMap<Integer,String>();
        } 
        useCache = b;
    }
    public void useLocalFiles(boolean b) {
        useLocalFiles = b;
    }
    public void setGenomePath(String genomePath){
    	this.genomePath = genomePath;
    	File test = new File(genomePath);
    	if(!test.exists()){
    		System.err.println("Genome sequence directory/file does not exist!");
    		System.exit(1);
    	}else{
    		// Skip the check for a file extension for galaxy
    		//if(test.isFile() && (genomePath.endsWith(".fa")||genomePath.endsWith(".fasta")||genomePath.endsWith(".seq")))
    		//	genomePathIsFullGenomeFile=true;
    		if(test.isFile())
    			genomePathIsFullGenomeFile=true;
    	}
    }
    
    /** cache the whole chromosome of this region */
    private void cache(X region) throws SQLException, IOException {
        int chromid = region.getGenome().getChromID(region.getChrom());
        synchronized(cache) {
            if (cache.containsKey(chromid)) {
                return;
            }
        }
        String chromseq = null;
        if (useLocalFiles) {
        	if (genomePath==null)
        		genomePath = "/scratch/" + region.getGenome().getVersion();
        	File f=null;
            if(genomePathIsFullGenomeFile){
            	f = new File(genomePath);
            }else{
	        	f = new File( genomePath + "/chr" + region.getChrom() + ".fa");
	            if (!f.exists())
	                f = new File( genomePath+ "/chr" + region.getChrom() + ".fasta");
	            if (!f.exists())
	                f = new File( genomePath+ "/" + region.getChrom() + ".fa");
	            if (!f.exists())
	                f = new File( genomePath+ "/" + region.getChrom() + ".fasta");
	            if (!f.exists())
	                f = new File( genomePath+ "/chromosome" + region.getChrom() + ".fa");
	            if (!f.exists())
	                f = new File( genomePath+ "/chromosome" + region.getChrom() + ".fasta");
	            if (!f.exists())
	                f = new File( genomePath+ "/chrom" + region.getChrom() + ".fa");
	            if (!f.exists())
	                f = new File( genomePath+ "/chrom" + region.getChrom() + ".fasta");
            }
            
            if (f.exists()) {
                FASTAStream stream = new FASTAStream(f);
                while (stream.hasNext()) {
                    Pair<String,String> pair = stream.next();
                    String pairchrom = pair.car().replaceFirst("^chromosome", "").replaceFirst("^chrom", "").replaceFirst("^chr","");
                    if (pairchrom.equals(region.getChrom())) {
                        chromseq = pair.cdr();
                        break;
                    }                    
                }
                stream.close();                
            }else{
            	System.out.print("FASTA file for chromosome "+region.getChrom() +" is not found at "+genomePath+". \n");
            	System.exit(-1);
            }
            
        } else {
            java.sql.Connection cxn = DatabaseConnectionManager.getConnection("core");
            PreparedStatement ps = cxn.prepareStatement("select sequence from chromsequence where id = ?");
            ps.setInt(1,chromid);
            ResultSet rs = ps.executeQuery();
            if (rs.next()) {
                chromseq = rs.getString(1);
            }   
            rs.close();
            ps.close();
            if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role core", ex); }
        }
        if (chromseq == null) {
            return;
        }
        synchronized(cache) {
            if (!cache.containsKey(chromid)) {
                cache.put(chromid, chromseq);
            }
        }
    }
    /**
     * get sequence of specified region (including start and end)
     */
    public String execute(X region) {
    	String result = null;
    	if(maxQuery==-1 || region.getWidth()<=maxQuery){
	    	if (regionIsCached)
	    		return getRegionCacheSequence(region);  
	    	
	    	String chromname = region.getChrom();
	        
	        try {
	            Genome genome = region.getGenome();
	            int chromid = genome.getChromID(chromname);
	            if (useCache) {
	                cache(region);
	                String chromString = null;
	                synchronized(cache) {
	                    if (!cache.containsKey(chromid)) {
	                        return null; 
	                    }
	                    chromString = cache.get(chromid);
	                }
	                //1-based version (string.substr is 0-based) 
	                result = chromString.substring(region.getStart()-1, region.getEnd());
	            }
	            if (result == null) {
	            	if(useLocalFiles) {
	            		System.out.println("Cannot find this region in provided genome files: "+region.getLocationString());
	    	        	System.exit(-1);
	            	}else {
		                java.sql.Connection cxn =
		                DatabaseConnectionManager.getConnection("core");
		                PreparedStatement ps;
		                //1-based version (mysql substr is 1-based)
		                int start = Math.max(region.getStart(), 0);
		                ps = cxn.prepareStatement("select substr(sequence,?,?) from chromsequence where id = ?");
		                ps.setInt(1,start);
		                ps.setInt(2,region.getEnd() - region.getStart() + 1);
		                ps.setInt(3,chromid);                   
		                ResultSet rs = ps.executeQuery();
		                if (rs.next()) {
		                    result = rs.getString(1);
		                } 
		                rs.close();
		                ps.close();
		                if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role core", ex); }
	            	}
	            }
	        } catch (SQLException ex) {
	            ex.printStackTrace();           
	        } catch (UnknownRoleException ex) {
	            ex.printStackTrace();
	            throw new DatabaseException("Couldn't connect to core",ex);
	        } catch (IOException ex) {
	            ex.printStackTrace();
	            throw new RuntimeException("Couldn't load file to cache " + ex.toString(), ex);
	        } catch (StringIndexOutOfBoundsException ex){
	        	System.err.println("Error in "+region.getLocationString());
	        	ex.printStackTrace();
	        }
	
	        if (result == null) {
	            throw new DatabaseException("Couldn't get any sequence for " + region);
	        }
	
	        if (result.length() != region.getWidth()) {
	            System.err.println("Wanted " + region + "(" + 
	            		region.getWidth() + ") but only got " + result.length());
	        }
    	}
        return result;
    }
    
    /**
     * Setup light-weight region cache of genome sequences, cover only the specified regions<br>
     * So that it does not cache the whole chromosome, save memory space. <br>
     * At the same time, retrieve some one-time sequences in rs.
     * @param regions sorted, non-overlapping regions for cache
     * @param rs regions for one-time sequence retrieval
     */
    public String[] setupRegionCache(List<Region> regions, List<Region> rs){
    	ArrayList<String> seqs = new ArrayList<String>();
    	if (regions==null||regions.isEmpty())
    		return null;    	
    	Collections.sort(regions);
    	
		// group one-time regions by chrom
		Map<String, ArrayList<Region>> chr2rs = new HashMap<String, ArrayList<Region>>();
		for(Region r:rs) {
			String chrom = r.getChrom();
			if(!chr2rs.containsKey(chrom))
				chr2rs.put(chrom, new ArrayList<Region>());
			chr2rs.get(chrom).add(r);
		}
    	
    	useCache(true);
    	regionCache = new HashMap<String, String[]>();
    	regionStarts = new HashMap<String, int[]>();
    	Genome g = regions.get(0).getGenome();
    	Region lastRegion = regions.get(regions.size()-1);
    	// setup the region space
    	String chrom = regions.get(0).getChrom();
    	int count = 0;
    	for(Region r: regions){
    		if (!r.getChrom().equals(chrom)){		// new chrom
    			regionCache.put(chrom, new String[count]);
    			regionStarts.put(chrom, new int[count]);
    			chrom = r.getChrom();
    			count = 1;
    		}
    		else		// same chrom
        		count ++;
    	}
		regionCache.put(chrom, new String[count]);
		regionStarts.put(chrom, new int[count]);
    	
		chrom = regions.get(0).getChrom();
    	count = 0;
    	for (int i=0;i<regions.size();i++){
    		Region r = regions.get(i);
    		if (!r.getChrom().equals(chrom)){	// new Chrom
    			if (cache!=null){							// for previous chrom
    				if (chr2rs.containsKey(chrom)){			// piggy-back to retrieve one-time sequences
    					for (Region r1:chr2rs.get(chrom))
    						seqs.add(execute((X)r1));
    				}
	    			synchronized(cache) {
	    				cache.put(g.getChromID(chrom), null);
	    				cache.remove(g.getChromID(chrom));	// clean cache
	    			}
	    	    	System.gc();
    			}
    			chrom = r.getChrom();
    			count = 0;
    		}		
    		// cache region using the current cache
    		synchronized(regionCache) {
    			regionStarts.get(chrom)[count]=r.getStart(); 
        		regionCache.get(chrom)[count]=execute((X)r);    			
    		}
    		count ++;
    	}
    	if (cache!=null){
			if (chr2rs.containsKey(lastRegion.getChrom())){			// piggy-back to retrieve one-time sequences
				for (Region r1:chr2rs.get(lastRegion.getChrom()))
					seqs.add(execute((X)r1));
			}
	    	synchronized(cache) {
	    		cache.put(g.getChromID(lastRegion.getChrom()), null);
	    		cache.remove(g.getChromID(lastRegion.getChrom()));
	    	}
	    	// retrieve those regions that are not in the chromosomes of cache regions
	    	for(String chr:chr2rs.keySet()){
	    		if (!regionCache.containsKey(chr)){
	    			for (Region r1:chr2rs.get(chr))
						seqs.add(execute((X)r1));
	    			synchronized(cache) {
	    				cache.put(g.getChromID(chrom), null);
	    				cache.remove(g.getChromID(chrom));	// clean cache for this chrom
	    			}
	    			System.gc();
	    		}
	    	}
	    	cache=null;
	    	System.gc();
    	}
    	
    	regionIsCached = true;
    	String[] result = new String[seqs.size()];
    	for (int i=0;i<result.length;i++)
    		result[i]=seqs.get(i);
    	return result;
    }
    
    private String getRegionCacheSequence(Region r){
    	int[] starts = regionStarts.get(r.getChrom());
    	int idx = Arrays.binarySearch(starts, r.getStart());
    	if( idx < 0 ) { idx = -idx - 2; }
    	if (!regionCache.containsKey(r.getChrom()))
    		return null;
	    synchronized(regionCache) {
	    	try{
	    		return regionCache.get(r.getChrom())[idx].substring(r.getStart()-starts[idx], r.getEnd()-starts[idx]+1);
	    	}
	    	catch(Exception e){
	    		e.printStackTrace(System.out);
	    		System.out.println(r.toString()+" idx="+idx+" starts[idx]="+starts[idx]);
	    	    return null;
	    	}
    	}
    }
    
    public static void clearCache() {
        synchronized(cache) {
            cache.clear();
        }
    }
    
    // the following method has bee add by akshay
    public static void setOffRegionCache(){regionIsCached = false;}
    
    private static final String[] inputNames = { "Regions" };
    private static final EchoType[] inputTypes = { new ClassType(Region.class) };
    private static final EchoType outputType = new ClassType(String.class);

    public EchoType[] getInputClasses() { return inputTypes; }

    public String[] getInputNames() { return inputNames; }

    public EchoType getOutputClass() { return outputType; }

    public EchoType[] getParameterClasses() {
        return null;
    }

    public String[] getParameterNames() {
        return null;
    }

    public void init(Map<String, Object> params) {
    }    
}
