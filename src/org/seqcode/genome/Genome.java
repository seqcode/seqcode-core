package org.seqcode.genome;


import java.util.*;
import java.io.*;
import java.sql.*;

import org.seqcode.data.connections.DatabaseConnectionManager;
import org.seqcode.data.connections.DatabaseException;
import org.seqcode.data.connections.Sequence;
import org.seqcode.data.connections.UnknownRoleException;
import org.seqcode.genome.location.ChromosomeInfo;
import org.seqcode.utils.*;

/**
 * Genome represents one version (or genome build) of some species.
 * <i>Note</i>: We assume 1-based, inclusive coordinate.
 */
public class Genome{
	//static cache of all Genomes
	private static Map<String,Genome> staticGenomes = new HashMap<String,Genome>();
    private static Map<Integer, Genome> genomeids = new HashMap<Integer,Genome>();
    
    private Species species;
    private String version;
    private int dbid;
    private Map<String,ChromosomeInfo> chromsByName;
    private Map<Integer,ChromosomeInfo> chromsByID;
    
    
    /**
     * Constructs a new Genome from a Species and a genome version.
     */
    public Genome(Species species, String version) throws NotFoundException {
        this.species = species;
        this.version = version;
        Connection cxn = null;
        Statement stmt = null;
        ResultSet rs = null;
        try {
            cxn = DatabaseConnectionManager.getConnection("core");       
            stmt = cxn.createStatement();
            
            rs = stmt.executeQuery("select id from genome where species = " + species.getDBID() + 
                                             " and version ='" + version + "'");
            
            if (rs.next()) {
                dbid = rs.getInt(1);
            } else {
                throw new NotFoundException("Couldn't find " + species.getName() + " version " + version);
            }
            
            fillChroms(cxn);
            
        } catch (SQLException ex) {
            ex.printStackTrace();
            throw new DatabaseException("Couldn't find " + species + ": "+ ex.toString(),ex);
        }  catch (UnknownRoleException ex) {
            ex.printStackTrace();
            throw new DatabaseException("Couldn't connect with role core");
        } finally {
        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
        	if (stmt != null) { try { stmt.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role core", ex); }
        }
    }
       
    /**
     * Constructs a new Genome from almost complete info
     * 	Only used when a connection is open already for populating chromosome lengths
     *  Only used within the static methods below to populate a table of all genomes
     */
    private Genome(Species s, int dbid, String version, Connection cxn) throws NotFoundException {
        this.species = s;
        this.version = version;
        this.dbid = dbid;
        
        try {
			fillChroms(cxn);
		} catch (SQLException e) {
			e.printStackTrace();
		}
    }
    
    /**
     * Construct a Genome from a file of chromosome lengths
     * @param tempName
     * @param chrLengths
     * @param inventids
     */
    public Genome(String tempName, File chrLengths, boolean inventids) {
    	species = new Species(-1, "FakeOrganism");
    	version = tempName;
    	dbid = -1;
    	chromsByName = new HashMap<String,ChromosomeInfo>();
    	chromsByID = new HashMap<Integer,ChromosomeInfo>();
    	if(!chrLengths.isFile()){System.err.println("Invalid genome info file name");System.exit(1);}
        BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(chrLengths));
		    String line;
	        int id=0;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(words.length>=2){
	            	String chr = words[0].replaceFirst("^chromosome", "");
	            	chr = chr.replaceFirst("^chrom", "");
	            	chr = chr.replaceFirst("^chr", "");
	            	ChromosomeInfo info;
	            	if (inventids) {
	            		info = new ChromosomeInfo(id++, Integer.parseInt(words[1]), chr);
	            	} else {
	            		info = new ChromosomeInfo(Integer.parseInt(words[2]), Integer.parseInt(words[1]), chr);
	            	}
	            	chromsByName.put(info.getName(), info);
	            	chromsByID.put(info.getDBID(), info);
	            }
	    	}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
    
    /**
     * Construct a genome from a Map of names and lengths
     *  (mostly used to merge fake genomes that are data generated)
     * @param tempName
     * @param chrLengthMap
     */
    public Genome(String tempName, Map<String, Integer> chrLengthMap) {
    	species = new Species(-1, "FakeOrganism");
    	version = tempName;
    	dbid = -1;
    	chromsByName = new HashMap<String,ChromosomeInfo>();
    	chromsByID = new HashMap<Integer,ChromosomeInfo>();
    	int id=0;
    	for(String s : chrLengthMap.keySet()){
    		ChromosomeInfo info = new ChromosomeInfo(id--, chrLengthMap.get(s), s);
        	chromsByName.put(info.getName(), info);
        	chromsByID.put(info.getDBID(), info);
    	}
    }
    
    /**
     * Retrieves the chromosomes for this Genome from the database and fills
     * the relevant data structures: chroms and chromsByID
     * @param cxn
     * @throws SQLException
     */
    private void fillChroms(Connection cxn) throws SQLException {
        chromsByName = new HashMap<String,ChromosomeInfo>();
        chromsByID = new HashMap<Integer,ChromosomeInfo>();

        Statement stmt = cxn.createStatement();
        ResultSet rs = stmt.executeQuery("select c.id, c.name, cs.len from chromosome c, chromsequence cs " +
                "where c.id=cs.id and c.genome=" + dbid);
        while(rs.next()) { 
            int dbid = rs.getInt(1);
            String name = rs.getString(2);
            int length = rs.getInt(3);
            
            ChromosomeInfo info = new ChromosomeInfo(dbid, length, name);
            if(chromsByName.containsKey(name) || chromsByID.containsKey(dbid)) { 
                throw new IllegalArgumentException("Duplicate name \"" + name + 
                        "\" seems to exist in genome " + version);
            }
            chromsByName.put(name, info);
            chromsByID.put(dbid, info);
        }
        
        if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
        if (stmt != null) { try { stmt.close();} catch (SQLException ex) { } }
    	
    }

    //Accessors
    public String getVersion() {return version;}
    public int getDBID() {return dbid;}
    public Species getSpecies(){return species;}
    public String getSpeciesName() {return species.getName();}    
    public int getSpeciesDBID() { return species.getDBID(); }
    public String toString() { return getSpeciesName() +","+getVersion(); }

    //Chromosome-related accessors
    public Collection<ChromosomeInfo> getChromInfo(){ if(chromsByName!=null){return chromsByName.values();}else{return null;}}
    public List<String> getChromList() { return new LinkedList<String>(chromsByName.keySet()); }
    public ChromosomeInfo getChrom(String name) { return chromsByName.get(name); }
    public boolean containsChromName(String chromName) { return chromsByName.containsKey(chromName); }
    public String getChromName(int chromID) { return chromsByID.get(chromID).getName(); }
    public int getChromID(String chromName) { 
        if (chromsByName.get(chromName) == null) {
            throw new NullPointerException("Null chromosome for " + chromName);
        }
        return chromsByName.get(chromName).getDBID(); 
    }
    public int getChromLength(String chromName) { return chromsByName.get(chromName).getLength(); }
    public Map<String,Integer> getChromLengthMap() { 
        Map<String,Integer> chromLengths = new HashMap<String,Integer>();
        for(String n : chromsByName.keySet()) { chromLengths.put(n, chromsByName.get(n).getLength()); }
        return chromLengths;
    }
    
    /** Returns the genome info string with chromosome name <tab> length format*/
    public String getGenomeInfo(){
    	StringBuilder sb = new StringBuilder();
    	for(String n : chromsByName.keySet()) { sb.append(n).append("\t").append(chromsByName.get(n).getLength()).append("\n"); }
        return sb.toString();
    }

    /**
     * Return total length of all chromosomes
     * @return
     */
    public double getGenomeLength() { 
        double totalLen=0;
        for(String n : chromsByName.keySet()) { totalLen+= (double)chromsByName.get(n).getLength();}
        return totalLen;
    }

    
    //Roman numeral to integer translation helpers
	private static int[] romvals;
    private static String[] intvals;
    public static String convertChromNameToRoman(String c) {
        return convertChromNameToRoman(Integer.parseInt(c));
    }
    public static String convertChromNameToRoman(int chrom) {
        if(intvals == null) { 
            intvals = new String[10];
            intvals[0] = "X";
            intvals[1] = "I";
            intvals[2] = "II";
            intvals[3] = "III";
            intvals[4] = "IV";
            intvals[5] = "V";
            intvals[6] = "VI";
            intvals[7] = "VII";
            intvals[8] = "VIII";
            intvals[9] = "IX";
        }
        StringBuilder sb = new StringBuilder();
        sb.append("chr");
        while(chrom >= 10) { 
            chrom -= 10;
            sb.append(intvals[0]);
        }
        if(chrom > 0)
            sb.append(intvals[chrom]);
        return sb.toString();
    }
    public static String convertChromNameFromRoman(String chrom) {
        if (romvals == null) {
            romvals = new int[Character.getNumericValue('Z')];
            romvals[Character.getNumericValue('X')] = 10;
            romvals[Character.getNumericValue('V')] = 5;
            romvals[Character.getNumericValue('I')] = 1;
        }
        String chr = chrom;
        chr = chr.replaceAll("\\.fa?s?t?a$","");
        if (chr.matches("^[cC][hH][rR].*")) {
            chr = chr.substring(3);
        } 

        if (chr.matches("^[1234567890MmtUnXY]+(_random)?[LRh]?$")) {
            return chr;
        } else {
            throw new NumberFormatException("Can't fix chrom name " + chrom + "," + chr);
        }
    }
    public static String convertYeastChromNameFromRoman(String chrom) {
        if (romvals == null) {
            romvals = new int[Character.getNumericValue('Z')];
            romvals[Character.getNumericValue('X')] = 10;
            romvals[Character.getNumericValue('V')] = 5;
            romvals[Character.getNumericValue('I')] = 1;
        }
        String chr = chrom;
        chr = chr.replaceAll("\\.fa?s?t?a$","");
        if (chr.matches("^[cC][hH][rR].*")) {
            chr = chr.substring(3);
        } 
        if (chr.matches("^[XVI]+$")) {
            int val = 0, pos = 1, curval, lastval, buffer; char cur, last;
            boolean random = false;
            if (chr.matches("_random$")) {
                random = true;
                chr.replaceFirst("_random$","");
            }            
            last = chr.charAt(0);
            lastval = romvals[Character.getNumericValue(last)];
            buffer = lastval;
            //            System.err.println("== " + buffer);
            while (pos < chr.length()) {
                cur = chr.charAt(pos);
                curval = romvals[Character.getNumericValue(cur)];
                if (curval > lastval) {
                    val += curval - lastval;
                    buffer = 0;
                } else if (cur != last) {
                    val += buffer;
                    buffer = curval;
                } else {
                    buffer += curval;
                }
                last = cur;
                lastval = curval;
                pos++;
            }
            val += buffer;
            if (random) {
                return Integer.toString(val) + "_random";
            } else {
                return Integer.toString(val);
            }
        } else 
        if (chr.matches("^[1234567890MUXY]+(_random)?[LRh]?$")) {
            return chr;
        } else if (chr.matches("Mito")) {
            return "mt";
        } else {
            throw new NumberFormatException("Can't fix chrom name " + chrom + "," + chr);
        }
    }

    
    
    /** 
     * Returns a read connection to the annotation database for this genome 
     */
    public Connection getAnnotationDBConnection() throws SQLException {
        try {
        	String v = this.getVersion().replaceAll("[^\\w\\-]","_");
        	//We should store these table names in core and load at runtime
            return DatabaseConnectionManager.getConnection("ucsc_" + v);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't create a database connection for genome " + 
                                        getVersion(),ex);
        }

    }
        
    
    /**
	 * Load all Genomes from database 
	 * @return
	 */
    public static Collection<Genome> getAllGenomes(boolean forceRefreshFromDB){
    	List<Genome> gens = new ArrayList<Genome>();
    	
    	if(staticGenomes.isEmpty() || forceRefreshFromDB){
    		staticGenomes.clear(); genomeids.clear();
	    	Connection cxn = null;
	        Statement stmt = null;
	        ResultSet rs = null;
	        try {
	            cxn = DatabaseConnectionManager.getConnection("core");
	            stmt = cxn.createStatement();
	            rs = stmt.executeQuery("select id, species, version from genome");
	            
	            while(rs.next()) { 
	            	Genome gen = new Genome(Species.getSpecies(rs.getInt(2)), rs.getInt(1), rs.getString(3), cxn);
	            	gens.add(gen);
	            	staticGenomes.put(gen.getVersion(), gen);
	            	genomeids.put(gen.getDBID(), gen);
	            }
	
	        } catch (SQLException ex) {
	            ex.printStackTrace();
	            throw new DatabaseException("mySQL error: " + ex.toString(), ex);
	        } catch (UnknownRoleException ex) {
	            ex.printStackTrace();
	            throw new DatabaseException("Couldn't connect with role core", ex);
	        } catch (NotFoundException e) {
				e.printStackTrace();
			} finally {
	        	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
		        if (stmt != null) { try { stmt.close();} catch (SQLException ex) { } }
	        	if (cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role core", ex); }
	        }
    	}else{
    		gens.addAll(staticGenomes.values());
    	}
    	return gens;
    }
    
    /**
     * Return all Genomes for a given Species
	 * @param genomeName
	 * @return
	 * @throws NotFoundException
	 */
	public static Collection<Genome> getAllGenomesBySpecies(Species s) throws NotFoundException {
		if(staticGenomes.isEmpty()){
        	getAllGenomes(true);
        }
		List<Genome> sGens = new ArrayList<Genome>();
		for(Genome g : staticGenomes.values()){
			if(g.getSpeciesName().equals(s.getName()))
				sGens.add(g);
		}
		return sGens;
	}
    
    
    /**
	 * @param gid
	 * @return
	 * @throws NotFoundException
	 */
	public static Genome findGenome(int gid) throws NotFoundException {
		if(staticGenomes.isEmpty()){
        	getAllGenomes(false);
        }
		if (genomeids.containsKey(gid)) {
	        return genomeids.get(gid);
	    }
	
	    Connection cxn=null;
	    Statement stmt = null;
	    ResultSet rs = null;
	    
	    try {
	        cxn = DatabaseConnectionManager.getConnection("core");
	        stmt = cxn.createStatement();
	        rs = stmt.executeQuery("select version, species from genome where id=" + gid);
	        Genome g = null;
	
	        if (rs.next()) {
	            String genomeName = rs.getString(1);
	            int orgID = rs.getInt(2);
	            Species org = Species.getSpecies(orgID);
		        g = new Genome(org, genomeName);
	        }
	
	        if (g == null) {
	            throw new NotFoundException("Couldn't find genome: " + gid);
	        }
	        genomeids.put(gid,g);
	        staticGenomes.put(g.getSpeciesName(), g);
	        return g;
	
	    } catch (SQLException se) {
	        throw new DatabaseException("SQLException: " + se.getMessage(), se);
	    } catch (UnknownRoleException ex) {
	        throw new DatabaseException("Couldn't connect with role core", ex);
	    } finally {
	    	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (stmt != null) { try { stmt.close();} catch (SQLException ex) { } }
	    	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role core", ex); }
	    }
	}

	/**
	 * @param genomeName
	 * @return
	 * @throws NotFoundException
	 */
	public static Genome findGenome(String genomeName) throws NotFoundException {
		if(staticGenomes.isEmpty()){
        	getAllGenomes(false);
        }
		if (staticGenomes.containsKey(genomeName)) {
	        return staticGenomes.get(genomeName);
	    }
	    Connection cxn=null;
	    Statement stmt = null;
	    ResultSet rs = null;
	    
	    try {
	        cxn = DatabaseConnectionManager.getConnection("core");
	        stmt = cxn.createStatement();
	        rs = stmt.executeQuery("select species from genome where version='" + genomeName + "'");
	        Genome g = null;
	
	        if (rs.next()) {
	            int orgID = rs.getInt(1);
	            Species org = Species.getSpecies(orgID);
		        g = new Genome(org, genomeName);
	        }
	
	        if (g == null) {
	            throw new NotFoundException("Couldn't find genome: " + genomeName);
	        }
	        staticGenomes.put(genomeName, g);
	        genomeids.put(g.getDBID(), g);
	        return g;
	
	    } catch (SQLException se) {
	        throw new DatabaseException("SQLException: " + se.getMessage(), se);
	    } catch (UnknownRoleException ex) {
	        throw new DatabaseException("Couldn't connect with role core", ex);
	    } finally {
	    	if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (stmt != null) { try { stmt.close();} catch (SQLException ex) { } }
	    	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role core", ex); }
	    }
	}
	
	/**
     * Returns all of the versions/builds for this species.
     * @return
     */
    public static Collection<String> getAllGenomeNames(boolean forceRefreshFromDB) {
    	if(staticGenomes.isEmpty() || forceRefreshFromDB){
        	getAllGenomes(forceRefreshFromDB);
        }
    	return staticGenomes.keySet();
    }

	/**
     * Insert a new Genome into the database 
     * @param version
     * @throws SQLException
     */
    public static void insertGenome(Species species, String version) throws SQLException {
        Connection cxn = null;
        Statement stmt = null;
        try {
            cxn = DatabaseConnectionManager.getConnection("core");
            stmt = cxn.createStatement();
            String nextIdString = Sequence.getInsertSQL(cxn, "genome_id");
            String insertSQL = String.format("insert into genome(id, species, version) values (%s, %d, '%s')", nextIdString,
            		species.getDBID(), version);
            stmt.executeUpdate(insertSQL);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role core", ex);
        } catch (SQLException se) {
            throw se;
        } finally {
            if (stmt != null) { try { stmt.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role core", ex); }
        }
    }
    
	public int hashCode() {
        return getSpeciesName().hashCode()*37 + getVersion().hashCode();
    }

    public boolean equals(Object o) {
        if (o instanceof Genome) {
            Genome other = (Genome)o;
            return (getSpeciesName().equals(other.getSpeciesName()) &&
                    getVersion().equals(other.getVersion()));
        } else {
            return false;
        }
    }

}

