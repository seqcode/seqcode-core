package org.seqcode.data.seqdata;

import java.util.*;
import java.sql.*;
import java.io.*;

import org.seqcode.data.connections.DatabaseConnectionManager;
import org.seqcode.data.connections.DatabaseException;
import org.seqcode.data.core.CellLine;
import org.seqcode.data.core.ExptCondition;
import org.seqcode.data.core.ExptTarget;
import org.seqcode.data.core.Lab;
import org.seqcode.data.core.MetadataLoader;
import org.seqcode.data.core.SeqDataUser;
import org.seqcode.data.readdb.Client;
import org.seqcode.data.readdb.ClientException;
import org.seqcode.data.readdb.PairedHit;
import org.seqcode.data.readdb.SingleHit;
import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.utils.NotFoundException;
import org.seqcode.utils.Pair;



/**
 * SeqDataLoader serves as a clearinghouse for query interactions with the seqdata database and
 * associated metadata from tables in the core database (via the MetadataLoader). 
 * 
 * Implements a simple access control by checking if the username is in the permissions field
 * in each SeqAlignment. Note that while this is access-control-lite for experiment metadata, 
 * access to the data stored in the underlying readdb entries has more robust access-control.
 * 
 * @author tdanford
 * @author mahony
 * 
 * Created as ChipSeqLoader on May 16, 2007
 */
public class SeqDataLoader implements org.seqcode.utils.Closeable {

	public static String role = "seqdata";


	public static void main(String[] args) throws Exception{
		try {
			SeqDataLoader loader = new SeqDataLoader();
			
			Collection<SeqAlignment> aligns = loader.loadAllAlignments();
			for (SeqAlignment align : aligns) {
				SeqExpt expt = align.getExpt();
				System.out.println(expt.getDBID() + "\t" + expt.getName() + ";"+ expt.getReplicate()+"\t"+align.getName()+"\t"+align.getDBID()+"\t"+align.getGenome());
			}
			System.out.println("DONE");
			loader.close();
			DatabaseConnectionManager.close();
		}
		catch (SQLException e) {
			e.printStackTrace();
		}
	}
	private String myusername = "";
	private SeqDataUser myUser=null;
    private Client client=null;
    private MetadataLoader mloader=null;
    private boolean closed=false;
    
    /**
     * Accessor for user
     * @return
     */
    public SeqDataUser getMyUser(){return myUser;}    
    /**
     * Accessor for ReadDB client
     * @return
     */
    public Client getClient(){return client;}
    /**
     * Accessor for MetadataLoader
     * @return
     */
    public MetadataLoader getMetadataLoader(){return mloader;}
    
	//Constructors
    public SeqDataLoader() throws SQLException, IOException{this(true, true);}
	public SeqDataLoader(boolean openClient, boolean cacheAllMetadata) throws SQLException, IOException {
		mloader = new MetadataLoader(cacheAllMetadata); 
		
		if(openClient){
	        try {
	            client = new Client();
	            
	        } catch (ClientException e) {
	            throw new IllegalArgumentException(e);
	        }
		}
		
		myusername = DatabaseConnectionManager.getUsername(role);
		myUser = mloader.loadSeqDataUser(myusername, false, false);
        
		closed=false;
	}
    
	/**
	 * Load the genomes that a SeqExpt is aligned to
	 * @param expt
	 * @return
	 * @throws SQLException
	 */
	public Collection<Genome> loadExperimentGenomes(SeqExpt expt) throws SQLException {
		Connection cxn = null;
		Statement stmt = null;
        ResultSet rs = null;
        LinkedList<Genome> genomes = new LinkedList<Genome>();
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            String query = String.format("select genome from seqalignment where expt=%d", expt.getDBID());
            //System.out.println(query);
			stmt = cxn.createStatement();
			rs = stmt.executeQuery(query);
			while (rs.next()) {
				int gid = rs.getInt(1);
				try {
					Genome g = Genome.findGenome(gid);
					genomes.add(g);
				}
				catch (NotFoundException e) {
					e.printStackTrace();
				}
			}
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (stmt != null) { try { stmt.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return genomes;
	}

	/**
	 * Load all SeqExpts from the database. 
	 * @return
	 * @throws SQLException
	 */
	public Collection<SeqExpt> loadAllExperiments() throws SQLException {
		LinkedList<SeqExpt> expts = new LinkedList<SeqExpt>();
		Connection cxn=null;
		PreparedStatement stmt = null;
        ResultSet rs = null;
        try {
	        //Loads experiments
	        cxn = DatabaseConnectionManager.getConnection(role);
	        stmt = SeqExpt.createLoadAll(cxn);
	        //System.out.println(stmt.toString());
	        stmt.setFetchSize(20000);
			rs = stmt.executeQuery();
			while (rs.next())
				expts.addLast(new SeqExpt(rs, this));
			Collections.sort(expts);
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (stmt != null) { try { stmt.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return expts;
	}

	/**
	 * Find a SeqExpt by name and replicate
	 * (returns null if not found)
	 * @param name
	 * @param rep
	 * @return
	 * @throws SQLException
	 */
	public SeqExpt findExperiment(String name, String rep) throws SQLException {
		SeqExpt expt = null;
		Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            cxn = DatabaseConnectionManager.getConnection(role);
		    ps = SeqExpt.createLoadByNameReplicate(cxn);
			ps.setString(1, name);
			ps.setString(2, rep);
			//System.out.println(ps.toString());
			rs = ps.executeQuery();
			if (rs.next()) {
				expt = new SeqExpt(rs, this);
			}
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try { ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return expt;
	}
	
	/**
	 * Load a SeqExpt by name and replicate
	 * (throws exception if not found)
	 * @param name
	 * @param rep
	 * @return
	 * @throws NotFoundException
	 * @throws SQLException
	 */
	public SeqExpt loadExperiment(String name, String rep) throws NotFoundException, SQLException {
		SeqExpt expt = null;
		Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            cxn = DatabaseConnectionManager.getConnection(role);
			ps = SeqExpt.createLoadByNameReplicate(cxn);
			ps.setString(1, name);
			ps.setString(2, rep);
			//System.out.println(ps.toString());
			rs = ps.executeQuery();
			if (rs.next()) {
				expt = new SeqExpt(rs, this);
			}
			if (expt == null) { throw new NotFoundException(name+";"+rep); }
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try { ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}	
		return expt;
	}

	/**
	 * Load a collection of SeqExpts by name
	 * @param name
	 * @return
	 * @throws SQLException
	 */
	public Collection<SeqExpt> loadExperiments(String name) throws SQLException {
		SeqExpt expt = null;
		LinkedList<SeqExpt> expts = new LinkedList<SeqExpt>();
		Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            cxn = DatabaseConnectionManager.getConnection(role);
			ps = SeqExpt.createLoadByName(cxn);
			ps.setString(1, name);
			//System.out.println(ps.toString());
			rs = ps.executeQuery();
			while (rs.next()) {
				expt = new SeqExpt(rs, this);
				expts.add(expt);
			}
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try { ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return expts;
	}
	
	/**
	 * Load a collection of SeqExpts by Lab
	 * @param lab
	 * @return
	 * @throws SQLException
	 */
	public Collection<SeqExpt> loadExperiments(Lab lab) throws SQLException {
		LinkedList<SeqExpt> expts = new LinkedList<SeqExpt>();
		Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = SeqExpt.createLoadByLab(cxn);
			ps.setInt(1, lab.getDBID());
			//System.out.println(ps.toString());
			rs = ps.executeQuery();
			SeqExpt expt = null;
			while (rs.next()) {
				expt = new SeqExpt(rs, this);
				expts.add(expt);
			}
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try { ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return expts;
	}
	
	/**
	 * Load a collection of SeqExpts by ExptCondition
	 * @param cond
	 * @return
	 * @throws SQLException
	 */
	public Collection<SeqExpt> loadExperiments(ExptCondition cond) throws SQLException {
		LinkedList<SeqExpt> expts = new LinkedList<SeqExpt>();
		Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            cxn = DatabaseConnectionManager.getConnection(role);
            ps = SeqExpt.createLoadByCondition(cxn);
			ps.setInt(1, cond.getDBID());
			//System.out.println(ps.toString());
			rs = ps.executeQuery();
			SeqExpt expt = null;
			while (rs.next()) {
				expt = new SeqExpt(rs, this);
				expts.add(expt);
			}
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try { ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}	
		return expts;
	}
	
	/**
	 * Load a collection of SeqExpts by ExptTarget
	 * @param target
	 * @return
	 * @throws SQLException
	 */
	public Collection<SeqExpt> loadExperiments(ExptTarget target) throws SQLException {
		LinkedList<SeqExpt> expts = new LinkedList<SeqExpt>();
		Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            cxn = DatabaseConnectionManager.getConnection(role);
	        ps = SeqExpt.createLoadByTarget(cxn);
			ps.setInt(1, target.getDBID());
			//System.out.println(ps.toString());
			rs = ps.executeQuery();
			SeqExpt expt = null;
			while (rs.next()) {
				expt = new SeqExpt(rs, this);
				expts.add(expt);
			}
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try { ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return expts;
	}

	/**
	 * Load a collection of SeqExpts by CellLine
	 * @param cell
	 * @return
	 * @throws SQLException
	 */
	public Collection<SeqExpt> loadExperiments(CellLine cell) throws SQLException {
		LinkedList<SeqExpt> expts = new LinkedList<SeqExpt>();
		Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            cxn = DatabaseConnectionManager.getConnection(role);
		    ps = SeqExpt.createLoadByCellline(cxn);
			ps.setInt(1, cell.getDBID());
			//System.out.println(ps.toString());
			rs = ps.executeQuery();
			SeqExpt expt = null;
			while (rs.next()) {
				expt = new SeqExpt(rs, this);
				expts.add(expt);
			}
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try { ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}	
		return expts;
	}

	/**
	 * Load a SeqExpt by ID
	 * @param dbid
	 * @return
	 * @throws NotFoundException
	 * @throws SQLException
	 */
	public SeqExpt loadExperiment(int dbid) throws NotFoundException, SQLException {
		SeqExpt expt = null;
		Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            cxn = DatabaseConnectionManager.getConnection(role);
			ps = SeqExpt.createLoadByDBID(cxn);
			ps.setInt(1, dbid);
			//System.out.println(ps.toString());
			rs = ps.executeQuery();
			if (rs.next()) {
				expt = new SeqExpt(rs, this);
			}
	
			if (expt == null) {
				String err = String.format("No such SeqExpt experiment with ID= %d", dbid);
				throw new NotFoundException(err);
			}
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return expt;
	}

	/**
	 * Load all SeqAlignments in the database
	 * @param g
	 * @return
	 * @throws SQLException
	 */
    public Collection<SeqAlignment> loadAllAlignments () throws SQLException {
    	Collection<SeqAlignment> aligns = new LinkedList<SeqAlignment>();
    	Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            Collection<SeqExpt> allexpts = loadAllExperiments();
	        Map<Integer,SeqExpt> exptmap = new HashMap<Integer,SeqExpt>();
	        for (SeqExpt e : allexpts) {
	            exptmap.put(e.getDBID(), e);
	        }
	        
	        cxn = DatabaseConnectionManager.getConnection(role);
			ps = SeqAlignment.createLoadAllStatement(cxn);
	        ps.setFetchSize(20000);
	        //System.out.println(ps.toString());
	        rs = ps.executeQuery();
			while (rs.next()) {
	            SeqAlignment align = new SeqAlignment(rs, exptmap.get(rs.getInt(2)), mloader.loadAlignType(rs.getInt(6), false));
	            aligns.add(align);
			}
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return filterAlignmentsByPermission(aligns);
    }
	/**
	 * Load all SeqAlignments for a given Genome
	 * @param g
	 * @return
	 * @throws SQLException
	 */
    public Collection<SeqAlignment> loadAllAlignments (Genome g) throws SQLException {
    	Collection<SeqAlignment> aligns = new LinkedList<SeqAlignment>();
    	Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            Collection<SeqExpt> allexpts = loadAllExperiments();
	        Map<Integer,SeqExpt> exptmap = new HashMap<Integer,SeqExpt>();
	        for (SeqExpt e : allexpts) {
	            exptmap.put(e.getDBID(), e);
	        }
	        
	        cxn = DatabaseConnectionManager.getConnection(role);
			ps = SeqAlignment.createLoadAllByGenomeStatement(cxn);
	        ps.setFetchSize(20000);
			ps.setInt(1, g.getDBID());
			//System.out.println(ps.toString());
	        rs = ps.executeQuery();
			while (rs.next()) {
	            SeqAlignment align = new SeqAlignment(rs, exptmap.get(rs.getInt(2)), mloader.loadAlignType(rs.getInt(6), false));
	            aligns.add(align);
			}
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return filterAlignmentsByPermission(aligns);
    }

    /**
     * Load all SeqAlignments for a given SeqExpt
     * @param expt
     * @return
     * @throws SQLException
     */
	public Collection<SeqAlignment> loadAlignmentsBySeqExpt(SeqExpt expt) throws SQLException {
		Collection<SeqAlignment> aligns = new LinkedList<SeqAlignment>();
		Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            cxn = DatabaseConnectionManager.getConnection(role);
			ps = SeqAlignment.createLoadAllByExptStatement(cxn);
			ps.setInt(1, expt.getDBID());
			//System.out.println(ps.toString());
			rs = ps.executeQuery();
			while (rs.next()) {
				SeqAlignment align = new SeqAlignment(rs, expt, mloader.loadAlignType(rs.getInt(6), false));
				aligns.add(align);
			}
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return filterAlignmentsByPermission(aligns);
	}
	
	/**
	 * Load SeqAlignment by SeqExpt, name, genome
	 * @param expt
	 * @param n
	 * @param g
	 * @return
	 * @throws SQLException
	 */
	public SeqAlignment loadAlignment(SeqExpt expt, String n, Genome g) throws SQLException {
		SeqAlignment align = null;
		Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            cxn = DatabaseConnectionManager.getConnection(role);
			ps = SeqAlignment.createLoadByNameAndExptStatement(cxn);
			ps.setString(1, n);
			ps.setInt(2, expt.getDBID());
			//System.out.println(ps.toString());
			rs = ps.executeQuery();        
			while (align == null && rs.next()) {
				align = new SeqAlignment(rs, expt, mloader.loadAlignType(rs.getInt(6), false));
	            if (!align.getGenome().equals(g))
	                align = null;
	            if(align!=null && !align.getPermissions().contains("public") && !align.getPermissions().contains(myusername))
					align=null;
			}
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try { ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		if (align == null) {
        	//Don't throw exception because sometimes we have replicates that don't match all alignment names
            //throw new NotFoundException("Couldn't find alignment " + n + " for " + expt + " in genome " + g);
        }
        return align;

	}
	
	/**
	 * Load SeqAlignment by ID
	 * @param dbid
	 * @return
	 * @throws NotFoundException
	 * @throws SQLException
	 */
	public SeqAlignment loadAlignment(int dbid) throws NotFoundException, SQLException {
		SeqAlignment align = null;
		Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
			ps = SeqAlignment.createLoadByIDStatement(cxn);
			ps.setInt(1, dbid);
			//System.out.println(ps.toString());
			rs = ps.executeQuery();
			if (rs.next()) {
				align = new SeqAlignment(rs, this);
				if(!align.getPermissions().contains("public") && !align.getPermissions().contains(myusername))
					align=null;
			}
			else {
				throw new NotFoundException("Couldn't find alignment with ID = " + dbid);
			}
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return align;
	}

	/**
	 * Load a collection of SeqAlignments by SeqLocator
	 * @param locator
	 * @param genome
	 * @return
	 * @throws SQLException
	 * @throws NotFoundException
	 */
	public Collection<SeqAlignment> loadAlignments(SeqLocator locator, Genome genome) throws SQLException, NotFoundException {
		List<SeqAlignment> output = new ArrayList<SeqAlignment>();
        if(locator.getAlignName()!=null){ //Alignment name provided
			if (locator.getReplicates().size() == 0) { //No replicate names provided
	            for (SeqExpt expt : loadExperiments(locator.getExptName())) { 
	                SeqAlignment align = loadAlignment(expt, locator.getAlignName(), genome);
	                if (align != null) {
	                    output.add(align);
	                }
	            }
	        } else {
	            for (String rep : locator.getReplicates()) { //Replicate names provided
	                try {
	                    SeqExpt expt = loadExperiment(locator.getExptName(), rep);
	                    SeqAlignment align = loadAlignment(expt, locator.getAlignName(), genome);
	                    if (align != null) {
	                        output.add(align);
	                    }
	                }
	                catch (IllegalArgumentException e) {
	                    throw new NotFoundException("Couldn't find alignment for " + locator);
	                }
	            }
	        }
        }else{ //No alignment name provided
        	if (locator.getReplicates().size() == 0) { //No alignment or replicate names provided
	            for (SeqExpt expt : loadExperiments(locator.getExptName())) { 
	                Collection<SeqAlignment> aligns = loadAlignmentsBySeqExpt(expt);
	                for(SeqAlignment a : aligns){
	                	if (a.getGenome().equals(genome)) { 
            				output.add(a);
    						break;
    					}
	                }
	            }
	        } else {
	            for (String rep : locator.getReplicates()) { //Replicate names provided, but not alignment name
	                try {
	                    SeqExpt expt = loadExperiment(locator.getExptName(), rep);
	                    Collection<SeqAlignment> aligns = loadAlignmentsBySeqExpt(expt);
		                for(SeqAlignment a : aligns){
		                	if (a.getGenome().equals(genome)) { 
	            				output.add(a);
	    						break;
	    					}
		                }
	                }
	                catch (IllegalArgumentException e) {
	                    throw new NotFoundException("Couldn't find alignment for " + locator);
	                }
	            }
	        }
        }
        Collection<SeqAlignment> filtered = filterAlignmentsByPermission(output);
        if (filtered.size() == 0)
            throw new NotFoundException("Locators were " + toString() + " but didn't get any alignments that you have permission to see.");
        
		return filtered;
	}
     
	/**
	 * Load a collection of SeqAlignments by querying multiple metadata labels 
	 * @param name
	 * @param replicate
	 * @param align
	 * @param exptType
	 * @param lab
	 * @param condition
	 * @param target
	 * @param cells
	 * @param readType
	 * @param genome
	 * @return
	 * @throws SQLException
	 */
	public Collection<SeqAlignment> loadAlignments(String name, String replicate, String align,
                                                       Integer exptType, Integer lab, Integer condition,
                                                       Integer target, Integer cells, Integer readType,
                                                       Genome genome) throws SQLException {
        
		Collection<SeqAlignment> output = new ArrayList<SeqAlignment>();
		Connection cxn=null;
		PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            cxn = DatabaseConnectionManager.getConnection(role);
			String query = "select id, expt, name, genome, permissions, aligntype, numhits, totalweight, numtype2hits, totaltype2weight, numpairs, totalpairweight, aligndir, alignfile, idxfile, collabalignid  from seqalignment";
	        if (name != null || replicate != null || align != null || exptType!=null || lab!=null || condition!=null || target != null || cells != null || readType != null || genome != null) {
	            query += " where ";
	        }
	        boolean and = false;
	        if (name != null || replicate != null || exptType!=null || lab!=null || condition!=null || target != null || cells != null || readType != null) {
	            query += " expt in ( select id from seqexpt where ";
	            if (name != null) { query += " name = ? "; and = true;}
	            if (replicate != null) { query += (and ? " and " : " ") + " replicate = ? "; and = true;}
	            if (exptType != null) { query += (and ? " and " : " ") + " expttype = " + exptType; and = true;}
	            if (lab != null) { query += (and ? " and " : " ") + " lab = " + lab; and = true;}
	            if (condition != null) { query += (and ? " and " : " ") + " exptcondition = " + condition; and = true;}
	            if (target != null) { query += (and ? " and " : " ") + " expttarget = " + target; and = true;}
	            if (cells != null) { query += (and ? " and " : " ") + " cellline = " + cells; and = true;}
	            if (readType != null) { query += (and ? " and " : " ") + " readtype = " + readType; and = true;}
	            query += ")";
	            and = true;
	        }
	        if (genome != null) {query += (and ? " and " : " ") + " genome = " + genome.getDBID(); and = true; }
	        if (align != null) {query += (and ? " and " : " ") + " name = ? "; and = true; }
	
	        ps = cxn.prepareStatement(query);
	        int index = 1;
	        if (name != null || replicate != null) {
	            if (name != null) { ps.setString(index++,name);}
	            if (replicate != null) { ps.setString(index++,replicate);}
	        }
	        if (align != null) {ps.setString(index++,align);}
	        //System.out.println(ps.toString());
	        rs = ps.executeQuery();
	        while (rs.next()) {
	            try {
	                output.add(new SeqAlignment(rs,this));
	            } catch (NotFoundException e) {
	                throw new DatabaseException(e.toString(),e);
	            }
	        }
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (ps != null) { try {ps.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
        return filterAlignmentsByPermission(output);
    }
	
	/**
	 * Filter out alignments that the user shouldn't see
	 * @param aligns
	 * @return
	 */
	private Collection<SeqAlignment> filterAlignmentsByPermission(Collection<SeqAlignment> aligns){
		Collection<SeqAlignment> output = new ArrayList<SeqAlignment>();
		for(SeqAlignment a : aligns)
			if(a.getPermissions().contains("public") || a.getPermissions().contains(myusername))
				output.add(a);
		return output;
	}

	/**
	 * Read a set of alignment parameters
	 * @param reader
	 * @return
	 * @throws IOException
	 */
	public static Map<String,String> readParameters(BufferedReader reader) throws IOException {
		Map<String, String> params = new HashMap<String, String>();
		String line = null;
		while ((line = reader.readLine()) != null) {
			int p = line.indexOf('=');
			String key = line.substring(0, p);
			String value = line.substring(p + 1);
			params.put(key, value);
		}
		reader.close();
        return params;
    }

	/**
	 * Add a set of alignment parameters to a SeqAlignment
	 * @param align
	 * @param paramsfile
	 * @throws SQLException
	 * @throws IOException
	 */
	public void addAlignmentParameters(SeqAlignment align, File paramsfile) throws SQLException, IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(paramsfile)));
		addAlignmentParameters(align, readParameters(reader));
	}


	public void addAlignmentParameters(SeqAlignment align, Map<String, ? extends Object> params) throws SQLException {
		Connection cxn=null;
		PreparedStatement insert = null;
        try {
            cxn = DatabaseConnectionManager.getConnection(role);
            insert = cxn.prepareStatement("insert into alignmentparameters(alignment,name,value) values(?,?,?)");
			insert.setInt(1, align.getDBID());
			for (String k : params.keySet()) {
				insert.setString(2, k);
				Object val = params.get(k);
				if (val == null) {
					val = "";
				} else {
					val = val.toString();
					if (val == null) {
						val = "";
					}
				}
	
	
				insert.setString(3, (String)val);
				insert.execute();
			}
			insert.close();
		} finally {
			if (insert != null) { try { insert.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
	}


	public Map<String, String> getAlignmentParameters(SeqAlignment align) throws SQLException {
		Map<String, String> output = new HashMap<String, String>();
		Connection cxn=null;
		Statement get = null;
        ResultSet rs = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
			get = cxn.createStatement();
			rs = get.executeQuery("select name, value from alignmentparameters where alignment = " + align.getDBID());
			//System.out.println("select name, value from alignmentparameters where alignment = " + align.getDBID());
			while (rs.next()) {
				output.put(rs.getString(1), rs.getString(2));
			}
			rs.close();
			get.close();
		} finally {
			if (rs != null) { try {rs.close(); } catch (SQLException ex) {  }}
	        if (get != null) { try {get.close();} catch (SQLException ex) { } }
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}		
		return output;
	}


	/*
	 * SeqHit loading  
	 */
	
	/**
	 * Convert SingleHits to SeqHits 
	 * @param input
	 * @param align
	 * @return
	 */
	public List<SeqHit> convert(Collection<SingleHit> input, SeqAlignment align) {
        Genome g = align.getGenome();
        ArrayList<SeqHit> output = new ArrayList<SeqHit>();
        for (SingleHit s : input) {
            int start = s.pos;
            int end = s.strand ? s.pos + s.length : s.pos - s.length;
            output.add(new SeqHit(g, g.getChromName(s.chrom), Math.min(start,end), Math.max(start,end),
                                      s.strand ? '+' : '-', s.weight));
        }
        return output;
    }
	
	/**
	 * Convert PairedHits to SeqHitsPairs 
	 * @param input
	 * @param align
	 * @return
	 */
	public List<SeqHitPair> convertPairs(Collection<PairedHit> input, SeqAlignment align) {
        Genome g = align.getGenome();
        ArrayList<SeqHitPair> output = new ArrayList<SeqHitPair>();
        for (PairedHit p : input) {
            int lstart = p.leftPos;
            int lend = p.leftStrand ? p.leftPos + p.leftLength : p.leftPos - p.leftLength;
            int rstart = p.rightPos;
            int rend = p.rightStrand ? p.rightPos + p.rightLength : p.rightPos - p.rightLength;
            SeqHit left = new SeqHit(g, g.getChromName(p.leftChrom), Math.min(lstart,lend), Math.max(lstart,lend),
            		p.leftStrand ? '+' : '-', p.weight);
            SeqHit right = new SeqHit(g, g.getChromName(p.rightChrom), Math.min(rstart,rend), Math.max(rstart,rend),
                    p.rightStrand ? '+' : '-', p.weight);
            output.add(new SeqHitPair(left, right, p.weight, p.pairCode));
        }
        return output;
    }

	/**
	 * Load a chromosome's worth of single hits 
	 * @param a
	 * @param chromid
	 * @return
	 * @throws IOException
	 */
	public List<SeqHit> loadByChrom(SeqAlignment a, int chromid, boolean loadR2) throws IOException {
		List<SeqHit> data = new ArrayList<SeqHit>();
        String alignid = Integer.toString(a.getDBID());
        try {
            data.addAll(convert(client.getSingleHits(alignid, chromid,loadR2, null,null,null,null),a));
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
		return data;
	}
	
	/**
	 * Load a chromosome's worth of paired hits 
	 * @param a
	 * @param chromid
	 * @return
	 * @throws IOException
	 */
	public List<SeqHitPair> loadPairsByChrom(SeqAlignment a, int chromid) throws IOException {
		List<SeqHitPair> data = new ArrayList<SeqHitPair>();
        String alignid = Integer.toString(a.getDBID());
        try {
            data.addAll(convertPairs(client.getPairedHits(alignid, chromid,true,null,null,null,null),a));
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
		return data;
	}
	
	/**
	 * Load single hits per region
	 * @param align
	 * @param r
	 * @return
	 * @throws IOException
	 */
	public List<SeqHit> loadByRegion(SeqAlignment align, Region r, boolean loadR2) throws IOException {
        try {
            return convert(client.getSingleHits(Integer.toString(align.getDBID()),
                                                r.getGenome().getChromID(r.getChrom()),
                                                loadR2,
                                                r.getStart(),
                                                r.getEnd(),
                                                null,
                                                null), align);
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
	}
	
	/**
	 * Load pairs per region
	 * @param align
	 * @param r
	 * @return
	 * @throws IOException
	 */
	public List<SeqHitPair> loadPairsByRegion(SeqAlignment align, Region r) throws IOException {
        try {
            return convertPairs(client.getPairedHits(Integer.toString(align.getDBID()),
                                                r.getGenome().getChromID(r.getChrom()),
                                                true,
                                                r.getStart(),
                                                r.getEnd(),
                                                null,
                                                null), align);
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
	}
	
	/**
	 * Load single hits by region
	 * @param alignments
	 * @param r
	 * @return
	 * @throws IOException
	 */
	public Collection<SeqHit> loadByRegion(List<SeqAlignment> alignments, Region r, boolean loadR2) throws IOException {
		if (alignments.size() < 1) {
			throw new IllegalArgumentException("Alignment List must not be empty.");
		}
        Collection<SeqHit> output = null;
        for (SeqAlignment a : alignments) {
            if (output == null) {
                output = loadByRegion(a,r, loadR2);
            } else {
                output.addAll(loadByRegion(a,r, loadR2));
            }
        }
		return output;
	}
	
	/**
	 * Load pairs by region
	 * @param alignments
	 * @param r
	 * @return
	 * @throws IOException
	 */
	public Collection<SeqHitPair> loadPairsByRegion(List<SeqAlignment> alignments, Region r) throws IOException {
		if (alignments.size() < 1) {
			throw new IllegalArgumentException("Alignment List must not be empty.");
		}
        Collection<SeqHitPair> output = null;
        for (SeqAlignment a : alignments) {
            if (output == null) {
                output = loadPairsByRegion(a,r);
            } else {
                output.addAll(loadPairsByRegion(a,r));
            }
        }
		return output;
	}
	
	/**
	 * Load integer positions of single hits by region.
	 * If Region is a StrandedRegion, then the positions returned are only for that strand
	 * @param alignments
	 * @param r
	 * @return
	 * @throws IOException
	 * @throws ClientException
	 */
    public List<Integer> positionsByRegion(List<SeqAlignment> alignments, Region r, boolean loadR2) throws IOException, ClientException {
		if (alignments.size() < 1) {
			throw new IllegalArgumentException("Alignment List must not be empty.");
		}
        List<Integer> output = new ArrayList<Integer>();
        for (SeqAlignment a : alignments) {
            int[] pos = client.getPositions(Integer.toString(a.getDBID()),
                                            r.getGenome().getChromID(r.getChrom()),
                                            loadR2,
                                            false,
                                            r.getStart(),
                                            r.getEnd(),
                                            null,
                                            null,
                                            r instanceof StrandedRegion ? (((StrandedRegion)r).getStrand() == '+') : null);
            for (int i = 0; i < pos.length; i++) {
                output.add(pos[i]);
            }                                            
        }
        return output;
    }
    /**
	 * Load integer positions of single hits by region.
	 * If Region is a StrandedRegion, then the positions returned are only for that strand
	 * @param alignment
	 * @param r
	 * @return
	 * @throws IOException
	 * @throws ClientException
	 */
    public List<Integer> positionsByRegion(SeqAlignment alignment, Region r, boolean loadR2) throws IOException {
        List<Integer> output = new ArrayList<Integer>();
        try {
            int[] pos = client.getPositions(Integer.toString(alignment.getDBID()),
                                            r.getGenome().getChromID(r.getChrom()),
                                            loadR2,
                                            false,
                                            r.getStart(),
                                            r.getEnd(),
                                            null,
                                            null,
                                            r instanceof StrandedRegion ? (((StrandedRegion)r).getStrand() == '+') : null);
            for (int i = 0; i < pos.length; i++) {
                output.add(pos[i]);
            }                                            
            return output;
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
    }

    /**
     * Count single hits in a region.
     * If Region is a StrandedRegion, then the counts returned are only for that strand
     * @param align
     * @param r
     * @return
     * @throws IOException
     */
	public int countByRegion(SeqAlignment align, Region r, boolean loadR2) throws IOException {
        try {
            return client.getCount(Integer.toString(align.getDBID()),
                                   r.getGenome().getChromID(r.getChrom()),
                                   loadR2,
                                   false,
                                   r.getStart(),
                                   r.getEnd(),
                                   null,
                                   null,
                                   r instanceof StrandedRegion ? (((StrandedRegion)r).getStrand() == '+') : null);
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
    }
    

	/**
	 * Count single hits in a region
	 * @param alignments
	 * @param r
	 * @return
	 * @throws IOException
	 */
	public int countByRegion(List<SeqAlignment> alignments, Region r, boolean loadR2) throws IOException {
		if (alignments.size() < 1) { 
			throw new IllegalArgumentException("Alignment List must not be empty."); 
		}
        int total = 0;
        for (SeqAlignment a : alignments) {
            total += countByRegion(a,r, loadR2);
        }
        return total;
	}

	/**
	 * Weigh single hits in a region.
	 * If Region is a StrandedRegion, then the counts returned are only for that strand
	 * @param alignments
	 * @param r
	 * @return
	 * @throws IOException
	 */
	public double weightByRegion(List<SeqAlignment> alignments, Region r, boolean loadR2) throws IOException {
		if (alignments.size() < 1) { 
			throw new IllegalArgumentException("Alignment List must not be empty."); 
		}
        double total = 0;
        for (SeqAlignment a : alignments) {
            try {
                total += client.getWeight(Integer.toString(a.getDBID()),
                                          r.getGenome().getChromID(r.getChrom()),
                                          loadR2,
                                          false,
                                          r.getStart(),
                                          r.getEnd(),
                                          null,
                                          null,
                                          r instanceof StrandedRegion ? (((StrandedRegion)r).getStrand() == '+') : null);
            } catch (ClientException e) {
                throw new IllegalArgumentException(e);
            }            
        }
        return total;
	}
	

	/**
	 * Count all single hits in the alignment
	 * @param align
	 * @return
	 * @throws SQLException
	 */
	public int countAllHits(SeqAlignment align, boolean loadR2) throws IOException {
        try {
            return client.getCount(Integer.toString(align.getDBID()),
                                   loadR2,false,false,null);
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
	}

	/**
	 * Weigh all single hits in the alignment
	 * @param align
	 * @return
	 * @throws IOException
	 */
	public double weighAllHits(SeqAlignment align, boolean loadR2) throws IOException {
        try {
            return client.getWeight(Integer.toString(align.getDBID()),
                                    loadR2,false,false,null);
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
	}

	/** Generates a histogram of the total weight of reads mapped to each bin.
     * Output maps bin center to weight centered around that bin.  Each read
     * is summarized by its start position.
     */
    public Map<Integer,Float> histogramWeight(SeqAlignment align, char strand, Region r, int binWidth, boolean loadR2) throws IOException {
        try {
            return client.getWeightHistogram(Integer.toString(align.getDBID()),
                                             r.getGenome().getChromID(r.getChrom()),
                                             loadR2,
                                             false,
                                             0,
                                             binWidth,
                                             r.getStart(),
                                             r.getEnd(),
                                             null,
                                             strand == '+');
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
    }
    /** Generates a histogram of the total number of reads mapped to each bin.
     * Output maps bin center to weight centered around that bin.  Each read
     * is summarized by its start position.
     */
    public Map<Integer,Integer> histogramCount(SeqAlignment align, char strand, Region r, int binWidth, boolean loadR2) throws IOException {        
        try {
            return client.getHistogram(Integer.toString(align.getDBID()),
                                       r.getGenome().getChromID(r.getChrom()),
                                       loadR2,
                                       false,
                                       0,
                                       binWidth,
                                       r.getStart(),
                                       r.getEnd(),
                                       null,
                                       strand == '+');
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
    }
    /**
     * Get the total # of hits and weight for an alignment but only include reads
     * on the specified strand.  
     */
    public Pair<Long,Double> getAlignmentStrandedCountWeight(SeqAlignment align, char strand, boolean loadR2) throws IOException {
        try {
            long count = client.getCount(Integer.toString(align.getDBID()), false, false, strand=='+', loadR2);
            double weight = client.getWeight(Integer.toString(align.getDBID()), false, false, strand=='+', loadR2);
            Pair<Long,Double> output = new Pair<Long,Double>(count,weight);
            return output;
        } catch (ClientException e) {
            throw new IllegalArgumentException(e);
        }
    }
    
	public void close() {
        if (!closed && client != null) {
            client.close();
            client = null;
        }
        closed=true;
	}
	public boolean isClosed() {
		return closed;
	}

}