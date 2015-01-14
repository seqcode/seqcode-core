/*
 * Created as ChipSeqLoader on May 16, 2007
 */
package edu.psu.compbio.seqcode.gse.datasets.seqdata;

import java.util.*;
import java.sql.*;
import java.io.*;


import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Organism;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.core.AlignType;
import edu.psu.compbio.seqcode.gse.datasets.core.CellLine;
import edu.psu.compbio.seqcode.gse.datasets.core.ExptCondition;
import edu.psu.compbio.seqcode.gse.datasets.core.ExptTarget;
import edu.psu.compbio.seqcode.gse.datasets.core.ExptType;
import edu.psu.compbio.seqcode.gse.datasets.core.Lab;
import edu.psu.compbio.seqcode.gse.datasets.core.MetadataLoader;
import edu.psu.compbio.seqcode.gse.datasets.core.ReadType;
import edu.psu.compbio.seqcode.gse.datasets.core.SeqDataUser;
import edu.psu.compbio.seqcode.gse.projects.readdb.Client;
import edu.psu.compbio.seqcode.gse.projects.readdb.ClientException;
import edu.psu.compbio.seqcode.gse.projects.readdb.PairedHit;
import edu.psu.compbio.seqcode.gse.projects.readdb.SingleHit;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseConnectionManager;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseException;

/**
 * @author tdanford
 * @author mahony
 * 
 * SeqDataLoader serves as a clearinghouse for query interactions with the seqdata database and
 * associated metadata from tables in the core database (via the MetadataLoader). 
 * 
 * Implements a simple access control by checking if the username is in the permissions field
 * in each SeqAlignment. Note that while this is access-control-lite for experiment metadata, 
 * access to the data stored in the underlying readdb entries has more robust access-control.
 */
public class SeqDataLoader implements edu.psu.compbio.seqcode.gse.utils.Closeable {

	public static String role = "seqdata";


	public static void main(String[] args) throws Exception{
		try {
			SeqDataLoader loader = new SeqDataLoader();
			Collection<SeqExpt> expts = loader.loadAllExperiments();
			for (SeqExpt expt : expts) {
				Collection<SeqAlignment> aligns = loader.loadAllAlignments(expt);
				for (SeqAlignment align : aligns) {
					System.out.println(expt.getDBID() + "\t" + expt.getName() + ";"+ expt.getReplicate()+"\t"+align.getName()+"\t"+align.getDBID()+"\t"+align.getGenome());
				}				
			}
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
    
    //Experiment descriptors
    private Collection<ExptType> exptTypes = null;
    private Collection<Lab> labs=null;
    private Collection<ExptTarget> exptTargets = null;
    private Collection<ExptCondition> exptConditions = null;
    private Collection<CellLine> celllines = null;
    private Collection<ReadType> readTypes = null;
    private Collection<AlignType> alignTypes = null;
    private Collection<SeqDataUser> seqDataUsers=null;
    public Collection<ExptType> getExptTypes() throws SQLException {if(exptTypes==null){exptTypes=mloader.loadAllExptTypes();} return exptTypes;}
    public Collection<Lab> getLabs() throws SQLException {if(labs==null){labs=mloader.loadAllLabs();} return labs;}
    public Collection<ExptTarget> getExptTargets() throws SQLException {if(exptTargets==null){exptTargets = mloader.loadAllExptTargets();} return exptTargets;}
    public Collection<ExptCondition> getExptConditions() throws SQLException {if(exptConditions==null){exptConditions = mloader.loadAllExptConditions();} return exptConditions;} 
    public Collection<CellLine> getCellLines() throws SQLException {if(celllines==null){celllines = mloader.loadAllCellLines();} return celllines;}
	public Collection<ReadType> getReadTypes() throws SQLException {if(readTypes==null){readTypes = mloader.loadAllReadTypes();} return readTypes;}
	public Collection<AlignType> getAlignTypes() throws SQLException {if(alignTypes==null){alignTypes = mloader.loadAllAlignTypes();} return alignTypes;}
	public SeqDataUser findSeqDataUser(String name) throws SQLException {return mloader.findSeqDataUser(name);}
	public Collection<SeqDataUser> getSeqDataUsers() throws SQLException {if(seqDataUsers==null){seqDataUsers = mloader.loadAllSeqDataUsers();} return seqDataUsers;}
    public SeqDataUser getMyUser(){return myUser;}    
	
    public SeqDataLoader() throws SQLException, IOException{this(true);}
	public SeqDataLoader(boolean openClient) throws SQLException, IOException {
		mloader = new MetadataLoader();
		if(openClient){
	        try {
	            client = new Client();
	            
	        } catch (ClientException e) {
	            throw new IllegalArgumentException(e);
	        }
		}
		Connection cxn = null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            myusername = DatabaseConnectionManager.getUsername(role);
            myUser = findSeqDataUser(myusername);
        } catch (SQLException e) {
            throw new DatabaseException(e.toString(),e);
        } finally {
        	if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
        }
		closed=false;
	}
    public Client getClient(){return client;}
    
	public Collection<Genome> loadExperimentGenomes(SeqExpt expt) throws SQLException {
		Connection cxn = null;
		LinkedList<Genome> genomes = new LinkedList<Genome>();
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            String query = String.format("select genome from seqalignment where expt=%d", expt.getDBID());
			Statement s = cxn.createStatement();
			ResultSet rs = s.executeQuery(query);
			while (rs.next()) {
				int gid = rs.getInt(1);
				try {
					Genome g = Organism.findGenome(gid);
					genomes.add(g);
				}
				catch (NotFoundException e) {
					e.printStackTrace();
				}
			}
			rs.close();
			s.close();
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return genomes;
	}


	public Collection<SeqExpt> loadAllExperiments() throws SQLException {
		LinkedList<SeqExpt> expts = new LinkedList<SeqExpt>();
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
	        labs=mloader.loadAllLabs();
			exptTypes=mloader.loadAllExptTypes();
			exptTargets = mloader.loadAllExptTargets();
			exptConditions = mloader.loadAllExptConditions();
	        celllines = mloader.loadAllCellLines();
	        readTypes = mloader.loadAllReadTypes();
	        alignTypes = mloader.loadAllAlignTypes();
	        seqDataUsers = mloader.loadAllSeqDataUsers();
	        PreparedStatement ps = SeqExpt.createLoadAll(cxn);
	        ps.setFetchSize(1000);
			ResultSet rs = ps.executeQuery();
			while (rs.next()) {
				expts.addLast(new SeqExpt(rs, this));
			}
			Collections.sort(expts);
			rs.close();
			ps.close();
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return expts;
	}


	public SeqExpt findExperiment(String name, String rep) throws SQLException {
		SeqExpt expt = null;
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
		    PreparedStatement ps = SeqExpt.createLoadByNameReplicate(cxn);
			ps.setString(1, name);
			ps.setString(2, rep);
			ResultSet rs = ps.executeQuery();
			if (rs.next()) {
				expt = new SeqExpt(rs, this);
			}
			rs.close();
			ps.close();
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return expt;
	}
	
	public SeqExpt loadExperiment(String name, String rep) throws NotFoundException, SQLException {
		SeqExpt expt = null;
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
			PreparedStatement ps = SeqExpt.createLoadByNameReplicate(cxn);
			ps.setString(1, name);
			ps.setString(2, rep);
			ResultSet rs = ps.executeQuery();
			if (rs.next()) {
				expt = new SeqExpt(rs, this);
			}
			rs.close();
			ps.close();
	
			if (expt == null) { throw new NotFoundException(name+";"+rep); }
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}	
		return expt;
	}


	public Collection<SeqExpt> loadExperiments(String name) throws SQLException {
		SeqExpt expt = null;
		LinkedList<SeqExpt> expts = new LinkedList<SeqExpt>();
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
			PreparedStatement ps = SeqExpt.createLoadByName(cxn);
			ps.setString(1, name);
			ResultSet rs = ps.executeQuery();
			while (rs.next()) {
				expt = new SeqExpt(rs, this);
				expts.add(expt);
			}
			rs.close();
			ps.close();
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return expts;
	}
	
	public Collection<SeqExpt> loadExperiments(Lab lab) throws SQLException {
		LinkedList<SeqExpt> expts = new LinkedList<SeqExpt>();
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement ps = SeqExpt.createLoadByLab(cxn);
			ps.setInt(1, lab.getDBID());
			ResultSet rs = ps.executeQuery();
			SeqExpt expt = null;
			while (rs.next()) {
				expt = new SeqExpt(rs, this);
				expts.add(expt);
			}
			rs.close();
			ps.close();
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return expts;
	}
	
	public Collection<SeqExpt> loadExperiments(ExptCondition cond) throws SQLException {
		LinkedList<SeqExpt> expts = new LinkedList<SeqExpt>();
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement ps = SeqExpt.createLoadByCondition(cxn);
			ps.setInt(1, cond.getDBID());
			ResultSet rs = ps.executeQuery();
			SeqExpt expt = null;
			while (rs.next()) {
				expt = new SeqExpt(rs, this);
				expts.add(expt);
			}
			rs.close();
			ps.close();
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}	
		return expts;
	}
	
	public Collection<SeqExpt> loadExperiments(ExptTarget target) throws SQLException {
		LinkedList<SeqExpt> expts = new LinkedList<SeqExpt>();
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
	        PreparedStatement ps = SeqExpt.createLoadByTarget(cxn);
			ps.setInt(1, target.getDBID());
			ResultSet rs = ps.executeQuery();
			SeqExpt expt = null;
			while (rs.next()) {
				expt = new SeqExpt(rs, this);
				expts.add(expt);
			}
			rs.close();
			ps.close();
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return expts;
	}

	public Collection<SeqExpt> loadExperiments(CellLine cell) throws SQLException {
		LinkedList<SeqExpt> expts = new LinkedList<SeqExpt>();
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
		    PreparedStatement ps = SeqExpt.createLoadByCellline(cxn);
			ps.setInt(1, cell.getDBID());
			ResultSet rs = ps.executeQuery();
			SeqExpt expt = null;
			while (rs.next()) {
				expt = new SeqExpt(rs, this);
				expts.add(expt);
			}
			rs.close();
			ps.close();
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}	
		return expts;
	}

	public SeqExpt loadExperiment(int dbid) throws NotFoundException, SQLException {
		SeqExpt expt = null;
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
			PreparedStatement ps = SeqExpt.createLoadByDBID(cxn);
			ps.setInt(1, dbid);
			ResultSet rs = ps.executeQuery();
			if (rs.next()) {
				expt = new SeqExpt(rs, this);
			}
			rs.close();
			ps.close();
	
			if (expt == null) {
				String err = String.format("No such SeqExpt experiment with ID= %d", dbid);
				throw new NotFoundException(err);
			}
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return expt;
	}

    public Collection<SeqAlignment> loadAlignments (Genome g) throws SQLException {
    	Collection<SeqAlignment> aligns = new LinkedList<SeqAlignment>();
    	Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
	    	Collection<SeqExpt> allexpts = loadAllExperiments();
	        Map<Integer,SeqExpt> exptmap = new HashMap<Integer,SeqExpt>();
	        for (SeqExpt e : allexpts) {
	            exptmap.put(e.getDBID(), e);
	        }
	
			PreparedStatement ps = SeqAlignment.createLoadAllByGenomeStatement(cxn);
	        ps.setFetchSize(1000);
			ps.setInt(1, g.getDBID());
	        ResultSet rs = ps.executeQuery();
			while (rs.next()) {
	            SeqAlignment align = new SeqAlignment(rs, exptmap.get(rs.getInt(2)), mloader.loadAlignType(rs.getInt(6)));
	            System.out.println(align.getName());
	            aligns.add(align);
			}
			rs.close();
			ps.close();
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return filterAlignmentsByPermission(aligns);
    }

	public Collection<SeqAlignment> loadAllAlignments(SeqExpt expt) throws SQLException {
		Collection<SeqAlignment> aligns = new LinkedList<SeqAlignment>();
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
			PreparedStatement ps = SeqAlignment.createLoadAllByExptStatement(cxn);
			ps.setInt(1, expt.getDBID());
	
			ResultSet rs = ps.executeQuery();
			while (rs.next()) {
				SeqAlignment align = new SeqAlignment(rs, expt, mloader.loadAlignType(rs.getInt(6)));
				aligns.add(align);
			}
			rs.close();
			ps.close();
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return filterAlignmentsByPermission(aligns);
	}
	

	public SeqAlignment loadAlignment(SeqExpt expt, String n, Genome g) throws SQLException {
		SeqAlignment align = null;
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
			PreparedStatement ps = SeqAlignment.createLoadByNameAndExptStatement(cxn);
			ps.setString(1, n);
			ps.setInt(2, expt.getDBID());
	
			ResultSet rs = ps.executeQuery();        
			while (align == null && rs.next()) {
				align = new SeqAlignment(rs, expt, mloader.loadAlignType(rs.getInt(6)));
	            if (!align.getGenome().equals(g))
	                align = null;
	            if(align!=null && !align.getPermissions().contains("public") && !align.getPermissions().contains(myusername))
					align=null;
			}
			rs.close();
			ps.close();
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		if (align == null) {
        	//Don't throw exception because sometimes we have replicates that don't match all alignment names
            //throw new NotFoundException("Couldn't find alignment " + n + " for " + expt + " in genome " + g);
        }
        return align;

	}
	public SeqAlignment loadAlignment(int dbid) throws NotFoundException, SQLException {
		SeqAlignment align = null;
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
			PreparedStatement ps = SeqAlignment.createLoadByIDStatement(cxn);
			ps.setInt(1, dbid);
	
			ResultSet rs = ps.executeQuery();
			if (rs.next()) {
				align = new SeqAlignment(rs, this);
				if(!align.getPermissions().contains("public") && !align.getPermissions().contains(myusername))
					align=null;
			}
			else {
				throw new NotFoundException("Couldn't find alignment with ID = " + dbid);
			}
			rs.close();
			ps.close();
		} finally {
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
		return align;
	}


	public Collection<SeqAlignment> loadAlignments(SeqLocator locator, Genome genome) throws SQLException, NotFoundException {
		List<SeqAlignment> output = new ArrayList<SeqAlignment>();
        if (locator.getReplicates().size() == 0) {
            for (SeqExpt expt : loadExperiments(locator.getExptName())) { 
                SeqAlignment align = loadAlignment(expt, locator.getAlignName(), genome);
                if (align != null) {
                    output.add(align);
                }
            }
        } else {
            for (String rep : locator.getReplicates()) {
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
		return filterAlignmentsByPermission(output);
	}
                                                             
	public Collection<SeqAlignment> loadAlignments(String name, String replicate, String align,
                                                       Integer exptType, Integer lab, Integer condition,
                                                       Integer target, Integer cells, Integer readType,
                                                       Genome genome) throws SQLException {
        
		Collection<SeqAlignment> output = new ArrayList<SeqAlignment>();
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
			String query = "select id, expt, name, genome from seqalignment";
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
	
	        PreparedStatement ps = cxn.prepareStatement(query);
	        int index = 1;
	        if (name != null || replicate != null) {
	            if (name != null) { ps.setString(index++,name);}
	            if (replicate != null) { ps.setString(index++,replicate);}
	        }
	        if (align != null) {ps.setString(index++,align);}
	        
	        ResultSet rs = ps.executeQuery();
	        while (rs.next()) {
	            try {
	                output.add(new SeqAlignment(rs,this));
	            } catch (NotFoundException e) {
	                throw new DatabaseException(e.toString(),e);
	            }
	        }
	        rs.close();
	        ps.close();
		} finally {
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

	public void addAlignmentParameters(SeqAlignment align, File paramsfile) throws SQLException, IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(paramsfile)));
		addAlignmentParameters(align, readParameters(reader));
	}


	public void addAlignmentParameters(SeqAlignment align, Map<String, ? extends Object> params) throws SQLException {
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
            PreparedStatement insert = cxn.prepareStatement("insert into alignmentparameters(alignment,name,value) values(?,?,?)");
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
			if(cxn!=null) try {cxn.close();}catch (Exception ex) {throw new DatabaseException("Couldn't close connection with role "+role, ex); }
		}
	}


	public Map<String, String> getAlignmentParameters(SeqAlignment align) throws SQLException {
		Map<String, String> output = new HashMap<String, String>();
		Connection cxn=null;
		try {
            cxn = DatabaseConnectionManager.getConnection(role);
			Statement get = cxn.createStatement();
			ResultSet rs = get.executeQuery("select name, value from alignmentparameters where alignment = " + align.getDBID());
			while (rs.next()) {
				output.put(rs.getString(1), rs.getString(2));
			}
			rs.close();
			get.close();
		} finally {
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
        if (client != null) {
            client.close();
            client = null;
        }
        closed=true;
	}
	public boolean isClosed() {
		return closed;
	}

}