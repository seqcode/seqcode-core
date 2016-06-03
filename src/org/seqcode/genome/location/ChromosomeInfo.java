package org.seqcode.genome.location;

/**
 * ChromosomeInfo: basic information on a Chromosome (name, length, database ID)
 * 
 * @author mahony
 *
 */
public class ChromosomeInfo {
	private int dbid;
    private int length;
    private String name;
    
    public ChromosomeInfo(int id, int len, String n) { 
        length = len;
        dbid = id;
        name = n;
    }
    
    public String getName() { return name; }
    public int getDBID() { return dbid; }
    public int getLength() { return length; }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ChromosomeInfo)) { return false; }
        ChromosomeInfo i = (ChromosomeInfo)o;
        if(dbid != i.dbid) { return false; }
        if(!name.equals(i.name)) { return false; }
        if(length != i.length) { return false; }
        return true;
    }
    
    public String toString() { return name  + " (" + length + " bp)"; }
    
    public int hashCode() { 
        int code = 17;
        code += dbid; code *= 37;
        code += name.hashCode(); code *= 37;
        code += length; code *= 37;
        return code;
    }
}
