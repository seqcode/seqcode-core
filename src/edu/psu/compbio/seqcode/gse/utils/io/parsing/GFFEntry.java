package edu.psu.compbio.seqcode.gse.utils.io.parsing;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;


public class GFFEntry {
    private String fChr, fSource, fFeature;
    private int fStart, fEnd;
    private double fScore;
    private char fStrand;
    private String fFrame;
    private String fAttribString;
    private Set<String> fSimpleAttribs;
    private Map<String,String> fComplexAttribs;
    
    public GFFEntry(String line) { 
        StringTokenizer st = new StringTokenizer(line, "\t");
        fChr = st.nextToken();
        if(fChr.startsWith("chr"))
        	fChr = fChr.replaceFirst("chr", "");
        fSource = st.nextToken();
        fFeature = st.nextToken();
        fStart = Integer.parseInt(st.nextToken());
        fEnd = Integer.parseInt(st.nextToken());
        try { 
            fScore = Double.parseDouble(st.nextToken());
        } catch(NumberFormatException nfe) { 
            fScore = 1.0;
        }
        
        fStrand = st.nextToken().charAt(0);
        fFrame = st.nextToken();
        String attribs = st.nextToken();
        fAttribString = attribs;
        st = new StringTokenizer(attribs, ";");
        fSimpleAttribs = new HashSet<String>();
        fComplexAttribs = new HashMap<String,String>();
        
        while(st.hasMoreTokens()) { 
            String aTok = st.nextToken();
            if(aTok.indexOf("=") != -1) { 
                String key = aTok.substring(0, aTok.indexOf("="));
                String value = aTok.substring(aTok.indexOf("=") + 1, aTok.length());
                fComplexAttribs.put(key, value);
			} else {
                fSimpleAttribs.add(aTok);
            }
        }
    }
    
    public void printLine() { printLine(System.out); }
    public void printLine(PrintStream ps) { 
        ps.println("[" + fFeature + "] (" + fChr + "," + fSource + ")");
        ps.println(String.valueOf(fScore) + ": " + fStart + " - " + fEnd);
        Iterator itr = fSimpleAttribs.iterator();
        while(itr.hasNext()) { 
            ps.println("\t" + itr.next());
        }
        
        itr = fComplexAttribs.entrySet().iterator();
        while(itr.hasNext()) { 
            Map.Entry e = (Map.Entry)itr.next();
            ps.println("\t" + e.getKey() + " ==> " + e.getValue());
        }
    }
    
    public boolean hasComplexAttrib(String key) { return fComplexAttribs.containsKey(key); }
    public String getComplexAttrib(String key) { 
        return (String)fComplexAttribs.get(key);
    }
    
    public Set<String> getPrefixedSimpleAttribs(String pref) { 
        Set<String> hashSet = new HashSet<String>();
        Iterator itr = fSimpleAttribs.iterator();
        while(itr.hasNext()) { 
            String attrib = (String)itr.next();
            if(attrib.startsWith(pref)) { hashSet.add(attrib); }
        }
        
        return hashSet;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof GFFEntry)) { return false; }
        GFFEntry l = (GFFEntry)o;
        if(!fChr.equals(l.fChr)) { return false; }
        if(!fSource.equals(l.fSource)) { return false; }
        if(!fFeature.equals(l.fFeature)) { return false; }
        if(fStart != l.fStart || fEnd != l.fEnd) { return false; }
        if(fStrand != l.fStrand) { return false; }
        if(!fFrame.equals(l.fFrame)) { return false; }
        if(fScore != l.fScore) { return false; }
        if(fSimpleAttribs.size() != l.fSimpleAttribs.size()) { return false; }
        if(fComplexAttribs.keySet().size() != l.fComplexAttribs.keySet().size()) { return false; }
        Iterator itr = fSimpleAttribs.iterator();
        while(itr.hasNext()) { 
            String key = (String)itr.next();
            if(!l.fSimpleAttribs.contains(key)) { return false; }
        }
        
        itr = fComplexAttribs.keySet().iterator();
        while(itr.hasNext()) { 
            String key = (String)itr.next();
            if(!l.fComplexAttribs.containsKey(key)) { return false; }
            String v1 = (String)fComplexAttribs.get(key);
            String v2 = (String)l.fComplexAttribs.get(key);
            if(!v1.equals(v2)) { return false; }
        }
        return true;
    }
    
    public String getChr() { return fChr; }
    public String getSource() { return fSource; }
    public String getFeature() { return fFeature; }
    public int getStart() { return fStart; }
    public int getEnd() { return fEnd; }
    public double getScore() { return fScore; }
    public char getStrand() { return fStrand; }
    public String getFrame() { return fFrame; }

    public int getMidPoint(){
    	int mid = (fStart+fEnd)/2;
    	if(fStrand == '-' && (fEnd-fStart)%2==1)
    		mid += 1;
    	return mid;
    }
    public String getAttribString() { 
		return fAttribString;
	}
    
    public boolean hasAttrib(String v) { 
        return fSimpleAttribs.contains(v) || fComplexAttribs.containsKey(v); 
    }
    
    public String getAttrib(String v) { 	
        if(!fComplexAttribs.containsKey(v)) { return null; }
        return (String)fComplexAttribs.get(v);
    }

}
