package edu.psu.compbio.seqcode.gse.utils.io.parsing;

import java.io.BufferedReader;
import java.io.*;
import java.util.*;
import java.net.URLEncoder;

public class ParseGFF implements Iterator { 
    
    public static void main(String[] args) { 
        try { 
            ParseGFF parser = new ParseGFF(new File(args[0]));
            parser.printLines(System.out);
        } catch(IOException ie) { 
            ie.printStackTrace();
        }
    }
    
    private BufferedReader br;    
    private int lineNum;
    private String line, filename;
    private boolean dirty;

    public ParseGFF(File f) throws IOException { 
        br = new BufferedReader(new FileReader(f));
        lineNum = 0;
        dirty = true;
        filename = f.getName();
    }

    public boolean hasNext() {
        if (dirty) {
            try {
                lineNum++;
                line = br.readLine();
            } catch (IOException ex) {
                throw new RuntimeException("Parsing Error, File \"" + filename + "\", line " + lineNum);
            }
            dirty = false;
        }
        if (line == null) {
            try {
                br.close();
            } catch (IOException ex) {
                throw new RuntimeException("Can't close " + filename);
            }
            return false;
        } else {
            if (line.startsWith("#")) {
                dirty = true;
                return hasNext();
            } else {
                return true;
            }
        }
    }
    public GFFEntry next() throws NoSuchElementException {
        if (line == null) {
            throw new NoSuchElementException("No more lines to parse");
        }
        dirty = true;
        line = line.trim();
        if(!line.startsWith("#")) { 
            try { 
            	GFFEntry gffLine = new GFFEntry(line);
                return gffLine;
            } catch(NoSuchElementException e) { 
                throw new RuntimeException("Parsing Error, File \"" + filename + "\", line " + lineNum);
            }
        } else {
            return next();
        }
    }
    public void remove() throws UnsupportedOperationException {
        throw new UnsupportedOperationException("Can't remove lines from GFF file");
    }
    
    public void printLines() { printLines(System.out); }
    public void printLines(PrintStream ps) { 
        while(hasNext()) { 
            next().printLine(ps);
        }
    }
    
}
