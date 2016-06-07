package org.seqcode.genome.sequence;

import java.util.Iterator;
import java.util.regex.*;
import java.io.*;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gsebricks.iterators.SingleIterator;
import org.seqcode.gsebricks.verbs.ExpanderIterator;
import org.seqcode.gsebricks.verbs.FileLineExpander;
import org.seqcode.gsebricks.verbs.Sink;
import org.seqcode.gseutils.Closeable;
import org.seqcode.gseutils.NotFoundException;


/**
 * Writes the sequence contained in the input Regions
 * to a FASTA file
 */
public class FastaWriter<X extends Region> implements Sink<X>, Closeable {
    
    public static void main(String[] args) { 
        FastaWriter writer = new FastaWriter(System.out);
		Pattern p = Pattern.compile("(.*):([\\d]+)-([\\d]+)");
        try {
			Genome genome = Genome.findGenome(args[0]);
			SingleIterator<File> fitr = new SingleIterator<File>(new File(args[1]));
			Iterator<String> lines = new ExpanderIterator<File,String>(new FileLineExpander(), fitr);
			while(lines.hasNext()) { 
				String line = lines.next().trim();
				Matcher m = p.matcher(line);
				if(m.matches()) { 
					String chrom = m.group(1);
					if(chrom.startsWith("chr")) { chrom = chrom.substring(3, chrom.length()); }
					int start = Integer.parseInt(m.group(2));
					int end = Integer.parseInt(m.group(3));
					Region r = new Region(genome, chrom, start, end);
					writer.consume(r);
				} else { 
					System.err.println("Can't match \"" + line + "\"");
				}
			}
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }
    
    private PrintStream ps;
    private FileOutputStream os;
    private SequenceGenerator seqgen;
    private int lineLength;
    
    public FastaWriter(File f) throws IOException { 
        os = new FileOutputStream(f);
        ps = new PrintStream(os);
        seqgen = new SequenceGenerator();
        lineLength = 100;        
    }

    public FastaWriter (String fname) throws IOException, FileNotFoundException {
        os = new FileOutputStream(fname);
        ps = new PrintStream(os);
        seqgen = new SequenceGenerator();
        lineLength = 100;
    }

    public FastaWriter(PrintStream ps) {
        seqgen = new SequenceGenerator();
        this.ps = ps;
        os = null;
        lineLength = 100;
    }
    
    public void setLineLength(int ll) { lineLength = ll; }

	public void init() {}
    public void useCache(boolean b) {
        seqgen.useCache(b);
    }

    public void consume(Iterator<X> iter) {
        while (iter.hasNext()) {
            consume(iter.next());
        }
    }

	public void finish() { close(); }

    public void consume(X r) {
        Genome g = r.getGenome();
        String s = seqgen.execute(r);
        if (r instanceof StrandedRegion) {
            if (((StrandedRegion)r).getStrand() == '-') {
                s = org.seqcode.genome.sequence.SequenceUtils.reverseComplement(s);
            }
        }
        int start = 0;
        int len = s.length();
        ps.println(">" + r.toString());
        while (start < len) {
            int end = start + lineLength;
            if (end > len) {
                end = len;
            }
            ps.println(s.substring(start,end));
            start+= lineLength;
        }
    }

    public void close() {
        try {
            ps.close();
            if (os != null) {
                os.close();                
            }
            ps.close();
            seqgen = null;
            os = null;
            ps = null;
        } catch (IOException ex) {
            throw new RuntimeException(ex.toString(),ex);
        }
    }

    public boolean isClosed() {return ps == null;}


}
