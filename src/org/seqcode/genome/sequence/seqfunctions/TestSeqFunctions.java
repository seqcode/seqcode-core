package org.seqcode.genome.sequence.seqfunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.NamedStrandedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.Pair;

public class TestSeqFunctions {

	List<Region> regions;
	
	public TestSeqFunctions(List<Region> regs){
		regions = regs;
		
	}
	
	public void execute(){
		
	}
	
	
	public static void main(String[] args) {
		try {
        	String infile = Args.parseString(args, "reg", "");
            Pair<Species,Genome> pair = Args.parseGenome(args);
            Genome genome = pair.cdr();
            List<Region> regs = new ArrayList<Region>();
            
            Region r=null;
            String line;
            BufferedReader reader;
            if (!infile.equals("")) {
            	reader = new BufferedReader(new FileReader(infile));
            } else {
            	reader = new BufferedReader(new InputStreamReader(System.in));
            }
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                r = StrandedRegion.fromString(genome, line);
                if (r != null) {
                    String pieces[] = line.split("\\t");
                    if (pieces.length > 1) {
                        r = new NamedStrandedRegion(r, pieces[1], r.getStrand());
                    }
                } else {
                    r = Region.fromString(genome,line);                       
                }
            }

            if(r!=null)
            	regs.add(r);

            TestSeqFunctions tester = new TestSeqFunctions(regs);
            tester.execute();
            
        } catch (Exception e) {
            e.printStackTrace();
        }

	}

}
