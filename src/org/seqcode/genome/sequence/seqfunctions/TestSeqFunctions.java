package org.seqcode.genome.sequence.seqfunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.NamedStrandedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.gseutils.Args;

public class TestSeqFunctions {
	
	Genome genome;
	List<Region> regions;
	SequenceGenerator seqgen;
	
	public TestSeqFunctions(Genome gen, List<Region> regs){
		genome = gen;
		regions = regs;
		seqgen = new SequenceGenerator(gen);
	}
	
	public void execute(){
		SeqFunction fn;
		
		try {
			for(Region r : regions){
				String seq = seqgen.execute(r);
				
				System.out.println("Region:\t"+r.getLocationString()+"\t"+seq);
				
				//BaseFrequency
				fn = new BaseFrequencyFunction();
				printFunction(fn.score(seq), fn.dimensionLabels());
				
				//GC
				fn = new GCContentFunction(3);
				printFunction(fn.score(seq), fn.dimensionLabels());
				
				//polyA
				fn = new PolyAFunction(3);
				printFunction(fn.score(seq), fn.dimensionLabels());
				
				//MGW
				fn = new MGWStructureFunction();
				printFunction(fn.score(seq), fn.dimensionLabels());
				
				//PropTwist
				fn = new PropTwistStructureFunction();
				printFunction(fn.score(seq), fn.dimensionLabels());
			
				//Roll
				fn = new RollStructureFunction();
				printFunction(fn.score(seq), fn.dimensionLabels());
				
				//HelixTwist
				fn = new HelixTwistStructureFunction();
				printFunction(fn.score(seq), fn.dimensionLabels());
			
			}
		} catch (SeqFunctionException e) {
			e.printStackTrace();
		}
	}
	
	public void printFunction(double[][] score, String[] labels){
		for(int x=0; x<score.length; x++){
			System.out.print(labels[x]);
			for(int y=0; y<score[x].length; y++){
				System.out.print("\t"+score[x][y]);
			}
			System.out.println("");
		}
	}
	
	public static void main(String[] args) {
		try {
			GenomeConfig gcon = new GenomeConfig(args);
        	String infile = Args.parseString(args, "reg", "");
            Genome genome = gcon.getGenome();
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

            TestSeqFunctions tester = new TestSeqFunctions(genome, regs);
            tester.execute();
            
        } catch (Exception e) {
            e.printStackTrace();
        }

	}

}
