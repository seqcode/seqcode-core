package edu.psu.compbio.seqcode.projects.shaun;

import java.io.IOException;
import java.text.ParseException;
import java.util.Random;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Organism;
import edu.psu.compbio.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.BackgroundModelIO;

public class RandomSequenceGenerator {

	private MarkovBackgroundModel markov;
	private int modelLen;
	private Random rand = new Random();
	
	
	public static void main(String[] args) throws IOException, ParseException, NotFoundException {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("len") || !ap.hasKey("num")||!ap.hasKey("markov") || !ap.hasKey("species")) { 
            System.err.println("Usage:\n" +
                               "RandomSequenceGenerator " +
                               "--len <length of sequences> " +
                               "--num <number of sequences> " +
                               "--markov <markov model file> " +
                               "--species <organism;genome> ");
            return;
        }
        int length = new Integer(ap.getKeyValue("len")).intValue();
        int number = new Integer(ap.getKeyValue("num")).intValue();
        String mFile = ap.getKeyValue("markov");
        Pair<Organism, Genome> pair = Args.parseGenome(args);        
        
        MarkovBackgroundModel model = BackgroundModelIO.parseMarkovBackgroundModel(mFile, pair.cdr());
        
        RandomSequenceGenerator gen = new RandomSequenceGenerator(model);
        for(int i=1; i<=number; i++){
        	String s = gen.execute(length);
        	System.out.println(">Sequence"+i);
        	System.out.println(s);
        }
	}
	
	public RandomSequenceGenerator(MarkovBackgroundModel m){
		markov=m;
		modelLen = markov.getMaxKmerLen();
	}
	
	public String execute(int len){
		String seq = new String();
		
		//Preliminary bases
		for(int i=1; i<modelLen && i<=len; i++){
			double prob = rand.nextDouble();
			double sum=0; int j=0;
			while(sum<prob){
				String test = seq.concat(int2base(j));
				sum += markov.getMarkovProb(test);
				if(sum>=prob){
					seq = test;
					break;
				}
				j++;
			}
		}
		//Remaining bases
		for(int i=modelLen; i<=len; i++){
			String lmer = seq.substring(seq.length()-(modelLen-1));
			double prob = rand.nextDouble();
			double sum=0; int j=0;
			while(sum<prob){
				String test = lmer.concat(int2base(j));
				sum += markov.getMarkovProb(test);
				if(sum>=prob){
					seq =seq.concat(int2base(j));
					break;
				}
				j++;
			}
		}
		
		return seq;
	}
	
	protected String int2base(int x){
		String c;
		switch(x){
			case 0:
				c="A"; break;
			case 1:
				c="C"; break;
			case 2:
				c="G"; break;
			case 3:
				c="T"; break;
			default:
				c="N";
		}
		return(c);
	}
}
