package edu.psu.compbio.seqcode.projects.shaun.teqseq.exptprops;

/**
 * SeqBiasModel: Estimate the sequencing bias in an experiment, as caused by hexamer priming or other effects. 
 * Classes of this type count the observed and expected occurrences of fixed-length k-mers at the start of reads. 
 * The way in which these counts are performed is implementation-specific.
 * The class produces a reweighting scheme based on the expected/observed or similar.   
 *  
 * @author	Shaun Mahony
 * @version	%I%, %G%
 * 
 */
import htsjdk.samtools.SAMRecord;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Organism;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.GenomeLoader;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.GenomeSegmenter;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.ReadLoader;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.SAMReadLoader;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.AGene;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.GTFAnnotationLoader;

public abstract class SeqBiasModel {

	protected GenomeLoader genLoader=null;
	protected Genome gen=null;
	protected int kmerLength=6;
	protected int numKmers;
	protected double[] kmerObsCounts=null;
	protected double[] kmerExpCounts=null;
	protected double[] kmerObsProbs=null;
	protected double[] kmerExpProbs=null;
	protected double[] kmerWeights=null;
	protected double contributingHits=0.0;
	protected byte[] kmer=null; //general purpose kmer
	private final double LOG2 = Math.log(2);
	protected String fileExt="";
	
	public SeqBiasModel(GenomeLoader g){this(g, 6);}
	/**
	 * Constructor
	 * @param k	K-mer length
	 * @param g Genome
	 */
	public SeqBiasModel(GenomeLoader g, int k){
		genLoader = g;
		gen=genLoader.getGenome();
		kmerLength = k;
		numKmers = (int)Math.pow(4, k);
		kmerObsCounts = new double[numKmers+1];
		kmerExpCounts = new double[numKmers+1];
		kmerObsProbs = new double[numKmers+1];
		kmerExpProbs = new double[numKmers+1];
		kmerWeights = new double[numKmers+1];
		for(int i=0; i<numKmers; i++){
			kmerObsCounts[i]=1.0;  //pseudo-count to avoid zero probabilities
			kmerExpCounts[i]=1.0;    //pseudo-count to avoid zero probabilities
			kmerObsProbs[i]=0.0;
			kmerExpProbs[i]=0.0;
			kmerWeights[i]=0.0;
		}
		kmer = new byte[kmerLength];
		fileExt = "."+kmerLength+"mer.unknown.seqbias";
	}
	
	/**
	 * Calculate the seq bias from a set of reads & annotation
	 * @param reads	ReadLoader for a given experiment
	 * @param genes A collection of genes that may be used to produce obs/exp k-mer counts
	 */
	public abstract void execute(ReadLoader reads, Collection<AGene> genes);
	
		
	/**
	 * Return the weight of a read (byte array)
	 * @param read A byte array containing the read sequence
	 * @return double representing the weight that should be associated with this read 
	 */
	public double getWeight(byte[] read){
		double weight = 1;
		for(int c=0; c<kmerLength; c++)
			kmer[c] = read[c];
		int index = seq2int(kmer);
		if(index!=-1)
			weight = kmerWeights[index];
		return(weight);
	}
	/**
	 * Return the weight of a read (byte array)
	 * @param read A String containing the read sequence
	 * @return double representing the weight that should be associated with this read
	 */
	public double getWeight(String read){
		double weight = 1;
		for(int c=0; c<kmerLength; c++)
			kmer[c] = (byte)read.charAt(c);
		int index = seq2int(kmer);
		if(index!=-1)
			weight = kmerWeights[index];
		return(weight);
	}
	
	/**
	 * Return the weight of a SAMRecord read
	 * @param s SAMRecord
	 * @return double representing the weight that should be associated with this read
	 */
	public double getWeight(SAMRecord s){
		if(s.getReadNegativeStrandFlag()){
			return(getWeight(SequenceUtils.reverseComplement(s.getReadBases())));
		}else{
			return(getWeight(s.getReadBases()));
		}
	}
	/**
	 * Print to screen
	 */
	public void print(){
		System.out.println("#K="+kmerLength+"\tContributingHits="+contributingHits);
		System.out.println("#i\tseq\tobsCounts\texpCounts\tobsProbs\texpProbs\tweights\tlog2(weight)");
		for(int i=0; i<numKmers; i++){
			System.out.println(i+"\t"+int2seq(i, kmerLength)+"\t"+kmerObsCounts[i]+"\t"+kmerExpCounts[i]+"\t"+kmerObsProbs[i]+"\t"+kmerExpProbs[i]+"\t"+kmerWeights[i]+"\t"+Math.log(kmerWeights[i])/LOG2);
		}
	}
	/**
	 * Save to file
	 * @param outFilename Output file name for k-mer seq bias file
	 */
	public void save(String outFilename){
		try {
			FileWriter outFW = new FileWriter(outFilename);
			outFW.write("#K="+kmerLength+"\tContributingHits="+contributingHits+"\n");
			outFW.write("#i\tseq\tobsCounts\texpCounts\tobsProbs\texpProbs\tweights\tlog2(weight)\n");
			for(int i=0; i<numKmers; i++){
				outFW.write(i+"\t"+int2seq(i, kmerLength)+"\t"+kmerObsCounts[i]+"\t"+kmerExpCounts[i]+"\t"+kmerObsProbs[i]+"\t"+kmerExpProbs[i]+"\t"+kmerWeights[i]+"\t"+Math.log(kmerWeights[i])/LOG2+"\n");
			}outFW.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
		
	/**
	 * Load a saved seq bias model from file
	 * @param filename
	 */
	public void load(String filename){
		try {
			BufferedReader reader = new BufferedReader(new FileReader(new File(filename)));
			String line;
			int pos=0;
			while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	        	if(line.charAt(0)!='#'){
		            String[] words = line.split("\\t");
		            
		            kmerObsCounts[pos] = new Double(words[2]);
		            kmerExpCounts[pos] = new Double(words[3]);
		            kmerObsProbs[pos] = new Double(words[4]);
		            kmerExpProbs[pos] = new Double(words[5]);
		            kmerWeights[pos] = new Double(words[6]);
		            
		            pos++;
	        	}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Try to load a saved model with the extension used for a given type of SeqBiasModel
	 * @param reads ReadLoader for a given experiment
	 * @return true for success, false for failure
	 */
	public boolean attemptLoad(ReadLoader reads) {
		String fName = reads.getSourcePath()+fileExt;
		File f = new File(fName);
		if(f.exists()){
			load(fName);
			return true;
		}
		return false;
	}
	
	//byte to int
	public int base2int(byte base) {
		int intVal = -1;
	    switch (base) {
	    	case 'A':
	    	case 'a':
	    		intVal = 0;
	    		break;
	    	case 'C':
	    	case 'c':
	    		intVal = 1;
	    		break;
	    	case 'G':
	    	case 'g':
	    		intVal = 2;
	    		break;
	    	case 'T':
	    	case 't':
	    		intVal = 3;
	    		break;
	    	default:
	    		intVal = -1;
	    }
	    return intVal;
	}
	
	//integer to character
	public char int2base(int x) {
		char base;
	    switch (x) {
	    	case 0:
	    		base = 'A';
	    		break;
	    	case 1:
	    		base = 'C';
	    		break;
	    	case 2:
	    		base = 'G';
	    		break;
	    	case 3:
	    		base = 'T';
	    		break;
	    	default:
	    		throw new IllegalArgumentException("Invalid int: " + x);
	    }
	    return (base);
	}
		
	//Sequence (in byte array) to integer
	public int seq2int(byte[] seq) {
		int intVal = 0;
	    int len = seq.length;
	
	    for (int i = 0; i < len; i++) {
	    	long currInt = base2int(seq[i]);
	    	if (currInt == -1) {
	    		return -1;
	    	}
	    	intVal = intVal << 2;
	    	intVal += currInt;
	    }
	    return intVal;
	}
	//Sequence (in String) to integer
	public int seq2int(String seq) {
		int intVal = 0;
	    int len = seq.length();
	
	    for (int i = 0; i < len; i++) {
	    	long currInt = base2int((byte)seq.charAt(i));
	    	if (currInt == -1) {
	    		return -1;
	    	}
	    	intVal = intVal << 2;
	    	intVal += currInt;
	    }
	    return intVal;
	}
	
	//integer back to sequence string
	public String int2seq(long x, int kmerLen) {
	    /**
	     * check that the x is valid for the specified maxKmerLen. Note: 4 << (2 *
	     * (kmerLen - 1)) = 4^kmerLen
	     */
		if (x > ((4 << (2 * (kmerLen - 1))) - 1)) {
			throw new IllegalArgumentException("Invalid int value, " + x + ", for kmerLen " + kmerLen);
	    }
		StringBuffer seq = new StringBuffer(kmerLen);
	    for (int i = 0; i < kmerLen; i++) {
	    	int baseVal = (int) (x % 4);
	    	seq.append(int2base(baseVal));
	    	x = x >> 2;
	    }
	    return seq.reverse().toString();
	}
	
	//byte array to String
	private String bytes2String(byte[] b){
		StringBuffer buf = new StringBuffer(b.length);
		for(int i=0; i<b.length; i++)
			buf.append((char)b[i]);
		return buf.toString();
	}
	
	
	
	/**
	 * Main method for testing
	 * @param args Command-line arguments
	 */
	public static void main(String[] args) {
		int k=6;
		ArgParser ap = new ArgParser(args);
		if((!ap.hasKey("species")&&!ap.hasKey("fa")) || !ap.hasKey("sam") || !ap.hasKey("gtf")) { 
            System.err.println("Usage:\n" +
                               "  --species <species;genome>\n" +
                               "    OR\n" +
                               "  --fa <FASTA file>\n" +
                               "  --sam <SAM file of reads>\n" +
                               "  --gtf <GTF file>\n"+
                               "  --k <kmer length>\n" +
                               "  --model [obsexp/hansen]"+
                               "");
            return;
        }
        try {
        	GenomeLoader gLoad=null;
        	if(ap.hasKey("species")){
				Pair<Organism, Genome> pair = Args.parseGenome(args);
				Genome currgen = pair.cdr();
				gLoad = new GenomeLoader(currgen);
        	}else{
        		String faFile = ap.getKeyValue("fa");
        		gLoad = new GenomeLoader(new File(faFile), true);
        	}
			List<File> samFiles = Args.parseFileHandles(args, "sam");
			k = Args.parseInteger(args, "k", 6);		
			String gtfFile = ap.getKeyValue("gtf");	
			String modelType = ap.getKeyValue("model");	
			
			//Load genes
			GTFAnnotationLoader reader = new GTFAnnotationLoader(new File(gtfFile), gLoad);
			List<AGene> geneSet = reader.loadGenes();
			GenomeSegmenter segmenter = new GenomeSegmenter(gLoad, 5000000, 100000);
			Map<String, List<Region>> subChrRegions = segmenter.segmentWithGenes(geneSet);

			//Load experiments, make bias models
			ArrayList<ReadLoader> expts = new ArrayList<ReadLoader>();
			HashMap<String, SeqBiasModel> biasModels = new HashMap<String, SeqBiasModel>();
			for(File f : samFiles){
				ReadLoader currLoader = new SAMReadLoader(f, "", gLoad.getGenome(), gLoad.getNameTranslator());
				expts.add(currLoader);
			
				String currName = f.getName();
				currName = currName.replaceAll(".bam", "");
				currName = currName.replaceAll(".sam", "");
				String outName = currName+"."+k+"mer."+modelType+".seqbias";
				
				if(modelType.equals("hansen") || modelType.equals("obsexp")){
					SeqBiasModel model=null;
					if(modelType.equals("hansen"))
						model = new HansenSeqBiasModel(gLoad, k);
					else
						model = new ObsExpSeqBiasModel(gLoad, k, subChrRegions);
					
					//Run model
					biasModels.put(currName, model);
					model.execute(currLoader, geneSet);					
				}else{
					System.err.println("Error: model type \""+modelType+"\" not recognized.");
				}
			}
	    } catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
