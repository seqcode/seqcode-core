package org.seqcode.motifs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.io.StreamGobbler;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gsebricks.verbs.motifs.WeightMatrixScoreProfile;
import org.seqcode.gsebricks.verbs.motifs.WeightMatrixScorer;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.Pair;

/**
 * Utility to run MEME within a java class and from a command line.  It also evaluates significance 
 * of the MEME discovered motifs against randomly picked genomic sequences and calculates ROC score. 
 * 
 * Input:
 * 		- Genome
 * 		- MEME path and MEME arguments
 * 		- Window around peaks in which to get sequences
 * Output:
 * 		- MEME motif files
 * 		- ROC score for each motif
 * 
 * @author akshay kakumanu
 */

public class MemeER {
	
	protected String MEMEpath;
	protected String MEMEargs;
	protected static String PWMfile = null;	// pwm output file name
	protected Float pseudo = (float)0.001;
	public static final int MOTIF_FINDING_NEGSEQ=10000;
	public static final double MOTIF_FINDING_ALLOWED_REPETITIVE = 0.2;
	protected double motifMinROC = 0.70;
	public MemeER(String path, String args) {
		this.MEMEpath = path;
		this.MEMEargs = args;		
	}
	
	// option to set motif minimum ROC
	public void setMotifMinROC(double minroc){motifMinROC = minroc;}
	public double getMotifMinROC(){return motifMinROC;}

	public Pair<List<WeightMatrix>,List<WeightMatrix>> execute(List<String> sequences, File memeOutDirFullName, boolean bestOnly){
		List<WeightMatrix> wm = new ArrayList<WeightMatrix>();
		List<WeightMatrix> fm = new ArrayList<WeightMatrix>();
		String memeOutDir = null;
		File workingDir = new File(System.getProperty("user.dir"));
		if(memeOutDirFullName == null){
			String wDir = System.getProperty("user.dir");
			memeOutDir = wDir+"/meme_out";
		}else{
			memeOutDir = memeOutDirFullName.getAbsolutePath();
		}
		
		try {
			//Set up the input file
			File seqFile= File.createTempFile("seq", ".fa", workingDir);
			String seqFilename = seqFile.getCanonicalPath();
			FileWriter fout = new FileWriter(seqFile);
			int sCount=1;
			for(String seq : sequences){
				fout.write(">Seq"+sCount+"\n"+seq+"\n");
				sCount++;
			}
			fout.close();
			
			//Test if meme directory exists. If it does, recursively delete contents
			File memeOutPath = new File(memeOutDir);
			if(memeOutPath.exists())
				deleteDirectory(memeOutPath);
			
			//Call MEME
			String MEMEcmd = MEMEpath+"/meme ";
			Process proc = Runtime.getRuntime().exec(MEMEcmd+" "+seqFilename+" "+MEMEargs +" -o "+memeOutDir);
			// any error message? 
			StreamGobbler errorGobbler = new 
			StreamGobbler(proc.getErrorStream(), "MEME_ERR", true);
			// any output? 
			StreamGobbler outputGobbler = new 
			StreamGobbler(proc.getInputStream(), "MEME_OUT", true);
			// kick them off 
			errorGobbler.start(); 
			outputGobbler.start(); 
			// any error??? 
			int exitVal = proc.waitFor(); 
			System.err.println("MEME ExitValue: " + exitVal);

            File memeOutFile = new File(memeOutDir+"/meme.txt");
            if (!memeOutFile.exists()) {
            	//Clean up intermediate files (fasta, etc)
                if(seqFile.exists())
                	seqFile.delete();
            	throw new FileNotFoundException("Can't find file " + memeOutFile.getName());
            }else{
	            BufferedReader memeReader = new BufferedReader(new FileReader(memeOutFile));
	            Map<String,Double> back = parseMEMEResultsForBackground(memeReader);
	            memeReader.close();
	            BufferedReader memeReader2 = new BufferedReader(new FileReader(memeOutFile));
	            List<Pair<WeightMatrix,Double>> currFM = parseMEMEResultsForFreqMatries(memeReader2);
	            memeReader2.close();
	            
	            //Extract best motif?
				if(bestOnly){
					double minScore = Double.MAX_VALUE; 
					WeightMatrix bestMotif=null; 
					for(Pair<WeightMatrix,Double> m : currFM)
						if(m.cdr()<minScore){
							minScore = m.cdr();
							bestMotif = m.car();
						}
					fm.add(bestMotif);
					WeightMatrix wMatrix = WeightMatrix.getLogOddsVersion(bestMotif, back);
					wm.add(wMatrix);
					//System.out.println(bestMotif.getName()+"\t"+minScore);
					//System.out.println(WeightMatrix.printMatrix(bestMotif));
				}else{
					for(Pair<WeightMatrix,Double> m : currFM){
						fm.add(m.car());
						WeightMatrix wMatrix = WeightMatrix.getLogOddsVersion(m.car(), back);
						wm.add(wMatrix);
					}
				}
			}
         
            //Clean up intermediate files (fasta, etc)
            if(seqFile.exists())
            	seqFile.delete();
            proc.destroy();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		return new Pair<List<WeightMatrix>,List<WeightMatrix>>(wm,fm);
	}
	
	protected Map<String,Double> parseMEMEResultsForBackground(BufferedReader memeOut){
		Map<String,Double> back = new HashMap<String,Double>();
		try {
			String line=memeOut.readLine();
			int lineno = 1; 
			while (line!=null && !line.matches(".*Background letter frequencies.*")) {
		      line = memeOut.readLine();
		      lineno++;
		    }
			line = memeOut.readLine();
			if(line!=null){
		    	try {
		    		String[] pieces = line.split("\\s+");
		    		double A = Double.parseDouble(pieces[1]);
		    		double C = Double.parseDouble(pieces[3]);
		    		double G = Double.parseDouble(pieces[5]);
		    		double T = Double.parseDouble(pieces[7]);
		    		back.put("A", A);
		    		back.put("C", C);
		    		back.put("G", G);
		    		back.put("T", T);
		    	}
		    	catch (NumberFormatException ex) {
		    		System.err.println("At line " + lineno + ": " + line);
		    		ex.printStackTrace();
		    		throw ex;
		    	}
		    	catch (ArrayIndexOutOfBoundsException ex) {
		    		System.err.println("At line " + lineno + ": " + line);
		    		ex.printStackTrace();
		    		throw ex;
		    	}
			}
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return back;
	}
	
	protected List<Pair<WeightMatrix,Double>> parseMEMEResultsForFreqMatries(BufferedReader memeOut){
		List<Pair<WeightMatrix,Double>> parsed = new ArrayList<Pair<WeightMatrix,Double>>();
		String line;
		int lineno = 1; 
		int motifCount=0;
		try {
			while((line = memeOut.readLine()) != null){
				while (line!=null && !line.matches(".*letter-probability matrix.*")) {
			      line = memeOut.readLine();
			      lineno++;
			    }
				if(line!=null){
					motifCount++;
				    String lenStr = line.replaceFirst("^.*w=\\s*", "");
				    lenStr = lenStr.replaceFirst("\\s*nsites=.*", "");
				    int length = Integer.parseInt(lenStr);
				    String EStr = line.replaceFirst("^.*E=\\s*", "");
				    double Eval = Double.parseDouble(EStr);
				    WeightMatrix matrix = new WeightMatrix(length);
				    matrix.setNameVerType("Motif"+motifCount, "freq", "MEME");
				    for (int i = 0; i < length; i++) {
				    	line = memeOut.readLine().replaceFirst("^\\s*", "");
				    	lineno++;
				    	try {
				    		String[] pieces = line.split("\\s+");
				    		float A = Float.parseFloat(pieces[0])+pseudo;
				    		float C = Float.parseFloat(pieces[1])+pseudo;
				    		float G = Float.parseFloat(pieces[2])+pseudo;
				    		float T = Float.parseFloat(pieces[3])+pseudo;
				    		float total = A+C+G+T;
				    		matrix.matrix[i]['A'] = A/total;
				    		matrix.matrix[i]['C'] = C/total;
				    		matrix.matrix[i]['G'] = G/total;
				    		matrix.matrix[i]['T'] = T/total;
				    	}
				    	catch (NumberFormatException ex) {
				    		System.err.println("At line " + lineno + ": " + line);
				    		ex.printStackTrace();
				    		throw ex;
				    	}
				    	catch (ArrayIndexOutOfBoundsException ex) {
				    		System.err.println("At line " + lineno + ": " + line);
				    		ex.printStackTrace();
				    		throw ex;
				    	}
				    }
				    matrix.setLogOdds();
				    parsed.add(new Pair<WeightMatrix,Double>(matrix, Eval));
				}
			}
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return parsed;
	}
	
	public static void main(String[] args) {
		try {
			ArgParser ap = new ArgParser(args);
			String memeargs = Args.parseString(args, "memeargs", " -dna -mod zoops -revcomp -nostatus ");
			int MEMEminw = Args.parseInteger(args, "mememinw", 6);
			//MEME maxw
			int MEMEmaxw = Args.parseInteger(args, "mememaxw", 18);
			//MEME nmotifs option
			int MEMEnmotifs = Args.parseInteger(args,"memenmotifs", 3);
			int WinSize = Args.parseInteger(args, "win", 200);
			
			if (!ap.hasKey("memepath")||!ap.hasKey("seq")||!ap.hasKey("locations")){
				System.err.println("Usage:\n " +
	                    "MemeER\n " +
	                    "--geninfo <genome info file> \n " +
	                    "--seq <fasta seq directory> \n " +
	                    "--memepath <path to the meme bin dir (default: meme is in $PATH)>\n " +
	                    "\nOPTIONS:\n " +
	                    "--win <window of sequence to take around peaks(default=200)>\n " +
	                	"--memenmotifs <number of motifs MEME should find for each condition (default=3)>\n " +
	                	"--mememinw <minw arg for MEME (default=6)>\n " +
	                	"--mememaxw <maxw arg for MEME (default=18)>\n " +
	                	"--memeargs <additional args for MEME (default=  -dna -mod zoops -revcomp -nostatus)>\n " +
	                	"--out <output file prefix>\n " +
	                	"--printPWM [flag to print PWM]\n " +
	                	"--minROC <min ROC required for pwm output (default=0.7)>\n " +
	                    "");
				System.exit(0);
				
			}
			
			String GenPath = ap.getKeyValue("seq");
			memeargs = memeargs + " -nmotifs "+MEMEnmotifs + " -minw "+MEMEminw+" -maxw "+MEMEmaxw;
			
			// Multiple input files
			String[] points = ap.getKeyValue("locations").split(";"); 
			List<WeightMatrix> selectedMotifs = new ArrayList<WeightMatrix>();
			List<Double> selectedMotifsRocs = new ArrayList<Double>();
			
			GenomeConfig gcon = new GenomeConfig(args);
			
			Genome gen = gcon.getGenome();
			
			
			MemeER meme = new MemeER(Args.parseString(args, "memepath", ""), memeargs);
			
			// if minROC is provided, set the new value
			if (ap.hasKey("minROC")){
				meme.setMotifMinROC(Args.parseDouble(args, "minROC", 0.7));
			}
			// specify output directory if provided
			String outFolderName = null;
			if (ap.hasKey("out")){
				outFolderName = ap.getKeyValue("out");
			}else{
				outFolderName = System.getProperty("user.dir");
			}	
			File outFolder = new File(outFolderName);
			outFolder.mkdirs();		
			PrintWriter writer = null;
			// print PWM if specified
			if (ap.hasKey("printPWM")){
				writer = new PrintWriter(new File(outFolderName+File.separator+"pwm.out"));
			}
			
			for(int p=0; p<points.length; p++){
				List<Region> search_regs = RegionFileUtilities.loadRegionsFromPeakFile(gen, points[p], WinSize);
			
			
				SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
				seqgen.useCache(true);
				seqgen.useLocalFiles(true);
				seqgen.setGenomePath(GenPath);
				
				SequenceGenerator.setOffRegionCache();
			
				String[] randomSequences=new String[5000];
			
				if (!seqgen.isRegionCached()){
					System.err.println("Caching sequences");
					List<Region> randomRegions = randomRegionPick(gen, null, MemeER.MOTIF_FINDING_NEGSEQ,WinSize);
					randomSequences = seqgen.setupRegionCache(search_regs, randomRegions);
					System.err.println("Caching completed");
				}
				
				List<String> seqs = new ArrayList<String>();
			
				for(int i=0; i<search_regs.size(); i++){
					String currSeq = seqgen.execute(search_regs.get(i));
					if(lowercaseFraction(currSeq)<=MOTIF_FINDING_ALLOWED_REPETITIVE){
						seqs.add(currSeq);
					}
				}
			
		//		Pair<List<WeightMatrix>,List<WeightMatrix>> matrices = meme.execute(seqs, null, false);
				// allowing to specify output directory
				Pair<List<WeightMatrix>,List<WeightMatrix>> matrices = meme.execute(seqs, new File(outFolderName+File.separator+"meme_out"), false);
				List<WeightMatrix> wm = matrices.car();
				List<WeightMatrix> fm = matrices.cdr();
			
				if(wm.size()>0){
					//Evaluate the significance of the discovered motifs
					double rocScores[] = meme.motifROCScores(wm,seqs,randomSequences);
					System.err.println("MEME results for:" );
					for(int w=0; w<fm.size(); w++){
						if(fm.get(w)!=null){
							System.err.println("\t"+fm.get(w).getName()+"\t"+ WeightMatrix.getConsensus(fm.get(w))+"\tROC:"+String.format("%.2f",rocScores[w]));
						}
						if(rocScores[w] > meme.getMotifMinROC()){
							selectedMotifs.add(fm.get(w));
							selectedMotifsRocs.add(rocScores[w]);
						}
					}
				}
			}
			
			
			//Printing the selected motifs
			
			for(int m =0; m<selectedMotifs.size(); m++){
				// make motif ID consistent between the ROC and PWM
				String out = WeightMatrix.printTransfacMatrix(selectedMotifs.get(m),selectedMotifs.get(m).getName());
//				String out = WeightMatrix.printTransfacMatrix(selectedMotifs.get(m),"Motif_"+Integer.toString(m));
				System.err.println(out);
				
				// optional : print selected motif ROC and PWM to a file
				if (writer != null){
					writer.println("\t"+selectedMotifs.get(m).getName()+"\t"+ WeightMatrix.getConsensus(selectedMotifs.get(m))+"\tROC:"+String.format("%.2f",selectedMotifsRocs.get(m)));
					writer.println(out);
				}
			}
			// close the pwm output file
			if (writer != null){writer.close();}
			
				
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	}
	
	public static double lowercaseFraction(String seq){
		double count = 0;
		for (char c:seq.toCharArray())
			if (Character.isLowerCase(c) || c=='N')
				count++;
		return count/(double)seq.length();
	}
	
	public static boolean deleteDirectory(File path) {
	    if( path.exists() ) {
	      File[] files = path.listFiles();
	      for(int i=0; i<files.length; i++) {
	         if(files[i].isDirectory()) {
	           deleteDirectory(files[i]);
	         }
	         else {
	           files[i].delete();
	         }
	      }
	    }
	    return( path.delete() );
	}
	
	public  double[] motifROCScores(List<WeightMatrix> matrices, List<String> posSeqs, String[] negSeqs){
		String[] testNegSeqs=negSeqs;
		double[] rocScores = new double[matrices.size()];
		
		//If lengths of posSeqs and negSeqs are not matched, 
		//we will get errors in our discriminative performance estimates
		int posAvgLen=0, negAvgLen=0;
		for(String posSeq : posSeqs){ posAvgLen+=posSeq.length(); }
		for(String negSeq : negSeqs){ negAvgLen+=negSeq.length(); }
		posAvgLen/=posSeqs.size();
		negAvgLen/=negSeqs.length;
		if(posAvgLen<negAvgLen && posAvgLen>0){
			testNegSeqs = new String[negSeqs.length];
			for(int s=0; s<negSeqs.length; s++)
				testNegSeqs[s] = negSeqs[s].substring(0, posAvgLen);
		}else if(posAvgLen > negAvgLen){
			System.err.println("WARNING: average length of positive sequences is "
					+ "greater than that of negatives. AUROCs will not be accurate");
		}
		
		int m=0;
		for(WeightMatrix motif : matrices){
			List<Double> posScores = new ArrayList<Double>();
			List<Double> negScores = new ArrayList<Double>();
			if(motif!=null){
				WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
				for(String posSeq : posSeqs){
					WeightMatrixScoreProfile profiler = scorer.execute(posSeq);
					posScores.add(profiler.getMaxScore());
				}
				for(int s=0; s<testNegSeqs.length; s++){
					WeightMatrixScoreProfile profiler = scorer.execute(testNegSeqs[s]);
					negScores.add(profiler.getMaxScore());
				}
			}
			rocScores[m] = calcROCAUC(posScores, negScores);
			m++;
		}
		return rocScores;
	}
	
	protected  double calcROCAUC(List<Double> posMaxScores, List<Double> negMaxScores) {
		double auc = 0;
		if(posMaxScores.size()==0)
			return 0;
		if(negMaxScores.size()==0)
			return 1;
		ArrayList<LabeledDouble> data = new ArrayList<LabeledDouble>();
		for(Double d : posMaxScores)
			data.add(new LabeledDouble(d, 1));
		for(Double d : negMaxScores)
			data.add(new LabeledDouble(d, 0));
		
		Collections.sort(data);
		double pCount = (double)posMaxScores.size();
		double nCount = (double)negMaxScores.size();
		int x=0;
		double possum=0;
		double lastsn=0;
		double lastfpr=0;
		double lastdval = 10000000;
		
		for(LabeledDouble d : data){
			possum+=d.label;
			if(d.dat!=lastdval){
				double sn = possum/pCount;
				double fp = (x+1)-possum;
				double sp = (nCount-fp)/nCount;
				double fpr=1-sp;
				if(x>0){
						    //Rectangle             //Triangle
					auc += ((fpr-lastfpr)*lastsn) + ((sn-lastsn)*(fpr-lastfpr)/2);
				}
				lastfpr=fpr;
				lastsn = sn;
			}
			lastdval = d.dat;
			x++;
		}
		return auc;
	}
	
	public  class LabeledDouble implements Comparable<LabeledDouble>{
		public Double dat;
		public Integer label;
		public LabeledDouble(Double d, Integer i){dat=d; label=i;}
		public int compareTo(LabeledDouble ld) {
			if(dat > ld.dat){return(-1);}
			else if(dat < ld.dat){return(1);}
			else{return 0;}
		}
	}
	
	public static List<Region> randomRegionPick(Genome gen, List<Region> blackList, int numSamples, int sampleSize){
		List<Region> regs = new ArrayList<Region>();
		Random rand = new Random();
		int validSamples=0;
		
		//First see how big the genome is:
		int numChroms=0;
		long genomeSize=0;
		long [] chromoSize = new long[gen.getChromList().size()];
		String [] chromoNames = new String[gen.getChromList().size()];
		Iterator<NamedRegion> chroms = new ChromRegionIterator(gen);
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			genomeSize += (double)currentChrom.getWidth();
			chromoSize[numChroms]=currentChrom.getWidth();
			chromoNames[numChroms]=currentChrom.getChrom();
			numChroms++;				
		}

		//Now, iteratively generate random positions and check if they are valid and not overlapping repeats. 
		while(validSamples<numSamples){
			Region potential;				
			long randPos = (long)(1+(rand.nextDouble()*genomeSize));
			//find the chr
			boolean found=false;
			long total=0;
			for(int c=0; c<numChroms && !found; c++){
				if(randPos<total+chromoSize[c]){
					found=true;
					if(randPos+sampleSize<total+chromoSize[c]){
						potential = new Region(gen, chromoNames[c], (int)(randPos-total), (int)(randPos+sampleSize-total));
						
						//is this region in the blacklist? 
						boolean valid=true;
						if(blackList!=null){
							for(Region r : blackList){
								if(potential.overlaps(r)){valid=false;}
							}
						}
						if(valid){
							validSamples++;
							regs.add(potential);
						}
					}
				}total+=chromoSize[c];
			}
		}
		return(regs);
	}
}