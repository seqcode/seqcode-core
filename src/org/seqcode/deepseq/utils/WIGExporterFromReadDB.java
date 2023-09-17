package org.seqcode.deepseq.utils;

import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import org.seqcode.data.readdb.Client;
import org.seqcode.data.readdb.ClientException;
import org.seqcode.data.seqdata.SeqAlignment;
import org.seqcode.data.seqdata.SeqDataLoader;
import org.seqcode.data.seqdata.SeqLocator;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.NotFoundException;
import org.seqcode.gseutils.Pair;


/**
 * Outputs a fixed-step WIG format file for a readDB experiment (uses Client.getWeightHistogram)
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class WIGExporterFromReadDB {
	private Genome gen;
//	private Sample sample;
	SeqDataLoader sdloader=null;
	protected int [] stackedHitCounts;
	private int winSize=1, winStep=1;
	private int readExt=100;
	private String outName="out";
	private String trackName="out";
	private String trackDesc="out";
	private String trackColor="0,0,255";
	private int trackYMax=-1;
	private int perBaseMax=2;
	private char strand='.';
	private boolean loadType1=true;
	private boolean loadType2=false;
	private boolean loadPairs=false;
	private Client client=null;
	private List<String> exptNames = new ArrayList<String>();
	private List<SeqAlignment> aligns = new ArrayList<SeqAlignment>();
	private Collection<String> alignIDs= new ArrayList<String>();
	private HashMap<SeqAlignment, Set<Integer>> availSingleChroms = new HashMap<SeqAlignment, Set<Integer>>();
	private HashMap<SeqAlignment, Set<Integer>> availSingleType2Chroms = new HashMap<SeqAlignment, Set<Integer>>();
	private HashMap<SeqAlignment, Set<Integer>> availPairedChroms = new HashMap<SeqAlignment, Set<Integer>>();
	private boolean hasPairedAligns=false;
	private List<Integer> fivePrimePosList = null;
	private List<Float> fivePrimeCountsList = null;
	
	private int conncount = 0; 

	
	public static void main(String[] args) throws SQLException, NotFoundException {
		
		WIGExporterFromReadDB wig = new WIGExporterFromReadDB(args);
		wig.execute();
		wig.close();
		
		System.exit(0);
	}
	
	
	public WIGExporterFromReadDB(String [] args) {
		if(args.length==0){
			System.err.println("WIGExporterFromReadDB usage:\n" +
					"\t--species <organism;genome>\n" +
					"\t--rdbexpt <experiment names>\n" +
					"\t--readext <extension length from 1st mapped base. if -1, plot actual read length>\n" +
					"\t--pbmax <max read count per base>\n" +
					"\t--strand <+/-/. limit to reads from one strand>\n" +
					"\t--winsize <window size/step in WIG file>\n" +
					"\t--name <string to use as track name>\n" +
					"\t--description <string to use as track description>\n" +
					"\t--ylimit <default track y max>\n" +
					"\t--color <R,G,B>\n" +
					"\t--out <output file name>");
			System.exit(1);
		}
		ArgParser ap = new ArgParser(args);
		GenomeConfig gcon = new GenomeConfig(args);
		
		// WIGExporter specific options
		outName = Args.parseString(args,"out",outName);
		trackName = Args.parseString(args,"name",trackName);
		trackDesc = Args.parseString(args,"description",trackDesc);
		trackColor = Args.parseString(args,"color",trackColor);
		readExt = Args.parseInteger(args,"readext",readExt);
		String strandStr = Args.parseString(args,"strand",".");
		strand = strandStr.charAt(0);
		winSize = Args.parseInteger(args,"winsize",winSize);
		
		//General options processed directly by ExptConfig
		
		perBaseMax = Args.parseInteger(args,"pbmax",perBaseMax);
		if(ap.hasKey("ylimit")){trackYMax=Args.parseInteger(args,"ylimit",-1);}
	    winStep=winSize;
	    gen = gcon.getGenome();
	    
	    try{
		    		
		    // Load the experiments
		    sdloader = new SeqDataLoader();
		    List<SeqLocator> rdbexpts = Args.parseSeqExpt(args,"rdbexpt");
		    if (rdbexpts.size()>0){
		    	//Start a new ReadDB client
				if(client==null)
					client = new Client();
	
				//Process the SeqLocators
				for(SeqLocator locator : rdbexpts){
					String exptName = locator.getExptName(); exptNames.add(exptName);
					aligns.addAll(sdloader.loadAlignments(locator, gen));
				}
				for(SeqAlignment alignment : aligns) {
		        	if(client.exists(Integer.toString(alignment.getDBID()))){
			            alignIDs.add(Integer.toString(alignment.getDBID()));
			            hasPairedAligns = hasPairedAligns || (alignment.getAlignType().getName().equals("PAIRED") || alignment.getAlignType().getName().equals("MIXED"));
			            
			            //Find the available chromosomes for each alignment
			            availSingleChroms.put(alignment, new HashSet<Integer>());
			            if(alignment.getNumHits()>0){
			            	availSingleChroms.get(alignment).addAll(client.getChroms(Integer.toString(alignment.getDBID()), false,false, null));
			            }
			        	availSingleType2Chroms.put(alignment, new HashSet<Integer>());
			            if(alignment.getNumType2Hits()>0){
			            	availSingleType2Chroms.get(alignment).addAll(client.getChroms(Integer.toString(alignment.getDBID()), true,false, null));
			            }
			            availPairedChroms.put(alignment, new HashSet<Integer>());
			            if(alignment.getNumPairs()>0){
			            	availPairedChroms.get(alignment).addAll(client.getChroms(Integer.toString(alignment.getDBID()), false,true, null));
			            }
		        	}else{
		        		System.err.println("ReadDBHitLoader: Error: "+alignment.getExpt().getName()+";"+alignment.getExpt().getReplicate()+";"+alignment.getName()+"\tRDBID:"+alignment.getDBID()+" does not exist in ReadDB.");
		        		System.exit(1);
		        	}
			    }
		        
	           
		    	//ReadDBHitLoader = new ReadDBHitLoader(sdloader, gen, rdbexpts, loadType1, loadType2, loadPairs);
		    }else {
		    	System.err.println("Must provide READDB identifier");
		    	System.exit(1);
		    }
	    } catch (IOException e) {
			e.printStackTrace();
		} catch (ClientException e) {
			//Do nothing here: ClientException could be thrown because chromosome doesn't contain any hist. 
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
	    	    
	}
	
	public void close(){
		if(client!=null)
			client.close();
		if(sdloader!=null)
			sdloader.close();
	}
	
	public void execute(){
		try {
			FileWriter fw = new FileWriter(outName+".wig");
			
			double basesDone=0, printStep=10000000,  numPrint=0;
			if(trackName.equals("out"))
				trackName=outName;
			if(trackDesc.equals("out"))
				trackDesc=outName;
			
			//Print the header
			fw.write("track type=wiggle_0 name=\""+trackName+"\" description=\""+trackDesc+" summary\""+" visibility=full color="+trackColor+" ");
			if(trackYMax >0)
				fw.write("autoScale=off viewLimits=0:"+trackYMax+" ");
			fw.write("\n");
			       			
			ChromRegionIterator chroms = new ChromRegionIterator(gen);
			while(chroms.hasNext()){
				NamedRegion currentRegion = chroms.next();
				
				//Split the job up into chunks of 10Mbp
				for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=10000000){
					int y = x+10000000; 
					if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
					Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
					
					//Load Hits
					fivePrimePosList = new ArrayList<Integer>();
					fivePrimeCountsList = new ArrayList<Float>();
					if(strand=='+' || strand=='.'){
						Pair<ArrayList<Integer>,ArrayList<Float>> coordhits = loadStrandedBaseCounts(currSubRegion, '+');
						fivePrimePosList.addAll(coordhits.car());
						fivePrimeCountsList.addAll(coordhits.cdr());
					}
					if(strand=='-' || strand=='.'){
						Pair<ArrayList<Integer>,ArrayList<Float>> coordhits = loadStrandedBaseCounts(currSubRegion, '-');
						fivePrimePosList.addAll(coordhits.car());
						fivePrimeCountsList.addAll(coordhits.cdr());
					}
					
					
                    double stackedHitCounts[] = makeLandscape(fivePrimePosList, fivePrimeCountsList, currSubRegion, perBaseMax, strand);
                    
                    boolean recording=false;
	                //Scan regions
					for(int i=currSubRegion.getStart(); i<currSubRegion.getEnd()-(int)winSize; i+=(int)winStep){
						Region currWin = new Region(gen, currentRegion.getChrom(), i, (int)(i+winSize-1));
						
						int binid = (int)Math.max(0, ((double)(currWin.getStart()-currSubRegion.getStart())/winStep));
						double winHits=(double)stackedHitCounts[binid];
						
						if(winHits>0){
							if(!recording){
								//WIG is 1-based, inclusive (although bigWig is 0 based)
								fw.write("fixedStep chrom=chr"+currSubRegion.getChrom()+" start="+i+" step="+winStep+" span="+winSize+"\n");
								recording=true;
							}
							fw.write(String.format("%.1f\n", winHits));
						}else{
							recording=false;
						}
						
						//Print out progress
						basesDone+=winStep;
						if(basesDone > numPrint*printStep){
							if(numPrint%10==0){System.out.print(String.format("(%.0f)", (numPrint*printStep)));}
							else{System.out.print(".");}
							if(numPrint%50==0 && numPrint!=0){System.out.print("\n");}
							numPrint++;
						}
					}
				}
			}
			System.out.print("\n");
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private Pair<ArrayList<Integer>,ArrayList<Float>> loadStrandedBaseCounts(Region r, char strand){    	
        TreeMap<Integer,Float> allHits = null;
        ArrayList<Integer> coords = new ArrayList<Integer>();
        ArrayList<Float> counts = new ArrayList<Float>();
        try {
        	//Start a new ReadDB client
    		if(client==null) {
    			client = new Client();
    		}
    		
    		allHits = client.getWeightHistogram(alignIDs,
                                                r.getGenome().getChromID(r.getChrom()),
                                                loadType2,
                                                loadPairs,
                                                readExt,
                                                winSize,
                                                perBaseMax,
                                                r.getStart(),
                                                r.getEnd(),
                                                null,
                                                strand == '+');
            
            if (allHits == null) {
                if (alignIDs.size() != 0) {
                    //throw new NullPointerException("ReadDBHitLoader: client.getWeightHistogram returned null");
                }
            } else {
                coords.addAll(allHits.keySet());
                counts.addAll(allHits.values());
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ClientException e) {
            //Do nothing here; ClientException could be thrown because a chromosome doesn't contain any hits
        }
        return new Pair<ArrayList<Integer>,ArrayList<Float>>(coords, counts);
    }

	protected double[] makeLandscape(List<Integer> posList, List<Float> countsList, Region currReg, int perBaseMax, char strand){
        int numBins = (int)(currReg.getWidth()/winStep);
        double[] landscape = new double[numBins+1];
        for(int i=0; i<=numBins; i++){landscape[i]=0;}
        for(int x=0; x<posList.size(); x++){
            int pos = posList.get(x);
            float weight = countsList.get(x);
            int bin = inBounds((int)((double)(pos-currReg.getStart())/winStep), 0, numBins);
        	landscape[bin]+=weight;
        }
        return(landscape);
    }
	//keep the number in bounds
	protected final double inBounds(double x, double min, double max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
	protected final int inBounds(int x, int min, int max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
}
