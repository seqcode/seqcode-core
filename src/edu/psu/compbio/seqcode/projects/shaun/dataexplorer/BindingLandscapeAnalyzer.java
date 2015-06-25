package edu.psu.compbio.seqcode.projects.shaun.dataexplorer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.PointParser;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RegionParser;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class BindingLandscapeAnalyzer {

	private Genome gen;
	private ArrayList<Region> towers = new ArrayList<Region>();
	//Settings
	private int winSize=200, winStep=100;
	private int readExt = 0;
	//Points
	private ArrayList<Region> regions = new ArrayList<Region>();
	private ArrayList<Point> points = new ArrayList<Point>();
	private int numRegions=0;
	private boolean entireGenome = false;
	//Data Sources
	private Collection<DataSource> motifSources;
	private Collection<DataSource> secondaryMotifSources;
	private Collection<DataSource> deepSeqSources;
	private Collection<DataSource> conservationSources;
	private Collection<DataSource> uniquenessSources; 
	//DataLoader here
	private String dataInfoFileName=null;
	private DataSourceLoader loader;
	//Combination rules
	private String motifComboRule;
	//Filter
	private DataCollectionFilter dataFilter;
	//Data Collection
	private DataCollection collection;
	
	//Main method
	public static void main(String[] args) throws IOException, ParseException {
		try {
			Pair<Species, Genome> pair = Args.parseGenome(args);
		
			if(pair==null){printError();return;}
		
			//Initialize
			BindingLandscapeAnalyzer analyzer = new BindingLandscapeAnalyzer(pair.cdr());
			
			//Options
			String workflow = Args.parseString(args, "workflow", "default");
			String outName = Args.parseString(args, "out", "out.data");
			//Towers
			String towerFileName = Args.parseString(args, "towers", null);
			if(towerFileName!=null)
				analyzer.loadTowers(towerFileName);
			
			if(Args.parseArgs(args).contains("win"))
				analyzer.setWinSize(Args.parseInteger(args, "win", 200));
			if(Args.parseArgs(args).contains("readext"))
				analyzer.setReadExt(Args.parseInteger(args, "readext", 0));
			//Points to load: entire genome, points, or random 
			analyzer.setEntireGenome(Args.parseFlags(args).contains("entiregenome"));
			if(Args.parseArgs(args).contains("points")){
			   analyzer.loadPointsFromFile(Args.parseString(args, "points",null));
			   analyzer.setEntireGenome(false);
			}if(Args.parseArgs(args).contains("random")){
				analyzer.loadRandomRegions(Args.parseInteger(args, "random", 10000));
				analyzer.setEntireGenome(false);
			}
			//Data Info file
			String infoFile = Args.parseString(args, "datainfo", null);
			analyzer.setDataInfoFile(infoFile);
			
			//Combination rules
			String mRule = Args.parseString(args, "motifcombo", "OR");
			analyzer.setMotifComboRule(mRule);
			
			
			//Generate the data
			if(workflow.equals("default"))
				analyzer.executeWorkflowDefault();
			else if (workflow.equals("peakdata"))
				analyzer.executeWorkflowPeakData();
			else if (workflow.equals("binned"))
				analyzer.executeWorkflowBinned();
			
			//Print
			analyzer.printData(outName);
			
			//Cleanup
			analyzer.close();
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
	}
	
	//Constructor
	public BindingLandscapeAnalyzer(Genome g){
		gen = g;
		loader = new DataSourceLoader(gen, readExt);
		System.out.println("BindingLandscapeAnalyzer initialized");		
	}
	
	//Accessors
	public void setEntireGenome(boolean e){entireGenome=e;}
	public void setWinSize(int w){winSize=w;}
	public void setReadExt(int e){readExt=e;}
	public void setDataInfoFile(String s){dataInfoFileName=s;}
	public void setMotifComboRule(String s){motifComboRule=s;}
	public boolean examineEntireGenome(){return(entireGenome);}
	
	//Load towers from a file
	public void loadTowers(String tFile){
		towers = loadRegionsFromFile(tFile);
	}
	
	//Execution
	/**
	 * Default workflow:
	 * Filter the input regions for those containing ANY motifs
	 * Get the DeepSeq data for the filtered set
	 * Combine all DeepSeq data with weights
	 */
	public void executeWorkflowDefault() throws IOException, ParseException {
		DataCollection dataCollection;
		motifSources=loader.getMotifDataSources(dataInfoFileName);
		secondaryMotifSources=loader.getMotifDataSources(dataInfoFileName, true);
		deepSeqSources =loader.getDeepSeqDataSource(dataInfoFileName); 
		conservationSources = loader.getConservationDataSource(dataInfoFileName); 
		uniquenessSources = loader.getUniquenessDataSource(dataInfoFileName); 
		dataFilter = new DataCollectionFilter(motifComboRule);

		//Primary motif collections first
		dataCollection = new DataCollection(motifSources, regions, -1);
		if(motifSources.size()>0){
			dataFilter.filter(dataCollection);
			DataSourceCombiner mCombo = new DataSourceCombiner("MAX", dataCollection, motifSources, "MAXMOTIF", 0, 1);
			dataCollection.addColumn(mCombo);
		}
		System.out.println(dataCollection.getRegions().size()+" regions pass the filters");
		
		//Secondary motifs 
		dataCollection.addColumns(secondaryMotifSources);
		
		//Now get the DeepSeqData, Conservation & Uniqueness sets
		System.out.println("Loading DeepSeqData");
		dataCollection.addColumns(deepSeqSources);
		System.out.println("Loading ConservationData");
		dataCollection.addColumns(conservationSources);
		System.out.println("Loading UniquenessData");
		dataCollection.addColumns(uniquenessSources);
		//Filter again?
		
		//A combined DeepSeq data source
		if(deepSeqSources.size()>0){
			System.out.println("Combining DeepSeqData");
			DataSourceCombiner dsCombo = new DataSourceCombiner("SUM", dataCollection, deepSeqSources, "DeepSeqCombination", 0, 1);
			dataCollection.addColumn(dsCombo);
		}
		
		collection=dataCollection;
		//dataCollection.print();
	}
	
	/**
	 * Workflow for simple analysis of read counts at peaks:
	 * Filter the input regions for those containing ANY motifs (if provided)
	 * Get the DeepSeq data for the filtered set
	 */
	public void executeWorkflowPeakData() throws IOException, ParseException {
		DataCollection dataCollection;
		motifSources=loader.getMotifDataSources(dataInfoFileName);
		secondaryMotifSources=loader.getMotifDataSources(dataInfoFileName, true);
		deepSeqSources =loader.getDeepSeqDataSource(dataInfoFileName); 
		conservationSources = loader.getConservationDataSource(dataInfoFileName); 
		uniquenessSources = loader.getUniquenessDataSource(dataInfoFileName); 
		dataFilter = new DataCollectionFilter(motifComboRule);

		//Primary motif collections first
		dataCollection = new DataCollection(motifSources, regions, -1);
		if(motifSources.size()>0){
			dataFilter.filter(dataCollection);
			DataSourceCombiner mCombo = new DataSourceCombiner("MAX", dataCollection, motifSources, "MAXMOTIF", 0, 1);
			dataCollection.addColumn(mCombo);
		}
		System.out.println(dataCollection.getRegions().size()+" regions pass the filters");
		
		//Secondary motifs 
		dataCollection.addColumns(secondaryMotifSources);
		
		//Now get the DeepSeqData, Conservation & Uniqueness sets
		System.out.println("Loading DeepSeqData");
		dataCollection.addColumns(deepSeqSources);
		System.out.println("Loading ConservationData");
		dataCollection.addColumns(conservationSources);
		System.out.println("Loading UniquenessData");
		dataCollection.addColumns(uniquenessSources);
		
		collection=dataCollection;
	}
	
	/**
	 * Workflow for simple analysis of binned read counts around peaks:
	 * Filter the input regions for those containing ANY motifs (if provided)
	 * Get the DeepSeq data for the filtered set
	 */
	public void executeWorkflowBinned() throws IOException, ParseException {
		DataCollection dataCollection;
		motifSources=loader.getMotifDataSources(dataInfoFileName);
		secondaryMotifSources=loader.getMotifDataSources(dataInfoFileName, true);
		deepSeqSources =loader.getDeepSeqDataSource(dataInfoFileName); 
		conservationSources = loader.getConservationDataSource(dataInfoFileName); 
		uniquenessSources = loader.getUniquenessDataSource(dataInfoFileName); 
		dataFilter = new DataCollectionFilter(motifComboRule);

		//Primary motif collections first
		dataCollection = new DataCollection(motifSources, regions, 50);
		if(motifSources.size()>0){
			dataFilter.filter(dataCollection);
			DataSourceCombiner mCombo = new DataSourceCombiner("MAX", dataCollection, motifSources, "MAXMOTIF", 0, 1);
			dataCollection.addColumn(mCombo);
		}
		System.out.println(dataCollection.getRegions().size()+" regions pass the filters");
		
		//Secondary motifs 
		dataCollection.addColumns(secondaryMotifSources);
		
		//Now get the DeepSeqData, Conservation & Uniqueness sets
		System.out.println("Loading DeepSeqData");
		dataCollection.addColumns(deepSeqSources);
		System.out.println("Loading ConservationData");
		dataCollection.addColumns(conservationSources);
		System.out.println("Loading UniquenessData");
		dataCollection.addColumns(uniquenessSources);
		
		collection=dataCollection;
	}
	
	////////////////////////Loaders/////////////////////////////
	//Points from a file
	public void loadPointsFromFile(String f){
		try {
			File pFile = new File(f);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
			
	        Pattern ptpatt = Pattern.compile("([^:\\s]+):(\\d+)");
    		Pattern strptpatt = Pattern.compile("([^:\\s]+):(\\d+):([^:\\s]+)");
	        reader.mark(500);
	        String line= reader.readLine(); //skip header line
	        if(!line.startsWith("Region")){reader.reset();} //don't skip header line if there is no header
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            Matcher m =null;
	            
	            if(words.length>=3){
	            	 m = ptpatt.matcher(words[2]);
	            }else if(words.length==1){
	            	m = ptpatt.matcher(words[0]);
	            }
	            if(m !=null && m.find()) { 
	            	String chrom = m.group(1);
	            	int location = Integer.parseInt(m.group(2));
	            	char strand = '?';
	                Matcher sm = strptpatt.matcher(line);
	                if(sm.find()){
	                	String strandstr = sm.group(3);
	                	if(strandstr.length() > 0) { strand = strandstr.charAt(0); }
	                }
	                Point p=null;

	                if(strand=='+' || strand == '-')
	            		p = new StrandedPoint(gen, chrom, location, strand);
	            	else
	            		p = new Point(gen, chrom, location);
		         
	                if(p!=null){
			         	int rstart = p.getLocation()-(winSize/2)<1 ? 1:p.getLocation()-(winSize/2);
			            int rend = p.getLocation()+(winSize/2)>gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+(winSize/2)-1;
			            Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
			            if(p instanceof StrandedPoint ){
			            	r = new StrandedRegion(r.getGenome(), r.getChrom(), r.getStart(), r.getEnd(), ((StrandedPoint) p).getStrand());
			            }
			            regions.add(r);
			            points.add(p);
			            numRegions++;
	                }
	            }
	        }System.out.println(points.size()+" regions loaded");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	//Regions from a file
	public ArrayList<Region> loadRegionsFromFile(String f){
		ArrayList<Region> res = new ArrayList<Region>();
		try {
			File pFile = new File(f);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
			
	        reader.mark(500);
	        String line= reader.readLine(); //skip header line
	        if(!line.startsWith("Region")){reader.reset();} //don't skip header line if there is no header
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            RegionParser rparser = new RegionParser(gen);
	            Region r=null;
	            if(words.length>=1)
		           	r = rparser.execute(words[0]);
	         	res.add(r);	            
	        }
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return(res);
	}
	//Points from random
	public void loadRandomRegions(int numSamples){
		int validSamples=0;
		int genomeSize=0, numChroms=0;
		long [] chromoSize = new long[gen.getChromList().size()];
		Random rand = new Random();
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
					if(randPos+winSize<total+chromoSize[c]){
						potential = new Region(gen, chromoNames[c], (int)(randPos-total), (int)(randPos+winSize-total));
						
						//is this overlapping an already sampled region?
						boolean over=false;
						for(Region r : regions){
							if(r.overlaps(potential)){
								over=true;
							}
						}
						if(!over){
							validSamples++;
							regions.add(potential);
							numRegions++;
							//System.out.println(potential.getChrom()+":"+potential.getStart()+"-"+potential.getEnd());
						}
					}
				}total+=chromoSize[c];
			}				
		}			
	}
	
	
	/////////////////////////////Extra methods///////////////////////////////////
	public void printData(){
		collection.print();
	}
	public void printData(String outFile){
		collection.print(outFile);
	}
	public void close(){
		collection.cleanup();
	}
	////////////////////////////////////////////////////////////////////////////
	public static void printError(){
		System.err.println("Usage:\n " +
                "BindingLandscapeAnalyzer \n " +
                "Required:\n  " +
                "--species <organism name;genome version> \n  " +
                "Options: \n  " +
                "--workflow <default/peakdata/binned>\n  " +
                "--entiregenome [make regions from whole genome> " +
                "--points <points file> --random <num random points>\n  " +
                "--datainfo <file listing motifs & experiments to load and their thresholds>\n" +
                "--motifcombo <AND/OR>\n" +
                "--towers <file with towers to screen out>\n  " +
                "--win <win size>\n  " +
                "--readext <read extension>\n  " +
                "--out <output file name>\n  "
                
                );
	}
}
