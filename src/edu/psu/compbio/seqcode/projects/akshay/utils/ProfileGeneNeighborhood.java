package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.motifs.WeightMatrixScoreProfile;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.motifs.WeightMatrixScorer;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.projects.shaun.MotifAnalysisSandbox;
import edu.psu.compbio.seqcode.projects.shaun.rnaseq.GTFReader;
import edu.psu.compbio.seqcode.projects.shaun.rnaseq.genemodels.GeneTUnit;

public class ProfileGeneNeighborhood {
	
	private GenomeConfig gcon;
	
	private List<Point> peaks;
	private Map<String,String> peaksSeqs;
	
	private Map<String,List<Point>> peaksAtgenes;
	private Map<String,List<PointCluster>> peakClustersAtgenes;
	
	private Map<String, StrandedPoint> genes;
	private Map<String, Region> geneDomains;
	
	private int radius; // Neighborhood radius
	
	private int cluster_distance;
	
	private List<WeightMatrix> motifs = new ArrayList<WeightMatrix>();
	
	private int win;
	
	@SuppressWarnings("unchecked")
	private SequenceGenerator<Region> seqgen = gcon.getSequenceGenerator();
	
	
	public ProfileGeneNeighborhood(GenomeConfig g) {
		gcon = g;
	}
	
	/**
	 * 
	 * @param args
	 * @throws ParseException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, ParseException{
		ArgParser ap = new ArgParser(args);
		GenomeConfig gc = new GenomeConfig(args);
		ProfileGeneNeighborhood profiler = new ProfileGeneNeighborhood(gc);
		
		int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():150;
		int rad = ap.hasKey("radius") ? new Integer(ap.getKeyValue("radius")).intValue() : 50000;
		int clusD = ap.hasKey("clusD") ? new Integer(ap.getKeyValue("clusD")).intValue() : 300;
		
		String motiffile = ap.getKeyValue("motiffile");
		String backfile = ap.getKeyValue("back");
		
		List<WeightMatrix> matrixList = MotifAnalysisSandbox.loadMotifFromFile(motiffile, backfile, gc.getGenome());
		
		String peaksFile = ap.getKeyValue("peaks");
		List<Point> peaks = RegionFileUtilities.loadPeaksFromPeakFile(gc.getGenome(), ap.getKeyValue("peaks"), win);
		String genefile = Args.parseString(args, "gtf", null);
		
		//Read selected genes into a list
		String geneListFile = Args.parseString(args, "geneList", null);
		File gFile = new File(geneListFile);
		if(!gFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
        BufferedReader reader = new BufferedReader(new FileReader(gFile));
        String line;
        List<String> geneList = new ArrayList<String>();
        while ((line = reader.readLine()) != null) {
        	line.trim();
        	geneList.add(line);
        }
		
		profiler.setWindow(win);
		profiler.setClusDistance(clusD);
		profiler.setRadius(rad);
		profiler.setPeaks(peaks);
		profiler.setGenes(genefile, geneList);
		profiler.setMotifs(matrixList);
		
		boolean clusterSyntax = ap.hasKey("ClusterSyntax");
		if(clusterSyntax){profiler.printPeakClusterSyntax();}
		boolean peakSyntax = ap.hasKey("PeakSyntax");
		if(peakSyntax){profiler.printPeakSyntax();}
		
		profiler.clear();
		
		
		
		
	}
	
	
	
	//Mutators
	public void setWindow(int w){win = w;}
	public void setClusDistance(int c){cluster_distance = c;}
	public void setRadius(int r){radius =r;}
	public void setPeaks(List<Point> pts){peaks =pts; loadPeakSeqs();}
	public void setGenes(String gfffile, List<String> geneList){loadGenes(gfffile,geneList);}
	public void setMotifs(List<WeightMatrix> m){motifs = m;}
	
	
	
	public void printPeakClusterSyntax(){
		loadpeakClustersAtgenes();
		StringBuilder sb = new StringBuilder();
		int nMotifs = motifs.size();
		for(String gname : peakClustersAtgenes.keySet()){
			if(peakClustersAtgenes.get(gname).size() == 0){
				sb.append(gname);sb.append("\t");
				sb.append("NA\tNA\t");
				for(int m=0; m<nMotifs; m++){sb.append("NA"); sb.append("\t");}
				sb.deleteCharAt(sb.length()-1);
				sb.append("\n");
			}else{
				for(PointCluster pc: peakClustersAtgenes.get(gname)){
					sb.append(gname);sb.append("\t");
					int strand = genes.get(gname).getStrand() == '+' ? 1: -1;
					int sign = strand*(genes.get(gname).getLocation()-pc.getLocation()) > 0 ? -1 : +1;
					int distance = sign*pc.distance(genes.get(gname));
					sb.append(pc.getPeakLocationsString()+"\t");
					sb.append(distance);sb.append("\t");
					String seq = seqgen.execute(pc.expand(win/2));
					for(int m=0; m<nMotifs; m++){
						Pair<Integer,Double> mScore = bestMotif(seq,motifs.get(m));
						sb.append(mScore.cdr());sb.append("\t");
					}
					sb.deleteCharAt(sb.length()-1);
					sb.append("\n");
				}
			}
		}
		
		StringBuilder header = new StringBuilder();
		header.append("Gene\tPeak\tDistance\t");
		for(WeightMatrix wm : motifs){
			header.append(wm.getName());header.append("\t");
		}
		header.deleteCharAt(header.length()-1);
		System.out.println(header.toString());
		System.out.println(sb.toString());
	}
	
	public void printPeakSyntax(){
		loadpeaksAtgenes();
		
		StringBuilder sb = new StringBuilder();
		int nMotifs = motifs.size();
		for(String gname : peaksAtgenes.keySet()){
			if(peaksAtgenes.get(gname).size() == 0){
				sb.append(gname);sb.append("\t");
				sb.append("NA\tNA\t");
				for(int m=0; m<nMotifs; m++){sb.append("NA"); sb.append("\t");}
				sb.deleteCharAt(sb.length()-1);
				sb.append("\n");
			}else{
				for(Point p: peaksAtgenes.get(gname)){
					sb.append(gname);sb.append("\t");
					int strand = genes.get(gname).getStrand() == '+' ? 1: -1;
					int sign = strand*(genes.get(gname).getLocation()-p.getLocation()) > 0 ? -1 : +1;
					int distance = sign*p.distance(genes.get(gname));
					sb.append(p.getLocationString()+"\t");
					sb.append(distance);sb.append("\t");
					String seq = peaksSeqs.get(p.getLocation()); 
					for(int m=0; m<nMotifs; m++){
						Pair<Integer,Double> mScore = bestMotif(seq,motifs.get(m));
						sb.append(mScore.cdr());sb.append("\t");
					}
					sb.deleteCharAt(sb.length()-1);
					sb.append("\n");
				}
			}
		}
		
		StringBuilder header = new StringBuilder();
		header.append("Gene\tPeak\tDistance\t");
		for(WeightMatrix wm : motifs){
			header.append(wm.getName());header.append("\t");
		}
		header.deleteCharAt(header.length()-1);
		System.out.println(header.toString());
		System.out.println(sb.toString());
	}
	
	
	
	
	//Loaders
	
	private void loadPeakSeqs(){
		peaksSeqs = new HashMap<String, String>();
		for(Point p : peaks){
			String seq = seqgen.execute(p.expand(win/2));
			peaksSeqs.put(p.getLocationString(), seq);
		}
	}
	
	private void loadGenes(String gtfFile, List<String> geneList){
		GTFReader gffreader = new GTFReader(new File(gtfFile), gcon.getGenome());
		List<GeneTUnit> geneObjects = gffreader.loadGenes();
		Map<String, StrandedPoint> allGenes = new HashMap<String, StrandedPoint>();
		genes = new HashMap<String, StrandedPoint>();
		for(GeneTUnit gu: geneObjects){
			allGenes.put(gu.getName().toUpperCase(), new StrandedPoint(gu.getTSS(),gu.getStrand()));
		}
		for(String s : geneList){
			String gene_name = s.toUpperCase();
			if(allGenes.containsKey(gene_name)){
				genes.put(gene_name, allGenes.get(gene_name));
			}
		}
		
		geneDomains = new HashMap<String,Region>();
		for(String gene_name : genes.keySet()){
			geneDomains.put(gene_name, genes.get(gene_name).expand(radius));
		}
	}
	
	
	private void loadpeakClustersAtgenes(){
		Map<String,List<Point>> peaksbyChr = hashbychrom(peaks);
		//Sort the points
		for(String chrom: peaksbyChr.keySet()){
			Collections.sort(peaksbyChr.get(chrom));
		}
		
		for(String gene_name : geneDomains.keySet() ){
			String geneChr = geneDomains.get(gene_name).getChrom();
			peakClustersAtgenes.put(gene_name, new ArrayList<PointCluster>());
			if(peaksbyChr.containsKey(geneChr)){
				List<Point> nearbyPeaks = new ArrayList<Point>();
				for(Point p: peaksbyChr.get(geneChr)){
					if(geneDomains.get(gene_name).contains(p)){
						nearbyPeaks.add(p);
					}
				}
				
				Collections.sort(nearbyPeaks); // They should be sorted already, but just to make sure
				PointCluster lastadded = null;
				for(Point p : nearbyPeaks){
					if(lastadded == null){
						peakClustersAtgenes.get(gene_name).add(new PointCluster(p));
						lastadded = peakClustersAtgenes.get(gene_name).get(peakClustersAtgenes.get(gene_name).size()-1);
					}else{
						if(lastadded.getLastPeak().distance(p) < cluster_distance){
							peakClustersAtgenes.get(gene_name).get(peakClustersAtgenes.get(gene_name).size()-1).addpeak(p);
							lastadded = peakClustersAtgenes.get(gene_name).get(peakClustersAtgenes.get(gene_name).size()-1);
						}else{
							peakClustersAtgenes.get(gene_name).add(new PointCluster(p));
							lastadded = peakClustersAtgenes.get(gene_name).get(peakClustersAtgenes.get(gene_name).size()-1);
						}
					}
				}
			}
		}
	}
	
	
	private void loadpeaksAtgenes(){
		Map<String,List<Point>> peaksbyChr = hashbychrom(peaks);
		for(String gene_name : geneDomains.keySet()){
			String geneChr = geneDomains.get(gene_name).getChrom();
			peaksAtgenes.put(gene_name, new ArrayList<Point>());
			if(peaksbyChr.containsKey(geneChr)){
				for(Point p: peaksbyChr.get(geneChr)){
					if(geneDomains.get(gene_name).contains(p)){
						peaksAtgenes.get(gene_name).add(p);
					}
				}
			}
		}
	}
	
	
	
	private Pair<Integer,Double> bestMotif(String seq, WeightMatrix motif){
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		WeightMatrixScoreProfile profiler = scorer.execute(seq);
		@SuppressWarnings("unchecked")
		Pair<Integer, Double> ret = new Pair(profiler.getMaxIndex(), profiler.getMaxScore());
		return ret;
	}
	
	
	private Map<String, List<Point>> hashbychrom(List<Point> pts){
		Map<String, List<Point>> byChr = new HashMap<String, List<Point>>();
		for(Point p : pts){
			if(!byChr.containsKey(p.getChrom()))
				byChr.put(p.getChrom(), new ArrayList<Point>());
			byChr.get(p.getChrom()).add(p);
		}
		
		return byChr;
	}
	
	
	public void clear(){SequenceGenerator.clearCache();}
	
	
	/**
	 * Cluster of points, usually used to represent a cluster of ChIPSeq peaks. PointCluster extends Point.
	 * The Point attributes you get from the PointCluster are those of the leftmost Point in the cluster
	 * 
	 * @author akshaykakumanu
	 *
	 */
	public class PointCluster extends Point {
		
		List<Point> peaks;

		public PointCluster(Point p) {
			super(p.getGenome(), p.getChrom(), p.getLocation());
			// TODO Auto-generated constructor stub
			peaks = new ArrayList<Point>();
			peaks.add(p);
		}
		
		public void addpeak(Point p){
			if(p.getChrom().equals(getChrom())){
				peaks.add(p);
				Collections.sort(peaks);
				if(p.compareTo(peaks.get(0))<0){
					location = p.getLocation();
				}
			}
		}
		
		public Region expand(int distance){
			int ns = Math.max(1, peaks.get(0).getLocation() - distance);
			int ne = Math.min(peaks.get(peaks.size()-1).getLocation() + distance, g.getChromLength(chrom));
			return new Region(g, chrom, ns, ne);
		}
		
		public String getPeakLocationsString(){
			StringBuilder retSB = new StringBuilder();
			for(Point p: peaks){
				retSB.append(p.getLocationString()+";");
			}
			retSB.deleteCharAt(retSB.length()-1);
			return retSB.toString();
		}
		
		public Point getLastPeak(){return peaks.get(peaks.size()-1);}
		
		public int getNumPeaks(){return peaks.size();}
		
	}
	
	
	
	
	

}
