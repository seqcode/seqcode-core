package edu.psu.compbio.seqcode.projects.shaun;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipData;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipDataset;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.ExptNameVersion;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.SQLData;
import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.species.Gene;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromosomeGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ExpanderIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.FileLineExpander;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RefGeneGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.iterators.SingleIterator;
import edu.psu.compbio.seqcode.gse.viz.scatter.Dataset2D;
import edu.psu.compbio.seqcode.gse.viz.scatter.ScatterPanel;

public class ArrayCorrelationAnalysis {

	private Genome gen;
	private boolean logScale =true;
	private static double valThres = 0.0;
	private final int geneTSSWindow=500;
	private HashMap<String, float[]> exprData;
	private boolean expressionLoaded=false;
	private boolean useAllArrays=false;
	private static boolean geneLevel=false;
	private static boolean t2tChange=true;
	private int inCommon=0, overThres=0;
	private ArrayList<Region> wellTiledRegions;
	
	
	private static String [] names = new String[] { "H3K79me2", "H3K4me3", "H3K27me3", "Suz12", "Ring1b" };
	private static String [] times = new String[] { "0", "2", "2+8hRA", "3", "4", "7"};
	
	public static void main(String[] args)  throws SQLException, NotFoundException{
		ArgParser ap = new ArgParser(args);
		if(!ap.hasKey("species") && !ap.hasKey("genome")){
			System.out.println("Usage: ArrayCorrelationAnalysis --species species --genome genomeVersion");
		}else{
			String [][][] experiments = new String[][][] { h3k79Versions, h3k4Versions, h3k27Versions, suz12Versions, ring1bVersions};
			Pair<Organism,Genome> pair = Args.parseGenome(args);
			
			ArrayCorrelationAnalysis aca = new ArrayCorrelationAnalysis(pair.cdr());
			
			aca.setGeneLevel(ap.hasKey("genes"));
			if(ap.hasKey("genes") && ap.hasKey("expression")){
				aca.loadExpression(ap.getKeyValue("expression"));
			}if(ap.hasKey("welltiled")){
				aca.loadWellTiledRegions(ap.getKeyValue("welltiled"));
			}
			
			String type;
			if(geneLevel){
				if(t2tChange)
					type = "T2T_Genes";
				else
					type = "Genes";
			}
			else{type = "Probes";}
			
			for(int i=0; i<experiments.length; i++){
				for(int j=-1; j<experiments.length; j++){
					if(i!=j){
					for(int t=1; t<times.length; t++){
						int it=t, jt=t;
						String [] exp1 = experiments[i][it];
						String [] exp1prev = t2tChange ? experiments[i][it-1]:null;
						String [] exp2=null, exp2prev=null;;
						if(j>=0){
							exp2 = experiments[j][jt];
							if(t2tChange)
								exp2prev = experiments[j][jt-1];
						}
							
						if(j==-1 && !exp1[0].equals("NONE") &&(!t2tChange || !exp1prev[0].equals("NONE"))){
							if(aca.isExpressionLoaded()){
								String currA = names[i]+"_"+times[it]+"_vs_"+"Expression_"+times[jt]+"_"+type;
								//System.err.println(currA);
								String scatterFile = currA+".jpg";
								Dataset2D data= aca.makeGeneExpressionDataset(exp1,exp1prev, jt);
								Dataset2D thresdata = aca.thresholdDataset(data, valThres);
								double correlation = aca.PearsonsCorrelation(thresdata);
								aca.MakeScatterPlot(data, scatterFile);
								System.out.println(type+"\t"+names[i]+"\t"+times[it]+"\tExpression\t"+times[jt]+"\t"+correlation+"\t"+scatterFile+"\t"+aca.getInCommon()+"\t"+aca.getOverThres());
							}
						}else if(j>=0 && !exp1[0].equals("NONE") && !exp2[0].equals("NONE") &&(!t2tChange || (!exp1prev[0].equals("NONE") && !exp2prev[0].equals("NONE")))){
							String currA = names[i]+"_"+times[it]+"_vs_"+names[j]+"_"+times[jt]+"_"+type;
							//System.err.println(currA);
							String scatterFile = currA+".jpg";
							Dataset2D data;
							if(geneLevel)
								data= aca.makeGeneDataset(exp1,exp1prev,exp2,exp2prev);
							else
								data= aca.makeProbeDataset(exp1, exp2);
							Dataset2D thresdata = aca.thresholdDataset(data, valThres);
							double correlation = aca.PearsonsCorrelation(thresdata);
							aca.MakeScatterPlot(data, scatterFile);
							System.out.println(type+"\t"+names[i]+"\t"+times[it]+"\t"+names[j]+"\t"+times[jt]+"\t"+correlation+"\t"+scatterFile+"\t"+aca.getInCommon()+"\t"+aca.getOverThres());
						}
					}
				}
				}
			}
		}	
	}

	public ArrayCorrelationAnalysis(Genome g){
		gen=g;
		//Initialize the allowed chips list (a hack, I know)
		//Unnormalized arrays
		allowedChips.add(new Integer(3104));
		allowedChips.add(new Integer(3344));
		allowedChips.add(new Integer(3324));
		allowedChips.add(new Integer(3364));
		allowedChips.add(new Integer(3103));
		allowedChips.add(new Integer(3102));
		allowedChips.add(new Integer(2122));
		allowedChips.add(new Integer(2882));
		allowedChips.add(new Integer(2462));
		allowedChips.add(new Integer(2463));
		allowedChips.add(new Integer(2702));
		allowedChips.add(new Integer(3203));
		allowedChips.add(new Integer(2782));
		allowedChips.add(new Integer(2124));
		allowedChips.add(new Integer(2125));
		allowedChips.add(new Integer(2602));
		allowedChips.add(new Integer(2603));
		allowedChips.add(new Integer(2123));
		allowedChips.add(new Integer(2883));
		allowedChips.add(new Integer(2464));
		allowedChips.add(new Integer(2465));
		//allowedChips.add(new Integer(2703));
		allowedChips.add(new Integer(3202));
		allowedChips.add(new Integer(2802));
		allowedChips.add(new Integer(2126));
		allowedChips.add(new Integer(2127));
		allowedChips.add(new Integer(2604));
		allowedChips.add(new Integer(3343));
		allowedChips.add(new Integer(3323));
		allowedChips.add(new Integer(3426));
		//allowedChips.add(new Integer(3363));
		allowedChips.add(new Integer(3425));
		allowedChips.add(new Integer(3424));
		allowedChips.add(new Integer(3423));
		allowedChips.add(new Integer(3422));
		//Normalized
		allowedChips.add(new Integer(3382));
		allowedChips.add(new Integer(3383));
		allowedChips.add(new Integer(3384));
		allowedChips.add(new Integer(3385));
		allowedChips.add(new Integer(3386));
		allowedChips.add(new Integer(3387));
		allowedChips.add(new Integer(3388));
		allowedChips.add(new Integer(3389));
		allowedChips.add(new Integer(3390));
		allowedChips.add(new Integer(3282));
		allowedChips.add(new Integer(3283));
		allowedChips.add(new Integer(3284));
		allowedChips.add(new Integer(3285));
		//allowedChips.add(new Integer(3286));
		allowedChips.add(new Integer(3287));
		allowedChips.add(new Integer(3288));
		allowedChips.add(new Integer(3289));
		allowedChips.add(new Integer(3290));
		allowedChips.add(new Integer(3291));
		allowedChips.add(new Integer(3292));
		allowedChips.add(new Integer(3262));
		allowedChips.add(new Integer(3263));
		allowedChips.add(new Integer(3264));
		allowedChips.add(new Integer(3265));
		allowedChips.add(new Integer(3266));
		allowedChips.add(new Integer(3267));
		allowedChips.add(new Integer(3268));
		allowedChips.add(new Integer(3269));
		allowedChips.add(new Integer(3270));
		allowedChips.add(new Integer(3271));
		allowedChips.add(new Integer(3443));
		allowedChips.add(new Integer(3444));
		allowedChips.add(new Integer(3442));		
	}
	
	public boolean isExpressionLoaded(){return(expressionLoaded);}
	public void setGeneLevel(boolean g){geneLevel = g;}
	public int getInCommon(){return inCommon;}
	public int getOverThres(){return overThres;}
	
	public void loadExpression(String exprFile){
		exprData = new HashMap<String, float[]>();
		
		try {
			File eFile = new File(exprFile);
			if(!eFile.isFile()){System.err.println("Invalid expression file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(eFile));
			String line = reader.readLine(); //Get rid of header line
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String [] tokens = line.split("\\s+");
	            String name = tokens[0];
	            float [] vals = new float[tokens.length-2];
	            for(int i=2; i<tokens.length; i++){
	            	vals[i-2] = new Double(tokens[i]).floatValue();
	            }
	            exprData.put(name, vals);
	        }
	        expressionLoaded=true;
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void MakeScatterPlot(Dataset2D data, String saveFile){
		JFreeChart chart = generateChart(data, 1);
		try {
			ChartUtilities.saveChartAsJPEG(new File(saveFile), chart, 1000, 1000);
			//System.err.println("Image saved to: "+saveFile);			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public Dataset2D makeProbeDataset(String [] e1, String [] e2){
		float[][] values = null;
		ChipChipDataset ccset = new ChipChipDataset(gen);
        ExptNameVersion expt1=null, expt2=null;
        if(e1.length==2)
        	expt1= new ExptNameVersion(e1[0], e1[1]);
        else if(e1.length==3)
        	expt1= new ExptNameVersion(e1[0], e1[1], e1[2]);
        if(e2.length==2)
        	expt2= new ExptNameVersion(e2[0], e2[1]);
        else if(e2.length==3)
        	expt2= new ExptNameVersion(e2[0], e2[1], e2[2]);
        
        ChipChipData data1, data2;
		try {
			data1 = ccset.getData(expt1);
			data2 = ccset.getData(expt2);
			
			if (data1 instanceof SQLData && data2 instanceof SQLData) {
				
				ChromosomeGenerator<Genome> chromgen = new ChromosomeGenerator<Genome>();
				Iterator<Region> chromiter = chromgen.execute(gen);
				Map<String,Double> values1 = new HashMap<String,Double>();
				Map<String,Double> values2 = new HashMap<String,Double>();
		
				int numCommon=0;
				while (chromiter.hasNext()) {
					Region chromregion = chromiter.next();
					data1.window(chromregion.getChrom(), chromregion.getStart(), chromregion.getEnd());
					data2.window(chromregion.getChrom(), chromregion.getStart(), chromregion.getEnd());
					
					for(int i = 0; i < data1.getCount(); i++) {
						values1.put(new String(chromregion.getChrom()+":"+data1.getPos(i)), getMeanRatio(data1, i));
					}
					for(int i = 0; i < data2.getCount(); i++) {
						String str = new String(chromregion.getChrom()+":"+data2.getPos(i));
						values2.put(str, getMeanRatio(data2, i));
						if(values1.containsKey(str)){
							numCommon++;
						}
					}
				}
				//System.err.println(numCommon+" values in common");
				inCommon = numCommon;
				values = new float[2][numCommon];
				int posind = 0;
				for(String pos : values1.keySet()) {
					if(values2.containsKey(pos)){
						//if(values1.get(pos)>=valThres || values2.get(pos)>=valThres){
						if(logScale){
							values[0][posind] = values1.get(pos)>0 ? (float)Math.log(values1.get(pos)) : 0;
							values[1][posind] = values2.get(pos)>0 ? (float)Math.log(values2.get(pos)) : 0;
						}else{
							values[0][posind]=values1.get(pos).floatValue();
							values[1][posind]=values2.get(pos).floatValue();
						}
						posind++;				
						//}
					}
				}
			}
		} catch (NotFoundException e3) {
			// TODO Auto-generated catch block
			e3.printStackTrace();
		}
		Dataset2D d = new Dataset2D(values,  expt1.getName(), expt2.getName());
		return(d);
	}
	
	public Dataset2D makeGeneDataset(String [] e1, String [] e1p,String [] e2,String [] e2p){
		float[][] values = null;
		ChipChipDataset ccset = new ChipChipDataset(gen);
		ExptNameVersion expt1=null, expt1prev=null, expt2=null, expt2prev=null;
        if(e1.length==2)
        	expt1= new ExptNameVersion(e1[0], e1[1]);
        else if(e1.length==3)
        	expt1= new ExptNameVersion(e1[0], e1[1], e1[2]);
        if(e1p !=null && e1p.length==2)
        	expt1prev= new ExptNameVersion(e1p[0], e1p[1]);
        else if(e1p !=null && e1p.length==3)
        	expt1prev= new ExptNameVersion(e1p[0], e1p[1], e1p[2]);
        if(e2.length==2)
        	expt2= new ExptNameVersion(e2[0], e2[1]);
        else if(e2.length==3)
        	expt2= new ExptNameVersion(e2[0], e2[1], e2[2]);
        if(e2p !=null && e2p.length==2)
        	expt2prev= new ExptNameVersion(e2p[0], e2p[1]);
        else if(e2p !=null && e2p.length==3)
        	expt2prev= new ExptNameVersion(e2p[0], e2p[1], e2p[2]);
        
        ChipChipData data1, data2,data1prev=null,data2prev=null;
		try {
			data1 = ccset.getData(expt1);
			data2 = ccset.getData(expt2);
			if(t2tChange){
				data1prev = ccset.getData(expt1prev);
				data2prev = ccset.getData(expt2prev);					
			}
			
			if (data1 instanceof SQLData && data2 instanceof SQLData) {
				Map<String,Double> values1 = new HashMap<String,Double>();
				Map<String,Double> values2 = new HashMap<String,Double>();
		
				RefGeneGenerator<NamedRegion> geneGen = new RefGeneGenerator<NamedRegion>(gen, "refGene");
				ChromRegionIterator chroms = new ChromRegionIterator(gen);
				Iterator<Gene> genes = new ExpanderIterator<NamedRegion,Gene>(geneGen, chroms);

				int numCommon=0;
				while(genes.hasNext()) { 
					Gene gene = genes.next();
					Region chromregion;
					if(gene.getStrand()=='+')
						chromregion= new Region(gen, gene.getChrom(), gene.getStart()-geneTSSWindow, gene.getStart()+geneTSSWindow);
					else 
						chromregion= new Region(gen, gene.getChrom(), gene.getEnd()-geneTSSWindow, gene.getEnd()+geneTSSWindow);
					
					boolean valid = false;
					if(useAllArrays){valid=true;}
					else{
						for(Region r: wellTiledRegions){
							if(chromregion.overlaps(r)){
								valid=true;
							}
						}
					}
					if(valid){
						data1.window(chromregion.getChrom(), chromregion.getStart(), chromregion.getEnd());
						data2.window(chromregion.getChrom(), chromregion.getStart(), chromregion.getEnd());
						if(t2tChange){
							data1prev.window(chromregion.getChrom(), chromregion.getStart(), chromregion.getEnd());
							data2prev.window(chromregion.getChrom(), chromregion.getStart(), chromregion.getEnd());								
						}
						String str = new String(chromregion.getChrom()+":"+chromregion.getStart()+"-"+chromregion.getEnd());
						if(data1.getCount()>0){
							if(t2tChange)
								values1.put(str, getMeanWindowRatio(data1)/getMeanWindowRatio(data1prev));
							else
								values1.put(str, getMeanWindowRatio(data1));
						}
						
						if(data2.getCount()>0){
							if(t2tChange)
								values2.put(str, getMeanWindowRatio(data2)/getMeanWindowRatio(data2prev));
							else
								values2.put(str, getMeanWindowRatio(data2));
							if(values1.containsKey(str)){
								numCommon++;
							}
						}
					}
				}
				//System.err.println(numCommon+" values in common");
				inCommon = numCommon;
				values = new float[2][numCommon];
				int posind = 0;
				for(String pos : values1.keySet()) {
					if(values2.containsKey(pos)){
						//if(values1.get(pos)>=valThres || values2.get(pos)>=valThres){
						if(logScale){
							values[0][posind] = values1.get(pos)>0 ? (float)Math.log(values1.get(pos)) : 0;
							values[1][posind] = values2.get(pos)>0 ? (float)Math.log(values2.get(pos)) : 0;
						}else{
							values[0][posind]=values1.get(pos).floatValue();
							values[1][posind]=values2.get(pos).floatValue();
						}
						posind++;				
						//}
					}
				}
			}
		} catch (NotFoundException e3) {
			// TODO Auto-generated catch block
			e3.printStackTrace();
		}
		Dataset2D d = new Dataset2D(values,  expt1.getName(), expt2.getName());
		return(d);
	}
	
	public Dataset2D makeGeneExpressionDataset(String [] e1,String [] e1p, int etime){
		if(!expressionLoaded){System.err.println("Expression data not loaded");System.exit(1);}
		float[][] values = null;
		ChipChipDataset ccset = new ChipChipDataset(gen);
        ExptNameVersion expt1=null, expt1prev=null;
        if(e1.length==2)
        	expt1= new ExptNameVersion(e1[0], e1[1]);
        else if(e1.length==3)
        	expt1= new ExptNameVersion(e1[0], e1[1], e1[2]);
        if(e1p !=null && e1p.length==2)
        	expt1prev= new ExptNameVersion(e1p[0], e1p[1]);
        else if(e1p !=null && e1p.length==3)
        	expt1prev= new ExptNameVersion(e1p[0], e1p[1], e1p[2]);
        
        ChipChipData data1, data1prev=null;
		try {
			data1 = ccset.getData(expt1);
			if(t2tChange)
				data1prev = ccset.getData(expt1prev);
			
			if (data1 instanceof SQLData) {
				Map<String,Double> values1 = new HashMap<String,Double>();
				Map<String,Double> values2 = new HashMap<String,Double>();
				
				RefGeneGenerator<NamedRegion> geneGen = new RefGeneGenerator<NamedRegion>(gen, "refGene");
				ChromRegionIterator chroms = new ChromRegionIterator(gen);
				Iterator<Gene> genes = new ExpanderIterator<NamedRegion,Gene>(geneGen, chroms);

				int numCommon=0;
				while(genes.hasNext()) { 
					Gene gene = genes.next();
					Region chromregion;
					if(gene.getStrand()=='+')
						chromregion= new Region(gen, gene.getChrom(), gene.getStart()-geneTSSWindow, gene.getStart()+geneTSSWindow);
					else 
						chromregion= new Region(gen, gene.getChrom(), gene.getEnd()-geneTSSWindow, gene.getEnd()+geneTSSWindow);
					
					boolean valid = false;
					if(useAllArrays){valid=true;}
					else{
						for(Region r: wellTiledRegions){
							if(chromregion.overlaps(r)){
								valid=true;
							}
						}
					}
					if(valid){
						data1.window(chromregion.getChrom(), chromregion.getStart(), chromregion.getEnd());
						if(t2tChange)
							data1prev.window(chromregion.getChrom(), chromregion.getStart(), chromregion.getEnd());
						
						String str = new String(chromregion.getChrom()+":"+chromregion.getStart()+"-"+chromregion.getEnd());
						if(data1.getCount()>0){
							if(t2tChange)
								values1.put(str, getMeanWindowRatio(data1)/getMeanWindowRatio(data1prev));
							else
								values1.put(str, getMeanWindowRatio(data1));
							
							if(exprData.containsKey(gene.getID())){
								float [] evals = exprData.get(gene.getID());
								if(t2tChange){
									values2.put(str, new Double(evals[etime]-evals[etime-1]));
								}else{
									values2.put(str, new Double(evals[etime]));
								}
								numCommon++;						
							}
						}
					}
				}
				//System.err.println(numCommon+" values in common");
				inCommon = numCommon;
				values = new float[2][numCommon];
				int posind = 0;
				for(String pos : values1.keySet()) {
					if(values2.containsKey(pos)){
						//if(values1.get(pos)>=valThres || values2.get(pos)>=valThres){
						if(logScale){
							values[0][posind] = values1.get(pos)>0 ? (float)Math.log(values1.get(pos)) : 0;
							values[1][posind]=values2.get(pos).floatValue();
						}else{
							values[0][posind]=values1.get(pos).floatValue();
							values[1][posind]=values2.get(pos).floatValue();
						}
						posind++;				
						//}
					}
				}
			}
		} catch (NotFoundException e3) {
			// TODO Auto-generated catch block
			e3.printStackTrace();
		}
		Dataset2D d = new Dataset2D(values,  expt1.getName(), new String("Expression_Day"+times[etime]));
		return(d);
	}
	
	public Dataset2D thresholdDataset(Dataset2D d, double threshold){
		float thres = (float)threshold;
		if(logScale)
			thres = (float)Math.log(threshold);
		
		ArrayList <Pair<Float, Float>> vals = new ArrayList <Pair<Float,Float>>();
		
		for(int i=0; i<d.getCount(); i++){
			if(d.getVal(0, i)>=thres || d.getVal(1, i)>=thres){
				vals.add(new Pair<Float,Float>(new Float(d.getVal(0,i)), new Float(d.getVal(1,i))));
			}
		}
		float[][] values = new float[2][vals.size()];
		int posind = 0;
		for(Pair<Float, Float> p : vals) {
			values[0][posind] = p.getFirst().floatValue();
			values[1][posind] = p.getLast().floatValue();
			posind++;
		}
		//System.err.println(posind+" values pass the threshold");
		overThres = posind;
		Dataset2D e = new Dataset2D(values, d.getLabelOne(), d.getLabelTwo());
		return(e);
	}
	
	private JFreeChart generateChart(Dataset2D data, double sample){
		
		String key = new String(data.getLabelOne()+"_vs_"+data.getLabelTwo());
		ScatterPanel panel;
		if(geneLevel)
			 panel = new ScatterPanel(key, data, 1.0, Color.blue);
		else
			panel = new ScatterPanel(key, data, 1.0);
			
		return(panel.getChart());
	}
	
	private double getMeanRatio(ChipChipData d, int i) { 
		double sum = 0.0;
		double count = 0.0;
		for(int j = 0; j < d.getReplicates(i); j++) { 
			if(useAllArrays || allowedChips.contains(new Integer(d.getExptID(i, j)))){
				if(!Double.isInfinite(d.getRatio(i, j)) && !Double.isNaN(d.getRatio(i, j))){
					sum += d.getRatio(i, j);
					count += 1;
				}
			}
		}
		return count > 0.0 ? sum/count : 1.0;
	}
	private double getMeanWindowRatio(ChipChipData d) { 
		double sum = 0.0;
		double count = 0.0;
		for(int i=0; i<d.getCount(); i++){
			for(int j = 0; j < d.getReplicates(i); j++) { 
				if(useAllArrays || allowedChips.contains(new Integer(d.getExptID(i, j)))){
					if(!Double.isInfinite(d.getRatio(i, j)) && !Double.isNaN(d.getRatio(i, j))){
						sum += d.getRatio(i, j);
						count += 1;
					}
				}
			}
		}
		return count > 0.0 ? sum/count : 1.0;
	}
	
	public double PearsonsCorrelation(Dataset2D data){
		double cor = 0;
		double count=0, sum1=0, sum2=0, mean1=0, mean2=0, stddev1=0, stddev2=0, crossmean=0;
		
		//Means
		for(int i=0; i<data.getCount(); i++){
			if(!Double.isInfinite(data.getVal(0, i))&&!Double.isNaN(data.getVal(0, i))&&!Double.isInfinite(data.getVal(1, i))&&!Double.isNaN(data.getVal(1, i))){
				sum1+= data.getVal(0, i);
				sum2+= data.getVal(1, i);
				count++;
			}
		}mean1 = sum1/count;
		mean2 = sum2/count;

		//Deviations & Std. Devs.
		double tmp1=0, tmp2=0, devhold1=0, devhold2=0, crosssum=0;
		for(int i=0; i<data.getCount(); i++){
			if(!Double.isInfinite(data.getVal(0, i))&&!Double.isNaN(data.getVal(0, i))&&!Double.isInfinite(data.getVal(1, i))&&!Double.isNaN(data.getVal(1, i))){
				tmp1 = (data.getVal(0, i)-mean1);
				devhold1+= tmp1*tmp1;
				tmp2 = (data.getVal(1, i)-mean2);
				devhold2+= tmp2*tmp2;
				
				crosssum+= tmp1*tmp2;
			}
		}
		stddev1 = Math.sqrt(devhold1/count);
		stddev2 = Math.sqrt(devhold2/count);
		crossmean = crosssum/count;
		
		cor = crossmean/(stddev1*stddev2);
		
		return(cor);
	}
	
	public void loadWellTiledRegions(String fname){
		try {
			File f = new File(fname);
			if(!f.isFile()){System.err.println("Invalid well-tiled file name");System.exit(1);}
			wellTiledRegions = new ArrayList<Region>();
			BufferedReader reader = new BufferedReader(new FileReader(f));
			String line;
			while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String [] tokens = line.split("\\s+"); 
	            if(tokens.length==3){
					String chr = tokens[0];
					int start = Integer.parseInt(tokens[1]);
					int end = Integer.parseInt(tokens[2]);
					Region r = new Region(gen, chr, start, end);
					wellTiledRegions.add(r);
	            }
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static String[][] h3k4Versions, h3k27Versions, h3k27subVersions, h3k79Versions, ring1bVersions, suz12Versions, rarChipSeq, backChipSeq;
	static { 
		h3k4Versions = new String[][] { 
                { "Mm H3K4me3:HBG3:ES Stage vs H3:HBG3:ES Stage", "median linefit, quantile norm" },
                { "Mm H3K4me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA", "median linefit, quantile norm" },
                { "Mm H3K4me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm" },
                { "Mm H3K4me3:HBG3:2+1 day vs H3:HBG3:2+1 day", "median linefit, quantile norm" },
                { "Mm H3K4me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm" },
                { "Mm H3K4me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm" }
        };

        h3k27Versions = new String[][] { 
                { "Mm H3K27me3:HBG3:ES Stage vs H3:HBG3:ES Stage", "median linefit, quantile norm" },
                { "Mm H3K27me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA", "median linefit, quantile norm" },
                { "Mm H3K27me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm", "2 (2/1/08)" },
                { "Mm H3K27me3:HBG3:2+1 day vs H3:HBG3:2+1 day", "median linefit, quantile norm" },
                { "Mm H3K27me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm" },
                { "Mm H3K27me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm" }
        }; 
        
        h3k79Versions = new String[][] { 
                { "Mm H3K79me2:HBG3:mES vs H3:HBG3:mES", "median linefit, quantile norm 6tp" },
                { "Mm H3K79me2:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage", "median linefit, quantile norm 6tp" },
                { "Mm H3K79me2:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA", "median linefit, quantile norm 6tp" },
                { "Mm H3K79me2:HBG3:2+1 day vs WCE:HBG3:2+1 day", "median linefit, quantile norm 6tp" },
                { "Mm H3K79me2:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage", "median linefit, quantile norm 6tp" },
                { "Mm H3K79me2:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage", "median linefit, quantile norm 6tp" }
        };
        ring1bVersions = new String[][] {
        		{"NONE","NONE"},
        		{"Mm Ring1b:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage", "median linefit, quantile norm"},	
        		{"Mm Ring1b:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA","median linefit, quantile norm"},	
        		{"Mm Ring1b:HBG3:2+1 day vs WCE:HBG3:2+1 day", "median linefit, quantile norm"},
        		{"NONE","NONE"},
        		{"NONE","NONE"}
        };
        suz12Versions = new String[][] {
        		{"NONE","NONE"},
        		{"Mm Suz12:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage","median linefit"},	
        		{"Mm Suz12:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA","median linefit"},	
        		{"Mm Suz12:HBG3:2+1 day vs WCE:HBG3:2+1 day","median linefit", "2 (7/24/08)"},
        		{"NONE","NONE"},
        		{"NONE","NONE"}
        };
        rarChipSeq = new String[][]{
        		{"PPG_Solexa_RAR_ES+2d", "ELAND_unique"},
        		{"PPG_Solexa_RAR_8hr", "ELAND_unique"}        		
        };
        backChipSeq = new String[][]{
        		{"PPG_Solexa_WCE_2+1", "ELAND_unique"},
        		{"PPG_Solexa_WCE_ES+2d", "ELAND_unique"}
        };
	}
	public static ArrayList<Integer> allowedChips= new ArrayList<Integer>();
	
}
