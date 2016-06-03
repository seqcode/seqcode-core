package org.seqcode.projects.shaun.viz;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import org.seqcode.gse.datasets.motifs.WeightMatrix;
import org.seqcode.gse.datasets.seqdata.SeqLocator;
import org.seqcode.gse.utils.Pair;
import org.seqcode.gse.utils.io.BackgroundModelIO;
import org.seqcode.projects.multigps.utilities.Utils;
import org.seqcode.projects.shaun.FreqMatrixImport;


public class FigureOptions {

	public Genome genome;
	public boolean reverseOrder=false;
	public boolean defaultFilledColumns=false;
	public boolean drawExptLabels=true;
	public boolean drawGeneLabels=true;
	public boolean drawTrackNames=false;
	public boolean useDBGenes=true;
	public int labelFontSize=20;
	public int geneFontSize=14;
	public int fontSize=12;
	public int screenSizeX=1000, screenSizeY=900;
	public int geneHeight=15;
	public int boxHeight=30, boxWidth = 18, siteWidth=10, motifWidth=5, motifHeight=10;
	public int defaultExptTrackHeight=120;
	public int readExt = 100;
	public Color geneColor = new Color(232,232,232);
	public Color specialGeneColor = new Color(232,232,232);
	public String specialGeneName = "";
	public Color motifColor = new Color(51,102,0);
	public Color loopColor = Color.lightGray;
	public Color interColor = Color.red;
	public Color defaultBgColor = Color.gray;
	public Color defaultPeakColor = Color.red;
	public int exptBgThick = 1;
	public int exptPeakThick=2;
	public int exptPeakWin=-1;
    public Color diffExptPosColor = new Color(100,0,0);
    public Color diffExptNegColor = new Color(0,100,0);
	public int topBorder=50, bottomBorder=50;
	public int leftBorder=25, rightBorder=25;
	public Region gRegion=null;
	public File transcriptGTF=null;
	public HashMap<String, SeqExperiment> experiments = new HashMap<String, SeqExperiment>();
	public ArrayList<String> exptNames=new ArrayList<String>();
	private HashMap<String, WeightMatrix> motifs = new HashMap<String, WeightMatrix>();
	private MarkovBackgroundModel motifBack;

	public FigureOptions(Genome g){
		genome = g;		            	
	}

	public void loadOptions(File f){
		try {
			if(f.isFile()){
				BufferedReader reader = new BufferedReader(new FileReader(f));
				String line;
				while((line= reader.readLine())!=null){
					String [] tokens = line.split("\\t");

					if(tokens.length>1 && !tokens[0].startsWith("#")){
						String oName = tokens[0];
						String oVal = tokens[1];
						
						//Colors
						if(tokens.length>3){
							if(oName.equals("geneColor")){
								geneColor = new Color(new Integer(tokens[1]), new Integer(tokens[2]), new Integer(tokens[3]));
							}if(oName.equals("specialGeneColor")){
								specialGeneColor = new Color(new Integer(tokens[2]), new Integer(tokens[3]), new Integer(tokens[4]));
								specialGeneName = tokens[1];
							}if(oName.equals("loopColor")){
								loopColor = new Color(new Integer(tokens[1]), new Integer(tokens[2]), new Integer(tokens[3]));
								if(tokens.length>4){
									loopColor = new Color(new Integer(tokens[1]), new Integer(tokens[2]), new Integer(tokens[3]), new Integer(tokens[4]));
								}
							}if(oName.equals("interColor")){
								interColor = new Color(new Integer(tokens[1]), new Integer(tokens[2]), new Integer(tokens[3]));
								if(tokens.length>4){
									interColor = new Color(new Integer(tokens[1]), new Integer(tokens[2]), new Integer(tokens[3]), new Integer(tokens[4]));
								}
							}if(oName.equals("motifColor")){
								motifColor = new Color(new Integer(tokens[1]), new Integer(tokens[2]), new Integer(tokens[3]));
							}if(oName.equals("exptBgColor") && tokens.length==4){
								defaultBgColor = new Color(new Integer(tokens[1]), new Integer(tokens[2]), new Integer(tokens[3]));
							}if(oName.equals("exptPeakColor") && tokens.length==4){
								defaultPeakColor = new Color(new Integer(tokens[1]), new Integer(tokens[2]), new Integer(tokens[3]));
							}if(oName.equals("diffExptPosColor") && tokens.length==4){
								diffExptPosColor = new Color(new Integer(tokens[1]), new Integer(tokens[2]), new Integer(tokens[3]));
							}if(oName.equals("diffExptNegColor") && tokens.length==4){
								diffExptNegColor = new Color(new Integer(tokens[1]), new Integer(tokens[2]), new Integer(tokens[3]));
							}
						}
						
						if(oName.equals("gRegion")){
							String[] gs1 =oVal.split(":");
							String[] gs2 =gs1[1].split("-");
							gRegion = new Region(genome, gs1[0], new Integer(gs2[0]), new Integer(gs2[1]));
						}if(oName.equals("reverse")){
							if(oVal.equals("true")){reverseOrder=true;
							}else{reverseOrder=false;}
						}if(oName.equals("filledColumns") && tokens.length==2){
							if(oVal.equals("true")){defaultFilledColumns=true;
							}else{defaultFilledColumns=false;}
						}if(oName.equals("drawExptLabels")){
							if(oVal.equals("true")){drawExptLabels=true;
							}else{drawExptLabels=false;}
						}if(oName.equals("drawGeneLabels")){
							if(oVal.equals("true")){drawGeneLabels=true;
							}else{drawGeneLabels=false;}
						}if(oName.equals("drawTrackNames")){
							if(oVal.equals("true")){drawTrackNames=true;
							}else{drawTrackNames=false;}
						}if(oName.equals("screenSizeX")){
							screenSizeX = new Integer(oVal);
						}if(oName.equals("screenSizeY")){
							screenSizeY = new Integer(oVal);
						}if(oName.equals("fontSize")){
							fontSize = new Integer(oVal);
						}if(oName.equals("labelFontSize")){
							labelFontSize = new Integer(oVal);
						}if(oName.equals("geneFontSize")){
							geneFontSize = new Integer(oVal);
						}if(oName.equals("geneHeight")){
							geneHeight = new Integer(oVal);
						}if(oName.equals("boxHeight")){
							boxHeight = new Integer(oVal);
						}if(oName.equals("boxWidth")){
							boxWidth = new Integer(oVal);
						}if(oName.equals("siteWidth")){
							siteWidth = new Integer(oVal);
						}if(oName.equals("motifWidth")){
							motifWidth = new Integer(oVal);
						}if(oName.equals("motifHeight")){
							motifHeight = new Integer(oVal);
						}if(oName.equals("exptTrackHeight") && tokens.length==2){
							defaultExptTrackHeight = new Integer(oVal);
						}if(oName.equals("exptBgThick")){
							exptBgThick = new Integer(oVal);
						}if(oName.equals("exptPeakThick")){
							exptPeakThick = new Integer(oVal);
						}if(oName.equals("exptPeakWin")){ //Can be -1 for data-sourced region
							exptPeakWin = new Integer(oVal);
						}if(oName.equals("topBorder")){
							topBorder = new Integer(oVal);
						}if(oName.equals("bottomBorder")){
							bottomBorder = new Integer(oVal);
						}if(oName.equals("rightBorder")){
							rightBorder = new Integer(oVal);
						}if(oName.equals("leftBorder")){
							leftBorder = new Integer(oVal);
						}if(oName.equals("readExt")){
							readExt = new Integer(oVal);
						}if(oName.equals("motifs")){
							loadMotifsFromFile(oVal);
						}if(oName.equals("motifBack")){
							loadBackgroundFromFile(oVal);
						}if(oName.equals("dbGenes")){
							if(oVal.equals("true")){useDBGenes=true;
							}else{useDBGenes=false;}
						}if(oName.equals("gtfFile")){
							transcriptGTF = new File(oVal);
						}if(oName.equals("rnaExpt")){
							if(!experiments.containsKey(oVal)){
								experiments.put(oVal, new SeqExperiment());
								System.err.println("RNA-seq Experiment: "+oVal);
								exptNames.add(oVal);
							}
						}


						//Experiments
						if(tokens.length>2){
							if(oName.equals("exptLoc") || oName.equals("diffExptLoc")){
                                if(!experiments.containsKey(oVal)){
                                     experiments.put(oVal, new SeqExperiment());
                                     System.err.println("Experiment: "+oVal);
                                }
                                if(oName.equals("diffExptLoc"))
                                	 experiments.get(oVal).isDiff=true;
                                int numExptTokens = oName.equals("exptLoc") ? 1 : 2;
                                for(int t=0; t<numExptTokens; t++){
                                     String[] pieces = tokens[2+t].split(";");
                                     Set<String> replist = new TreeSet<String>();

                                     if(pieces.length == 2){
                                          if(t==0){
                                        	  experiments.get(oVal).expt = new SeqLocator(pieces[0], replist, pieces[1]);
                                        	  exptNames.add(oVal);
                                          }else{
                                               experiments.get(oVal).baseExpt = new SeqLocator(pieces[0], replist, pieces[1]);
                                          }
                                     }else if (pieces.length == 3){
                                          System.err.println(pieces[1]);
                                          String[] reps = pieces[1].split(",");
                                          if (reps.length > 1) {
                                               for (int i=0; i<reps.length; i++) {
                                                    replist.add(reps[i]);
                                               }
                                               if(t==0){
                                                    experiments.get(oVal).expt = new SeqLocator(pieces[0],replist,pieces[2]);
                                                    exptNames.add(oVal);
                                               }else{
                                                    experiments.get(oVal).baseExpt = new SeqLocator(pieces[0],replist,pieces[2]);
                                               }
                                          } else {
                                        	  replist.add(pieces[1]);
                                               if(t==0){
                                                    experiments.get(oVal).expt = new SeqLocator(pieces[0], replist, pieces[2]);
                                                    exptNames.add(oVal);
                                               }else{
                                                    experiments.get(oVal).baseExpt = new SeqLocator(pieces[0], replist, pieces[2]);
                                               }
                                          }
                                     }else{
                                          System.err.println("Couldn't parse a SeqLocator from " + tokens[2]);
                                     }
                                }
                           }
							//Hack option to allow loading of old ChIP-chip data files
							if(oName.equals("exptDataFile") || oName.equals("diffExptDataFile")){
                                if(!experiments.containsKey(oVal)){
                                     experiments.put(oVal, new SeqExperiment());
                                     System.err.println("Experiment: "+oVal);
                                }
                                if(oName.equals("diffExptDataFile"))
                                	experiments.get(oVal).isDiff=true;
                                exptNames.add(oVal);
                                experiments.get(oVal).preFormattedDataFile = new File(tokens[2]);
							}
							
							if(oName.equals("exptYMax")){
								experiments.get(oVal).yMax = new Integer(tokens[2]);
							}if(oName.equals("exptPBMax")){
								experiments.get(oVal).pbMax = new Integer(tokens[2]);
							}if(oName.equals("exptPeaks")){
								List<Region> regs = Utils.loadRegionsFromFile(tokens[2], genome,  exptPeakWin);
								for(Region r : regs){
									if(gRegion.overlaps(r)){
										experiments.get(oVal).peaks.add(r);
									}
								}
							}if(oName.equals("exptInters")){//ccr
								List<Pair<Point,Point>> inters = Utils.loadIntersFromFile(tokens[2], genome);
								for (Pair<Point,Point> inter : inters) {
									if (gRegion.contains(inter.car()) && gRegion.contains(inter.cdr())) {
										experiments.get(oVal).inters.add(inter);
									}
								}
							}if(oName.equals("exptMotif")){
								if(motifs.containsKey(tokens[2])){
									experiments.get(oVal).motif = motifs.get(tokens[2]);
								}else{
									System.err.println("Cannot find motif: " + tokens[2]);
								}
							}if(oName.equals("exptMotifThres")){
								experiments.get(oVal).motifThres = new Double(tokens[2]);
							}
							if(oName.equals("exptBgColor") && tokens.length==5){
								experiments.get(oVal).exptBgColor = new Color(new Integer(tokens[2]), new Integer(tokens[3]), new Integer(tokens[4]));
							}if(oName.equals("exptPeakColor") && tokens.length==5){
								experiments.get(oVal).exptBgColor = new Color(new Integer(tokens[2]), new Integer(tokens[3]), new Integer(tokens[4]));
							}
							if(oName.equals("filledColumns") && tokens.length==3){
								if(oVal.equals("true")){experiments.get(oVal).filledColumns=true;
								}else{experiments.get(oVal).filledColumns=false;}
							}
							if(oName.equals("exptTrackHeight") && tokens.length==3){
								experiments.get(oVal).exptTrackHeight = new Integer(tokens[2]);
							}
							//Differential binding stuff
							if(oName.equals("scalingFactor") && tokens.length==3){
								experiments.get(oVal).scaling = new Double(tokens[2]);
							}if(oName.equals("diffWinWidth") && tokens.length==3){
								experiments.get(oVal).diffWinWidth = new Integer(tokens[2]);
							}if(oName.equals("diffWinStep") && tokens.length==3){
								experiments.get(oVal).diffWinStep = new Integer(tokens[2]);
							}
							
							//RNA-seq stuff
							if(oName.equals("junctionsFile")){
								experiments.get(oVal).junctionsFile = new File(tokens[2]);
							}if(oName.equals("samFile")){
								experiments.get(oVal).samFile = new File(tokens[2]);
							}if(oName.equals("readDepth")){
								if(tokens[2].equals("true")){
									experiments.get(oVal).readDepth=true;
								}else{
									experiments.get(oVal).readDepth=false;
								}
							}if(oName.equals("pairedReads")){
								if(tokens[2].equals("true")){
									experiments.get(oVal).pairedReads=true;
								}else{
									experiments.get(oVal).pairedReads=false;
								}
							}
						}
					}
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

	//Load freq matrices
	public void loadMotifsFromFile(String filename){
		FreqMatrixImport motifImport = new FreqMatrixImport();
		motifImport.setBackground(motifBack);
		for(WeightMatrix wm : motifImport.readTransfacMatrices(filename)){
			motifs.put(wm.name, wm);
		}
	}
	//Load background model
	public void loadBackgroundFromFile(String backFile){		
		try {
			motifBack = BackgroundModelIO.parseMarkovBackgroundModel(backFile, genome);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	protected class SeqExperiment{
		SeqLocator expt = null;
		int exptTrackHeight = defaultExptTrackHeight;
		int yMax=100;
		int pbMax=100;
		boolean filledColumns = defaultFilledColumns;
		//Colors
		Color exptBgColor = defaultBgColor;
		Color exptPeakColor = defaultPeakColor;
		//ChIP-seq 
		ArrayList<Region> peaks = new ArrayList<Region>();
		WeightMatrix motif = null;
		double motifThres = 0;
		//ChIA-PET
		List<Pair<Point,Point>> inters = new ArrayList<Pair<Point,Point>>();
		//RNA-seq
		File samFile=null;
		File junctionsFile = null;
		boolean readDepth=false;
		boolean pairedReads=false;
        //Differential experiments
		boolean isDiff=false;
        double scaling=1;
        SeqLocator baseExpt = null;
    	File preFormattedDataFile = null;
        int diffWinWidth=500;
        int diffWinStep=100;
	}
}

