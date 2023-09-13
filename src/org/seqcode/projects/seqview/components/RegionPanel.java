package org.seqcode.projects.seqview.components;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

import java.util.*;
import java.util.regex.Pattern;
import java.util.regex.Matcher;


/* this is all for saveImage() */
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.dom.GenericDOMImplementation;
import org.seqcode.data.connections.DatabaseException;
import org.seqcode.data.motifdb.*;
import org.seqcode.data.seqdata.*;
import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Gene;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.NamedStrandedRegion;
import org.seqcode.genome.location.NamedTypedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.ScoredRegion;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.gsebricks.*;
import org.seqcode.gsebricks.verbs.*;
import org.seqcode.gsebricks.verbs.location.RefGeneGenerator;
import org.seqcode.gsebricks.verbs.motifs.PerBaseMotifMatch;
import org.seqcode.gseutils.*;
import org.seqcode.gseutils.models.Model;
import org.seqcode.projects.seqview.*;
import org.seqcode.projects.seqview.model.*;
import org.seqcode.projects.seqview.paintable.*;
import org.w3c.dom.Document;
import org.w3c.dom.DOMImplementation;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Writer;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;

/**
 *  RegionFrame encapsulates a set of painters that display a particular genomic region
 */
public class RegionPanel extends JPanel 
implements ActionListener, KeyListener, 
Listener<EventObject>, PainterContainer, MouseListener {

	//status lives in the parent frame
	private SeqViewStatus status;
	// controls at the bottom of the panel
	private JPanel buttonPanel;
	private JScrollPane scrollPane=null;
	private RegionContentPanel mainPanel;
	// scrolling controls
	private JButton leftButton, rightButton, zoomInButton, zoomOutButton, farLeftButton, farRightButton, refreshButton;
	private JTextField locationField;
	// maps a track name to its coordinates in the current layout
	private Hashtable<String,Integer> ulx, uly, lrx, lry;
	// maps a track name to the set of painters in that track
	private Hashtable<String,ArrayList<RegionPaintable>> painters;
	// keeps track of the order in which painters are to be drawn
	private ArrayList<RegionPaintable> allPainters= new ArrayList<RegionPaintable>();
	private Hashtable<String,Integer> trackPaintOrder;
	private Hashtable<String,Double> trackSizeFactor;
	// set of all the models and painters
	private HashSet<RegionModel> allModels;
	private Genome genome;
	private boolean closed=false;
	// painterCount is the total number of painters.
	// readyCount is reset to 0 when the current region changes and keep track
	// of how many painters have reported in as ready.  When readyCount ==
	// painterCount, we can try to draw.
	private int readyCount, painterCount;
	private Region currentRegion;
	private SeqViewOptions currentOptions;

	// mgp is here because there's at most one instance per RegionPanel
	// and this is an easy way to keep track of it.
	//private MultiGenePainter mgp = null;
	private ExonGenePainter egp = null;
	private CDSGenePainter cgp = null;
	private static Color transparentWhite = new Color(255,255,255,100);
	private static boolean paintedTransparent=false, firstPaint=true;
	private boolean forceupdate = false, firstconfig = true;
	private File currDirectory = new File(System.getProperty("user.home"));
	private Hashtable<RegionPaintable, ArrayList<RegionModel>> painterModelMap = new Hashtable<RegionPaintable, ArrayList<RegionModel>>();
	
	public RegionPanel(Genome g, SeqViewStatus s, File currDir) {
		super();
		status = s;
		currDirectory = currDir;
		init(g);        
		currentOptions = new SeqViewOptions(g);
		addPaintersFromOpts(currentOptions);
		setVisible(true);
	}

	public RegionPanel(SeqViewOptions opts, SeqViewStatus s, File currDir) {
		super();
		status = s;
		currDirectory = currDir;
		Genome g = opts.getGenome();
		init(g);
		currentOptions = opts;
		addPaintersFromOpts(currentOptions);
		setVisible(true);
		//Find our initial region.
		Region startingRegion = null;
		if (opts.gene != null && opts.gene.matches("...*")) {
			startingRegion = regionFromString(genome,opts.gene);
		}
		if (startingRegion != null) {
			setRegion(startingRegion);
		} else if (opts.start >= 0 && opts.chrom != null) {
			setRegion(new Region(g,opts.chrom,opts.start,opts.stop));
		} else if (opts.position != null && opts.position.length() > 0) {
			Region r = regionFromString(genome,opts.position);
			if (r != null) {
				setRegion(r);
			} else {
				r = regionFromString(genome,opts.gene);
				if (r != null) {
					setRegion(r);
				} else {
					throw new NullPointerException("Need a valid starting position in either chrom or gene");
				}
			}
		} else {
			throw new NullPointerException("Need a starting position in either chrom or gene");
		}
		if (opts.regionListFile != null) {
			java.util.List<Region> regions = readRegionsFromFile(g,opts.regionListFile,false);
			RegionListPanel p = new RegionListPanel(this,
					regions);
			RegionListPanel.makeFrame(p);
		}
	}

	public void handleWindowClosing() { 
		for(RegionPaintable rp : allPainters) { 
			rp.cleanup();
		}
		synchronized(allModels){
			for (RegionModel m : allModels) {
				synchronized(m) {
					m.stopRunning();
			   	}
	 	   	}
		}
		closed=true;
	}

	public void init(Genome g) {
		allModels = new HashSet<RegionModel>();
		painters = new Hashtable<String,ArrayList<RegionPaintable>>();

		genome = g;
		trackPaintOrder = new Hashtable<String,Integer>();

		trackSizeFactor = new Hashtable<String,Double>();
		allPainters = new ArrayList<RegionPaintable>();
		ulx = new Hashtable<String,Integer>();
		uly = new Hashtable<String,Integer>();
		lrx = new Hashtable<String,Integer>();
		lry = new Hashtable<String,Integer>();
		painterCount = 0;
		readyCount = 0;

		currentRegion = new Region(g,"1",1,1000);

		buttonPanel = new JPanel();      
		buttonPanel.setLayout(new GridBagLayout());  
		leftButton = new JButton("<-");
		leftButton.setToolTipText("step left");
		rightButton = new JButton("->");
		rightButton.setToolTipText("step right");
		zoomInButton = new JButton("++");
		zoomInButton.setToolTipText("zoom in");
		zoomOutButton = new JButton("--");        
		zoomOutButton.setToolTipText("zoom out");
		farLeftButton = new JButton("<<<-");
		farLeftButton.setToolTipText("jump left");
		farRightButton = new JButton("->>>");
		farRightButton.setToolTipText("jump right");
		refreshButton = new JButton("R");
		refreshButton.setToolTipText("refresh");
		locationField = new JTextField();
		Dimension buttonSize = new Dimension(30,20);
		leftButton.setMaximumSize(buttonSize);
		rightButton.setMaximumSize(buttonSize);
		zoomInButton.setMaximumSize(buttonSize);
		zoomOutButton.setMaximumSize(buttonSize);
		farLeftButton.setMaximumSize(buttonSize);
		farRightButton.setMaximumSize(buttonSize);
		refreshButton.setMaximumSize(buttonSize);
		locationField.setMinimumSize(new Dimension(160,20));
		locationField.setPreferredSize(new Dimension(300,20));

		buttonPanel.add(farLeftButton);
		buttonPanel.add(leftButton);
		buttonPanel.add(zoomOutButton);
		buttonPanel.add(locationField);
		buttonPanel.add(zoomInButton);
		buttonPanel.add(rightButton);
		buttonPanel.add(farRightButton);
		buttonPanel.add(refreshButton);

		mainPanel = new RegionContentPanel();
		mainPanel.addMouseListener(this);
		
		scrollPane = new JScrollPane();
		scrollPane.setViewportView(mainPanel);
		scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
		scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		scrollPane.getVerticalScrollBar().setUnitIncrement(16);
		
		mainPanel.setSize(new Dimension(getWidth(),getHeight()-buttonPanel.getHeight()));
		mainPanel.setPreferredSize(new Dimension(getWidth(),getHeight()-buttonPanel.getHeight()));
		
		leftButton.addActionListener(this);
		rightButton.addActionListener(this);
		locationField.addActionListener(this);
		zoomInButton.addActionListener(this);
		zoomOutButton.addActionListener(this);     
		farLeftButton.addActionListener(this);
		farRightButton.addActionListener(this);
		refreshButton.addActionListener(this);
		buttonPanel.addKeyListener(this);
		mainPanel.addKeyListener(this);    
		scrollPane.addKeyListener(this);
		
		setLayout(new BorderLayout());
		setBackground(Color.WHITE);
		add(scrollPane,BorderLayout.CENTER);
		add(buttonPanel, BorderLayout.SOUTH);	
	}
	
	/**
	 * reinit is like the constructor, but assumes the superclass has already been called,
	 * and performs some cleanup first.  
	 */
	public void reinit(SeqViewOptions opts) {
		//Cleanup
		synchronized(allModels) {
			for (RegionModel m : allModels) {
				synchronized(m){
					m.stopRunning();
					m.notifyAll();
				}
 		   	}
 	   	}
		try {
			Thread.sleep(400);
 	   	} catch (Exception e) { }
		//Re-initialize (but don't make all the panels again)
		genome = opts.getGenome();
		allModels = new HashSet<RegionModel>();
		painters = new Hashtable<String,ArrayList<RegionPaintable>>();
		trackPaintOrder = new Hashtable<String,Integer>();
		trackSizeFactor = new Hashtable<String,Double>();
		allPainters = new ArrayList<RegionPaintable>();
		ulx = new Hashtable<String,Integer>();
		uly = new Hashtable<String,Integer>();
		lrx = new Hashtable<String,Integer>();
		lry = new Hashtable<String,Integer>();
		painterCount = 0;
		readyCount = 0;
		currentOptions = opts;
		
		//Set our initial region.
		Region startingRegion = null;
		if (opts.gene != null && opts.gene.matches("...*")) {
			startingRegion = regionFromString(genome,opts.gene);
		}
		if (startingRegion != null) {
			setRegion(startingRegion);
		} else if (opts.start >= 0 && opts.chrom != null) {
			setRegion(new Region(genome,opts.chrom,opts.start,opts.stop));
		} else if (opts.position != null && opts.position.length() > 0) {
			Region r = regionFromString(genome,opts.position);
			if (r != null) {
				setRegion(r);
			} else {
				r = regionFromString(genome,opts.gene);
				if (r != null) {
					setRegion(r);
				} else {
					throw new NullPointerException("Need a valid starting position in either chrom or gene");
				}
			}
		} else {
			throw new NullPointerException("Need a starting position in either chrom or gene");
		}
		
		//Add painters
		addPaintersFromOpts(opts);
		setVisible(true);
		
		//Load region list if available
		if (opts.regionListFile != null) {
			java.util.List<Region> regions = readRegionsFromFile(genome,opts.regionListFile,false);
			RegionListPanel p = new RegionListPanel(this,
					regions);
			RegionListPanel.makeFrame(p);
		}
	}

	public void addPaintersFromOpts(SeqViewOptions opts) {        

		if (!(opts.getGenome().equals(genome))){
			System.err.println("Should be the same");
		}
		opts.mergeInto(currentOptions);

		if (opts.hash) {
			HashMarkPaintable p = new HashMarkPaintable();
			p.setLabel("Chromosomal position");
			p.addEventListener(this);
			addPainter(p);
		}
		RegionMapperModel seqmodel = null;
		if (opts.gccontent || opts.cpg || opts.seqletters || opts.regexmatcher || opts.pyrpurcontent) {
			seqmodel = new RegionMapperModel(new SequenceGenerator(genome, opts.maxSequenceQuery));
			addModel(seqmodel);
			Thread t = new Thread(seqmodel);
			t.start();
		}

		if (opts.gccontent) {
			GCContentPainter p = new GCContentPainter(seqmodel);
			p.setLabel("GC content");
			p.addEventListener(this);
			p.setOption(SeqViewOptions.GCCONTENT,null);
			addPainter(p);
			addModelToPaintable(p,seqmodel);
		}
		if (opts.pyrpurcontent) {
			GCContentPainter p = new GCContentPainter(seqmodel);
			p.setLabel("Pyr (red) Pur (blue)");
			p.addEventListener(this);
			p.setOption(SeqViewOptions.GCCONTENT,null);
			GCContentProperties props = p.getProperties();
			props.BlueBases = "AG";
			props.RedBases = "CT";
			addPainter(p);
			addModelToPaintable(p,seqmodel);
		}        
		if (opts.cpg) {
			CpGPainter p = new CpGPainter(seqmodel);
			p.setLabel("CpG");
			p.addEventListener(this);
			p.setOption(SeqViewOptions.CPG,null);
			addPainter(p);
			addModelToPaintable(p,seqmodel);
		}
		if (opts.regexmatcher) {
			RegexMatchPainter p = new RegexMatchPainter(seqmodel);
			p.setLabel("regexes");
			p.addEventListener(this);
			p.setOption(SeqViewOptions.REGEXMATCHER,null);
			addPainter(p);
			addModelToPaintable(p,seqmodel);
			for (String r : opts.regexes.keySet()) {
				p.addRegex(r,opts.regexes.get(r));
			}
		}

		if (opts.seqletters) {
			BasePairPainter p = new BasePairPainter(seqmodel);
			p.setLabel("Sequence");
			p.addEventListener(this);
			p.setOption(SeqViewOptions.SEQLETTERS,null);
			addPainter(p);
			addModelToPaintable(p,seqmodel);
		}

		//Loading experiment painters
		if (opts.seqExpts.size() > 0) {
			try {
				SeqDataLoader loader = new SeqDataLoader();

				for(int i = 0; i < opts.seqExpts.size(); i++) { 
					Collection<SeqAlignment> alignments = opts.seqExpts.get(i).aligns;

					boolean isChIP=true;
					boolean allPaired = true;
			    	for(SeqAlignment a : alignments){
			    		allPaired = allPaired && (a.getAlignType().getName().equals("PAIRED") || a.getAlignType().getName().equals("MIXED"));
			    		isChIP = isChIP && (a.getExpt().getExptType().equals("CHIPSEQ") || a.getExpt().getExptType().equals("CHIPEXO")); 
			    	}
			    	
					RegionModel histomod, arcmod=null;
					RegionPaintable p;
					if (opts.seqHistogramPainter) {
						histomod = new SeqHistogramModel(alignments);
						if(allPaired){
							arcmod = new InteractionArcModel(alignments);
							p = new SeqHistogramPainter((SeqHistogramModel)histomod, (InteractionArcModel)arcmod, !isChIP);
						}else{
							p = new SeqHistogramPainter((SeqHistogramModel)histomod, null, false);
						}
					} else {
						System.err.println("Using old ChipSeq painters");
						histomod = new SeqDataModel(new org.seqcode.data.readdb.Client(),
								alignments);
						p = new SeqAboveBelowStrandPainter((SeqDataModel)histomod);
					}
					addModel(histomod);
					Thread t1 = new Thread((Runnable)histomod); t1.start();
					if(arcmod!=null){
						addModel(arcmod);
						Thread t2 = new Thread((Runnable)arcmod); t2.start();
					}
					p.setLabel(opts.seqExpts.get(i).toString());

					p.addEventListener(this);
					p.setOption(SeqViewOptions.SEQDATA,opts.seqExpts.get(i));
					addPainter(p);
					addModelToPaintable(p,histomod);
					if(arcmod!=null){
						addModelToPaintable(p,arcmod);
					}
				}
				loader.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		/*
		if (opts.pairedChipseqExpts.size() > 0) {
			try {
				SeqDataLoader loader = new SeqDataLoader(true);
				for(int i = 0; i < opts.pairedChipseqExpts.size(); i++) { 

					Collection<SeqAlignment> alignments = loader.loadAlignments(opts.pairedChipseqExpts.get(i), genome);
					PairedEndModel m = new PairedEndModel(alignments);
					PairedEndPainter p = new PairedEndPainter(m);
					addModel(m);
					Thread t = new Thread((Runnable)m); t.start();
					p.setLabel("Paired " + opts.pairedChipseqExpts.get(i).toString());

					p.addEventListener(this);
					addPainter(p);
					addModelToPaintable(p,m);
				}
				loader.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		if (opts.chipseqAnalyses.size() > 0) {
			try {
				SeqDataLoader loader = new SeqDataLoader(true);
				for (int i = 0; i < opts.chipseqAnalyses.size(); i++) {
					SeqAnalysis a = opts.chipseqAnalyses.get(i);
					ChipSeqAnalysisModel m = new ChipSeqAnalysisModel(a);
					SeqAnalysisPainter p = new SeqAnalysisPainter(a,m);
					addModel(m);
					Thread t = new Thread((Runnable)m); t.start();
					p.setLabel(a.toString());

					p.addEventListener(this);
					addPainter(p);
					addModelToPaintable(p,m);
				}
				loader.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}*/

		RegionExpanderFactoryLoader<Gene> gfLoader;
		RegionExpanderFactoryLoader<NamedTypedRegion> annotLoader;
		gfLoader = new RegionExpanderFactoryLoader<Gene>("gene");
		annotLoader = new RegionExpanderFactoryLoader<NamedTypedRegion>("annots");

		if(opts.genes.size() > 0 && (egp == null || cgp == null)) {
			int gCount=0, cCount=0;
			GeneModel geneModel = new GeneModel();
			GeneModel cdsModel = new GeneModel();
			
			for(int i = 0; i < opts.genes.size(); i++) {
				RegionExpanderFactory<Gene> genefactory = gfLoader.getFactory(genome,
						opts.genes.get(i).toString());
				Expander<Region,Gene> expander = genefactory.getExpander(genome);
				if (genefactory.getProduct().equals("Gene")) {
					geneModel.addExpander(expander); gCount++;
				}else if (genefactory.getProduct().equals("CDS")) {
					cdsModel.addExpander(expander); cCount++;
				}
			}        

			if(gCount>0){
				addModel(geneModel);
				Thread t = new Thread(geneModel); t.start();
				egp = new ExonGenePainter(geneModel);
				egp.setLabel("genes");
				egp.addEventListener(this);
				addPainter(egp);
				addModelToPaintable(egp, geneModel);
			}
			if(cCount>0){
				addModel(cdsModel);
				Thread t = new Thread(cdsModel); t.start();
				cgp = new CDSGenePainter(cdsModel);
				cgp.setLabel("CDS");
				cgp.addEventListener(this);
				addPainter(cgp);
				addModelToPaintable(cgp, cdsModel);
			}
		}

		for (int i = 0; i < opts.otherannots.size(); i++) {

			RegionExpanderFactory factory = annotLoader.getFactory(genome,
					opts.otherannots.get(i).toString());
			Expander expander = factory.getExpander(genome);
			RegionExpanderModel m;
			RegionPaintable p;
			if (factory.getProduct().equals("NamedTypedRegion")) {                
				m = new RegionExpanderModel<NamedTypedRegion>((Expander<Region,NamedTypedRegion>)expander);
				p = new NamedTypedPainter((RegionExpanderModel<NamedTypedRegion>)m);
			} else if (factory.getProduct().equals("ScoredRegion")) {
				double max = 1.0;
				if (opts.otherannots.get(i).toString().equals("readcoverage")) {
					max = 40;
				}
				m = new RegionExpanderModel<ScoredRegion>((Expander<Region,ScoredRegion>)expander);
				p = new HeightScoredPainter((RegionExpanderModel<ScoredRegion>)m,max);
			} else if (factory.getProduct().equals("NamedStrandedRegion")) {
				m = new RegionExpanderModel<NamedStrandedRegion>((Expander<Region,NamedStrandedRegion>)expander);
				p = new NamedStrandedPainter((RegionExpanderModel<NamedStrandedRegion>)m);
			} else if (factory.getProduct().equals("StrandedRegion")) {
				m = new RegionExpanderModel<StrandedRegion>((Expander<Region,StrandedRegion>)expander);
				p = new NamedStrandedPainter((RegionExpanderModel<StrandedRegion>)m);
			} else if (factory.getProduct().equals("NamedRegion")) {
				m = new RegionExpanderModel<NamedRegion>((Expander<Region,NamedRegion>)expander);
				p = new NamedStrandedPainter((RegionExpanderModel<NamedRegion>)m); // yes, this works.  NamedStrandedPainter can handle non-stranded Regions
			} else {
				throw new RuntimeException("Don't understand product type " + factory.getProduct());
			}

			addModel(m);
			Thread t = new Thread(m); t.start();               
			p.setLabel(opts.otherannots.get(i).toString());
			p.addEventListener(this);
			p.setOption(SeqViewOptions.OTHERANNOTS,opts.otherannots.get(i));
			addPainter(p);
			addModelToPaintable(p, m);
		}

		for (int i = 0; i < opts.motifs.size(); i++) {
			WeightMatrix matrix = opts.motifs.get(i);
			MarkovBackgroundModel bgModel = null;
			try {
				String bgmodelname = "whole genome zero order";
				BackgroundModelMetadata md = BackgroundModelLoader.getBackgroundModel(bgmodelname,
						1,
						"MARKOV",
						genome.getDBID());
				if (md != null) {
					bgModel = BackgroundModelLoader.getMarkovModel(md);
				} else {
					System.err.println("Couldn't get metadata for " + bgmodelname);
				}
			} catch (Exception e) {
				e.printStackTrace();
			}

			if (bgModel != null) {
				matrix.toLogOdds(bgModel);
			} else {
				matrix.toLogOdds();
			}

			PerBaseMotifMatch match = new PerBaseMotifMatch(matrix);
			RegionMapperModel<Double[]> m = new RegionMapperModel<Double[]>(new Mapper.Compose<Region,String,Double[]>(new SequenceGenerator(genome),
					match));
			addModel(m);
			Thread t = new Thread(m);
			t.start();
			PerBaseScorePainter p = new PerBaseScorePainter<Double>(m,matrix.getMinScore(),0.0,matrix.getMaxScore());
			p.setLabel(opts.motifs.get(i).toString());
			p.addEventListener(this);
			p.setOption(SeqViewOptions.MOTIFS,opts.motifs.get(i));            
			addPainter(p);
			addModelToPaintable(p, m);
		}

		
		for (String k : opts.regionTracks.keySet()) {
			if(k.endsWith(".gff"))
				addTrackFromFile(k,opts.regionTracks.get(k),true);
			else
				addTrackFromFile(k,opts.regionTracks.get(k),false);
		}

		if (firstconfig) {
			firstconfig = false;
		} else {
			forceupdate = true;
			setRegion(getRegion());
		}
	}

	/* when we take a paintable out of the display, we also
       need to take it out of currentOptions.  Failing to
       remove the paintable from currentOptions makes it
       impossible to add the paintable back into the display 
       later (when we add paintables later, we take the
       difference of the new options and the currentOptions.  If
       a removed paintable is still in currentOptions, the difference
       will remove it from the new options and it won't be re-added) */
	public void removePainterFromOpts(RegionPaintable p) {
		switch (p.getOptionKey()) {
			case SeqViewOptions.SEQDATA:
				currentOptions.seqExpts.remove(p.getOptionInfo());
			case SeqViewOptions.GENES:
				currentOptions.genes.remove(p.getOptionInfo());
			case SeqViewOptions.NCRNAS:
				currentOptions.ncrnas.remove(p.getOptionInfo());
			case SeqViewOptions.OTHERANNOTS:
				currentOptions.otherannots.remove(p.getOptionInfo());
			case SeqViewOptions.SEQLETTERS:
				currentOptions.seqletters = false;
			case SeqViewOptions.GCCONTENT:
				currentOptions.gccontent = false;            
		}
	}

	/* removes all the painters in a track from the visualizer.  
       This removes the track from the visualizer and unregisters 
       the RegionPanel as a listener from the Paintable.  The Paintable
       should keep track of how many listeners it has and should 
       unregister itself from its models when it has no listeners.  The Models, in turn,
       should call their stopRunning() method when they have no listeners */    
	public void removeTrack(String trackname) {
		if (!painters.containsKey(trackname)) {return;}
		ArrayList<RegionPaintable> plist = painters.get(trackname);            
		for (int i = 0; i < plist.size(); i++) {
			plist.get(i).removeEventListener(this);
			if (!plist.get(i).hasListeners()) {
				removePainterFromOpts(plist.get(i));
				allPainters.remove(plist.get(i));
				painterModelMap.remove(plist.get(i));
				painterCount--;
			}            
		}
		if(trackPaintOrder.containsKey(trackname))  
			removeTrackOrder(trackname,trackPaintOrder);
		painters.remove(trackname);
		for (RegionModel m : (HashSet<RegionModel>)allModels.clone()) {
			if (!m.hasListeners()) {
				allModels.remove(m);
			}
		}        
		repaint();
	}
	/* Removes a track with specified name from the track order.  Shifts all
       the other tracks up one to fill in the empty space */
	private void removeTrackOrder(String trackName, Hashtable<String,Integer> table) {
		if (!table.containsKey(trackName)) {return;}
		int value = table.get(trackName);
		table.remove(trackName);
		for (String k : table.keySet()) {
			if (table.get(k) > value) {
				table.put(k,table.get(k) - 1);
			}
		}
	}

	public SeqViewOptions getCurrentOptions() {
		return currentOptions;
	}

	public void actionPerformed(ActionEvent e) { 
		if(e.getSource() == leftButton) { 
			int increment = (int)(currentRegion.getWidth() * .25);
			mainPanel.paintTransBlock(this.getGraphics());
			setRegion(new Region(currentRegion.getGenome(),
					currentRegion.getChrom(),
					currentRegion.getStart() - increment,
					currentRegion.getEnd() - increment ));
		}

		if(e.getSource() == rightButton) { 
			int increment = (int)(currentRegion.getWidth() * .25);
			mainPanel.paintTransBlock(this.getGraphics());
			setRegion(new Region(currentRegion.getGenome(),
					currentRegion.getChrom(),
					currentRegion.getStart() + increment,
					currentRegion.getEnd() + increment ));
		}

		if (e.getSource() == farLeftButton) {
			int increment = (int)(currentRegion.getWidth() * .85);
			mainPanel.paintTransBlock(this.getGraphics());
			setRegion(new Region(currentRegion.getGenome(),
					currentRegion.getChrom(),
					currentRegion.getStart() - increment,
					currentRegion.getEnd() - increment ));
		}

		if (e.getSource() == farRightButton) {
			int increment = (int)(currentRegion.getWidth() * .85);
			mainPanel.paintTransBlock(this.getGraphics());
			setRegion(new Region(currentRegion.getGenome(),
					currentRegion.getChrom(),
					currentRegion.getStart() + increment,
					currentRegion.getEnd() + increment ));
		}

		if(e.getSource() == zoomInButton) { 
			int increment = (int)(currentRegion.getWidth() * .25);
			mainPanel.paintTransBlock(this.getGraphics());
			setRegion(new Region(currentRegion.getGenome(),
					currentRegion.getChrom(),
					currentRegion.getStart() + increment,
					currentRegion.getEnd() - increment));
		}

		if(e.getSource() == zoomOutButton) { 
			int increment = (int)(currentRegion.getWidth() * .5);
			mainPanel.paintTransBlock(this.getGraphics());
			setRegion(new Region(currentRegion.getGenome(),
					currentRegion.getChrom(),
					currentRegion.getStart() - increment,
					currentRegion.getEnd() + increment));
		}        
		if (e.getSource() == locationField) {
			Region r = regionFromString(genome,locationField.getText().trim());
			if (r != null) {
				mainPanel.paintTransBlock(this.getGraphics());
				setRegion(r);
			}
		}
		if (e.getSource() == refreshButton) {
			repaint();
		}
	}
	/* need to finish grabbing the key stuff */
	public void keyPressed(KeyEvent e) {}
	public void keyReleased(KeyEvent e) {}
	public void keyTyped(KeyEvent e) {}
	
	// this is a custom String -> int routine that handles suffixes such as k and m 
	public void setRegion (Region newRegion) {
		if (newRegion.getChrom().matches("^chr.*")) {            
			newRegion = new Region(newRegion.getGenome(),
					newRegion.getChrom().replaceAll("^chr",""),
					newRegion.getStart(),
					newRegion.getEnd());
		}
		if (newRegion.getEnd() - newRegion.getStart() < 30) {
			newRegion = new Region(newRegion.getGenome(),
					newRegion.getChrom(),
					newRegion.getStart() - 15,
					newRegion.getEnd() + 15);
		}
		if (newRegion.getStart() < 1) {
			newRegion = new Region(newRegion.getGenome(),
					newRegion.getChrom(),
					1,
					newRegion.getEnd());
		}
		if (newRegion.getEnd() > newRegion.getGenome().getChromLength(newRegion.getChrom())) {
			newRegion = new Region(newRegion.getGenome(),
					newRegion.getChrom(),
					newRegion.getStart(),
					newRegion.getGenome().getChromLength(newRegion.getChrom()));
		}

		if (!newRegion.equals(currentRegion) || forceupdate) {
			currentRegion = newRegion;
			locationField.setText(currentRegion.getLocationString());
			
			/* kick the painters here to give them a little extra time
               to update their data before we try to paint them */
    		   readyCount = 0;
    		   // set the new region in the paintables
    		   for (RegionPaintable p : allPainters) {
    			   p.setRegion(newRegion);
    		   	}
    		   // set the new region in the models
    		   // and call notify on them.  This gives a RegionModel
    		   // the option of wait()ing in a separate thread
    		   // if it so desires.
    		   synchronized(allModels){
	    		   for (RegionModel m : allModels) {
	    			   synchronized(m) {
	    				   m.setRegion(newRegion);
	    				   m.notifyAll();
	    			   }
	    		   }
    		   }
    		   //repaint();
		}
	}
       
	public boolean isClosed(){return closed;}
       
	public void close() { 
		synchronized(allModels){
			for (RegionModel m : allModels) {
				synchronized(m) {
					m.stopRunning();
					m.notifyAll();
			   	}
	 	   	}
		}
		for(RegionPaintable rp : allPainters) { 
			rp.cleanup();
		}
		closed=true;
		try {
			Thread.sleep(400);
		} catch (Exception e) {
			
		}
	}

	public static Region regionFromGFFString(Genome genome, String input){
    	   String [] pieces = input.split("\t");
    	   if(!pieces[0].startsWith("#")){
    		   String chromStr = pieces[0];
    		   int start = Integer.parseInt(pieces[3].replaceAll(",", ""));
    		   int end = Integer.parseInt(pieces[4].replaceAll(",", ""));
    		   if (chromStr.startsWith("chr")) {
    			   chromStr = chromStr.substring(3, chromStr.length());
    		   }
    		   char strand = pieces[6].charAt(0);
    		   if(strand=='.')
    			   return new Region(genome, chromStr, start, end);
    		   else
    			   return new StrandedRegion(genome, chromStr, start, end, strand);
    	   }else{
    		   return null;
    	   }
	}
       
       /* This parses input from the region/location text area and turns it into a region.
       	  If the text is parseable as a Region, this is easy.  Otherwise,
       	  try to turn it into a gene and use that.
        */
       public static Region regionFromString(Genome genome, String input) {
    	   String trimmed = input.trim().replaceAll("\\s+","");        
    	   Region r = Region.fromString(genome,input);
    	   if (r == null) {
    		   try {
    			   int upstream = 0;
    			   int downstream = 0;
    			   Pattern pattern = Pattern.compile("(.*)\\+(\\d+)\\-(\\d+)$");
    			   Matcher matcher = pattern.matcher(trimmed);
    			   if (matcher.matches()) {
    				   trimmed = matcher.group(1);
    				   upstream = Integer.parseInt(matcher.group(2));
    				   downstream = Integer.parseInt(matcher.group(3));
    			   }
    			   RefGeneGenerator generator = new RefGeneGenerator(genome);
    			   Iterator<Gene> iter = generator.byName(trimmed);
    			   if (iter.hasNext()) {
    				   Gene gene = iter.next();
    				   return new Region(gene.getGenome(),
							   gene.getChrom(),
							   gene.getStart(), gene.getEnd());
    			   }
    			   RegionExpanderFactoryLoader<Gene> gfLoader = new RegionExpanderFactoryLoader<Gene>("gene");
    			   for(String type : gfLoader.getTypes(genome)) {
    				   RegionExpanderFactory<Gene> genefactory = gfLoader.getFactory(genome,
    						   type);
    				   Expander<Region,Gene> expander = genefactory.getExpander(genome);
    				   if (expander instanceof RefGeneGenerator) {
    					   iter = ((RefGeneGenerator)expander).byName(trimmed);
    					   if (iter.hasNext()) {
    						   Gene gene = iter.next();
    						   return new Region(gene.getGenome(),
    								   gene.getChrom(),
    								   gene.getStart() - (gene.getStrand() == '+' ? upstream : downstream),
    								   gene.getEnd() + (gene.getStrand() == '+' ? downstream : upstream));
    					   }   
    				   }
    			   }
    		   } catch (DatabaseException ex) {
    			   // this means we couldn't get a ref gene table for this species
    			   ex.printStackTrace();
    		   }
    	   }
    	   return r;
       }
       public static java.util.List<Region> readRegionsFromFile(Genome g, String filename, boolean isGFF) {
    	   ArrayList<Region> regions = new ArrayList<Region>();
    	   try {
    		   BufferedReader r = new BufferedReader(new FileReader(filename));
    		   String s;
    		   while ((s = r.readLine()) != null) {
    			   Region region = isGFF ? regionFromGFFString(g,s) : regionFromString(g,s);
    			   if (region != null) {
    				   regions.add(region);
    			   } else {
    				   System.err.println("Couldn't parse " + s);
    			   }

    		   }
    		   r.close();
    	   } catch (IOException ex) {
    		   throw new RuntimeException("Can't read " + filename,ex);
    	   }
    	   return regions;
       }

       /** add a new painter to this panel and have it load its default
        * properties
        */
       public void addPainter(RegionPaintable p) {
    	   if (p == null) {return;}
    	   painterCount++;
    	   String pk = p.getLabel();
    	   if (painters.get(pk) == null) {
    		   painters.put(pk,new ArrayList<RegionPaintable>());
    		   trackSizeFactor.put(pk, 1.0);
    		   ulx.put(pk,0); uly.put(pk,0);
    		   lrx.put(pk,0); lry.put(pk,0);            
    	   }
    	   painters.get(pk).add(p);
    	   p.getProperties().loadDefaults();
    	   allPainters.add(p);
    	   
       }

       public void changePainter(RegionPaintable oldPainter, RegionPaintable newPainter) {
    	   String pk = oldPainter.getLabel();
    	   if ((oldPainter == null) || (newPainter == null) || (!painters.containsKey(pk))) {return;}
    	   newPainter.addEventListener(this);
    	   newPainter.setLabel(oldPainter.getLabel());

    	   ArrayList<RegionPaintable> painterList = painters.get(pk);
    	   painterList.set(painterList.indexOf(oldPainter), newPainter);
    	   allPainters.set(allPainters.indexOf(oldPainter), newPainter);
    	   painterModelMap.put(newPainter,painterModelMap.get(oldPainter));
    	   painterModelMap.remove(oldPainter);
    	   oldPainter.cleanup();
    	   for (RegionModel m : painterModelMap.get(newPainter)) {
    		   m.notifyListeners();
    	   }
    	   repaint();
       }     

       public void addModel(RegionModel m) {
    	   if (m == null) {
    		   throw new NullPointerException("Don't you give me a null model");
    	   }
    	   synchronized(allModels){
    		   allModels.add(m);
    	   }
       }
       public void addModelToPaintable(RegionPaintable p, RegionModel m) {
    	   if (painterModelMap.get(p) == null) {
    		   painterModelMap.put(p,new ArrayList<RegionModel>());
    	   }
    	   painterModelMap.get(p).add(m);
    	   m.getProperties().loadDefaults();
       }
       public void removeModel(RegionModel m) {
    	   for (RegionPaintable k : painterModelMap.keySet()) {
    		   ArrayList<RegionModel> l = painterModelMap.get(k);
    		   l.remove(m);
    	   }
    	   synchronized(allModels){
    		   allModels.remove(m);
    	   }
       }
       /* recompute the layout for this panel.  Newly added painters
       will not be visible until you call this method */
       public int computeLayout(int x, int y, int width, int height) {
    	   height-=1; //otherwise you always get a scrollbar
    	   Set<String> keyset = painters.keySet();
    	   String[] keys = new String[keyset.size()];
    	   int i = 0;
    	   for (String s : painters.keySet()) { keys[i++] = s; }
    	   Arrays.sort(keys,new AddedOrderComparator());

    	   //Find the minimum space required by all tracks
    	   int minTotalHeight = 0;
    	   for (i = 0; i < keys.length; i++) {
    		   String s = keys[i];
    		   int maxspace = 0;
    		   ArrayList<RegionPaintable> plist = painters.get(s);
    		   for (int j = 0; j < plist.size(); j++) {
    			   double requestd = plist.get(j).getMinVertSpace()*trackSizeFactor.get(s);
    			   int request = (int)requestd;
    			   if (maxspace < request)
    				   maxspace = request;
    		   }
    		   minTotalHeight+=maxspace;
    	   }
    	   boolean fillSpace=false;
    	   if(minTotalHeight<height)
    		   fillSpace=true;
    	   else
    		   height = minTotalHeight;
    	   int ypos = height;
    	   int thickSpace = height;
    	   int thickCount = 0;
    	   
    	   //Find height requests
    	   HashMap<String,Boolean> thickMap = new HashMap<String,Boolean>();
    	   Hashtable<String,Integer> requests = new Hashtable<String,Integer>();
    	   for (i = 0; i < keys.length; i++) {
    		   String s = keys[i];
    		   boolean isthick = false;
    		   int maxspace = 0;
    		   ArrayList<RegionPaintable> plist = painters.get(s);
    		   for (int j = 0; j < plist.size(); j++) {
    			   int request = fillSpace ? plist.get(j).getMaxVertSpace() : plist.get(j).getMinVertSpace();
    			   if (request == -1) {
    				   isthick = true;
    				   continue;
    			   }
    			   if (request>0){
    				   double tmpr = request*trackSizeFactor.get(s);
    				   request = (int)tmpr;
    				   if(maxspace < request)
    					   maxspace = request;
    			   }
    		   }
    		   requests.put(s,maxspace);
    		   thickMap.put(s, isthick);
    		   if(!trackPaintOrder.containsKey(s)) { 
    			   trackPaintOrder.put(s, trackPaintOrder.size() + 1);
    		   }
    	   }

    	   //Allocate non-autofill track space first
    	   for(String s : thickMap.keySet()) {
    		   if(!thickMap.get(s)) {
    			   int allocated = requests.get(s);
    			   thickSpace -= allocated;
    		   } else { 
    			   thickCount += 1;
    		   }
    	   }
    	   //Split the remaining space between autofill tracks
    	   double totalSizeFactor = 0;
    	   for(String s : thickMap.keySet()) {
    		   if(thickMap.get(s))
    			   totalSizeFactor+=trackSizeFactor.get(s);
    	   }
    	   //int thickAlloc = (thickSpace - 5*thickCount) / (Math.max(1, thickCount));
    	   double thickAllocd = (thickSpace) / (Math.max(1, totalSizeFactor));
    	   int thickAlloc = (int)thickAllocd;

    	   keys = new String[trackPaintOrder.size()];
    	   i = 0;
    	   for(String s : trackPaintOrder.keySet()) { keys[i++] = s; }
    	   Arrays.sort(keys, new HashtableComparator(trackPaintOrder));

    	   // layout from the bottom up
    	   for (i = keys.length - 1; i >=0; i--) {
    		   String s = keys[i]; 
    		   int allocated;
    		   if(thickMap.get(s)) {
    			   double a = thickAlloc*trackSizeFactor.get(s);
    			   allocated = (int) a;
    		   } else {
    			   allocated = requests.get(s);
    		   }
    		   ulx.put(s,x);
    		   lrx.put(s,width + x-5);

    		   if(thickMap.get(s)) { 
    			   lry.put(s,ypos-5);
    			   uly.put(s,ypos-allocated+5);            	
    		   } else { 
    			   lry.put(s,ypos);
    			   uly.put(s,ypos-allocated);
    		   }

    		   ypos -= allocated;  
    	   }
    	   return height;
       }
       
       
       /* this is a callback from a painter to the RegionPanel saying that
       the painter is ready to be painted.  
       This implementation is just a heuristic, but the goal is to avoid
       calling repaint() until all the painters are done */
       public synchronized void eventRegistered(EventObject e) {        
    	   if (e.getSource() instanceof VizPaintable) {            
    		   readyCount++;
    		   if (readyCount >= painterCount) {
    			   repaint();
    			   status.setStatus("Ready", Color.black);
    		   }else{
    			   status.setStatus("Waiting for data (received "+readyCount+"/"+painterCount+")", Color.red);
    		   }
    	   }
       }

       public void paintComponent(Graphics g) {
    	   int oheight = scrollPane.getViewport().getHeight();
    	   int mheight = computeLayout(getX(),getY(),getWidth(),oheight);
    	   mainPanel.setPreferredSize(new Dimension(getWidth(), mheight));
    	   mainPanel.revalidate();
       }

       public boolean allCanPaint() {
    	   boolean canpaint = true;
    	   for (String s : painters.keySet()) {
    		   ArrayList<RegionPaintable> plist = painters.get(s);
    		   for (int i = 0; i < plist.size(); i++) {
    			   canpaint = canpaint && plist.get(i).canPaint();
    		   }
    	   }
    	   return canpaint;
       }

       class RegionContentPanel extends JPanel {
    	   private Hashtable painted = new Hashtable<Object,Boolean>();

    	   //Transparent block when transitioning
    	   public void paintTransBlock(Graphics g){
    		   if(!paintedTransparent){
				   g.setColor(transparentWhite);
				   paintedTransparent=true;
				   g.fillRect(0,0,mainPanel.getWidth(),mainPanel.getHeight());
			   }
    	   }
    	   
    	   public void paintComponent(Graphics g) {
    		   paintComponent(g,mainPanel.getX(),mainPanel.getY(),
    				   mainPanel.getWidth(),mainPanel.getHeight());
    	   }

    	   public void paintComponent(Graphics g, int x, int y, int width, int height) {
    		   /* two passes: first make sure everyone can paint.  
               then do the painting */
    		   Graphics2D graphics = (Graphics2D) g;
    		   boolean canpaint = true;
    		   for (String s : painters.keySet()) {
    			   ArrayList<RegionPaintable> plist = painters.get(s);
    			   for (int i = 0; i < plist.size(); i++) {
    				   canpaint = canpaint && plist.get(i).canPaint();
    			   }
    		   }
    		   if (!canpaint) {
    			   if(!paintedTransparent){
    				   g.setColor(transparentWhite);
    				   paintedTransparent=true;
    				   g.fillRect(0,0,width,height);
    			   }
    			   return;
    		   }
    		   //If we're here, everyone can paint
    		   paintedTransparent=false;
    		   g.setColor(Color.WHITE);
    		   g.fillRect(0,0,width,height);
    		   for (String s : painters.keySet()) {
    			   ArrayList<RegionPaintable> plist = painters.get(s);
    			   for (int i = 0; i < plist.size(); i++) {
    				   try {
    					   plist.get(i).paintItem(graphics,
    							   ulx.get(s),
    							   uly.get(s),
    							   lrx.get(s),
    							   lry.get(s));
    				   } catch (Exception e) {
    					   e.printStackTrace();
    					   graphics.setColor(Color.RED);
    					   graphics.drawString("Error: " + e.toString(),ulx.get(s),lry.get(s));
    				   }
    			   }
    		   }
    	   }
       }
       public void saveImage(File f, int w, int h, boolean raster) throws IOException { 
    	   if (raster) {
    		   BufferedImage im = 
    				   new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
    		   Graphics g = im.getGraphics();
    		   Graphics2D g2 = (Graphics2D)g;
    		   g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
    		   mainPanel.paintComponent(g,0,0,w,h);
    		   ImageIO.write(im, "png", f);
    		   g.dispose();
    	   } else {
    		   DOMImplementation domImpl =
    				   GenericDOMImplementation.getDOMImplementation();
    		   // Create an instance of org.w3c.dom.Document
    		   Document document = domImpl.createDocument(null, "svg", null);
    		   // Create an instance of the SVG Generator
    		   SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
    		   svgGenerator.setSVGCanvasSize(new Dimension(w,h));
    		   // Ask the test to render into the SVG Graphics2D implementation
    		   mainPanel.paintComponent(svgGenerator);

    		   // Finally, stream out SVG to the standard output using UTF-8
    		   // character to byte encoding
    		   boolean useCSS = true; // we want to use CSS style attribute
    		   Writer out = new OutputStreamWriter(new FileOutputStream(f), "UTF-8");
    		   svgGenerator.stream(out, useCSS);
    	   }
       }
       public Genome getGenome() {return genome;}
       public Region getRegion() {return currentRegion;}

       public boolean equals(Object o) {
    	   if (o instanceof RegionPanel) {
    		   return (o == this);
    	   } else {
    		   return false;
    	   }
       }

       public void addTrackFromFile(boolean isGFF) {
    	   JFileChooser chooser;
    	   chooser = new JFileChooser(currDirectory);
    	   int v = chooser.showOpenDialog(null);
    	   if(v == JFileChooser.APPROVE_OPTION) { 
    		   File f = chooser.getSelectedFile();
    		   currDirectory = chooser.getCurrentDirectory();
    		   addTrackFromFile(f.getAbsolutePath(), f.getName(), isGFF);
    	   }
       }
       public void addTrackFromFile(String fname, String trackname, boolean isGFF) {
    	   System.err.println("Adding " + fname + " with name " + trackname);
    	   status.setStatus("Adding track from "+fname, Color.red);
    	   try {
    		   java.util.List<Region> regions = readRegionsFromFile(genome,fname,isGFF);
    		   StaticExpander<Region,Region> expander = new StaticExpander<Region,Region>(regions);
    		   RegionExpanderModel<Region> model = new RegionExpanderModel<Region>(expander);
    		   addModel(model);
    		   Thread t = new Thread(model);
    		   t.start();
    		   NamedStrandedPainter p = new NamedStrandedPainter(model);
    		   p.setLabel(trackname);
    		   p.addEventListener(this);
    		   addPainter(p);
    		   setRegion(currentRegion);
    		   this.forceModelUpdate();
    	   } catch (Exception e) {
    		   e.printStackTrace();
    		   status.setStatus("Error adding track from "+fname, Color.red);
    	   }
       }

       private class HashtableComparator implements Comparator<String> {
    	   private Hashtable<String,Integer> t;
    	   public HashtableComparator(Hashtable<String,Integer> t) {
    		   this.t = t;
    	   }
    	   public int compare(String a, String b) {
    		   return t.get(a) - t.get(b);
    	   }
    	   public boolean equals(Object o) {
    		   return (o == this);
    	   }
       }

       private class AddedOrderComparator implements Comparator<String> {
    	   public int compare(String a, String b) {
    		   int mina, minb, i;
    		   mina = 1000000;
    		   minb = 1000000;
    		   ArrayList<RegionPaintable> plist = painters.get(a);
    		   for (i = 0; i < plist.size(); i++) {
    			   int added = allPainters.lastIndexOf(plist.get(i));
    			   if (added < mina) {
    				   mina = added;
    			   }
    		   }
    		   plist = painters.get(b);
    		   for (i = 0; i < plist.size(); i++) {
    			   int added = allPainters.lastIndexOf(plist.get(i));
    			   if (added < minb) {
    				   minb = added;
    			   }
    		   }
    		   return mina - minb;
    	   }
    	   public boolean equals(Object o) {
    		   return (o == this);
    	   }
       }

       class ConfigureActionListener implements ActionListener {
    	   private String k;
    	   public ConfigureActionListener(String k) {this.k = k;}
    	   public void actionPerformed(ActionEvent e) {configureTrack(k);}
       }
       class RemoveActionListener implements ActionListener {
    	   private String k;
    	   private RegionPanel panel;
    	   public RemoveActionListener(String k, RegionPanel p) {this.k = k; this.panel = p;}
    	   public void actionPerformed(ActionEvent e) {removeTrack(k); panel.repaint();}
       }
       class SavePrefsActionListener implements ActionListener {
    	   private ArrayList<RegionPaintable> p;
    	   public SavePrefsActionListener (ArrayList<RegionPaintable> p) { this.p = p;}
    	   public void actionPerformed(ActionEvent e) {
    		   for (RegionPaintable paintable : p) {
    			   paintable.getProperties().saveToFile();
    			   if (painterModelMap.get(paintable) != null) {
    				   for (RegionModel m : painterModelMap.get(paintable)) {
    					   if (m instanceof SeqViewModel) {
    						   ((SeqViewModel)m).getProperties().saveToFile();
    					   }
    				   }
    			   }
    		   }
    	   }
       } 
       class LoadPrefsActionListener implements ActionListener {
    	   private ArrayList<RegionPaintable> p;
    	   public LoadPrefsActionListener (ArrayList<RegionPaintable> p) { this.p = p;}
    	   public void actionPerformed(ActionEvent e) {
    		   for (RegionPaintable paintable : p) {
    			   // save TrackLabel so that it is not overwritten
    			   String trackLabel = paintable.getProperties().TrackLabel;
    			   paintable.getProperties().loadFromFile();
    			   paintable.getProperties().TrackLabel = trackLabel;
    			   if (painterModelMap.get(paintable) != null) {
    				   for (RegionModel m : painterModelMap.get(paintable)) {
    					   if (m instanceof SeqViewModel) {
    						   ((SeqViewModel)m).getProperties().loadFromFile();
    					   }
    				   }
    			   }
    		   }
    	   }
       } 
       class SaveAllPrefsActionListener implements ActionListener {
    	   private RegionPanel panel;
    	   public SaveAllPrefsActionListener (RegionPanel p) {this.panel = p;}
    	   public void actionPerformed(ActionEvent e) {
    		   JFileChooser chooser = new JFileChooser(currDirectory);
    		   chooser.setDialogTitle("Save preferences to...");
    		   chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
    		   int returnVal = chooser.showSaveDialog(null);
    		   if (returnVal == JFileChooser.APPROVE_OPTION) {
    			   File f = chooser.getSelectedFile();
    			   currDirectory = chooser.getCurrentDirectory();
    			   for (String k : painters.keySet()) {
    				   for (RegionPaintable p : painters.get(k)) {
    					   p.savePropsInDir(f);
    				   }
    			   }
    		   }
    	   }
       } 
       class LoadAllPrefsActionListener implements ActionListener {
    	   private RegionPanel panel;
    	   public LoadAllPrefsActionListener (RegionPanel p) {this.panel = p;}
    	   public void actionPerformed(ActionEvent e) {
    		   JFileChooser chooser = new JFileChooser(currDirectory);
    		   chooser.setDialogTitle("Load preferences from...");
    		   chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
    		   int returnVal = chooser.showOpenDialog(null);
    		   if (returnVal == JFileChooser.APPROVE_OPTION) {
    			   File f = chooser.getSelectedFile();
    			   currDirectory = chooser.getCurrentDirectory();
    			   for (String k : painters.keySet()) {
    				   for (RegionPaintable p : painters.get(k)) {
    					   p.loadPropsInDir(f);
    				   }
    			   }
    		   }
    		   panel.repaint();
    	   }
       } 


       class MoveInLayoutActionListener implements ActionListener {
    	   private String key;
    	   private int dir;
    	   private RegionPanel panel;
    	   public MoveInLayoutActionListener(String k, int dir, RegionPanel p) {
    		   this.key = k;
    		   this.dir = dir;
    		   this.panel = p;
    	   }
    	   public void actionPerformed(ActionEvent e) {
    		   Hashtable<String,Integer> t = null;

    		   /*
	            if (trackPaintOrderThick.containsKey(key)) {
	                t = trackPaintOrderThick;
	            } else if (trackPaintOrderThin.containsKey(key)) {
	                t = trackPaintOrderThin;
	            } else {
	                return;
	            }
    		    */
    		   if(trackPaintOrder.containsKey(key)) { 
    			   t = trackPaintOrder;
    		   } else { 
    			   return;
    		   }

    		   int old = t.get(key);
    		   int desired = old + dir;
    		   for (String k : t.keySet()) {
    			   if (t.get(k) == desired) {
    				   t.put(k,old);
    				   t.put(key,desired);
    				   break;
    			   }
    		   }
    		   panel.repaint();
    	   }
       }
       class ChangeTrackSize implements ActionListener {
    	   private String key;
    	   private double factor;
    	   private RegionPanel panel;
    	   public ChangeTrackSize(String key, double factor, RegionPanel p) {
    		   this.key = key;
    		   this.factor = factor;
    		   this.panel = p;
    	   }
    	   public void actionPerformed(ActionEvent e) {
    		   trackSizeFactor.put(key ,
    				   (double)(trackSizeFactor.get(key) * factor));
    		   panel.repaint();
    	   }
       }

       /* respond to a right-click on a track.  
       Since tracks may overlap, keys is a set of Strings
       that describes the propertykeys for all of the tracks that
       we need to be able to handle */
       public JPopupMenu trackRightClick(ArrayList<String> keys) {
    	   JPopupMenu menu = new JPopupMenu("Track Setup");
    	   JMenuItem item;
    	   for (String k : keys) {
    		   item = new JMenuItem("Configure " + k.substring(0,Math.min(60,k.length())));
    		   item.addActionListener(new ConfigureActionListener(k));
    		   menu.add(item);
    		   item = new JMenuItem("Remove " + k.substring(0,Math.min(60,k.length())));
    		   item.addActionListener(new RemoveActionListener(k,this));
    		   menu.add(item);

    		   item = new JMenuItem("Save Prefs");
    		   item.addActionListener(new SavePrefsActionListener(painters.get(k)));
    		   menu.add(item);
    		   item = new JMenuItem("Load Prefs");
    		   item.addActionListener(new LoadPrefsActionListener(painters.get(k)));
    		   menu.add(item);
    		   item = new JMenuItem("Save All Prefs");
    		   item.addActionListener(new SaveAllPrefsActionListener(this));
    		   menu.add(item);
    		   item = new JMenuItem("Load All Prefs");
    		   item.addActionListener(new LoadAllPrefsActionListener(this));
    		   menu.add(item);

    		   item = new JMenuItem("Move Up in Layout");
    		   item.addActionListener(new MoveInLayoutActionListener(k,-1,this));
    		   menu.add(item);
    		   item = new JMenuItem("Move Down in Layout");
    		   item.addActionListener(new MoveInLayoutActionListener(k,1,this));
    		   menu.add(item);
    		   item = new JMenuItem("Increase track size");
    		   item.addActionListener(new ChangeTrackSize(k,1.3,this));
    		   menu.add(item);
    		   item = new JMenuItem("Decrease track size");
    		   item.addActionListener(new ChangeTrackSize(k,.75,this));
    		   menu.add(item);
    	   }
    	   return menu;
       }

       public void configureTrack(String track) {
    	   if (painters.containsKey(track)) {
    		   ArrayList<SeqViewProperties> props = new ArrayList<SeqViewProperties>();
    		   ArrayList<RegionPaintable> plist = painters.get(track);            
    		   for (RegionPaintable p : plist) {                
    			   props.add(p.getProperties());
    			   if (painterModelMap.get(p) == null) {
    				   continue;
    			   }
    			   for (RegionModel m : painterModelMap.get(p)) {
    				   props.add(m.getProperties());
    			   }
    		   }
    		   SeqViewProperties.configure(props,this, false);
    	   }
       }

       /**
        * Batch configure all SeqHistogram tracks
        */
       public void configSeqDataTracksBatch(){
    	   ArrayList<SeqViewProperties> props = new ArrayList<SeqViewProperties>();
    	   SeqDataBatchProperties batchProp = new SeqDataBatchProperties();
    	   batchProp.loadDefaults();
    	   props.add(batchProp);
    	   SeqViewProperties.configure(props,this, true);
       }
       
       public void batchUpdateModels(Collection<? extends Model> models){
		   for(RegionPaintable p : allPainters){
    		   PaintableProperties currProps = p.getProperties();
    		   
    		   for(Model i : models){
    			   currProps.setFromModel(i);
    			   if (painterModelMap.get(p) != null) {
        			   for (RegionModel m : painterModelMap.get(p)) {
        				   ModelProperties modProps = m.getProperties();
        				   modProps.setFromModel(i);
        			   }
        		   }

    		   }
		   }
       }
       
       public void forceModelUpdate(){
    	   synchronized(allModels){
	    	   for (RegionModel m : allModels) {
				   synchronized(m) {
					   m.resetRegion(currentRegion);
					   m.notifyAll();
				   }
			   }
    	   }
		   repaint(); 
       }
       
       public void mouseClicked(MouseEvent e) {
    	   if (e.getButton() == MouseEvent.BUTTON3 || e.isPopupTrigger()) {
    		   int xpos = e.getX();
    		   int ypos = e.getY();
    		   ArrayList<String> keys = new ArrayList<String>();
    		   for (String pk : painters.keySet()) {
    			   if (ulx.get(pk) <= xpos && 
    					   lrx.get(pk) >= xpos &&
    					   uly.get(pk) <= ypos &&
    					   lry.get(pk) >= ypos) {
    				   keys.add(pk);
    			   }
    		   }
    		   JPopupMenu m = trackRightClick(keys);
    		   m.show(mainPanel,xpos,ypos);
    	   } else {
    		   int xpos = e.getX();
    		   int ypos = e.getY();
    		   int totalitems = 0;
    		   JPopupMenu m = new JPopupMenu("Stuff");
    		   for (String pk : painters.keySet()) {
    			   if (ulx.get(pk) <= xpos && 
    					   lrx.get(pk) >= xpos &&
    					   uly.get(pk) <= ypos &&
    					   lry.get(pk) >= ypos) {
    				   for (RegionPaintable p : painters.get(pk)) {
    					   p.mouseClicked(e);
    					   ArrayList<JMenuItem> items = p.mouseClickedMenu(e);
    					   if (items == null) {continue;}
    					   String k = p.getLabel();
    					   if (totalitems > 0) {
    						   m.addSeparator();
    					   }
    					   JMenuItem label = new JMenuItem(k);
    					   label.setEnabled(false);
    					   m.add(label);
    					   for (JMenuItem item : items) {
    						   m.add(item);
    					   }
    					   totalitems++;
    				   }
    			   }
    		   }
    		   if (totalitems > 0){
    			   m.show(e.getComponent(),xpos,ypos);
    		   }
    	   }
       }

       public void mouseEntered(MouseEvent e) {}

       public void mouseExited(MouseEvent e) {}

       public void mousePressed(MouseEvent e) {
    	   int xpos = e.getX();
    	   int ypos = e.getY();
    	   for (String pk : painters.keySet()) {
    		   if (ulx.get(pk) <= xpos && 
    				   lrx.get(pk) >= xpos &&
    				   uly.get(pk) <= ypos &&
    				   lry.get(pk) >= ypos) {
    			   for (RegionPaintable p : painters.get(pk)) {
    				   p.mousePressed(e);
    			   }
    		   }
    	   }
       }

       public void mouseReleased(MouseEvent e) {
    	   int xpos = e.getX();
    	   int ypos = e.getY();
    	   for (String pk : painters.keySet()) {
    		   if (ulx.get(pk) <= xpos && 
    				   lrx.get(pk) >= xpos &&
    				   uly.get(pk) <= ypos &&
    				   lry.get(pk) >= ypos) {
    			   for (RegionPaintable p : painters.get(pk)) {
    				   p.mouseReleased(e);
    			   }
    		   }
    	   }
       }
}
