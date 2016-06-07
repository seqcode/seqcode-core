package org.seqcode.projects.seqview.components;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import org.seqcode.data.seqdata.SeqLocatorMatchedAligns;
import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.Gene;
import org.seqcode.genome.location.NamedTypedRegion;
import org.seqcode.gsebricks.RegionExpanderFactoryLoader;
import org.seqcode.gseutils.Closeable;
import org.seqcode.gseutils.NotFoundException;
import org.seqcode.projects.seqview.SeqViewOptions;

import java.util.Collection;
import java.util.List;
import java.util.ResourceBundle;
import java.util.HashMap;
import java.io.*;


public class SeqViewOptionsPane 
	extends JTabbedPane 
	implements ItemListener, ActionListener, Closeable {
	//Current genome settings
	Genome currGen = null;
	Species currSpecies = null;
    // regexes get special handling at the moment because there's no gui component for them,
    // so cache them if neccessary
    private HashMap<String,String> regexes;
    private boolean handlingChange, closed;
    private RegionExpanderFactoryLoader<Gene> gfLoader;
    private RegionExpanderFactoryLoader<NamedTypedRegion> annotLoader;

    private JPanel speciesLocationPanel,
        annotationsPanel,
        seqPanel,
        optionsPanel;
    private SeqViewOptions createdFrom;
    
    // species/location tab
    private JComboBox speciesCBox, genomeCBox;
    private JTextField position, gene;
    private JLabel specieslabel, genomelabel, positionlabel, genelabel;

    // options tab
    private JCheckBox hash, oldchipseq;
    
    // seqdata tab
    private SeqAlignmentSelectPanel seqSelect;

    // annotations tab
    private JList genes, ncrnas, otherfeats;
    private DefaultListModel genesmodel, ncrnasmodel, otherfeatsmodel;
    private JCheckBox seqletters, polyA, gccontent, pyrpurcontent, cpg, regexmatcher;
    private JLabel geneslabel, ncrnaslabel, otherfeatslabel, seqlimitlabel;
        
    // file-based tracks
    private FileBasedTracksPanel filetracks;

    public SeqViewOptionsPane () throws NotFoundException {
        super();
        gfLoader = new RegionExpanderFactoryLoader<Gene>("gene");
        annotLoader = new RegionExpanderFactoryLoader<NamedTypedRegion>("annots");
        handlingChange = true;
        init(null);
        handlingChange = false;
        closed = false;
        setSpeciesGenomeDefaults();
    }

    
    public SeqViewOptionsPane(SeqViewOptions opts) throws NotFoundException {
        super();
        gfLoader = new RegionExpanderFactoryLoader<Gene>("gene");
        annotLoader = new RegionExpanderFactoryLoader<NamedTypedRegion>("annots");
        closed = false;
        handlingChange=true;
        init(opts);
        handlingChange=false;
        if (opts.getGenome() == null) {
            setSpeciesGenomeDefaults();
        }        
        regexes = opts.regexes;
    }
    
    public boolean isClosed() { return closed; }
    
    public void close() { 
        seqSelect.close();
        closed = true; 
    }
    
    private void reloadGenome() { 
        String orgName = (String)speciesCBox.getSelectedItem();
        String genName = (String)genomeCBox.getSelectedItem();
        if(orgName != null && genName != null && !orgName.equals(currSpecies.getName()) && !genName.equals(currGen.getVersion())) { 
            Species org;
            try {
                org = Species.getSpecies(orgName);
                currSpecies = org;
                Genome gen = new Genome(org, genName);
                currGen = gen;
            } catch (NotFoundException e) {
                e.printStackTrace();
                throw new IllegalArgumentException(orgName + ":" + genName);
            }
        }
    }

    public void init(SeqViewOptions opts) throws NotFoundException {
    	speciesLocationPanel = new JPanel();
        annotationsPanel = new JPanel();
        seqPanel = new JPanel();
        optionsPanel = new JPanel();
        createdFrom = opts;
        
        // First tab lets the user select the species, genome version,
        // and the genomic coordinates they want to view
        specieslabel = new JLabel("Species");
        genomelabel = new JLabel("Genome Version");
        positionlabel = new JLabel("Genome position\nto view");
        genelabel = new JLabel("Gene to view");
        speciesCBox = new JComboBox();
        genomeCBox = new JComboBox();
        speciesCBox.addItemListener(this);
        genomeCBox.addItemListener(this);                
        position = new JTextField();
        position.setText("1:100-200");
        gene = new JTextField();

        /* need to fill the species and genome boxes here */
        Collection<String> organisms = Species.getAllSpeciesNames(false);
        for (String o : organisms) {
            speciesCBox.addItem(o);
        }
        //Set selected species & genome
        if (opts!=null && opts.getGenome() != null) {
            speciesCBox.removeItemListener(this);
            genomeCBox.removeItemListener(this);
            speciesCBox.setSelectedItem(opts.getGenome().getSpeciesName());
            Collection<String> genomes = opts.getGenome().getSpecies().getGenomeNames();
            genomeCBox.removeAllItems();
            for (String o : genomes) {
            	genomeCBox.addItem(o);
            }
            genomeCBox.setSelectedItem(opts.getGenome().getVersion());            
            currSpecies = opts.getGenome().getSpecies();
            currGen = opts.getGenome();
            genomeCBox.addItemListener(this);
            speciesCBox.addItemListener(this);
        } else {
            updateGenomeSelection();
            setSpeciesGenomeDefaults();
        }
        
        //Location panel
        speciesLocationPanel.setLayout(new GridLayout(4,2));
        speciesLocationPanel.add(specieslabel);
        speciesLocationPanel.add(speciesCBox);
        speciesLocationPanel.add(genomelabel);
        speciesLocationPanel.add(genomeCBox);
        speciesLocationPanel.add(positionlabel);
        speciesLocationPanel.add(position);
        speciesLocationPanel.add(genelabel);
        speciesLocationPanel.add(gene);
        
        // seqdata tab
        seqSelect = new SeqAlignmentSelectPanel();        
        seqPanel.setLayout(new BorderLayout());
        seqPanel.add(seqSelect, BorderLayout.CENTER);

        // Options tab
        optionsPanel.setLayout(new GridLayout(4,1));
        hash = new JCheckBox("Show chromosome coordinates");
        oldchipseq = new JCheckBox("Use old ChipSeq painter");
        hash.setSelected(true);
        optionsPanel.add(hash);
        optionsPanel.add(oldchipseq);
        
        // Annotations tab
        JPanel lists = new JPanel();
        lists.setLayout(new GridLayout(3,2));
        JPanel boxes = new JPanel();
        boxes.setLayout(new GridLayout(4,2));
        genesmodel = new DefaultListModel();
        ncrnasmodel = new DefaultListModel();
        otherfeatsmodel = new DefaultListModel();
        genes = new JList(genesmodel);
        ncrnas = new JList(ncrnasmodel);
        otherfeats = new JList(otherfeatsmodel);
        genes.setVisibleRowCount(7);genes.setLayoutOrientation(JList.VERTICAL);
        otherfeats.setVisibleRowCount(7); otherfeats.setLayoutOrientation(JList.VERTICAL);
        ncrnas.setVisibleRowCount(7); ncrnas.setLayoutOrientation(JList.VERTICAL);
        
        seqlimitlabel=new JLabel("Sequence features (visible only below "+SeqViewOptions.DEFAULT_MAX_SEQUENCE_VIEW+"bp):");
        seqletters = new JCheckBox("DNA sequence");
        seqletters.setSelected(true);
        polyA = new JCheckBox("PolyA sequences");
        gccontent = new JCheckBox("GC content");
        pyrpurcontent = new JCheckBox("Pyr/Pur content");
        cpg = new JCheckBox("CpG");
        regexmatcher = new JCheckBox("Regex Matcher");
        geneslabel = new JLabel("Genes");
        ncrnaslabel = new JLabel("ncRNAs");
        otherfeatslabel = new JLabel("Other annotations");

        lists.add(geneslabel);
        lists.add(new JScrollPane(genes));
        lists.add(ncrnaslabel);
        lists.add(new JScrollPane(ncrnas));
        lists.add(otherfeatslabel);
        lists.add(new JScrollPane(otherfeats));
        boxes.add(seqlimitlabel);
        boxes.add(new JLabel(""));
        boxes.add(seqletters);
        boxes.add(polyA);
        boxes.add(gccontent);
        boxes.add(pyrpurcontent);
        boxes.add(cpg);
        boxes.add(regexmatcher);

        annotationsPanel.setLayout(new BorderLayout());
        annotationsPanel.add(lists,BorderLayout.CENTER);
        annotationsPanel.add(boxes,BorderLayout.SOUTH);

        // file tracks tab
        filetracks = new FileBasedTracksPanel();        

        // use this to make the spacing be not-stupid
        JPanel dummy = new JPanel();        
        dummy.add(speciesLocationPanel);
        dummy.add(new JPanel());
        addTab("Species & Location",new JScrollPane(dummy));

        dummy = new JPanel();  dummy.add(annotationsPanel); dummy.add(new JPanel());
        addTab("Annotations",new JScrollPane(dummy));

        addTab("Seq Data", seqPanel);

        dummy = new JPanel();  dummy.add(filetracks); dummy.add(new JPanel());
        addTab("File Tracks",new JScrollPane(dummy));
                
        dummy = new JPanel();  
        dummy.setLayout(new BorderLayout());
        
        dummy = new JPanel();  dummy.add(optionsPanel); dummy.add(new JPanel());
        addTab("Display Options",new JScrollPane(dummy));
        
        if(opts !=null){
        	if (opts.gene != null &&
                !opts.gene.equals("")) {
                gene.setText(opts.gene);
            }
            hash.setSelected(opts.hash);
            gccontent.setSelected(opts.gccontent);        
            pyrpurcontent.setSelected(opts.pyrpurcontent);
            polyA.setSelected(opts.polya);
            cpg.setSelected(opts.cpg);
            regexmatcher.setSelected(opts.regexmatcher);
            seqletters.setSelected(opts.seqletters);
            oldchipseq.setSelected(!opts.seqHistogramPainter);
            seqlimitlabel.setText("Sequence features (visible only below "+opts.maxSequenceQuery+"bp):");
            int[] selected = new int[opts.genes.size()];
            for (int i = 0; i < opts.genes.size(); i++) {
                selected[i] = genesmodel.indexOf(opts.genes.get(i));            
            }
            genes.setSelectedIndices(selected);
            selected = new int[opts.ncrnas.size()];
            for (int i = 0; i < opts.ncrnas.size(); i++) {
                selected[i] = ncrnasmodel.indexOf(opts.ncrnas.get(i));            
            }
            ncrnas.setSelectedIndices(selected);
            selected = new int[opts.otherannots.size()];
            for (int i = 0; i < opts.otherannots.size(); i++) {
                selected[i] = otherfeatsmodel.indexOf(opts.otherannots.get(i));
            }
            otherfeats.setSelectedIndices(selected);
            
            if (opts.position != null &&
                !opts.position.equals("")) {
                position.setText(opts.position);
            } else {
                if (opts.chrom != null &&
                    !opts.chrom.equals("")) {
                    position.setText(opts.chrom + ":" + opts.start + "-" + opts.stop);
                }
            }
            seqSelect.addToSelected(opts.seqExpts);
            filetracks.fill(opts.regionTracks);
        }
        
        //Update the experiment list
        updateExptSelection();
    }

    
    public void setSpeciesGenomeDefaults() {
        String species = null, genome = null;
        try {
            ResourceBundle res = ResourceBundle.getBundle("defaultgenome");
            species = res.getString("species");
            genome = res.getString("genome");
        } catch (Exception e) {
            // don't do anything.  If it fails, then we just fall back to the defaults.
        }
        if (species == null || genome == null) {
            try {
                String homedir = System.getenv("HOME");
                String basename = "defaultSpeciesGenome";
                String fname = homedir + "/." + basename;
                File propfile = new File(fname);
                if (!(propfile.exists() && propfile.canRead())) {
                    homedir = System.getProperty("user.dir");
                    fname = homedir + "/" + basename;
                    propfile = new File(fname);
                }
                if (propfile.exists() && propfile.canRead()) {
                    InputStream is = new FileInputStream(propfile);
                    BufferedReader reader = new BufferedReader(new InputStreamReader(is));        
                    String line;
                    while ((line = reader.readLine()) != null) {
                        int p = line.indexOf('=');
                        String key = line.substring(0,p);
                        String value = line.substring(p+1);
                        if (key.equals("species")) {
                            species = value;
                        }
                        if (key.equals("genome")) {
                            genome = value;
                        }
                    }
                }
            } catch (Exception e) {
                // don't do anything.  If it fails, then we just fall back to the defaults.
            }

        }
        if (species == null || genome == null) {
            species = (String)this.speciesCBox.getSelectedItem();
            genome = (String)this.genomeCBox.getSelectedItem();
        }
        this.speciesCBox.setSelectedItem(species);
        this.genomeCBox.setSelectedItem(genome);
    }

    /* fills in and returns a SeqViewOptions object based on the current selections 
     */
    public SeqViewOptions parseOptions() {
        SeqViewOptions these = new SeqViewOptions();
        // parse the species and location tab
        String speciesStr = speciesCBox.getSelectedItem().toString();
        String genomeStr = genomeCBox.getSelectedItem().toString();
        try {
        	Species org = new Species(speciesStr);
        	these.setGenome(Genome.findGenome(genomeStr));
        } catch (NotFoundException e) {
            e.printStackTrace();
            throw new IllegalArgumentException(genomeStr);
        }
        these.position = position.getText();
        these.gene = gene.getText();

        // parse the options tab
        these.hash = hash.isSelected();
        these.seqletters = seqletters.isSelected();
        these.seqHistogramPainter = !oldchipseq.isSelected();

        // parse the annotations tab
        for (Object o : genes.getSelectedValues()) {
            these.genes.add(o.toString());
        }
        for (Object o : ncrnas.getSelectedValues()) {
            these.ncrnas.add(o.toString());
        }
        for (Object o : otherfeats.getSelectedValues()) {
            these.otherannots.add(o.toString());
        }
        these.gccontent = gccontent.isSelected();
        these.pyrpurcontent = pyrpurcontent.isSelected();
        these.polya = polyA.isSelected();
        these.cpg = cpg.isSelected();
        these.regexmatcher = regexmatcher.isSelected();

        // parse the expression tab
        //Collection<Experiment> expts = exprSelect.getSelected();
        //these.exprExperiments.addAll(expts);
        
        //Parse sequencing experiments
        for(SeqLocatorMatchedAligns lme : seqSelect.getSelected()) {
        	these.seqExpts.add(lme);
        }
        
       
        
        filetracks.parse(these.regionTracks);
        these.regexes = regexes;
        return these;
    }
    
    public SeqViewOptions parseAndDiff() {
        SeqViewOptions these = parseOptions();
        // need to see if we have existing options and if they're compatible.
        // if they are, return the difference.  Otherwise, return the complete
        // options.
        if (createdFrom != null &&
            these.getGenome().equals(createdFrom.getGenome())) {
            these.differenceOf(createdFrom);
        }
        return these;
    }
    
    
    /* updates the choice of experiments based on the
       currently selected genome and species */
    private void updateExptSelection() {
        reloadGenome();       
        seqSelect.setGenome(currGen);

        // update the set of Gene annotations
        genesmodel.clear();
        otherfeatsmodel.clear();

        if(currGen != null) { 
            for(String type : gfLoader.getTypes(currGen)) {
                genesmodel.addElement(type);
            }
            for(String type : annotLoader.getTypes(currGen)) {
                otherfeatsmodel.addElement(type);
            }
            List<String> chroms = currGen.getChromList();
            if (chroms.size() == 0) {
                throw new RuntimeException("Empty chromosome list for " + currGen);
            }
            java.util.Collections.sort(chroms);

            position.setText(chroms.get(0) + ":10000-20000");
        }
    }

    private void updateGenomeSelection () {
        try {
            Species org = new Species(speciesCBox.getSelectedItem().toString());
            Collection<String> genomes = org.getGenomeNames();
            genomeCBox.removeAllItems();
            for (String o : genomes) {
                genomeCBox.addItem(o);
            }
            genomeCBox.setSelectedIndex(0);                
        } catch (NotFoundException ex) {
            System.err.println("Couldn't find species " + speciesCBox.getSelectedItem());
            ex.printStackTrace();
        }
    }

    public void itemStateChanged(ItemEvent e) {
    	if (e.getStateChange() == ItemEvent.SELECTED) {    	
	        if (handlingChange) {return;}
	        Object source = e.getItemSelectable();
	        if (source == speciesCBox) {
	            updateGenomeSelection();
	        }
	        if (source == genomeCBox) {
	            synchronized(this) {
	                if (!handlingChange) {
	                    handlingChange = true;
	                    updateExptSelection();
	                    handlingChange = false;
	                }
	            }
	        }
    	}
    }
    
    public void actionPerformed (ActionEvent e) {

    }
}
