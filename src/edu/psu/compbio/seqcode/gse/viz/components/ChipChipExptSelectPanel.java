/*
 * Created on Nov 6, 2006
 */
package edu.psu.compbio.seqcode.gse.viz.components;

import java.util.*;
import java.util.regex.*;
import java.sql.SQLException;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.*;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.*;
import edu.psu.compbio.seqcode.gse.datasets.general.CellLine;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptCondition;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptTarget;
import edu.psu.compbio.seqcode.gse.datasets.general.MetadataLoader;
import edu.psu.compbio.seqcode.gse.datasets.locators.*;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;


/**
 * @author tdanford
 */
public class ChipChipExptSelectPanel extends GenericSelectPanel<ExptLocator> {
    
    private MetadataLoader loader;
    private ChipChipMetadataLoader chipLoader;
    
    private ExptLocatorFilter filter;
    private DefaultComboBoxModel cellsModel, condModel, factorModel;
    private JComboBox cellsBox, condBox, factorBox;
    private JTextField regex;
    private boolean hasAgilent, hasReps, hasMSP, hasBayes;
    private JRadioButton agilentButton, agilentReplicateButton, mspButton, bayesButton;
    private ButtonGroup typeGroup;
    private Wrapper<CellLine> noCells;
    private Wrapper<ExptCondition> noCond;
    private Wrapper<ExptTarget> noFactor;
    private ExptTableModel selectedModel, filteredModel;
    private Collection<CellLine> allCells;
    private Collection<ExptCondition> allConds;
    private Collection<ExptTarget> allFactors;

    private Collection<ExptLocator> scans;    
    private JButton addDiffButton;

    public ChipChipExptSelectPanel(Genome g) {
        super();
        hasAgilent = true;
        hasReps = true;
        hasMSP = true;
        hasBayes = true;
        init(g);       
    }

    /* this constructor lets you determine which experiment type selection buttons are shown.  The choices
       are experiment, experiment w/ replicate, MSP, and Bayes
    */
    public ChipChipExptSelectPanel(Genome g,
                           boolean hasExperimentButton,
                           boolean hasReplicateButton,
                           boolean hasMSPButton,
                           boolean hasBayesButton) {
        super();
        hasAgilent = hasExperimentButton;
        hasReps = hasReplicateButton;
        hasMSP = hasMSPButton;
        hasBayes = hasBayesButton;
        init(g);       
    }

    public void init (Genome g) {
        scans = new TreeSet<ExptLocator>();
        try { 
            loader = new MetadataLoader();
            chipLoader = new ChipChipMetadataLoader(loader);
        } catch(Exception e) { 
            e.printStackTrace(System.err);
            throw new RuntimeException(e.getMessage(), e);
        }
        filter = new ExptLocatorFilter();
        if (g != null) {filter.setGenome(g);}
        noCells = new Wrapper<CellLine>("<NONE>", null);
        noCond = new Wrapper<ExptCondition>("<NONE>", null);
        noFactor = new Wrapper<ExptTarget>("<NONE>", null);		
        setBorder(new TitledBorder("Selected Experiments:"));
        selectedModel = new ExptTableModel();
        filteredModel = new ExptTableModel();
        
        super.init(filteredModel, selectedModel);

        buttonPanel.add(addDiffButton = new JButton("AddDiff"));
        addDiffButton.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) { 
                addDiff();
            }
        });
        setGenome(g);

    }
    
    public void addDiff() { 
        int[] inds = filteredList.getSelectedRows();
        if(inds.length==2) { 
            ExptLocator e1 = filteredModel.getObject(inds[0]);
            ExptLocator e2 = filteredModel.getObject(inds[1]);
            
            if(e1 instanceof ChipChipLocator && e2 instanceof ChipChipLocator) {
                ChipChipLocator al1 = (ChipChipLocator)e1;
                ChipChipLocator al2 = (ChipChipLocator)e2;

                ChipChipLocator dl = new ChipChipDifferenceLocator(getGenome(), al1, al2);
                if(!selectedModel.contains(dl)) { 
                    selectedModel.addObject(dl);
                }
            }
        }
    }
    
    public JPanel getInputsPanel() {
        JPanel inputPanel = new JPanel(); inputPanel.setLayout(new BorderLayout());

        JPanel boxPanel = new JPanel(); boxPanel.setLayout(new GridLayout(3, 1));
        
        JPanel buttonPanel = new JPanel(); buttonPanel.setLayout(new GridLayout(1, 3));
        JPanel cellsPanel = new JPanel(); cellsPanel.setLayout(new BorderLayout());
        JPanel condPanel = new JPanel(); condPanel.setLayout(new BorderLayout());
        JPanel regexPanel = new JPanel(); regexPanel.setLayout(new BorderLayout());
        JPanel factorPanel = new JPanel(); factorPanel.setLayout(new BorderLayout());
        
        JPanel comboPanel = new JPanel(); comboPanel.setLayout(new GridLayout(1, 3));
        comboPanel.add(cellsPanel);
        comboPanel.add(condPanel);
        comboPanel.add(factorPanel);
		
        boxPanel.add(buttonPanel);
        boxPanel.add(comboPanel);
        boxPanel.add(regexPanel);
        
        agilentButton = new JRadioButton("ChipChip");
        agilentReplicateButton = new JRadioButton("ChipChip w/ Replicate");
        mspButton = new JRadioButton("Rosetta");
        bayesButton = new JRadioButton("Bayes");
        typeGroup = new ButtonGroup();
        boolean selected = false;
        if (hasAgilent) {
            typeGroup.add(agilentButton);
            agilentButton.setSelected(true);
            selected = true;
            buttonPanel.add(agilentButton);
        }
        if (hasReps) {
            typeGroup.add(agilentReplicateButton);
            buttonPanel.add(agilentReplicateButton);
            if (!selected) {
                agilentReplicateButton.setSelected(true);
                selected = true;
            }
        }
        if (hasMSP) {
            typeGroup.add(mspButton);
            buttonPanel.add(mspButton);
            if (!selected) {
                mspButton.setSelected(true);
                selected = true;
            }
        }
        if (hasBayes) {
            typeGroup.add(bayesButton);
            buttonPanel.add(bayesButton);
            if (!selected) {
                bayesButton.setSelected(true);
            }

        }
		
        cellsModel = new DefaultComboBoxModel();
        condModel = new DefaultComboBoxModel();
        factorModel = new DefaultComboBoxModel();
        regex = new JTextField();
		
        cellsBox = new JComboBox(cellsModel);
        condBox = new JComboBox(condModel);
        factorBox = new JComboBox(factorModel);       
		
        cellsPanel.add(cellsBox, BorderLayout.CENTER);
        condPanel.add(condBox, BorderLayout.CENTER);
        factorPanel.add(factorBox, BorderLayout.CENTER);
        regexPanel.add(regex, BorderLayout.CENTER);
		
        buttonPanel.setBorder(new TitledBorder("Experiment Type"));
        cellsPanel.setBorder(new TitledBorder("Cells"));
        condPanel.setBorder(new TitledBorder("Condition"));
        factorPanel.setBorder(new TitledBorder("Factor"));
        regexPanel.setBorder(new TitledBorder("Match"));
		
        inputPanel.add(boxPanel, BorderLayout.CENTER);
        return inputPanel;
    }
    /* this filter() only gets called from the swing thread */
    public void filter() {
        CellLine cells = ((Wrapper<CellLine>)(cellsModel.getSelectedItem())).value;
        ExptCondition cond = ((Wrapper<ExptCondition>)(condModel.getSelectedItem())).value;
        ExptTarget factor = ((Wrapper<ExptTarget>)(factorModel.getSelectedItem())).value;
        filter(cells,cond,factor);

        filteredModel.clear();
        synchronized(scans) {
            for(ExptLocator bs : scans) {
                filteredModel.addObject(bs);
            }
        }
    }
    /* this gets called from anywhere and just updates scans */
    public void filter(CellLine cells, ExptCondition cond, ExptTarget factor) {
        synchronized(scans) {
            scans.clear();
            try { 
                ExptLocatorFilterOptions opts = new ExptLocatorFilterOptions();
                // use this as the default in case nothing is set yet
                opts.setChipChipType();
                if(mspButton.isSelected()) { opts.setMSPType(); }
                if(bayesButton.isSelected()) { opts.setBayesType(); }
                if(cells != null) { opts.setCells(cells); } 
                if(cond != null) { opts.setCondition(cond); }
                if(factor != null) { opts.setFactor(factor); }
            
                Collection<ExptLocator> scanstemp = filter.findLocators(opts);
                String reg = regex.getText().trim();
                Pattern patt = null;
                if(reg != null && reg.length() > 0) {
                    patt = Pattern.compile(reg);
                }
                for(ExptLocator bs : scanstemp) {
                    if (agilentButton.isSelected() && bs instanceof ChipChipLocator) {
                        ((ChipChipLocator)bs).clearReplicate();
                    }
                    String str = bs.toString();
                    Matcher m = patt != null ? patt.matcher(str) : null;
                    if(m == null || m.find()) {                  
                        scans.add(bs);
                    }
                }
            } catch(SQLException se) { 
                se.printStackTrace(System.err);
            }        
        }
    }
    
    public void retrieveData() {
        try {
            filter.setGenome(getGenome());
            allCells = new TreeSet<CellLine>();
            allConds = new TreeSet<ExptCondition>();
            allFactors = new TreeSet<ExptTarget>();
            allCells.addAll(chipLoader.loadAllCells(getGenome()));
            allConds.addAll(chipLoader.loadAllConditions(getGenome()));
            allFactors.addAll(chipLoader.loadAllFactors(getGenome()));
            filter(null,null,null);
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    public void updateComponents() {
        synchronized(scans) {
            selectedModel.clear();
            if(cellsModel.getSize() > 0) { cellsModel.removeAllElements(); } 
            if(condModel.getSize() > 0) { condModel.removeAllElements(); } 
            if(factorModel.getSize() > 0) { factorModel.removeAllElements(); }
		
            cellsModel.addElement(noCells);
            condModel.addElement(noCond);
            factorModel.addElement(noFactor);
            for(CellLine cells : allCells) { 
                Wrapper<CellLine> wrapper = new Wrapper<CellLine>(cells.getName(), cells);
                cellsModel.addElement(wrapper);
            }
		
            for(ExptCondition cond : allConds) { 
                Wrapper<ExptCondition> wrapper = 
                    new Wrapper<ExptCondition>(cond.getName(), cond);
                condModel.addElement(wrapper);
            }

            for(ExptTarget factor : allFactors){ 
                Wrapper<ExptTarget> wrapper = 
                    new Wrapper<ExptTarget>(factor.getName(), factor);
                factorModel.addElement(wrapper);
            }
            cellsModel.setSelectedItem(noCells);
            condModel.setSelectedItem(noCond);
            factorModel.setSelectedItem(noFactor);
		
            for(ExptLocator bs : scans) {
                filteredModel.addObject(bs);
            }
            if (agilentButton.isSelected()) {
                for(ExptLocator bs : scans) {
                    if (bs instanceof ChipChipLocator) {
                        ((ChipChipLocator)bs).clearReplicate();
                    }
                }
            }
        }
    }
        
    public Collection<ExptLocator> getSelected() { 
        LinkedList<ExptLocator> scans = new LinkedList<ExptLocator>();
        for(int i = 0; i < selectedModel.getSize(); i++) { 
            scans.addLast(selectedModel.getObject(i));
        }
        return scans;
    }
    
    public void close() {
        if (!loader.isClosed()) {
            loader.close();
        }
        if (!filter.isClosed()) {
            filter.close();
        }
    }
    
    public boolean isClosed() { 
        return loader.isClosed() && filter.isClosed();
    }
    private static class Wrapper<X> { 
        public X value;
        public String name;
        public Wrapper(String n, X v) { name = n; value = v; }
        public String toString() { return name; }
        public int hashCode() { return name.hashCode(); }
        public boolean equals(Object o) { 
            if(!(o instanceof Wrapper)) { return false; }
            Wrapper w = (Wrapper)o;
            return w.name.equals(name);
        }
    }

}
