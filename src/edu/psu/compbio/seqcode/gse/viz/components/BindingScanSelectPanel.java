/*
 * Created on Nov 6, 2006
 */
package edu.psu.compbio.seqcode.gse.viz.components;

import java.util.*;
import java.sql.SQLException;
import java.util.regex.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.*;

import edu.psu.compbio.seqcode.gse.datasets.binding.*;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.*;
import edu.psu.compbio.seqcode.gse.datasets.general.CellLine;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptCondition;
import edu.psu.compbio.seqcode.gse.datasets.general.ExptTarget;
import edu.psu.compbio.seqcode.gse.datasets.general.MetadataLoader;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseException;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;

/**
 * @author tdanford
 */
public class BindingScanSelectPanel extends GenericSelectPanel<BindingScan> {
	
    private BindingScanLoader bindingloader;
    private MetadataLoader metaloader;
    private ChipChipMetadataLoader chipchiploader;
    
    private BindingScanTableModel selectedModel, filteredModel;
    private DefaultComboBoxModel cellsModel, condModel, factorModel;
    private JComboBox cellsBox, condBox, factorBox;
    private JTextField regex;
    private BindingScanFilter filter;
    private Wrapper<CellLine> noCells;
    private Wrapper<ExptCondition> noCond;
    private Wrapper<ExptTarget> noFactor;
    private Collection<CellLine> allCells;
    private Collection<ExptCondition> allConds;
    private Collection<ExptTarget> allFactors;
    private Collection<BindingScan> scans;

    public BindingScanSelectPanel() 
    	throws SQLException, UnknownRoleException {
        super();
        bindingloader = new BindingScanLoader();
        metaloader = new MetadataLoader();
        chipchiploader = new ChipChipMetadataLoader(metaloader);

        selectedModel = new BindingScanTableModel();        
        filteredModel = new BindingScanTableModel();
        noCells = new Wrapper<CellLine>("<NONE>", null);
        noCond = new Wrapper<ExptCondition>("<NONE>", null);
        noFactor = new Wrapper<ExptTarget>("<NONE>", null);	
        scans = new ArrayList<BindingScan>();
        init(filteredModel,selectedModel);
    }
    
    public JPanel getInputsPanel() {
        JPanel inputPanel = new JPanel(); inputPanel.setLayout(new BorderLayout());		
        JPanel boxPanel = new JPanel(); boxPanel.setLayout(new GridLayout(2, 1));
        JPanel cellsPanel = new JPanel(); cellsPanel.setLayout(new BorderLayout());
        JPanel condPanel = new JPanel(); condPanel.setLayout(new BorderLayout());
        JPanel regexPanel = new JPanel(); regexPanel.setLayout(new BorderLayout());
        JPanel factorPanel = new JPanel(); factorPanel.setLayout(new BorderLayout());
		
        
        JPanel comboPanel = new JPanel(); comboPanel.setLayout(new GridLayout(1, 3));
        comboPanel.add(cellsPanel);
        comboPanel.add(condPanel);
        comboPanel.add(factorPanel);
        
        boxPanel.add(comboPanel);
        boxPanel.add(regexPanel);
		
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
		
        cellsPanel.setBorder(new TitledBorder("Cells"));
        condPanel.setBorder(new TitledBorder("Condition"));
        factorPanel.setBorder(new TitledBorder("Factor"));
        regexPanel.setBorder(new TitledBorder("Match"));
		
        inputPanel.add(boxPanel, BorderLayout.CENTER);
        setBorder(new TitledBorder("Binding Scans:"));
        return inputPanel;
    }

    public BindingScanLoader getBindingLoader() { return bindingloader; }
    
    public void clearSelected() { selectedModel.clear(); }

    public void filter() {
        CellLine cells = ((Wrapper<CellLine>)(cellsModel.getSelectedItem())).value;
        ExptCondition cond = ((Wrapper<ExptCondition>)(condModel.getSelectedItem())).value;
        ExptTarget factor = ((Wrapper<ExptTarget>)(factorModel.getSelectedItem())).value;

        filter(cells,cond,factor);
        filteredModel.clear();
        for (BindingScan bs : scans) {
            filteredModel.addObject(bs);
        }
    }
    public void filter(CellLine cells, ExptCondition cond, ExptTarget factor) {
        String reg = regex.getText().trim();
        Pattern patt = null;
        if(reg.length() > 0) {
            patt = Pattern.compile(reg);
        }        
        try { 
            synchronized(scans) {
                filter = new BindingScanFilter(getGenome(),bindingloader);
                Collection<BindingScan> tempscans = filter.findScans(cells, cond, factor);        
                scans.clear();
                for(BindingScan bs : tempscans) {
                    String str = bs.toString();
                    Matcher m = patt != null ? patt.matcher(str) : null;
                    if(m == null || m.find()) {             
                        scans.add(bs);
                    }
                }
            }
        } catch(Exception se) { 
            se.printStackTrace(System.err);
        }        
    }

    public void retrieveData() {
        try {
            allCells = new TreeSet<CellLine>(chipchiploader.loadAllCells(getGenome()));
            allConds = new TreeSet<ExptCondition>(chipchiploader.loadAllConditions(getGenome()));
            allFactors = new TreeSet<ExptTarget>(chipchiploader.loadAllFactors(getGenome()));
            filter(null,null,null);
        } catch (SQLException e) {
            e.printStackTrace();
        }
        
    }
    public void updateComponents() {
        synchronized(scans) {
            filteredModel.clear();
            for (BindingScan bs : scans) {
                filteredModel.addObject(bs);
            }
        }
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
    }
    
    public void close() {
        if (!bindingloader.isClosed()) {
            bindingloader.close();
        }
        if (!metaloader.isClosed()) {
            metaloader.close();
        }
    }
    
    public boolean isClosed() { 
    	return bindingloader.isClosed() && metaloader.isClosed();
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
