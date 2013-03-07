/*
 * Created on Jan 11, 2008
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.psu.compbio.seqcode.gse.seqview.components;

import java.util.*;
import java.util.regex.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.sql.*;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.viz.components.GenericSelectPanel;

public class SeqSelectPanel extends GenericSelectPanel<SeqLocator> {
    
    private SeqDataLoader seqLoader;
    private TreeSet<SeqLocator> locators;
    private ArrayList<SeqAlignment> alignments;
    private JTextField regex;
    private ChipSeqTableModel selectedModel, filteredModel;

    public SeqSelectPanel() { 
        try {
            seqLoader = new SeqDataLoader(true);
        } catch (Exception e) {
            e.printStackTrace();
            seqLoader = null;
        }
        locators = new TreeSet<SeqLocator>();
        alignments = new ArrayList<SeqAlignment>();
        selectedModel = new ChipSeqTableModel();
        filteredModel = new ChipSeqTableModel();
        init(filteredModel,selectedModel);
    }
    public JPanel getInputsPanel() {
        JPanel inputPanel = new JPanel(); inputPanel.setLayout(new BorderLayout());
        inputPanel.setLayout(new BorderLayout());
        regex = new JTextField();
        inputPanel.add(new JLabel("pattern to filter alignments"), BorderLayout.WEST);
        inputPanel.add(regex, BorderLayout.CENTER);        
        return inputPanel;
    }

    /* this is different than CollapseLocatorsByName because it
       keys on the name and alignment
    */
    public Collection<SeqLocator> getFilteredForSelected() {
        Map<Pair<String,String>, Set<String>> experiments = new HashMap<Pair<String,String>, Set<String>>();
        for (SeqLocator l : super.getFilteredForSelected()) {
            Pair<String,String> key = new Pair<String,String>(l.getExptName(), l.getAlignName());
            if (!experiments.containsKey(key)) {
                experiments.put(key, new HashSet<String>());
            }
            Set<String> reps = experiments.get(key);
            reps.addAll(l.getReplicates());
        }
        ArrayList<SeqLocator> output = new ArrayList<SeqLocator>();
        for (Pair<String,String> nv : experiments.keySet()) {
            output.add(new SeqLocator(nv.car(),
                                          experiments.get(nv),
                                          nv.cdr()));
        }
        return output;
    }

    public void retrieveData() {
        try {
            synchronized(locators) {
                locators.clear();
                alignments.clear();
                alignments.addAll(seqLoader.loadAlignments(getGenome()));
                for(SeqAlignment align : alignments) { 
                    locators.add(new SeqLocator(align.getExpt().getName(),
                                                       align.getExpt().getReplicate(),
                                                       align.getName()));
                }
            }
        } catch (SQLException e) {
            throw new RuntimeException(e.toString(), e);
        }
    }
    public void updateComponents() {
        selectedModel.clear();
        filteredModel.clear();
        synchronized(locators) {
            for (SeqLocator l : locators) {
                filteredModel.addObject(l);
            }
        }
    }
    public void filter() {
        String reg = regex.getText().trim();
        Pattern patt = null;
        if(reg != null && reg.length() > 0) {
            patt = Pattern.compile(reg);
        }
        synchronized(locators) {
            locators.clear();
            for (SeqAlignment align : alignments){
                if (patt == null || patt.matcher(align.toString()).find()) {
                    locators.add(new SeqLocator(align.getExpt().getName(),
                                                    align.getExpt().getReplicate(),
                                                    align.getName()));
                }
            }
            filteredModel.clear();
            for (SeqLocator l : locators) {
                filteredModel.addObject(l);
            }
        }
    }

    
    public static Collection<SeqLocator> collapseLocatorsByName(Collection<SeqLocator> locs) { 
        LinkedHashMap<String,Map<String,Set<String>>> map = 
            new LinkedHashMap<String,Map<String,Set<String>>>();
        
        for(SeqLocator loc : locs) { 
            String exptName = loc.getExptName();
            String alignName = loc.getAlignName();
            if(!map.containsKey(exptName)) { map.put(exptName, new LinkedHashMap<String,Set<String>>()); }
            if(!map.get(exptName).containsKey(alignName)) { map.get(exptName).put(alignName, new TreeSet<String>()); }
            map.get(exptName).get(alignName).addAll(loc.getReplicates());
        }
        
        LinkedList<SeqLocator> collapsed = new LinkedList<SeqLocator>();
        
        for(String exptName : map.keySet()) { 
            for(String alignName : map.get(exptName).keySet()) { 
                SeqLocator newloc = new SeqLocator(exptName, map.get(exptName).get(alignName), alignName);
                collapsed.add(newloc);
            }
        }
        
        return collapsed;
    }

    public void close() {
        super.close();
        if (seqLoader != null) {
            seqLoader.close();
        }
    }
}
