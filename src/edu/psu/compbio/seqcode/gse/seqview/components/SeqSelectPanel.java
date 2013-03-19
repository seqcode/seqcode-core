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
import java.util.List;
import java.util.regex.*;
import java.awt.*;
import javax.swing.*;
import javax.swing.border.Border;

import java.sql.*;

import edu.psu.compbio.seqcode.gse.datasets.general.ExptType;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.viz.components.GenericSelectPanel;

public class SeqSelectPanel extends GenericSelectPanel<SeqLocatorMatchedExpt> {
    
	//Pass this in from SeqView later
    private SeqDataLoader seqLoader;
    
    private TreeSet<SeqLocatorMatchedExpt> lme;
    private ArrayList<SeqAlignment> alignments;
    private JComboBox jcbType;
    private JTextField regexLab, regexCond, regexTarget, regexCell, regexAlign, regexRep;
    private SeqTableModel selectedModel, filteredModel;

    public SeqSelectPanel() { 
        try {
            seqLoader = new SeqDataLoader(true);
        } catch (Exception e) {
            e.printStackTrace();
            seqLoader = null;
        }
        lme = new TreeSet<SeqLocatorMatchedExpt>();
        alignments = new ArrayList<SeqAlignment>();
        selectedModel = new SeqTableModel();
        filteredModel = new SeqTableModel();
        init(filteredModel,selectedModel);
    }
    public JPanel getInputsPanel() {
    	JPanel inputPanel = new JPanel();
		inputPanel.setLayout(new GridLayout(2,7));
		
		try {
    		ArrayList<String> types = new ArrayList<String>();
    		types.add("");
        	for(ExptType e : seqLoader.getExptTypes())
				types.add(e.getName());
			Collections.sort(types);
			jcbType = new JComboBox(types.toArray()); jcbType.setEditable(true);
			
	        //regexType = new JTextField();
	        regexLab = new JTextField();
	        regexCond = new JTextField(); 
	        regexTarget = new JTextField(); 
	        regexCell = new JTextField(); 
	        regexAlign = new JTextField(); 
	        regexRep = new JTextField();
	        
	        Font labelFont = new Font("SansSerif", Font.PLAIN, 10);
	        Border paddingBorder = BorderFactory.createEmptyBorder(5,5,5,5);
	        JLabel labelExpt = new JLabel("ExptType"); labelExpt.setFont(labelFont); labelExpt.setBorder(paddingBorder);
	        JLabel labelLab = new JLabel("Lab"); labelLab.setFont(labelFont); labelLab.setBorder(paddingBorder);
	        JLabel labelCond = new JLabel("ExptCondition"); labelCond.setFont(labelFont); labelCond.setBorder(paddingBorder); 
	        JLabel labelTarget = new JLabel("ExptTarget"); labelTarget.setFont(labelFont); labelTarget.setBorder(paddingBorder);
	        JLabel labelCell = new JLabel("CellLine"); labelCell.setFont(labelFont); labelCell.setBorder(paddingBorder);
	        JLabel labelAlign = new JLabel("Alignment"); labelAlign.setFont(labelFont); labelAlign.setBorder(paddingBorder);
	        JLabel labelRep = new JLabel("Replicate"); labelRep.setFont(labelFont); labelRep.setBorder(paddingBorder);
	        inputPanel.add(labelExpt); 
	        inputPanel.add(labelLab);
	        inputPanel.add(labelCond);
	        inputPanel.add(labelTarget);
	        inputPanel.add(labelCell);
	        inputPanel.add(labelAlign);
	        inputPanel.add(labelRep);
	        //inputPanel.add(regexType);
	        inputPanel.add(jcbType);
	        inputPanel.add(regexLab);
	        inputPanel.add(regexCond);
	        inputPanel.add(regexTarget);
	        inputPanel.add(regexCell);
	        inputPanel.add(regexAlign);
	        inputPanel.add(regexRep);
        
        } catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
        return inputPanel;
    }

    /* this is different than CollapseLocatorsByName because it
       keys on the SeqExpt
    */
    public Collection<SeqLocatorMatchedExpt> getFilteredForSelected() {
        Map<Pair<String,String>, Set<String>> reps = new HashMap<Pair<String,String>, Set<String>>();
        Map<Pair<String,String>, List<SeqExpt>> expts = new HashMap<Pair<String,String>, List<SeqExpt>>();
        for (SeqLocatorMatchedExpt le : super.getFilteredForSelected()) {
        	List<SeqExpt> e = le.expts;
        	SeqLocator l = le.locator;
            Pair<String,String> key = new Pair<String,String>(l.getExptName(), l.getAlignName());
            if (!reps.containsKey(key)) {
                reps.put(key, new HashSet<String>());
                expts.put(key, new ArrayList<SeqExpt>());
            }
            if(reps.containsKey(key)){
            	//TODO: we should really check that all the SeqExpts are consistent, since
            	//only the first one is used to print expt information. 
            	reps.get(key).addAll(l.getReplicates());
            	expts.get(key).addAll(e);
            }
        }
        ArrayList<SeqLocatorMatchedExpt> output = new ArrayList<SeqLocatorMatchedExpt>();
        for (Pair<String,String> nv : reps.keySet()) {
            output.add(new SeqLocatorMatchedExpt(expts.get(nv),
            		new SeqLocator(nv.car(), reps.get(nv), nv.cdr())));
        }
        return output;
    }

    public void retrieveData() {
        try {
            synchronized(lme) {
                lme.clear();
                alignments.clear();
                alignments.addAll(seqLoader.loadAlignments(getGenome()));
                for(SeqAlignment align : alignments) { 
                	List<SeqExpt> e = new ArrayList<SeqExpt>();
                	e.add(align.getExpt());
                    lme.add(new SeqLocatorMatchedExpt(e,
                    		new SeqLocator(align.getExpt().getName(), align.getExpt().getReplicate(), align.getName())));
                }
            }
        } catch (SQLException e) {
            throw new RuntimeException(e.toString(), e);
        }
    }
    public void updateComponents() {
        selectedModel.clear();
        filteredModel.clear();
        synchronized(lme) {
            for (SeqLocatorMatchedExpt l : lme) {
                filteredModel.addObject(l);
            }
        }
    }
    public void filter() {
    	Pattern pattType=null, pattLab=null, pattCond=null, pattTarget=null, pattCell=null, pattAlign=null, pattRep=null;
    	//String regType = regexType.getText().trim();
    	String regType = (String) jcbType.getSelectedItem();
        if(regType != null && regType.length() > 0)
        	pattType = Pattern.compile(regType);
        String regLab = regexLab.getText().trim();
        if(regLab != null && regLab.length() > 0)
        	pattLab = Pattern.compile(regLab);
        String regCond = regexCond.getText().trim();
        if(regCond != null && regCond.length() > 0)
        	pattCond = Pattern.compile(regCond);
        String regTarget = regexTarget.getText().trim();
        if(regTarget != null && regTarget.length() > 0)
        	pattTarget = Pattern.compile(regTarget);
        String regCell = regexCell.getText().trim();
        if(regCell != null && regCell.length() > 0)
        	pattCell = Pattern.compile(regCell);
        String regAlign = regexAlign.getText().trim();
        if(regAlign != null && regAlign.length() > 0)
        	pattAlign = Pattern.compile(regAlign);
        String regRep = regexRep.getText().trim();
        if(regRep != null && regRep.length() > 0)
        	pattRep = Pattern.compile(regRep);
        
        synchronized(lme) {
            lme.clear();
            for (SeqAlignment align : alignments){
                if( (pattType == null || pattType.matcher(align.getExpt().getExptType().getName()).find()) &&
                	(pattLab == null || pattLab.matcher(align.getExpt().getLab().getName()).find()) &&
                	(pattCond == null || pattCond.matcher(align.getExpt().getExptCondition().getName()).find()) &&
                	(pattTarget == null || pattTarget.matcher(align.getExpt().getExptTarget().getName()).find()) &&
                	(pattCell == null || pattCell.matcher(align.getExpt().getCellLine().getName()).find()) &&
                	(pattAlign == null || pattAlign.matcher(align.getName()).find()) &&
                	(pattRep == null || pattRep.matcher(align.getExpt().getReplicate()).find())
                		) {
                	List<SeqExpt> e = new ArrayList<SeqExpt>();
                	e.add(align.getExpt());
                    lme.add(new SeqLocatorMatchedExpt(e,
                    		new SeqLocator(align.getExpt().getName(),
                                                    align.getExpt().getReplicate(),
                                                    align.getName())));
                }
            }
            filteredModel.clear();
            for (SeqLocatorMatchedExpt l : lme) {
                filteredModel.addObject(l);
            }
        }
    }

    /*
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
    }*/

    public void close() {
        super.close();
        if (seqLoader != null) {
            seqLoader.close();
        }
    }
}
