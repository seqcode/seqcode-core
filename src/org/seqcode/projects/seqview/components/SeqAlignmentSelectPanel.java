package org.seqcode.projects.seqview.components;

import java.util.*;
import java.util.List;
import java.util.regex.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.*;
import javax.swing.border.Border;

import org.seqcode.data.core.ExptType;
import org.seqcode.data.seqdata.*;
import org.seqcode.data.seqdata.SeqLocatorMatchedAligns.SeqLMENotReplicatesException;
import org.seqcode.gseutils.Pair;
import org.seqcode.viz.components.GenericSelectPanel;

import java.sql.*;


public class SeqAlignmentSelectPanel extends GenericSelectPanel<SeqLocatorMatchedAligns> {
    
	private SeqDataLoader seqLoader;
    
    private TreeSet<SeqLocatorMatchedAligns> lmes;
    private ArrayList<SeqAlignment> alignments;
    private JComboBox jcbType;
    private JTextField regexLab, regexCond, regexTarget, regexCell, regexAlign, regexRep;
    private SeqAlignmentTableModel selectedModel, filteredModel;

    public SeqAlignmentSelectPanel() { 
        try {
            seqLoader = new SeqDataLoader(false, true);
        } catch (Exception e) {
            e.printStackTrace();
            seqLoader = null;
        }
        lmes = new TreeSet<SeqLocatorMatchedAligns>();
        alignments = new ArrayList<SeqAlignment>();
        selectedModel = new SeqAlignmentTableModel();
        filteredModel = new SeqAlignmentTableModel();
        init(filteredModel,selectedModel);
    }
    public JPanel getInputsPanel() {
    	JPanel inputPanel = new JPanel();
		inputPanel.setLayout(new GridLayout(2,7));
		
		try {
    		ArrayList<String> types = new ArrayList<String>();
    		types.add("");
        	for(ExptType e : seqLoader.getMetadataLoader().loadAllExptTypes(false))
				types.add(e.getName());
			Collections.sort(types);
			jcbType = new JComboBox(types.toArray()); jcbType.setEditable(true);
			jcbType.addActionListener(new ActionListener(){public void actionPerformed(ActionEvent e){filter();}});
			
	        regexLab = new JTextField(); regexLab.addActionListener(new ActionListener(){public void actionPerformed(ActionEvent e){filter();}});
	        regexCond = new JTextField(); regexCond.addActionListener(new ActionListener(){public void actionPerformed(ActionEvent e){filter();}});
	        regexTarget = new JTextField(); regexTarget.addActionListener(new ActionListener(){public void actionPerformed(ActionEvent e){filter();}});
	        regexCell = new JTextField(); regexCell.addActionListener(new ActionListener(){public void actionPerformed(ActionEvent e){filter();}});
	        regexAlign = new JTextField(); regexAlign.addActionListener(new ActionListener(){public void actionPerformed(ActionEvent e){filter();}});
	        regexRep = new JTextField(); regexRep.addActionListener(new ActionListener(){public void actionPerformed(ActionEvent e){filter();}});
	        
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
	        inputPanel.add(jcbType);
	        inputPanel.add(regexLab);
	        inputPanel.add(regexCond);
	        inputPanel.add(regexTarget);
	        inputPanel.add(regexCell);
	        inputPanel.add(regexAlign);
	        inputPanel.add(regexRep);
        
        } catch (SQLException e) {
			e.printStackTrace();
		}
        
        return inputPanel;
    }

    public Collection<SeqLocatorMatchedAligns> getFilteredForSelected() {
        Map<Pair<String,String>, Set<String>> reps = new HashMap<Pair<String,String>, Set<String>>();
        Map<Pair<String,String>, List<SeqAlignment>> aligns = new HashMap<Pair<String,String>, List<SeqAlignment>>();
        for (SeqLocatorMatchedAligns lme : super.getFilteredForSelected()) {
        	Collection<SeqAlignment> a = lme.aligns;
            SeqLocator l = lme.locator;
            Pair<String,String> key = new Pair<String,String>(l.getExptName(), l.getAlignName());
            if (!reps.containsKey(key)) {
	            reps.put(key, new HashSet<String>());
	            aligns.put(key, new ArrayList<SeqAlignment>());
	        }
	        if(reps.containsKey(key)){
	            reps.get(key).addAll(l.getReplicates());
	            aligns.get(key).addAll(a);
	        }
        }
        ArrayList<SeqLocatorMatchedAligns> output = new ArrayList<SeqLocatorMatchedAligns>();
        for (Pair<String,String> nv : reps.keySet()) {
            try {
				output.add(new SeqLocatorMatchedAligns(aligns.get(nv),
				            new SeqLocator(nv.car(), reps.get(nv), nv.cdr())));
			} catch (SeqLMENotReplicatesException e) {
				System.err.println(e.getMessage());
				e.printStackTrace();
			}
        }
        return output;
    }

    public void retrieveData() {
    	try {
            synchronized(lmes) {
                lmes.clear();
                alignments.clear();
                alignments.addAll(seqLoader.loadAllAlignments(getGenome()));
                for(SeqAlignment align : alignments) {
                	List<SeqAlignment> a = new ArrayList<SeqAlignment>();
                	a.add(align);
                    lmes.add(new SeqLocatorMatchedAligns(a,
                                new SeqLocator(align.getExpt().getName(), align.getExpt().getReplicate(), align.getName())));
                }
            }
        } catch (SQLException e) {
            throw new RuntimeException(e.toString(), e);
        } catch (SeqLMENotReplicatesException e1) {
        	System.err.println(e1.getMessage());
			e1.printStackTrace();
		}
    }
    public void updateComponents() {
        selectedModel.clear();
        filteredModel.clear();
        synchronized(lmes) {
            for (SeqLocatorMatchedAligns l : lmes) {
                filteredModel.addObject(l);
            }
        }
    }
    public void info() {
    	int[] inds = filteredList.getSelectedRows();
        for (int i = 0; i < inds.length; i++) {
        	SeqLocatorMatchedAligns lma = filteredModel.getObject(inds[i]);
        	
			Collection<SeqAlignment> aligns = lma.aligns;
			for(SeqAlignment align : aligns){
				SeqExpt expt = align.getExpt();
				String info=" "+expt.getName()+";"+expt.getReplicate()+";"+align.getName();
				info = info+ "\n\n SeqExpt info:\n";
				info = info+ "\tID: "+expt.getDBID()+"\n";
				info = info+ "\tName: "+expt.getName()+"\n";
				info = info+ "\tRep: "+expt.getReplicate()+"\n";
				info = info+ "\tSpecies: "+expt.getOrganism().getName()+"\n";
				info = info+ "\tLab: "+expt.getLab().getName()+"\n";
				info = info+ "\tExptType: "+expt.getExptType().getName()+"\n";
				info = info+ "\tExptCondition: "+expt.getExptCondition().getName()+"\n";
				info = info+ "\tExptTarget: "+expt.getExptTarget().getName()+"\n";
				info = info+ "\tExptCellLine: "+expt.getCellLine().getName()+"\n";
				info = info+ "\tReadType: "+expt.getReadType().getName()+"\n";
				info = info+ "\tReadLen: "+expt.getReadLength()+"\n";
				info = info+ "\tNumRead: "+expt.getNumRead()+"\n";
				info = info+ "\tCollabID: "+expt.getCollabID()+"\n";
				info = info+ "\tPublicSource: "+expt.getPublicSource()+"\n";
				info = info+ "\tPublicDBID: "+expt.getPublicDBID()+"\n";
				info = info+ "\tFQFile: "+expt.getFQFile()+"\n";
				info = info+ "\tExptNote: "+expt.getExptNote()+"\n";
	            info = info+ " SeqAlignment info:\n";
				info = info+ "\tID: "+align.getDBID()+"\n";
				info = info+ "\tName: "+align.getName()+"\n";
				info = info+ "\tGenome: "+align.getGenome().getVersion()+"\n";
				info = info+ "\tAlignType: "+align.getAlignType().getName()+"\n";
				info = info+ "\tPermissions: "+align.getPermissions()+"\n";
				info = info+ "\tNumHits: "+align.getNumHits()+"\n";
				info = info+ "\tTotalWeight: "+align.getTotalWeight()+"\n";
				info = info+ "\tAlignDir: "+align.getAlignDir()+"\n";
				info = info+ "\tAlignFile: "+align.getAlignFile()+"\n";
				info = info+ "\tIDXFile: "+align.getIDXFile()+"\n";
				info = info+ "\tCollabAlignID: "+align.getCollabAlignID()+"\n";
				
				JFrame frame = new JFrame();
				JTextArea textArea = new JTextArea(info);
				frame.add(textArea);
				frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
				frame.pack();
				frame.setVisible(true);
			}
        }
    }
    public void filter() {
    	Pattern pattType=null, pattLab=null, pattCond=null, pattTarget=null, pattCell=null, pattAlign=null, pattRep=null;
    	//String regType = regexType.getText().trim();
    	String regType = (String) jcbType.getSelectedItem(); 
    	regType=regType.toLowerCase();
        if(regType != null && regType.length() > 0)
        	pattType = Pattern.compile(regType);
        String regLab = regexLab.getText().trim();
        regLab = regLab.toLowerCase();
        if(regLab != null && regLab.length() > 0)
        	pattLab = Pattern.compile(regLab);
        String regCond = regexCond.getText().trim();
        regCond = regCond.toLowerCase();
        if(regCond != null && regCond.length() > 0)
        	pattCond = Pattern.compile(regCond);
        String regTarget = regexTarget.getText().trim();
        regTarget = regTarget.toLowerCase();
        if(regTarget != null && regTarget.length() > 0)
        	pattTarget = Pattern.compile(regTarget);
        String regCell = regexCell.getText().trim();
        regCell = regCell.toLowerCase();
        if(regCell != null && regCell.length() > 0)
        	pattCell = Pattern.compile(regCell);
        String regAlign = regexAlign.getText().trim();
        regAlign=regAlign.toLowerCase();
        if(regAlign != null && regAlign.length() > 0)
        	pattAlign = Pattern.compile(regAlign);
        String regRep = regexRep.getText().trim();
        regRep = regRep.toLowerCase();
        if(regRep != null && regRep.length() > 0)
        	pattRep = Pattern.compile(regRep);
        
        synchronized(lmes) {
            lmes.clear();
            for (SeqAlignment align : alignments){
                if( (pattType == null || pattType.matcher(align.getExpt().getExptType().getName().toLowerCase()).find()) &&
                        (pattLab == null || pattLab.matcher(align.getExpt().getLab().getName().toLowerCase()).find()) &&
                        (pattCond == null || pattCond.matcher(align.getExpt().getExptCondition().getName().toLowerCase()).find()) &&
                        (pattTarget == null || pattTarget.matcher(align.getExpt().getExptTarget().getName().toLowerCase()).find()) &&
                        (pattCell == null || pattCell.matcher(align.getExpt().getCellLine().getName().toLowerCase()).find()) &&
                        (pattAlign == null || pattAlign.matcher(align.getName().toLowerCase()).find()) &&
                        (pattRep == null || pattRep.matcher(align.getExpt().getReplicate().toLowerCase()).find())
                                ) {
                	List<SeqAlignment> a = new ArrayList<SeqAlignment>();
                	a.add(align);
                    try {
						lmes.add(new SeqLocatorMatchedAligns(a,
						            new SeqLocator(align.getExpt().getName(),
						                                align.getExpt().getReplicate(),
						                                align.getName())));
					} catch (SeqLMENotReplicatesException e1) {
						System.err.println(e1.getMessage());
						e1.printStackTrace();
					}
                }
            }
            filteredModel.clear();
            for (SeqLocatorMatchedAligns l : lmes) {
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
