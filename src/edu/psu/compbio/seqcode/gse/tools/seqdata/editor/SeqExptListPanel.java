package edu.psu.compbio.seqcode.gse.tools.seqdata.editor;

import java.awt.Font;
import java.awt.GridLayout;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

import javax.swing.BorderFactory;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.border.Border;

import edu.psu.compbio.seqcode.gse.datasets.general.ExptType;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExpt;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocatorMatchedExpt;
import edu.psu.compbio.seqcode.gse.seqview.components.SeqTableModel;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class SeqExptListPanel {
private SeqDataLoader seqLoader;
    
    private TreeSet<SeqLocatorMatchedExpt> lme;
    private ArrayList<SeqAlignment> alignments;
    private JComboBox jcbType;
    private JTextField regexLab, regexCond, regexTarget, regexCell, regexAlign, regexRep;
    private SeqTableModel filteredModel;

    public SeqExptListPanel() { 
        try {
            seqLoader = new SeqDataLoader(true);
        } catch (Exception e) {
            e.printStackTrace();
            seqLoader = null;
        }
        lme = new TreeSet<SeqLocatorMatchedExpt>();
        alignments = new ArrayList<SeqAlignment>();
        filteredModel = new SeqTableModel();
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
    public void info() {
    	ArrayList<SeqLocator> locs = new ArrayList<SeqLocator>();
    	int[] inds = filteredList.getSelectedRows();
        for (int i = 0; i < inds.length; i++) {
        	SeqLocatorMatchedExpt o = filteredModel.getObject(inds[i]);
        	locs.add(o.locator);
        }
        for(SeqLocator loc : locs){
        	try {
				Collection<SeqAlignment> aligns = seqLoader.loadAlignments(loc, getGenome());
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
			} catch (SQLException e) {
				e.printStackTrace();
			} catch (NotFoundException e) {
				e.printStackTrace();
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
        
        synchronized(lme) {
            lme.clear();
            for (SeqAlignment align : alignments){
                if( (pattType == null || pattType.matcher(align.getExpt().getExptType().getName().toLowerCase()).find()) &&
                	(pattLab == null || pattLab.matcher(align.getExpt().getLab().getName().toLowerCase()).find()) &&
                	(pattCond == null || pattCond.matcher(align.getExpt().getExptCondition().getName().toLowerCase()).find()) &&
                	(pattTarget == null || pattTarget.matcher(align.getExpt().getExptTarget().getName().toLowerCase()).find()) &&
                	(pattCell == null || pattCell.matcher(align.getExpt().getCellLine().getName().toLowerCase()).find()) &&
                	(pattAlign == null || pattAlign.matcher(align.getName().toLowerCase()).find()) &&
                	(pattRep == null || pattRep.matcher(align.getExpt().getReplicate().toLowerCase()).find())
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

    public void close() {
        if (seqLoader != null) {
            seqLoader.close();
        }
    }
}
