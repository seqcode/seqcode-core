package edu.psu.compbio.seqcode.gse.tools.seqdata.editor;

import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
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
import edu.psu.compbio.seqcode.gse.seqview.components.SeqLMETableModel;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.viz.components.GenericEditTablePanel;

public class SeqExptEditTablePanel  extends GenericEditTablePanel<SeqExpt> {
	private SeqDataLoader seqLoader;
    private TreeSet<SeqExpt> expts;
    private JComboBox jcbType;
    private JTextField regexLab, regexCond, regexTarget, regexCell, regexRep;
    private SeqExptTableModel filteredModel;

    public SeqExptEditTablePanel() { 
        try {
            seqLoader = new SeqDataLoader(true);
        } catch (Exception e) {
            e.printStackTrace();
            seqLoader = null;
        }
        expts = new TreeSet<SeqExpt>();
        filteredModel = new SeqExptTableModel();
        init(filteredModel);
    }
    public JPanel getInputsPanel() {
    	JPanel inputPanel = new JPanel();
		inputPanel.setLayout(new GridLayout(2,6));
		
		try {
    		ArrayList<String> types = new ArrayList<String>();
    		types.add("");
        	for(ExptType e : seqLoader.getExptTypes())
				types.add(e.getName());
			Collections.sort(types);
			jcbType = new JComboBox(types.toArray()); jcbType.setEditable(true);
			jcbType.addActionListener(new ActionListener(){public void actionPerformed(ActionEvent e){filter();}});
			
	        regexLab = new JTextField(); regexLab.addActionListener(new ActionListener(){public void actionPerformed(ActionEvent e){filter();}});
	        regexCond = new JTextField(); regexCond.addActionListener(new ActionListener(){public void actionPerformed(ActionEvent e){filter();}});
	        regexTarget = new JTextField(); regexTarget.addActionListener(new ActionListener(){public void actionPerformed(ActionEvent e){filter();}});
	        regexCell = new JTextField(); regexCell.addActionListener(new ActionListener(){public void actionPerformed(ActionEvent e){filter();}});
	        regexRep = new JTextField(); regexRep.addActionListener(new ActionListener(){public void actionPerformed(ActionEvent e){filter();}});
	        
	        Font labelFont = new Font("SansSerif", Font.PLAIN, 10);
	        Border paddingBorder = BorderFactory.createEmptyBorder(5,5,5,5);
	        JLabel labelExpt = new JLabel("ExptType"); labelExpt.setFont(labelFont); labelExpt.setBorder(paddingBorder);
	        JLabel labelLab = new JLabel("Lab"); labelLab.setFont(labelFont); labelLab.setBorder(paddingBorder);
	        JLabel labelCond = new JLabel("ExptCondition"); labelCond.setFont(labelFont); labelCond.setBorder(paddingBorder); 
	        JLabel labelTarget = new JLabel("ExptTarget"); labelTarget.setFont(labelFont); labelTarget.setBorder(paddingBorder);
	        JLabel labelCell = new JLabel("CellLine"); labelCell.setFont(labelFont); labelCell.setBorder(paddingBorder);
	        JLabel labelRep = new JLabel("Replicate"); labelRep.setFont(labelFont); labelRep.setBorder(paddingBorder);
	        inputPanel.add(labelExpt); 
	        inputPanel.add(labelLab);
	        inputPanel.add(labelCond);
	        inputPanel.add(labelTarget);
	        inputPanel.add(labelCell);
	        inputPanel.add(labelRep);
	        inputPanel.add(jcbType);
	        inputPanel.add(regexLab);
	        inputPanel.add(regexCond);
	        inputPanel.add(regexTarget);
	        inputPanel.add(regexCell);
	        inputPanel.add(regexRep);
        
        } catch (SQLException e) {
			e.printStackTrace();
		}
        
        return inputPanel;
    }

    public void retrieveData() {
        try {
            synchronized(expts) {
                expts.clear();
                expts.addAll(seqLoader.loadAllExperiments());
            }
        } catch (SQLException e) {
            throw new RuntimeException(e.toString(), e);
        }
    }
    public void updateComponents() {
        filteredModel.clear();
        synchronized(expts) {
            for (SeqExpt l : expts) {
                filteredModel.addObject(l);
            }
        }
    }
    

    public void filter() {
    	Pattern pattType=null, pattLab=null, pattCond=null, pattTarget=null, pattCell=null, pattRep=null;
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
        String regRep = regexRep.getText().trim();
        regRep = regRep.toLowerCase();
        if(regRep != null && regRep.length() > 0)
        	pattRep = Pattern.compile(regRep);
        
        synchronized(expts) {
        	List<SeqExpt> filteredExpts = new ArrayList<SeqExpt>();
            for (SeqExpt expt : expts){
                if( (pattType == null || pattType.matcher(expt.getExptType().getName().toLowerCase()).find()) &&
                	(pattLab == null || pattLab.matcher(expt.getLab().getName().toLowerCase()).find()) &&
                	(pattCond == null || pattCond.matcher(expt.getExptCondition().getName().toLowerCase()).find()) &&
                	(pattTarget == null || pattTarget.matcher(expt.getExptTarget().getName().toLowerCase()).find()) &&
                	(pattCell == null || pattCell.matcher(expt.getCellLine().getName().toLowerCase()).find()) &&
                	(pattRep == null || pattRep.matcher(expt.getReplicate().toLowerCase()).find()) &&
                	(((String)this.speciesCBox.getSelectedItem()).equals("") || ((String)this.speciesCBox.getSelectedItem()).equals(expt.getOrganism().getName()))
                		) {
                		filteredExpts.add(expt);
                }
            }
            filteredModel.clear();
            for (SeqExpt l : filteredExpts) {
                filteredModel.addObject(l);
            }
        }
    }

    public void edit(Collection<SeqExpt> toEdit){
    	SeqDataEditEntryForm editForm = new SeqDataEditEntryForm(seqLoader, toEdit);
    }
    
    public void close() {
        if (seqLoader != null) {
            seqLoader.close();
        }
    }
}
