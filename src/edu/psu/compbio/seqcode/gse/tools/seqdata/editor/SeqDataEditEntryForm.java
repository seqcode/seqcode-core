package edu.psu.compbio.seqcode.gse.tools.seqdata.editor;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.Border;

import edu.psu.compbio.seqcode.gse.datasets.general.ExptType;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExpt;

public class SeqDataEditEntryForm extends JFrame implements ActionListener {
	/**
	 * Constants for default values and limits on settings
	 */
	public static final int DEFAULT_WINDOW_WIDTH = 900;
	public static final int DEFAULT_WINDOW_HEIGHT = 650;
	public static final int DEFAULT_TOP_LEFT_X = 50;
	public static final int DEFAULT_TOP_LEFT_Y = 50;
	
	private SeqDataLoader seqLoader;
	private List<SeqExpt> editExpts;
	private JComboBox jcbType;
	private JTextField editLab, editCond, editTarget, editCell, editRep;
	private boolean sharedExptType=true, sharedLab=true, sharedCond=true, sharedTarget=true, sharedCell=true, sharedRep=true;
	private String exptType, lab, cond, target, cell, rep;
	private JButton okButton, cancelButton;
	private JPanel buttonPanel;
	
	public SeqDataEditEntryForm(SeqDataLoader loader, Collection<SeqExpt> expts){
		super();
		seqLoader = loader;
		editExpts = new ArrayList(expts);
		
		if(editExpts!=null && editExpts.size()>0){
			exptType = editExpts.get(0).getExptType().getName();
			lab = editExpts.get(0).getLab().getName();
			cond = editExpts.get(0).getExptCondition().getName();
			target = editExpts.get(0).getExptTarget().getName();
			cell = editExpts.get(0).getCellLine().getName();
			rep = editExpts.get(0).getReplicate();
		}
		
		init();
	}
	
	public void init(){
		//First off, check which fields are shared by all experiments
		for(SeqExpt e : editExpts){
			if(!e.getExptType().getName().equals(exptType)){sharedExptType=false; exptType="";}
			if(!e.getLab().getName().equals(lab)){sharedLab=false; lab="";}
			if(!e.getExptCondition().getName().equals(cond)){sharedCond=false; cond="";}
			if(!e.getExptTarget().getName().equals(target)){sharedTarget=false; target="";}
			if(!e.getCellLine().getName().equals(cell)){sharedCell=false; cell="";}
			if(!e.getReplicate().equals(rep)){sharedRep=false; rep="";}
		}
		
		if(!sharedExptType && !sharedLab && !sharedCond && !sharedTarget && !sharedCell && !sharedRep){
			//Error handle here - nothing shared between selected experiments
		}else{
			
			JPanel actPanel = new JPanel(); actPanel.setLayout(new BorderLayout());
			okButton = new JButton("OK");
			cancelButton = new JButton("Cancel");
			buttonPanel = new JPanel();
	        buttonPanel.setLayout(new GridBagLayout());
	        buttonPanel.add(okButton);
	        buttonPanel.add(cancelButton);
	        JPanel inputsPanel = getInputsPanel();
	        JPanel allInputsPanel = new JPanel(); allInputsPanel.setLayout(new BorderLayout());
	        allInputsPanel.add(inputsPanel,BorderLayout.CENTER);
	        allInputsPanel.add(buttonPanel,BorderLayout.SOUTH);
	        actPanel.add(allInputsPanel,BorderLayout.SOUTH);
	        add(actPanel);
		}
		this.setVisible(true);
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
			jcbType = new JComboBox(types.toArray());
			if(!sharedExptType){
				jcbType.setOpaque(true);
				jcbType.setEnabled(false);
			}else{
				jcbType.setSelectedItem(exptType);
			}
			
			editLab = new JTextField(); 
			if(!sharedLab){
				editLab.setEnabled(false); editLab.setOpaque(true);
			}else{ editLab.setText(lab); }
	        
			editCond = new JTextField();
	        if(!sharedCond){
				editCond.setEnabled(false); editCond.setOpaque(true);
			}else{ editCond.setText(cond); }
	        
	        editTarget = new JTextField(); 
	        if(!sharedTarget){
				editTarget.setEnabled(false); editTarget.setOpaque(true);
			}else{ editTarget.setText(target); }
	        
	        editCell = new JTextField(); 
	        if(!sharedCell){
				editCell.setEnabled(false); editCell.setOpaque(true);
			}else{ editCell.setText(cell); }
	        
	        editRep = new JTextField(); 
	        if(!sharedRep){
				editRep.setEnabled(false); editRep.setOpaque(true);
			}else{ editRep.setText(rep); }
	        
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
	        inputPanel.add(editLab);
	        inputPanel.add(editCond);
	        inputPanel.add(editTarget);
	        inputPanel.add(editCell);
	        inputPanel.add(editRep);
        
        } catch (SQLException e) {
			e.printStackTrace();
		}
        
        return inputPanel;
    }
	
	public void actionPerformed(ActionEvent e) {
		
	}
}
