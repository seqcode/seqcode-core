package edu.psu.compbio.seqcode.gse.tools.seqdata.editor;

import java.awt.BorderLayout;
import java.awt.Checkbox;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.security.AccessControlException;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.Border;

import edu.psu.compbio.seqcode.gse.datasets.core.ExptType;
import edu.psu.compbio.seqcode.gse.datasets.core.SeqDataUser;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataModifier;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataModifier.DuplicateDatabaseEntryException;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExpt;
import edu.psu.compbio.seqcode.gse.projects.readdb.ACLChangeEntry;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseException;

/**
 * Data entry form for editing one or more related SeqExpts
 * 
 * @author mahony
 *
 */
public class SeqDataEditEntryForm extends JFrame implements ActionListener {
	/**
	 * Constants for default values and limits on settings
	 */
	public static final int DEFAULT_WINDOW_WIDTH = 900;
	public static final int DEFAULT_WINDOW_HEIGHT = 350;
	public static final int DEFAULT_TOP_LEFT_X = 100;
	public static final int DEFAULT_TOP_LEFT_Y = 100;
	public static final String FINISHED = "SDEEFinished";
	
	private SeqDataLoader seqLoader;
	private SeqDataModifier seqModifier;
	private List<SeqExpt> editExpts;
	private List<SeqAlignment> editAligns=new ArrayList<SeqAlignment>();
	private JComboBox jcbType;
	private JTextField editLab, editCond, editTarget, editCell, editRep, editPubSrc, editPubID, editCollabExptID;
	private boolean sharedExptType=true, sharedLab=true, sharedCond=true, sharedTarget=true, sharedCell=true, sharedRep=true, sharedPubSrc=true, sharedPubID=true, sharedCollabExptID=true;
	private String exptType, lab, cond, target, cell, rep, pubSrc, pubID, collabExptID;
	private List<SeqDataUser> users;
	private HashMap<String, Boolean> sharedUsers = new HashMap<String, Boolean>();
	private HashMap<SeqDataUser, Checkbox> userCheckboxes = new HashMap<SeqDataUser, Checkbox>();
	private JButton okButton, cancelButton;
	private JPanel buttonPanel;
	
	public SeqDataEditEntryForm(SeqDataLoader loader, Collection<SeqExpt> expts){
		super("Edit SeqExperiment Metadata");
		seqLoader = loader;
		try {
			seqModifier = new SeqDataModifier(seqLoader);
		} catch (AccessControlException e1) {
			e1.printStackTrace();
		} catch (SQLException e1) {
			e1.printStackTrace();
		}
		editExpts = new ArrayList(expts);
		
		if(editExpts!=null && editExpts.size()>0){
			exptType = editExpts.get(0).getExptType().getName();
			lab = editExpts.get(0).getLab().getName();
			cond = editExpts.get(0).getExptCondition().getName();
			target = editExpts.get(0).getExptTarget().getName();
			cell = editExpts.get(0).getCellLine().getName();
			rep = editExpts.get(0).getReplicate();
			pubSrc = editExpts.get(0).getPublicSource();
			pubID = editExpts.get(0).getPublicDBID();
			collabExptID = editExpts.get(0).getCollabID();
			
			try{
				users = new ArrayList(seqLoader.getSeqDataUsers());
				for(SeqDataUser u : users)
					sharedUsers.put(u.getName(), false);
				for(SeqExpt e : editExpts)
					editAligns.addAll(seqLoader.loadAllAlignments(e));
				for(String us : editAligns.get(0).getPermissions())
					sharedUsers.put(us, true);
				
			} catch (SQLException e) {
				e.printStackTrace();
			}
		}
		
		init();
	}
	
	public void init(){
		setSize(DEFAULT_WINDOW_WIDTH,DEFAULT_WINDOW_HEIGHT);
        setLocation(DEFAULT_TOP_LEFT_X,DEFAULT_TOP_LEFT_Y);
    	
		//First off, check which fields are shared by all experiments
		for(SeqExpt e : editExpts){
			if(!e.getExptType().getName().equals(exptType)){sharedExptType=false; exptType="";}
			if(!e.getLab().getName().equals(lab)){sharedLab=false; lab="";}
			if(!e.getExptCondition().getName().equals(cond)){sharedCond=false; cond="";}
			if(!e.getExptTarget().getName().equals(target)){sharedTarget=false; target="";}
			if(!e.getCellLine().getName().equals(cell)){sharedCell=false; cell="";}
			if(!e.getReplicate().equals(rep)){sharedRep=false; rep="";}
			if(!e.getPublicSource().equals(pubSrc)){sharedPubSrc=false; pubSrc="";}
			if(!e.getPublicDBID().equals(pubID)){sharedPubID=false; pubID="";}
			if(!e.getCollabID().equals(collabExptID)){sharedCollabExptID=false; collabExptID="";}
		}
		for(SeqAlignment a : editAligns){
			for(String currPerm : a.getPermissions())
				sharedUsers.put(currPerm, sharedUsers.get(currPerm)&&true);
		}
		
		if(!sharedExptType && !sharedLab && !sharedCond && !sharedTarget && !sharedCell && !sharedRep && !sharedPubSrc && !sharedPubID && !sharedCollabExptID){
			JOptionPane.showMessageDialog(null, "Cannot update this selection: no information is shared by all experiments", "Error", JOptionPane.INFORMATION_MESSAGE);
			dispose();
		}else{
			JPanel actPanel = new JPanel(); actPanel.setLayout(new BorderLayout());
			
			okButton = new JButton("OK");
			cancelButton = new JButton("Cancel");
			okButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                	try {
						updateSeqExpts();
						firePropertyChange(FINISHED, false, true);
						//Flash a message confirming update?
						dispose();
					} catch (SQLException e1) {
						e1.printStackTrace();
					}
                }
            });
			cancelButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                	firePropertyChange(FINISHED, true, false);
                	dispose();
                }
            });
			
			//Layout
			buttonPanel = new JPanel();
	        buttonPanel.setLayout(new GridBagLayout());
	        buttonPanel.add(okButton);
	        buttonPanel.add(cancelButton);
	        JPanel inputsPanel = getInputsPanel(); 
	        JPanel permPanel = getPermissionsPanel();
	        JPanel allInputsPanel = new JPanel(); allInputsPanel.setLayout(new BoxLayout(allInputsPanel, BoxLayout.PAGE_AXIS));
	        allInputsPanel.add(inputsPanel);
	        allInputsPanel.add(Box.createHorizontalGlue());
	        allInputsPanel.add(permPanel);
	        allInputsPanel.add(Box.createVerticalGlue());
	        allInputsPanel.add(buttonPanel);
	        actPanel.add(allInputsPanel,BorderLayout.CENTER);
	        add(actPanel);
	        this.setVisible(true);
	        this.repaint();
		}
	}
	
	public JPanel getInputsPanel() {
    	JPanel inputPanel = new JPanel();
		inputPanel.setLayout(new GridLayout(2,1));
		JPanel inputPanelUpper = new JPanel();
		inputPanelUpper.setLayout(new GridLayout(2,6));
		JPanel inputPanelLower = new JPanel();
		inputPanelLower.setLayout(new GridLayout(2,6));
		
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
			}else{ editLab.setText(lab); editLab.setCaretPosition(0);}
	        
			editCond = new JTextField(); 
	        if(!sharedCond){
				editCond.setEnabled(false); editCond.setOpaque(true);
			}else{ editCond.setText(cond); editCond.setCaretPosition(0);}
	        
	        editTarget = new JTextField(); 
	        if(!sharedTarget){
				editTarget.setEnabled(false); editTarget.setOpaque(true);
			}else{ editTarget.setText(target); editTarget.setCaretPosition(0);}
	        
	        editCell = new JTextField(); 
	        if(!sharedCell){
				editCell.setEnabled(false); editCell.setOpaque(true);
			}else{ editCell.setText(cell); editCell.setCaretPosition(0);}
	        
	        editRep = new JTextField(); 
	        if(!sharedRep){
				editRep.setEnabled(false); editRep.setOpaque(true);
			}else{ editRep.setText(rep); editRep.setCaretPosition(0);}
	        
	        editPubSrc = new JTextField(); 
	        if(!sharedPubSrc){
				editPubSrc.setEnabled(false); editPubSrc.setOpaque(true);
			}else{ editPubSrc.setText(pubSrc); editPubSrc.setCaretPosition(0);}
	        
	        editPubID = new JTextField(); 
	        if(!sharedPubID){
				editPubID.setEnabled(false); editPubID.setOpaque(true);
			}else{ editPubID.setText(pubID); editPubID.setCaretPosition(0);}
	        
	        editCollabExptID = new JTextField(); 
	        if(!sharedCollabExptID){
				editCollabExptID.setEnabled(false); editCollabExptID.setOpaque(true);
			}else{ editCollabExptID.setText(collabExptID); editCollabExptID.setCaretPosition(0);}
	        
	        Font labelFont = new Font("SansSerif", Font.BOLD, 12);
	        Border paddingBorder = BorderFactory.createEmptyBorder(5,5,5,5);
	        JLabel labelExpt = new JLabel("ExptType"); labelExpt.setFont(labelFont); labelExpt.setBorder(paddingBorder);
	        JLabel labelLab = new JLabel("Lab"); labelLab.setFont(labelFont); labelLab.setBorder(paddingBorder);
	        JLabel labelCond = new JLabel("ExptCondition"); labelCond.setFont(labelFont); labelCond.setBorder(paddingBorder); 
	        JLabel labelTarget = new JLabel("ExptTarget"); labelTarget.setFont(labelFont); labelTarget.setBorder(paddingBorder);
	        JLabel labelCell = new JLabel("CellLine"); labelCell.setFont(labelFont); labelCell.setBorder(paddingBorder);
	        JLabel labelRep = new JLabel("Replicate"); labelRep.setFont(labelFont); labelRep.setBorder(paddingBorder);
	        JLabel labelPubSrc = new JLabel("PublicSource"); labelPubSrc.setFont(labelFont); labelPubSrc.setBorder(paddingBorder);
	        JLabel labelPubID = new JLabel("PublicDBID"); labelPubID.setFont(labelFont); labelPubID.setBorder(paddingBorder);
	        JLabel labelCollabExptID = new JLabel("CollabID"); labelCollabExptID.setFont(labelFont); labelCollabExptID.setBorder(paddingBorder);
	        inputPanelUpper.add(labelExpt); 
	        inputPanelUpper.add(labelLab);
	        inputPanelUpper.add(labelCond);
	        inputPanelUpper.add(labelTarget);
	        inputPanelUpper.add(labelCell);
	        inputPanelUpper.add(labelRep);
	        inputPanelUpper.add(jcbType);
	        inputPanelUpper.add(editLab);
	        inputPanelUpper.add(editCond);
	        inputPanelUpper.add(editTarget);
	        inputPanelUpper.add(editCell);
	        inputPanelUpper.add(editRep);
	        inputPanelLower.add(new Label()); inputPanelLower.add(new Label()); inputPanelLower.add(new Label());
	        inputPanelLower.add(labelPubSrc);
	        inputPanelLower.add(labelPubID);
	        inputPanelLower.add(labelCollabExptID);
	        inputPanelLower.add(new Label()); inputPanelLower.add(new Label()); inputPanelLower.add(new Label());
	        inputPanelLower.add(editPubSrc);
	        inputPanelLower.add(editPubID);
	        inputPanelLower.add(editCollabExptID);
	        inputPanel.add(inputPanelUpper);
	        inputPanel.add(inputPanelLower);
        } catch (SQLException e) {
			e.printStackTrace();
		}
        
        return inputPanel;
    }
	
	public JPanel getPermissionsPanel(){
		JPanel permissionsPanel = new JPanel();
		Font labelFont = new Font("SansSerif", Font.BOLD, 12);
		Border paddingBorder = BorderFactory.createEmptyBorder(5,5,5,5);
		JLabel labelPerm = new JLabel("Permissions"); labelPerm.setFont(labelFont); labelPerm.setBorder(paddingBorder); labelPerm.setAlignmentX(Component.LEFT_ALIGNMENT);
		permissionsPanel.setLayout(new BoxLayout(permissionsPanel, BoxLayout.Y_AXIS));
		permissionsPanel.add(labelPerm);
		
		JPanel permCBPanel = new JPanel();
		permCBPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		permCBPanel.setLayout(new GridLayout(0,2));
		for(SeqDataUser u : users){
			Checkbox cb = new Checkbox(u.getName(), sharedUsers.containsKey(u.getName()) && sharedUsers.get(u.getName()));
			if(u.isAdmin())
				cb.setEnabled(false);
			userCheckboxes.put(u, cb);
			permCBPanel.add(cb);
		}
		permissionsPanel.add(permCBPanel);
		permissionsPanel.add(Box.createVerticalGlue());
		return permissionsPanel;
	}
	
	public void updateSeqExpts() throws SQLException{
		//Get the new entries first
		String newExptType = jcbType.isEnabled() ? (String)jcbType.getSelectedItem() : "";
		String newLab = editLab.isEnabled() ? (String)editLab.getText(): "";
		String newCond = editCond.isEnabled() ? (String)editCond.getText(): "";
		String newTarget = editTarget.isEnabled() ? (String)editTarget.getText(): "";
		String newCell = editCell.isEnabled() ? (String)editCell.getText(): "";
		String newRep = editRep.isEnabled() ? (String)editRep.getText(): "";
		String newPubSrc = editPubSrc.isEnabled() ? (String)editPubSrc.getText(): "";
		String newPubID = editPubID.isEnabled() ? (String)editPubID.getText(): "";
		String newCollabExptID = editCollabExptID.isEnabled() ? (String)editCollabExptID.getText(): "";
		
		String newPermissions ="";
		for(SeqDataUser u : users){
			if(userCheckboxes.get(u).getState())
				newPermissions = newPermissions+u.getName()+";";
		}
		
		//Update each each of the SeqExpts
		for(SeqExpt expt : editExpts){
			String updateExptType = newExptType.equals("") ? expt.getExptType().getName() : newExptType;
			String updateLab = newLab.equals("") ? expt.getLab().getName() : newLab;
			String updateCond = newCond.equals("") ? expt.getExptCondition().getName() : newCond;
			String updateTarget = newTarget.equals("") ? expt.getExptTarget().getName() : newTarget;
			String updateCell = newCell.equals("") ? expt.getCellLine().getName() : newCell;
			String updateRep = newRep.equals("") ? expt.getReplicate() : newRep;
			String updatePubSrc = newPubSrc.equals("") ? expt.getPublicSource() : newPubSrc;
			String updatePubID = newPubID.equals("") ? expt.getPublicDBID() : newPubID;
			String updateCollabExptID = newCollabExptID.equals("") ? expt.getCollabID() : newCollabExptID;
			String updateName = updateLab+" "+updateCond+" "+updateTarget+" "+updateCell;
			
			System.err.println("UPDATING:\t"+updateName+"\t"+updateExptType+"\t"+updateLab+"\t"+updateCond+"\t"+updateTarget+"\t"+updateCell+"\t"+updateRep+"\t"+updatePubSrc+"\t"+updatePubID+"\t"+updateCollabExptID+"\t"+expt.getDBID());
			try {
				seqModifier.updateSeqExpt(expt, updateExptType, updateLab, updateCond, updateTarget, updateCell, updateRep, updatePubSrc, updatePubID, updateCollabExptID);
			} catch (DuplicateDatabaseEntryException e) {
				System.err.println("Cannot update "+expt.getName()+";"+expt.getReplicate()+" to "+updateName+";"+updateRep+" as this will result in a duplicate SeqExpt key.");
				JOptionPane.showMessageDialog(null, "Cannot update "+expt.getName()+";"+expt.getReplicate()+" to "+updateName+";"+updateRep+" as this will result in a duplicate SeqExpt key.", "Error: attempted duplicate", JOptionPane.INFORMATION_MESSAGE);
			}
			
			//Update the permissions of each alignment
            for(SeqAlignment align : seqLoader.loadAllAlignments(expt)){
            	seqModifier.updateSeqAlignmentPermissions(align, newPermissions);
            }
			
            for(SeqAlignment align : seqLoader.loadAllAlignments(expt)){
				//System.out.println(align.getDBID()+"\t"+newPermissions);
				String[] princs = new String[users.size()];
				String[] ops = new String[users.size()];
				String[] acls = new String[users.size()];
				int i=0;
				for(SeqDataUser u : users){
					princs[i] = u.getName();
					ops[i] = userCheckboxes.get(u).getState() ? "add" : "delete";
					acls[i] = "read";                            
					i++;
				}
				seqModifier.changeAlignmentACLmulti(align, princs, ops, acls);
			}
		}
		//This will delete any now-redundant Lab, ExptTargets, etc from the db
		seqModifier.coreCleanup();
	}
	
	public void actionPerformed(ActionEvent e) {
		
	}
}
