package edu.psu.compbio.seqcode.gse.tools.seqdata.editor;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import edu.psu.compbio.seqcode.gse.utils.Closeable;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class SDEEPane extends JTabbedPane implements ItemListener, ActionListener, Closeable {

	private JPanel seqPanel;
	private SeqExptEditTablePanel seqSelect;
	private boolean handlingChange, closed;

	public SDEEPane ()throws NotFoundException {
        super();
        handlingChange=false;
        init();
	}
	public boolean isClosed() { return closed; }
    
    public void close() { 
        closed = true; 
    }
    
    private void init() throws NotFoundException {
    	seqPanel = new JPanel();
    	seqSelect = new SeqExptEditTablePanel();
    	seqPanel.setLayout(new BorderLayout());
        seqPanel.add(seqSelect, BorderLayout.CENTER);
        seqSelect.updateData();
        addTab("Seq Data", seqPanel);
        setVisible(true);
    }
    
    //updates the choice of experiments
    private void updateExptSelection() {
    	seqSelect.updateData();    	
    }
	
	
    public void itemStateChanged(ItemEvent e) {
    	if (handlingChange) {return;}
    	Object source = e.getItemSelectable();
    	synchronized(this) {
    		if (!handlingChange) {
    			handlingChange = true;
    			updateExptSelection();
    			handlingChange = false;
    		}
    	}
	     
    }
    public void actionPerformed (ActionEvent e) {

    }

}
