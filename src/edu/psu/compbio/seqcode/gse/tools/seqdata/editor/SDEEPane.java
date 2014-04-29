package edu.psu.compbio.seqcode.gse.tools.seqdata.editor;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.Collection;
import java.util.List;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.utils.Closeable;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class SDEEPane extends JTabbedPane implements ItemListener, ActionListener, Closeable {

	private JPanel seqPanel;
	private boolean handlingChange, closed;

	public SDEEPane ()throws NotFoundException {
        super();
        handlingChange=false;
	}
	public boolean isClosed() { return closed; }
    
    public void close() { 
        closed = true; 
    }
    
    private void init() throws NotFoundException {
    	seqPanel = new JPanel();
    	
        addTab("Seq Data", seqPanel);

    }
    
    //updates the choice of experiments
    private void updateExptSelection() {
    	//seqSelect.setGenome(lg);
    	
    	
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
