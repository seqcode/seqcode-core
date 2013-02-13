/*
 * Created on Nov 20, 2006
 */
package edu.psu.compbio.seqcode.gse.seqview.components;

import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.*;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.expression.*;
import edu.psu.compbio.seqcode.gse.viz.components.ObjectTableModel;

public class ExprExperimentTableModel extends ObjectTableModel<Experiment> {
	
	private boolean sortByName;
	
	public ExprExperimentTableModel() { 
		sortByName = true;
	}
    
    public int getColumnCount() {
        return 1;
    }
    
    public Class getColumnClass(int i) {
        if(i==0) { return String.class; }
        if(i==1) { return String.class; }
        return null;
    }
    
    public String getColumnName(int i) {
        if(i==0) { return "Name"; }
        return null;
    }

    public Object getValueAt(int rowIndex, int c) {
        if(c==0) { return getObject(rowIndex).getName(); }
        return null;
    }

    public void sortByName() {
        sort(new ExperimentNameComparator());
        sortByName = true;
    }
    
}

class ExperimentNameComparator implements Comparator<Experiment> {
    public int compare(Experiment a, Experiment b) {
        return a.getName().compareTo(b.getName());
    }
}

