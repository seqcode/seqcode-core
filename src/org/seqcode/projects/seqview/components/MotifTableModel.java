package org.seqcode.projects.seqview.components;

import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.*;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.viz.components.ObjectTableModel;

import java.util.*;


public class MotifTableModel extends ObjectTableModel<WeightMatrix> {

    private boolean sortByVersion = false;

    public int getColumnCount() {
        return 3;
    }

    public Class getColumnClass(int i) {
        if(i==0) { return String.class; }
        if(i==1) { return String.class; }
        if(i==2) { return String.class; }
        return null;
    }

    public String getColumnName(int i) {
        if(i==0) { return "Name";}
        if(i==1) { return "Version"; }
        if(i==2) { return "Type"; }
        return null;
    }

    public Object getValueAt(int rowIndex, int c) {
        if(c==0) { return getObject(rowIndex).name; }
        if(c==1) { return getObject(rowIndex).version; }
        if(c==2) { return getObject(rowIndex).type; }
        return null;
    }

    public void sortByVersion() {
        sort(new WeightMatrixVersionComparator());
        sortByVersion = true;
    }
    
    public void sortByName() {
        sort(new WeightMatrixNameComparator());
        sortByVersion = false;
    }
}

class WeightMatrixVersionComparator implements Comparator<WeightMatrix> {
    public int compare(WeightMatrix a, WeightMatrix b) {
        return a.version.compareTo(b.version);
    }
}
class WeightMatrixNameComparator implements Comparator<WeightMatrix> {
    public int compare(WeightMatrix a, WeightMatrix b) {
        return a.name.compareTo(b.name);
    }
}


