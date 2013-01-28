/*
 * Created on Dec 5, 2006
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.assignment;

import java.awt.BorderLayout;
import java.awt.Container;
import java.sql.SQLException;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.binding.BindingScan;
import edu.psu.compbio.seqcode.gse.datasets.binding.BindingScanLoader;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;
import edu.psu.compbio.seqcode.gse.ewok.verbs.binding.BindingScanGenerator;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;

public class AnnotatedRegionTableModel implements TableModel {

    // columns:
    // 0: chrom
    // 1: start
    // 2: stop
    // 3: bit-vector

    private Annotations<Region, BindingEvent> annotations;
    private LinkedList<TableModelListener> listeners;
    
    public static final int CHROMCOL = 0;
    public static final int STARTCOL = 1;
    public static final int ENDCOL = 2;
    public static final int BITSCOL = 3;
    
    public AnnotatedRegionTableModel(Annotations<Region, BindingEvent> annots) { 
        annotations = annots;
        listeners = new LinkedList<TableModelListener>();
    }

    public void addTableModelListener(TableModelListener tml) {
        listeners.add(tml);
    }

    public void removeTableModelListener(TableModelListener tml) {
        listeners.remove(tml);
    }

    public Class< ? > getColumnClass(int c) {
        switch(c) { 
        case CHROMCOL: return String.class;
        case STARTCOL: return Integer.class;
        case ENDCOL: return Integer.class;
        case BITSCOL: return String.class;
        }
        return null;
    }

    public int getColumnCount() {
        return 4;
    }

    public String getColumnName(int c) {
        switch(c) { 
        case CHROMCOL: return "Chrom:";
        case STARTCOL: return "Start:";
        case ENDCOL: return "End:";
        case BITSCOL: return "Binding:";
        }
        return null;
    }

    public int getRowCount() {
        return annotations.getNumItems();
    }

    public Object getValueAt(int r, int c) {
        switch(c) { 
        case CHROMCOL: return annotations.getItem(r).getChrom();
        case STARTCOL: return annotations.getItem(r).getStart();
        case ENDCOL: return annotations.getItem(r).getEnd();
        case BITSCOL: return getBits(r);
        }
        return null;
    }
    
    private String getBits(int r) {
        return annotations.getAnnotationBitVector(annotations.getItem(r));
    }

    public boolean isCellEditable(int arg0, int arg1) {
        return false;
    }

    public void setValueAt(Object arg0, int arg1, int arg2) {
        throw new UnsupportedOperationException("Can't edit cell.");
    }
}
