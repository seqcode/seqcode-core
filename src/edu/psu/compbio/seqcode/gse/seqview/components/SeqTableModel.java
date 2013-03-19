package edu.psu.compbio.seqcode.gse.seqview.components;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocatorMatchedExpt;
import edu.psu.compbio.seqcode.gse.viz.components.ObjectTableModel;

public class SeqTableModel extends ObjectTableModel<SeqLocatorMatchedExpt> {
    
    private int compareNameVersions(String n1, String v1, String n2, String v2) { 
        int c = v1.compareTo(v2);
        if(c != 0) { return c; }
        return n1.compareTo(n2);
    }
    
    private int findNewIndex(SeqLocatorMatchedExpt bs) {
        String n = bs.locator.getExptName(), v = bs.locator.getAlignName();
        for(int i = 0; i < getSize(); i++) { 
            SeqLocator os = getObject(i).locator;
            String on = os.getExptName(), ov = os.getAlignName();
            int c = compareNameVersions(n, v, on, ov);
            if(c <= 0) { return i; }
        }
        return getSize();
    }

    public int getColumnCount() {
        return 7;
    }
    
    public Class getColumnClass(int i) {
        if(i==0) { return String.class; }
        if(i==1) { return String.class; }
        if(i==2) { return String.class; }
        if(i==3) { return String.class; }
        if(i==4) { return String.class; }
        if(i==5) { return String.class; }
        if(i==6) { return String.class; }
        return null;
    }


    public String getColumnName(int i) {
        if(i==0) { return "ExptType"; }
        if(i==1) { return "Lab"; }
        if(i==2) { return "ExptCondition"; }
        if(i==3) { return "ExptTarget"; }
        if(i==4) { return "CellLine"; }
        if(i==5) { return "Alignment"; }
        if(i==6) { return "Replicate"; }
        return null;
    }

    public Object getValueAt(int rowIndex, int c) {
        if(c==0) { return getObject(rowIndex).expts.get(0).getExptType().getName(); }
        if(c==1) { return getObject(rowIndex).expts.get(0).getLab().getName(); }
        if(c==2) { return getObject(rowIndex).expts.get(0).getExptCondition().getName(); }
        if(c==3) { return getObject(rowIndex).expts.get(0).getExptTarget().getName(); }
        if(c==4) { return getObject(rowIndex).expts.get(0).getCellLine().getName(); }
        if(c==5) { return getObject(rowIndex).locator.getAlignName(); }
        if(c==6) { return getObject(rowIndex).locator.getReplicateString(); }
        return null;
    }

}


