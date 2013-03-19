package edu.psu.compbio.seqcode.gse.seqview.components;

import java.util.*;
import java.awt.*;
import javax.swing.*;

import edu.psu.compbio.seqcode.gse.seqview.paintable.SeqViewPaintable;

public class SeqViewPaintableHTConfigurationPanel extends HashtableConfigurationPanel {

    private SeqViewPaintable wp;
    private Hashtable<String,Object> table;
    public SeqViewPaintableHTConfigurationPanel (Hashtable<String,Object> table, SeqViewPaintable w) {
        super(table,w.configurationKeyOrder());
        this.table = table;
        wp = w;
        System.err.println("WPHTCP for " + w);
        table.put("TrackLabel",w.getLabel());
    }
    public void parse() {
        super.parse();
        if (table.containsKey("TrackLabel") &&
            ((String)table.get("TrackLabel")).length() > 0) {
            wp.setLabel((String)table.get("TrackLabel"));
            table.remove("TrackLabel");
        }
    }
}
