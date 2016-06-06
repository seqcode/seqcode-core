package org.seqcode.projects.seqview.components;

import javax.swing.*;

import org.seqcode.genome.sequence.*;

import java.awt.event.*;


public class SeqViewToolsMenu extends JMenu {

    private final RegionPanel panel;

    public SeqViewToolsMenu(RegionPanel p) {
        super("Tools");
        panel = p;
        final RegionPanel thispanel = panel;
        /*JMenuItem item = new JMenuItem("Genome Browser");
        add(item);
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    SeqViewOptionsFrame.main(new String[0]);
                }
            });
        */
        JMenuItem item = new JMenuItem("Browse Weight Matrices");
        add(item);
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    MotifDisplayPane.main(new String[0]);
                }
            });
        
        item = new JMenuItem("Sequence Browser");
        add(item);
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    if (panel == null) {
                        SequenceViewer seqview = new SequenceViewer();
                    }else {
                        SequenceViewer seqview = new SequenceViewer(panel.getRegion());
                    }
                }
            });
    }


}
