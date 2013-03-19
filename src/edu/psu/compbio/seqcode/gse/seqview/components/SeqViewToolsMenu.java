package edu.psu.compbio.seqcode.gse.seqview.components;

import javax.swing.*;
import java.awt.event.*;

import edu.psu.compbio.seqcode.gse.datasets.function.*;
import edu.psu.compbio.seqcode.gse.tools.sequence.*;

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
        item = new JMenuItem("Gene Annotations");
        add(item);
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    try {
                        DatabaseFunctionLoader loader = new DatabaseFunctionLoader();
                        GOAnnotationPanel panel = new GOAnnotationPanel(loader, thispanel.getRegion().getGenome());
                        JFrame f = new JFrame();
                        f.add(panel);
                        f.pack();
                        f.setVisible(true);
                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }
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
        if (panel != null) {
            item = new JMenuItem("Run CGH Analysis");
            add(item);
            item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        if (panel != null) {
                            new CGHAnalysisPanel(panel);
                        } 
                    }
                }
                );
        }
    }


}
