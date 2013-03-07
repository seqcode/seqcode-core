package edu.psu.compbio.seqcode.gse.seqview.components;

import java.awt.*;
import java.sql.SQLException;
import java.util.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.Collection;
import java.util.Iterator;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.binding.BindingScan;
import edu.psu.compbio.seqcode.gse.datasets.binding.BindingScanLoader;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.binding.BindingExpander;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;
import edu.psu.compbio.seqcode.gse.viz.components.BindingScanSelectPanel;


public class BindingScanSelectFrame extends JFrame implements ActionListener {

    private Genome genome;
    private BindingScanSelectPanel bssp;
    private ChipSeqAnalysisSelectPanel csasp;
    private JButton ok, cancel;    
    private RegionList list;

    public static void main(String args[]) throws Exception {
        Genome genome = new Genome(args[0],args[1]);
        BindingEventAnnotationPanel beap = new BindingEventAnnotationPanel(null, new ArrayList<BindingEvent>());        
        JFrame f = new BindingEventAnnotationPanel.Frame(beap);
        f.pack();
        new BindingScanSelectFrame(genome,beap);
    }
    
    public BindingScanSelectFrame (Genome g, RegionList list) {
        super();

        try {
			bssp = new BindingScanSelectPanel();
            csasp = new ChipSeqAnalysisSelectPanel(g);
		} catch (UnknownRoleException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		} catch(SQLException e) { 
    		e.printStackTrace();
    		this.dispose();			
		}

        if(g != null ) { 
        	try {
        		bssp.setGenome(g);
        	} catch (Exception ex) {
        		ex.printStackTrace();
        		this.dispose();
        	}
        }
        
        this.list = list;
        this.genome = g;
        
        JPanel buttonPanel = new JPanel();
        JTabbedPane tabbed = new JTabbedPane();

        buttonPanel.setLayout(new GridBagLayout());
        Dimension buttonSize = new Dimension(30,20);
        ok = new JButton("OK");
        cancel = new JButton("Cancel");
        ok.setMaximumSize(buttonSize);
        cancel.setMaximumSize(buttonSize);
        ok.addActionListener(this);
        cancel.addActionListener(this);
        buttonPanel.add(ok);
        buttonPanel.add(cancel);

        tabbed.addTab("ChIP-Chip", bssp);
        tabbed.addTab("ChIP-Seq",csasp);

        Container content = getContentPane();
        content.setLayout(new BorderLayout());
        content.add(buttonPanel,BorderLayout.SOUTH);
        content.add(tabbed,BorderLayout.CENTER);
        setSize(500,400);
        setLocation(50,50);
        setVisible(true);
    }

    public void actionPerformed (ActionEvent e) {
        if (e.getSource() == ok) {
            Collection<BindingScan> scans = bssp.getObjects();
            Collection<SeqAnalysis> analyses = csasp.getObjects();
            Populator p = new Populator(scans,analyses,genome);
            Thread t = new Thread(p);
            t.start();
            this.dispose();
        } else if (e.getSource() == cancel) {
            this.dispose();
        }
    }

    private class Populator implements Runnable {
        
        private Collection<BindingScan> scans;
        private Collection<SeqAnalysis> analyses;
        private Genome genome;
        
        public Populator(Collection<BindingScan> scans, Collection<SeqAnalysis> analyses, Genome g) {
            this.scans = scans;
            this.analyses = analyses;
            genome = g;
        }
        
        public void run() {
            try {
                BindingScanLoader loader = new BindingScanLoader();

                Collection<BindingExpander> expanders = new ArrayList<BindingExpander>();
                for (BindingScan scan : scans) {
                    BindingExpander expander = new BindingExpander(loader,scan);
                    expanders.add(expander);
                }

                ChromRegionIterator iter = new ChromRegionIterator(genome);                
                while (iter.hasNext()) {
                    Region chrom = iter.next();
                    for (BindingExpander expander : expanders) {
                        Iterator<BindingEvent> events = expander.execute(chrom);
                        while (events.hasNext()) {
                            list.addRegion(events.next());
                        }
                    }
                    for (SeqAnalysis analysis : analyses) {
                        for (SeqAnalysisResult r : analysis.getResults(genome, chrom)) {
                            list.addRegion(new BindingEvent(r,
                                                            r.foldEnrichment,
                                                            r.pvalue,
                                                            analysis.getProgramName()));
                        }
                    }
                }                               
                loader.close();

            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }

    }
}

