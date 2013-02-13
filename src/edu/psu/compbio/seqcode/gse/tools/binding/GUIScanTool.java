package edu.psu.compbio.seqcode.gse.tools.binding;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.sql.SQLException;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.TitledBorder;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingParameters;
import edu.psu.compbio.seqcode.gse.datasets.binding.BindingScanLoader;
import edu.psu.compbio.seqcode.gse.datasets.locators.BayesLocator;
import edu.psu.compbio.seqcode.gse.datasets.locators.ChipChipLocator;
import edu.psu.compbio.seqcode.gse.datasets.locators.ExptLocator;
import edu.psu.compbio.seqcode.gse.datasets.locators.MSPLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.binding.CallerMapper;
import edu.psu.compbio.seqcode.gse.utils.Listener;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;
import edu.psu.compbio.seqcode.gse.viz.components.ChipChipExptSelectPanel;
import edu.psu.compbio.seqcode.gse.viz.utils.GenomeSelectPanel;

public class GUIScanTool extends JFrame implements Listener<ActionEvent> {

	private static final long serialVersionUID = 1L;

	public static void main(String[] args) { 
		try {
			GUIScanTool tool = new GUIScanTool();
			tool.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		} catch (SQLException e) {
			e.printStackTrace(System.err);
		} catch (UnknownRoleException e) {
			e.printStackTrace(System.err);
		}
	}
	
	private BindingScanLoader loader;
	private JButton scanButton;	
	private GenomeSelectPanel genomeSelector;
	private ChipChipExptSelectPanel selectPanel;
    private JTextField paramsField;

	public GUIScanTool() throws SQLException, UnknownRoleException { 
		super("Scan Tool");
		
		loader = new BindingScanLoader();
		
		Container c = (Container)getContentPane();
		c.setLayout(new BorderLayout());
		c.add(genomeSelector = new GenomeSelectPanel("Homo sapiens", "hg17"), 
				BorderLayout.NORTH);
		c.add(selectPanel = new ChipChipExptSelectPanel(getGenome()), 
				BorderLayout.CENTER);
		
		JPanel scanPanel = new JPanel(); 
		scanPanel.setLayout(new GridLayout(1, 1));
		scanPanel.add(scanButton = new JButton("Scan"));
        
        paramsField = new JTextField();
        JPanel paramsPanel = new JPanel(); paramsPanel.setLayout(new BorderLayout());
        paramsPanel.add(paramsField, BorderLayout.NORTH);
        paramsPanel.setBorder(new TitledBorder("Parameters"));
        
        JPanel bottomPanel = new JPanel();
        bottomPanel.setLayout(new BorderLayout());
        bottomPanel.add(scanPanel, BorderLayout.SOUTH);
        bottomPanel.add(paramsPanel, BorderLayout.NORTH);
		
        c.add(bottomPanel, BorderLayout.SOUTH);
        
        scanButton.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) { 
				try {
					scan();
				} catch (SQLException e1) {
					e1.printStackTrace(System.err);
				} catch (UnknownRoleException e1) {
					e1.printStackTrace(System.err);
				}
			}
		});
        
        selectPanel.setGenome(getGenome());
        genomeSelector.addEventListener(this);
		
		setVisible(true);
		pack();
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
	}
	
	public void scan() throws SQLException, UnknownRoleException { 
		while(selectPanel.getNumSelected() > 0) {
            
            BindingParameters params = new BindingParameters(paramsField.getText());
			ExptLocator loc = selectPanel.removeFirstSelected();
			CallerMapper cm = null;

			if(loc instanceof ChipChipLocator) {
			    if(params.containsKey("type") && params.get("type").equals("hmm")) { 
			        cm = new ScanTool.HMMCallerMapper(params);
			    } else { 
			        cm = new ScanTool.ChipChipCallerMapper(params);
			    }
			} else if (loc instanceof MSPLocator) {
			    if(params.containsKey("type") && params.get("type").equals("simple")) { 
			        cm = new ScanTool.SimpleMSPCallerMapper(params);
			    } else { 
			        cm = new ScanTool.MSPCallerMapper(params);
			    }
			} else if (loc instanceof BayesLocator) { 
			    cm = new ScanTool.JBDCallerMapper(params);
			}

			if(cm != null && params.containsKey("type") && params.get("type").equals("domain")) { 
			    cm = new ScanTool.DomainCallerMapper(params, cm);
			}

			if(cm != null) { 
			    ScanTool scanner = new ScanTool(getGenome(), loader, cm);
			    scanner.runScan(loc);
			}
		}
	}
	
	public Genome getGenome() {
		try { 
			return Organism.findGenome(genomeSelector.getGenome());
		} catch(NotFoundException nfe) { 
			throw new RuntimeException("Unknown genome", nfe);
		}
	}

	public void eventRegistered(ActionEvent e) {
	    selectPanel.setGenome(getGenome());
	    System.out.println("New Genome: " + getGenome().getVersion());
	}
}
