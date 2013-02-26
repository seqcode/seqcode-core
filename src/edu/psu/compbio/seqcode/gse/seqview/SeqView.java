package edu.psu.compbio.seqcode.gse.seqview;

import java.awt.*;
import java.util.*;
import java.awt.event.*;

import javax.swing.*;

import edu.psu.compbio.seqcode.gse.seqview.components.PainterContainer;
import edu.psu.compbio.seqcode.gse.seqview.components.RegionFrame;
import edu.psu.compbio.seqcode.gse.seqview.components.SeqViewOptionsDialog;
import edu.psu.compbio.seqcode.gse.seqview.components.SeqViewOptionsPane;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;


public class SeqView extends JFrame implements ActionListener {

	//static Logger logger = Logger.getLogger(WarpOptionsFrame.class);
	
    private ArrayList<PainterContainer> pcs;
    private SeqViewOptionsPane pane;
    private JButton ok, cancel;  
    
    //variables for menus
    private JMenuBar menuBar;
    private JMenu fileMenu;
    private JMenuItem openSessionItem;
    private JMenuItem saveSessionItem;
    private JMenuItem exitItem;
    private JMenu toolsMenu;
    private JMenuItem optionsItem;
    
    public SeqView() throws NotFoundException {
        super();
        setTitle("SeqView");
        pcs = new ArrayList<PainterContainer>();
        pane = new SeqViewOptionsPane();
        init();
    }

    public SeqView(String species, String genome) throws NotFoundException {
        super();
        setTitle("SeqView for " + species + ", " + genome);
        pcs = new ArrayList<PainterContainer>();
        pane = new SeqViewOptionsPane(species,genome);
        init();
    }
    
    public SeqView(SeqViewOptions opts) throws NotFoundException {
        super();
        setTitle("SeqView"); 
        pcs = new ArrayList<PainterContainer>();
        pane = new SeqViewOptionsPane(opts);
        init();
    }

    private void init() {
        JPanel buttonPanel = new JPanel();
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
        Container content = getContentPane();
        content.setLayout(new BorderLayout());
        content.add(buttonPanel,BorderLayout.SOUTH);
        content.add(pane,BorderLayout.CENTER);
        
        SeqViewOptions options = pane.parseOptions();
        this.setSize(options.getPreferredWindowWidth(), options.getPreferredWindowHeight());
        if (options.isWindowCentered()) {
        	this.setLocationRelativeTo(null);	
        }
        else {
        	this.setLocation(options.getPreferredTopLeftX(), options.getPreferredTopLeftY());
        }                            
        
        this.createMenu();
        
        setVisible(true);
    }

    
    /**
     * Create a JMenuBar for this GUI
     */
    private void createMenu() {
    	menuBar = new JMenuBar();
    	
    	//build the file menu
    	fileMenu = new JMenu("File");
    	fileMenu.setMnemonic(KeyEvent.VK_F);
    	menuBar.add(fileMenu);
    	
    	openSessionItem = new JMenuItem("Open Session", KeyEvent.VK_O);
    	openSessionItem.setToolTipText("Open a saved Warp Drive session");
    	openSessionItem.addActionListener(new ActionListener() {
    		public void actionPerformed(ActionEvent e) {
    			openSession_actionPerformed(e);
    		}
    	});
    	
    	saveSessionItem = new JMenuItem("Save Session", KeyEvent.VK_S);
    	saveSessionItem.setToolTipText("Save this Warp Drive session");
    	saveSessionItem.addActionListener(new ActionListener() {
    		public void actionPerformed(ActionEvent e) {
    			saveSession_actionPerformed(e);
    		}
    	});

    	exitItem = new JMenuItem("Exit", KeyEvent.VK_X);
    	exitItem.addActionListener(new ActionListener() {
    		public void actionPerformed(ActionEvent e) {
    			exit_actionPerformed(e);
    		}
    	});

    	
    	//TODO: delete these lines once the functionality for these buttons is
    	//implemented
    	openSessionItem.setEnabled(false);
    	saveSessionItem.setEnabled(false);
    	
    	fileMenu.add(openSessionItem);
    	fileMenu.add(saveSessionItem);
    	fileMenu.addSeparator();
    	fileMenu.add(exitItem);
    	//end building file menu
    	
    	
    	//build the tools menu
    	toolsMenu = new JMenu("Tools");
    	toolsMenu.setMnemonic(KeyEvent.VK_T);
    	menuBar.add(toolsMenu);
    	
    	optionsItem = new JMenuItem("Options...", KeyEvent.VK_O);
    	optionsItem.addActionListener(new ActionListener() {
    		public void actionPerformed(ActionEvent e) {
    			options_actionPerformed(e);
    		}
    	});

    	
    	toolsMenu.add(optionsItem);
    	//end building edit menu
    	
    	this.setJMenuBar(menuBar);
    }
    
    public void addPainterContainer(PainterContainer pc) {
        pcs.add(pc);
    }

    public void actionPerformed (ActionEvent e) {
        if (e.getSource() == ok) {
            SeqViewOptions opts = pane.parseOptions();
            PainterContainer pc = null;
            for (int i = 0; i < pcs.size(); i++) {
                if (pcs.get(i).getGenome().getVersion().equals(opts.genome)) {
                    pc = pcs.get(i);
                    break;
                }
            }
            if (pc == null) {
                /* this is a bit of a hack to let us create an 
                   initial RegionFrame.  In theory, we should
                   1) be able to handle different types of Frames 
                   2) know whether it's even appropriate to create a RegionFrame 
                */
                new RegionFrame(opts);
            } else {
                opts = pane.parseAndDiff();
                pc.addPaintersFromOpts(opts);
            }
            pane.close();
            this.dispose();
        } else if (e.getSource() == cancel) {
            pane.close();
            this.dispose();
        }
    }


    /**
     * 
     * @param e
     */
    void openSession_actionPerformed(ActionEvent e) {
    	//TODO Implement this method
    }
    
    
    /**
     * 
     * @param e
     */
    void saveSession_actionPerformed(ActionEvent e) {
    	//TODO Implement this method
    }

    
    /**
     * Exit the program
     * @param e
     */
    void exit_actionPerformed(ActionEvent e) {
    	//TODO Check if any data is being viewed, if so prompt for confirmation
    	//but if only the WarpOptionFrame is open then just exit
    	boolean dataWindowsOpen = true;
    	
    	if (dataWindowsOpen) {
    		int confirmResult = 
    			JOptionPane.showConfirmDialog(this, "Are you sure you want to exit Warp Drive?", 
    					"Confirm Exit", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);

    		if (confirmResult == JOptionPane.NO_OPTION) {
    			return;
    		}
    	}
        pane.close();
        try {
            Thread.sleep(400);
        } catch (Exception ex) {

        }
    	System.exit(0);
    }

    
    /**
     * Open the dialog to set preferences
     * @param e
     */
    void options_actionPerformed(ActionEvent e) {
    	new SeqViewOptionsDialog(this, this, pane.parseOptions());
    }

    
    /**
     * Configure log4j
     */
    public static void configureLogging() {
    	ClassLoader loader = SeqView.class.getClassLoader();
    	//PropertyConfigurator.configure(loader.getResource("edu/psu/compbio/seqcode/gse/utils/config/log4j.properties"));    	
    }
    
    
    /**
     * 
     * @param args
     */
    public static void main(String args[]) {        
        try {
        	SeqView.configureLogging();
        	
            SeqViewOptions opts = SeqViewOptions.parseCL(args);
			new SeqView(opts);
		} catch (Exception e) {
			e.printStackTrace();
		}
    }
}
