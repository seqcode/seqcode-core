package edu.psu.compbio.seqcode.gse.seqview;


import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JTextField;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.seqview.components.BindingScanSelectFrame;
import edu.psu.compbio.seqcode.gse.seqview.components.RegionListPanel;
import edu.psu.compbio.seqcode.gse.seqview.components.RegionPanel;
import edu.psu.compbio.seqcode.gse.seqview.components.SaveRegionsAsFasta;
import edu.psu.compbio.seqcode.gse.seqview.components.SeqViewOptionsFrame2;
import edu.psu.compbio.seqcode.gse.seqview.components.SeqViewOptionsPane;
import edu.psu.compbio.seqcode.gse.seqview.components.SeqViewStatusBar;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.viz.DynamicAttribute;

/**
 * SeqView is a driver class that derives a little from WarpDrive's RegionFrame. 
 * 
 * It aims to be the top-level class, running the others, whereas WarpOptionsFrame
 * was the top-level class in WarpDrive.
 * This class assumes you will view one genome version at a time (RegionFrame did not).
 * If you switch genomes, the current session should be wiped out.
 * 
 * @author mahony
 *
 */
public class SeqView extends JFrame {
	protected SeqViewOptions options=new SeqViewOptions();
	protected SeqViewOptionsFrame2 optionsFrame;
	protected SeqViewOptionsPane optionsPane;
	protected RegionPanel regPanel=null;
	protected SeqViewStatusBar statusBar;
	protected boolean imageraster;
	protected int imageheight = 1200, imagewidth = 1600;
	
    public SeqView(String[] args) throws NotFoundException, SQLException, IOException {
    	//Set up the browser window
    	setSize(600,400);
        setLocation(50,50);
    	setJMenuBar(createDefaultJMenuBar());
    	
    	this.setLayout(new BorderLayout(0,0));
    	statusBar = new SeqViewStatusBar();
        statusBar.setPreferredSize(new Dimension(getWidth(), SeqViewOptions.STATUS_BAR_HEIGHT));
        statusBar.updateStatus("Loading genome", Color.orange);
        setTitle("Loading genome information...");
        add(statusBar, BorderLayout.PAGE_END);
    	
        //Load command-line options
        options = SeqViewOptions.parseCL(args);
        setTitle(options.genome.getSpecies() + " " + options.genome.getVersion());

        regPanel = new RegionPanel(options); 
    	add(regPanel, BorderLayout.CENTER);
    	
        setVisible(true);

        imageraster = true;
        statusBar.updateStatus("Genome loaded", Color.green);
        
        //Initiate options pane & frame
        optionsPane = new SeqViewOptionsPane(options);
        optionsFrame = new SeqViewOptionsFrame2(optionsPane, this);
        
        //Close operations
        this.setDefaultCloseOperation(EXIT_ON_CLOSE);
        this.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent arg0) {
                if(regPanel!=null)
                	regPanel.handleWindowClosing();
                if(!optionsFrame.isClosed()){
                	optionsFrame.close();
                }
                System.err.println("WindowClosing: " + arg0.toString());
            }
        });
    }
    
    public RegionPanel getRegionPanel() { return regPanel; }
    
    public SeqViewOptions getOptions(){return options;}
    
    /**
     * Update the viewer's options, changing genomes if required. 
     * @param opts
     */
    public void updateOptions(SeqViewOptions opts){
    	if(regPanel ==null){
    		//regPanel should have been created in the constructor, so if this is true, something is wrong. 
    		//Nevertheless:
    		regPanel = new RegionPanel(opts); 
        	add(regPanel, BorderLayout.PAGE_START);
    	}else if ( options==null  || !options.genome.equals(opts.genome)) {
            //Overwrite panel for different genome
        	setTitle(opts.genome.getSpecies() + " " + opts.genome.getVersion());
            regPanel.reinit(opts);
        } else {
        	if (options != null && opts.genome.equals(options.genome)) {
        		SeqViewOptions diffopts = opts.clone();
        		diffopts.differenceOf(options);
        		regPanel.addPaintersFromOpts(diffopts);
            }
        }
        options=opts;
    	this.repaint();
    }
    
    private JMenuBar createDefaultJMenuBar() { 
        JMenuBar jmb = new JMenuBar();
        JMenu filemenu, imagemenu, navigationmenu, displaymenu, toolsmenu; 
        JMenuItem item;
        final SeqView thisviewer = this;
        final RegionPanel thispanel = regPanel;
        jmb.add((filemenu = new JMenu("File")));
        filemenu.add((item = new JMenuItem("Configure Tracks")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {                    
                    try {
                    	optionsPane = new SeqViewOptionsPane(options);
                    	optionsFrame = new SeqViewOptionsFrame2(optionsPane, thisviewer);
                    } catch (NotFoundException e1) {
                        e1.printStackTrace();
                    }
                }
            });
        filemenu.add((item = new JMenuItem("Add a track from a file")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {    
                	if(thispanel!=null && !thispanel.isClosed())
                		thispanel.addTrackFromFile();
                }
            });
        filemenu.add((item = new JMenuItem("Save Region as FASTA")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    new SaveRegionsAsFasta(thispanel.getRegion());
                }
            });
        filemenu.addSeparator();
        filemenu.add((item = new JMenuItem("Exit")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) { 
                	if(thispanel!=null && !thispanel.isClosed())
                		thispanel.close();
                    System.exit(0);
                }
            });        

        jmb.add((imagemenu = new JMenu("Image")));
        imagemenu.add((item = new JMenuItem("Save Image")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    JFileChooser chooser;
                    chooser = new JFileChooser(new File(System.getProperty("user.dir")));
                    
                    int v = chooser.showSaveDialog(null);
                    if(v == JFileChooser.APPROVE_OPTION) { 
                        File f = chooser.getSelectedFile();
                        boolean raster = imageraster;
                        if (f.getAbsolutePath().matches(".svg")) {
                            raster = false;
                        }
                        try {
                            regPanel.saveImage(f,imagewidth,imageheight,raster);
                        } catch (IOException ex) {
                            ex.printStackTrace();
                        }
                    }
                } 
            });
        imagemenu.add((item = new JMenuItem("Image Settings")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    new ImageConfigurationFrame(thisviewer);
                }
            });
        jmb.add((navigationmenu = new JMenu("Navigation")));
        
        
        navigationmenu.add((item = new JMenuItem("Binding Event List")));
        item.addActionListener(new ActionListener()  {
                public void actionPerformed(ActionEvent e) {
                    RegionListPanel rlp = new RegionListPanel(thispanel,null);
                    RegionListPanel.makeFrame(rlp);
                    new BindingScanSelectFrame(thispanel.getGenome(),rlp);
                }
            });               

        navigationmenu.add((item = new JMenuItem("Open Region List")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    JFileChooser chooser;
                    chooser = new JFileChooser(new File(System.getProperty("user.dir")));
                    
                    int v = chooser.showOpenDialog(null);
                    if(v == JFileChooser.APPROVE_OPTION) { 
                        File f = chooser.getSelectedFile();
                        java.util.List<Region> regions = RegionPanel.readRegionsFromFile(regPanel.getGenome(),f.getAbsolutePath());
                        RegionListPanel p = new RegionListPanel(regPanel,
                                                                regions);
                        RegionListPanel.makeFrame(p);
                    }
                    
                }
            });
        navigationmenu.add((item = new JMenuItem("New Region List")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    RegionListPanel p = new RegionListPanel(regPanel,
                                                            new ArrayList<Region>());
                    RegionListPanel.makeFrame(p);
                }
            });
        ButtonGroup group = new ButtonGroup();
        jmb.add(displaymenu = new JMenu("Display"));
        displaymenu.add((item = new JRadioButtonMenuItem("Screen")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    DynamicAttribute.getGlobalAttributes().setType(DynamicAttribute.SCREEN);
                    thisviewer.repaint();
                }
            });
        item.setSelected(true);
        group.add(item);
        displaymenu.add((item = new JRadioButtonMenuItem("Display")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    DynamicAttribute.getGlobalAttributes().setType(DynamicAttribute.DISPLAY);
                    thisviewer.repaint();
                }
            });
        group.add(item);
        displaymenu.add((item = new JRadioButtonMenuItem("Print")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    DynamicAttribute.getGlobalAttributes().setType(DynamicAttribute.PRINT);
                    thisviewer.repaint();
                }
            });
        group.add(item);
        
        //TODO: uncomment this to bring back the WeightMatrix browser option
        //jmb.add(new SeqViewToolsMenu(panel));

        return jmb;
    }
    
    class ImageConfigurationFrame extends JFrame implements ActionListener {
        private SeqView parent;
        private JCheckBox rasterbox;
        private JTextField widthfield, heightfield;
        private JButton okbutton, cancelbutton;

        public ImageConfigurationFrame(SeqView p) {
            parent = p;
            JLabel boxlabel = new JLabel("Configure Save-as-image");
            rasterbox = new JCheckBox("Raster Image?",parent.imageraster);
            JLabel widthlabel = new JLabel("Width");
            JLabel heightlabel = new JLabel("Height");
            widthfield = new JTextField(Integer.toString(parent.imagewidth));
            heightfield = new JTextField(Integer.toString(parent.imageheight));
            okbutton = new JButton("OK");
            cancelbutton = new JButton("Cancel");
            okbutton.addActionListener(this);
            cancelbutton.addActionListener(this);
            
            JPanel toppanel = new JPanel();
            toppanel.setLayout(new BorderLayout());
            
            JPanel buttonpanel = new JPanel();
            buttonpanel.add(okbutton);
            buttonpanel.add(cancelbutton);
            toppanel.add(buttonpanel,BorderLayout.SOUTH);

            JPanel infopanel = new JPanel();
            infopanel.setLayout(new BorderLayout());
            infopanel.add(rasterbox,BorderLayout.NORTH);
            JPanel textpanel = new JPanel();
            textpanel.setLayout(new GridLayout(2,2));
            textpanel.add(widthlabel);
            textpanel.add(widthfield);
            textpanel.add(heightlabel);
            textpanel.add(heightfield);
            infopanel.add(textpanel,BorderLayout.CENTER);
            
            toppanel.add(infopanel,BorderLayout.CENTER);

            getContentPane().add(toppanel);
            setMinimumSize(new Dimension(150,150));
            setSize(getPreferredSize());
            pack();
            setVisible(true);
        }

        public void actionPerformed (ActionEvent e) {
            if (e.getSource() == okbutton) {
                parent.imageraster = rasterbox.isSelected();
                try {
                    parent.imageheight = Integer.parseInt(heightfield.getText());
                } catch (NumberFormatException ex) {
                }
                try {
                    parent.imagewidth = Integer.parseInt(widthfield.getText());
                } catch (NumberFormatException ex) {
                }


                this.dispose();
            } else if (e.getSource() == cancelbutton) {
                this.dispose();
            }
        }     
   
    }

	
    /** 
     * The main driver for SeqView
     * @param args
     * @throws NotFoundException
     * @throws SQLException
     * @throws IOException
     */
	public static void main(String args[]) throws NotFoundException, SQLException, IOException {
        try{
        	SeqView viewer = new SeqView(args);
        } catch (Exception e) {
			e.printStackTrace();
		}
    }
}
