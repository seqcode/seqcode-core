package edu.psu.compbio.seqcode.gse.seqview;


import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.sql.SQLException;
import java.util.ArrayList;

import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JTextField;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.seqview.components.RegionListPanel;
import edu.psu.compbio.seqcode.gse.seqview.components.RegionPanel;
import edu.psu.compbio.seqcode.gse.seqview.components.SaveRegionsAsFasta;
import edu.psu.compbio.seqcode.gse.seqview.components.SeqViewOptionsFrame2;
import edu.psu.compbio.seqcode.gse.seqview.components.SeqViewOptionsPane;
import edu.psu.compbio.seqcode.gse.seqview.components.SeqViewStatus;
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
	protected SeqViewStatus status;
	protected SeqViewStatusBar statusBar;
	protected boolean imageraster;
	private File currDirectory = new File(System.getProperty("user.home"));
	
    public SeqView(String[] args) throws NotFoundException, SQLException, IOException {
    	//Set up the browser window
    	setSize(700,500);
        setLocation(50,50);
    	setJMenuBar(createDefaultJMenuBar());
    	imageraster = true;
        
    	this.setLayout(new BorderLayout(0,0));
    	status = new SeqViewStatus();
    	status.setStatus("Loading genome", Color.red);
    	statusBar = new SeqViewStatusBar(status);
        statusBar.setPreferredSize(new Dimension(getWidth(), SeqViewOptions.STATUS_BAR_HEIGHT));
        setTitle("Loading genome information...");
        this.setBackground(Color.black);
        add(statusBar, BorderLayout.PAGE_END);
        setVisible(true);
        
        //Load command-line options
        options = SeqViewOptions.parseCL(args);
        setTitle(options.genome.getSpecies() + " " + options.genome.getVersion());
        regPanel = new RegionPanel(options, status, currDirectory); 
    	add(regPanel, BorderLayout.CENTER);
    	status.setStatus("Genome loaded... loading experiment list", Color.red);
        
        //Initiate options pane & frame
        optionsPane = new SeqViewOptionsPane(options);
        optionsFrame = new SeqViewOptionsFrame2(optionsPane, this);
        status.setStatus("Genome loaded", Color.black);
        
        //Close operations
        this.setDefaultCloseOperation(EXIT_ON_CLOSE);
        final SeqView thisviewer = this;
        this.addWindowListener(new WindowAdapter() {
        	public void windowClosing(WindowEvent arg0) {
                if(regPanel!=null && !regPanel.isClosed())
                	regPanel.handleWindowClosing();
                if(!optionsFrame.isClosed()){
                	optionsFrame.close();
                }
                thisviewer.dispose();
            }
        });
    }
    
    public String getVersion(){
    	return "Version 0.1, Aug-2013";
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
    		regPanel = new RegionPanel(opts, status, currDirectory); 
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
        JMenu filemenu, editmenu, imagemenu, navigationmenu, displaymenu, aboutmenu, toolsmenu; 
        JMenuItem item;
        final SeqView thisviewer = this;
        jmb.add((filemenu = new JMenu("File")));
        filemenu.add((item = new JMenuItem("Add Data Tracks")));
        item.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {                    
        		try {
        			optionsPane = new SeqViewOptionsPane(regPanel.getCurrentOptions());
        			optionsFrame = new SeqViewOptionsFrame2(optionsPane, thisviewer);
        		} catch (NotFoundException e1) {
        			e1.printStackTrace();
        		}
        	}
        });
        filemenu.add((item = new JMenuItem("Add Track From Region File")));
        item.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {    
        		if(thisviewer.getRegionPanel()!=null && !thisviewer.getRegionPanel().isClosed())
        			thisviewer.getRegionPanel().addTrackFromFile(false);
        	}
        });
        filemenu.add((item = new JMenuItem("Add Track From GFF File")));
        item.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {    
        		if(thisviewer.getRegionPanel()!=null && !thisviewer.getRegionPanel().isClosed())
        			thisviewer.getRegionPanel().addTrackFromFile(true);
        	}
        });
        filemenu.add((item = new JMenuItem("Save Region as FASTA")));
        item.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {
        		new SaveRegionsAsFasta(thisviewer.getRegionPanel().getRegion());
        	}
        });
        filemenu.addSeparator();
        filemenu.add((item = new JMenuItem("Reconnect To Database")));
        item.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {
        		if(thisviewer.getRegionPanel()!=null && !thisviewer.getRegionPanel().isClosed())
        			thisviewer.getRegionPanel().reconnectModels();
        	}
        });
        filemenu.addSeparator();
        filemenu.add((item = new JMenuItem("Exit")));
        item.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) { 
        		if(thisviewer.getRegionPanel()!=null && !thisviewer.getRegionPanel().isClosed())
        			thisviewer.getRegionPanel().close();
        		System.exit(0);
        	}
        });        

        jmb.add((editmenu = new JMenu("Edit")));
        editmenu.add((item = new JMenuItem("Configure All SeqData Tracks")));
        item.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {
        		thisviewer.getRegionPanel().configSeqDataTracksBatch();
        	}
        });
        
        jmb.add((imagemenu = new JMenu("Image")));
        imagemenu.add((item = new JMenuItem("Save Image")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    new ImageConfigurationFrame(thisviewer);
                }
            });
        jmb.add((navigationmenu = new JMenu("Navigation")));
        
        
        /*navigationmenu.add((item = new JMenuItem("Binding Event List")));
        item.addActionListener(new ActionListener()  {
                public void actionPerformed(ActionEvent e) {
                    RegionListPanel rlp = new RegionListPanel(thisviewer.getRegionPanel(),null);
                    RegionListPanel.makeFrame(rlp);
                    new BindingScanSelectFrame(thisviewer.getRegionPanel().getGenome(),rlp);
                }
            });               
		*/
        navigationmenu.add((item = new JMenuItem("Open Region List")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    JFileChooser chooser;
                    chooser = new JFileChooser(currDirectory);
                    int v = chooser.showOpenDialog(null);
                    if(v == JFileChooser.APPROVE_OPTION) { 
                        File f = chooser.getSelectedFile();
                        currDirectory = chooser.getCurrentDirectory();
                        java.util.List<Region> regions = RegionPanel.readRegionsFromFile(regPanel.getGenome(),f.getAbsolutePath(),false);
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
        //jmb.add(toolsmenu = new SeqViewToolsMenu(panel));

        jmb.add((aboutmenu = new JMenu("Help")));
        aboutmenu.add((item = new JMenuItem("About")));
        item.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {
        		new AboutMessage(thisviewer.getVersion());
        	}
        });
        
        return jmb;
    }
    
    class ImageConfigurationFrame extends JFrame implements ActionListener {
        private SeqView parent;
        private ButtonGroup rasterGroup;
        private JRadioButton pngButton, svgButton;
        private JTextField widthfield, heightfield;
        private JButton savebutton, cancelbutton;
        private int imagewidth, imageheight;

        public ImageConfigurationFrame(SeqView p) {
            parent = p;
            rasterGroup = new ButtonGroup();
            pngButton = new JRadioButton("PNG");
            pngButton.setMnemonic(KeyEvent.VK_P);
            pngButton.setSelected(parent.imageraster);
            svgButton = new JRadioButton("SVG");
            svgButton.setMnemonic(KeyEvent.VK_S);
            svgButton.setSelected(!parent.imageraster);
            rasterGroup.add(pngButton);
            rasterGroup.add(svgButton);
            JLabel widthlabel = new JLabel("Width");
            JLabel heightlabel = new JLabel("Height");
            imagewidth = parent.getRegionPanel().getWidth();
            imageheight = parent.getRegionPanel().getHeight();
            widthfield = new JTextField(Integer.toString(imagewidth));
            heightfield = new JTextField(Integer.toString(imageheight));
            savebutton = new JButton("Save");
            cancelbutton = new JButton("Cancel");
            savebutton.addActionListener(this);
            cancelbutton.addActionListener(this);
            
            JPanel toppanel = new JPanel();
            toppanel.setLayout(new BorderLayout());
            
            JPanel buttonpanel = new JPanel();
            buttonpanel.add(savebutton);
            buttonpanel.add(cancelbutton);
            toppanel.add(buttonpanel,BorderLayout.SOUTH);

            JPanel infopanel = new JPanel();
            infopanel.setLayout(new BorderLayout());
            JPanel formatpanel = new JPanel(new GridLayout(1, 0));
            formatpanel.add(pngButton);
            formatpanel.add(svgButton);
            infopanel.add(formatpanel,BorderLayout.NORTH);
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
            if (e.getSource() == savebutton) {
                parent.imageraster = pngButton.isSelected();
                
                JFileChooser chooser;
                chooser = new JFileChooser(currDirectory);
                
                int v = chooser.showSaveDialog(null);
                if(v == JFileChooser.APPROVE_OPTION) { 
                    File f = chooser.getSelectedFile();
                    currDirectory = chooser.getCurrentDirectory();
                    boolean raster = imageraster;
                    if (f.getAbsolutePath().matches(".svg")) {
                        raster = false;
                    }else if (f.getAbsolutePath().matches(".png")) {
                        raster = true;
                    }
                    
                    try {
                        regPanel.saveImage(f,imagewidth,imageheight,raster);
                    } catch (IOException ex) {
                        ex.printStackTrace();
                    }
                }
            
                this.dispose();
            } else if (e.getSource() == cancelbutton) {
                this.dispose();
            }
        }     
   
    }

    class AboutMessage extends JFrame implements ActionListener {
        private JButton okbutton;

        public AboutMessage(String version) {
        	JLabel iconlabel=null;
        	try{
        		URL picURL = getClass().getResource("/edu/psu/compbio/seqcode/gse/seqview/SeqViewLogo.png");
        		if(picURL!=null){
	        		ImageIcon icon = new ImageIcon(picURL, "SeqView");
	                iconlabel = new JLabel(icon);
	                iconlabel.setPreferredSize(new Dimension(249,161));
        		}
                String text = "<html><center>"+
            			"<p>"+version+
            			"<p>"+
            			"<p>Written by Shaun Mahony"+
            			"<p>"+
            			"<p>Based on the GSE library,"+
            			"<p>by Alex Rolfe & Tim Danford"+
            			"<p>(Gifford Lab, CSAIL, MIT)"+
            			"</html>";
            			
            	javax.swing.JLabel messagelabel = new javax.swing.JLabel(text);
                okbutton = new JButton("OK");
                okbutton.addActionListener(this);
                
                JPanel toppanel = new JPanel();
                toppanel.setLayout(new GridBagLayout());
                
                if(iconlabel!=null)
                	toppanel.add(iconlabel);
                else{
                	javax.swing.JLabel svlabel = new javax.swing.JLabel(new String("<html><center><font size=\"+2\" color=\"blue\">SeqView</font></html>"));
                	toppanel.add(svlabel);
                }
                toppanel.add(messagelabel);
                
                JPanel buttonpanel = new JPanel();
                buttonpanel.add(okbutton);
                toppanel.add(buttonpanel);

                getContentPane().add(toppanel);
                pack();
                setMinimumSize(new Dimension(300,250));
                setPreferredSize(new Dimension(300,250));
                setLocation(100,100);
                setVisible(true);
        	}  catch (NullPointerException e) {
        		e.printStackTrace();
        	}
        }

        public void actionPerformed (ActionEvent e) {
            if (e.getSource() == okbutton) {
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
