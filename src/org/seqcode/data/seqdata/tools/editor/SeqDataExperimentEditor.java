package org.seqcode.data.seqdata.tools.editor;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.IOException;
import java.net.URL;
import java.sql.SQLException;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;

import org.seqcode.gseutils.NotFoundException;


public class SeqDataExperimentEditor extends JFrame implements ActionListener {
	/**
	 * Constants for default values and limits on settings
	 */
	public static final int DEFAULT_WINDOW_WIDTH = 900;
	public static final int DEFAULT_WINDOW_HEIGHT = 650;
	public static final int DEFAULT_TOP_LEFT_X = 50;
	public static final int DEFAULT_TOP_LEFT_Y = 50;
		
	//Program variables
	private SeqDataExptEditPane pane=null; 
    private JMenuBar menuBar;
    private JMenu fileMenu;
    private JMenuItem exitItem;
    private JMenuItem aboutItem;

	public SeqDataExperimentEditor() throws NotFoundException, SQLException, IOException {
		setTitle("SeqDataExperimentEditor");
        pane = new SeqDataExptEditPane();
        init();
        
        //Close operations
        this.setDefaultCloseOperation(EXIT_ON_CLOSE);
        final SeqDataExperimentEditor thiseditor = this;
        this.addWindowListener(new WindowAdapter() {
        	public void windowClosing(WindowEvent arg0) {
                if(pane!=null && !pane.isClosed())
                	pane.close();
                thiseditor.dispose();
            }
        });
    }

    private void init() {
    	setSize(DEFAULT_WINDOW_WIDTH,DEFAULT_WINDOW_HEIGHT);
        setLocation(DEFAULT_TOP_LEFT_X,DEFAULT_TOP_LEFT_Y);
    	
    	Container content = getContentPane();
        content.setLayout(new BorderLayout());
        content.add(pane,BorderLayout.CENTER);
        
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
    	
    	exitItem = new JMenuItem("Exit", KeyEvent.VK_X);
    	exitItem.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {
        		if(pane!=null && !pane.isClosed())
                	pane.close();
                dispose();
        		System.exit(0);
        	}
        });        
    	
    	aboutItem = new JMenuItem("About");
    	aboutItem.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {
        		new AboutMessage(getVersion());
        	}
        });

    	fileMenu.add(aboutItem);
    	fileMenu.add(exitItem);
    	//end building file menu
    	
    	this.setJMenuBar(menuBar);
    }
	
    public String getVersion(){
    	return "Version 0.1, May-2014";
    }

   
    
	class AboutMessage extends JFrame implements ActionListener {
        private JButton okbutton;

        public AboutMessage(String version) {
        	JLabel iconlabel=null;
        	try{
        		URL picURL = getClass().getResource("/edu/psu/compbio/seqcode/gse/seqview/SeqViewLogo.png");
        		if(picURL!=null){
	        		ImageIcon icon = new ImageIcon(picURL, "SeqDataExperimentEditor");
	                iconlabel = new JLabel(icon);
	                iconlabel.setPreferredSize(new Dimension(249,161));
        		}
                String text = "<html><center>"+
            			"<p>"+version+
            			"<p>"+
            			"<p>Written by Shaun Mahony"+
            			"<p>"+
            			"<p>An editor for SeqData DB experiment entries. "+
            			"</html>";
            			
            	javax.swing.JLabel messagelabel = new javax.swing.JLabel(text);
                okbutton = new JButton("OK");
                okbutton.addActionListener(this);
                
                JPanel toppanel = new JPanel();
                toppanel.setLayout(new GridBagLayout());
                
                if(iconlabel!=null)
                	toppanel.add(iconlabel);
                else{
                	javax.swing.JLabel svlabel = new javax.swing.JLabel(new String("<html><center><font size=\"+2\" color=\"blue\">SeqDataExperimentEditor</font></html>"));
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
     * The main driver for SeqDataExperimentEditor
     * @param args
     * @throws NotFoundException
     * @throws SQLException
     * @throws IOException
     */
	public static void main(String args[]) throws NotFoundException, SQLException, IOException {
        try{
        	SeqDataExperimentEditor editor = new SeqDataExperimentEditor();
        } catch (Exception e) {
			e.printStackTrace();
		}
    }

	public void actionPerformed(ActionEvent e) {
		repaint();
	}
}
