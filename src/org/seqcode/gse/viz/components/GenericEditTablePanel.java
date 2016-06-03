package org.seqcode.gse.viz.components;

import java.util.*;
import java.awt.*;

import javax.swing.*;

import java.awt.event.*;
import javax.swing.event.*;
import javax.swing.table.*;

import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.gse.datasets.core.ExptType;
import org.seqcode.gse.utils.*;


/**
 * GenericTablePanel is a single pane panel that holds a table list
 * 
 * @author mahony
 *
 * @param <X>
 */
public abstract class GenericEditTablePanel<X> extends JPanel implements Closeable, Runnable {

    private JButton editButton, filterButton;
    protected JComboBox speciesCBox;
    
    protected JTable filteredList;
    protected ObjectTableModel<X> filteredModel;
    protected JPanel buttonPanel;
    private boolean handlingUpdate, dataReady;
    private Thread thread;

    public GenericEditTablePanel() {
        super();
        handlingUpdate = false;
        dataReady = false;
    }

    public GenericEditTablePanel(Genome g) {
        super();
        handlingUpdate = false;
        dataReady = false;
    }
    /**
     * UpdateData
     *  Starts a new thread to update data.
       Database work should happen in that thread and the results should
       be cached.  Swing components should be updated with a call
       to SwingUtilities.invokeLater().  Since this also invokes run(),
       the run() method must determine the current operation by querying:
       handlingNewGenome() and dataReady()
    */
    public void updateData() {
        if (thread != null && thread.isAlive()) {
            thread.interrupt();
            thread = null;
        }
        thread = new Thread(this);
        handlingUpdate = true;
        dataReady = false;
        thread.start();        
    }
    public boolean handlingUpdate() {return handlingUpdate;}
    public boolean dataReady() {return dataReady;}
    public void run() {
        if (handlingUpdate()) {
            retrieveData();
            handlingUpdate = false;
            dataReady = true;
            SwingUtilities.invokeLater(this);
        } else if (dataReady()) {
            updateComponents();
            dataReady = false;
        }
    }
    /* retrieve data from the database; must cache somewhere
     */
    public abstract void retrieveData();
    /* called from swing thread to update the swing components
       for the cached data */
    public abstract void updateComponents();

    /* sub-class constructors must call one of the two init() methods */
    public void init(JTable filteredList,
                     ObjectTableModel<X> filteredModel) {
        final JTableHeader header = filteredList.getTableHeader();
        final ObjectTableModel filteredmodel = filteredModel;
        header.addMouseListener(new MouseAdapter() {
                public void mouseClicked(MouseEvent e) {
                    filteredmodel.sortByColumn(header.columnAtPoint(e.getPoint()));
                }
            });

        this.filteredList = filteredList;
        this.filteredModel = filteredModel;
        editButton = new JButton("Edit");
        filterButton = new JButton("Filter");
        
        Font tableFont = new Font("SansSerif", Font.PLAIN, 11);
        filteredList.setFont(tableFont);
        
        ArrayList<String> orgs = new ArrayList<String>();
		orgs.add("");
		Collection<String> organisms = Species.getAllSpeciesNames(false);
		for (String o : organisms) { orgs.add(o); }
		Collections.sort(orgs);
        speciesCBox = new JComboBox(orgs.toArray());
        speciesCBox.setEditable(true);
        
        final GenericEditTablePanel panel = this;
        editButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    panel.edit(getFilteredForEdit());
                }
            });
        filterButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    panel.filter();
                }
            });
        speciesCBox.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                panel.filter();
            }
        });
        filteredList.addMouseListener(new MouseAdapter() {
                public void mouseClicked(MouseEvent e) {
                    panel.filterToEdit(e);
                }
            });
        JPanel actPanel = new JPanel(); actPanel.setLayout(new BorderLayout());
        actPanel.add(new JScrollPane(filteredList), BorderLayout.CENTER);
        buttonPanel = new JPanel();
        buttonPanel.setLayout(new GridBagLayout());
        buttonPanel.add(editButton);
        buttonPanel.add(filterButton);
        buttonPanel.add(new JLabel("Species:"));
        buttonPanel.add(speciesCBox);
        JPanel inputsPanel = getInputsPanel();
        JPanel allInputsPanel = new JPanel(); allInputsPanel.setLayout(new BorderLayout());
        allInputsPanel.add(inputsPanel,BorderLayout.CENTER);
        allInputsPanel.add(buttonPanel,BorderLayout.SOUTH);
        actPanel.add(allInputsPanel,BorderLayout.SOUTH);
        add(actPanel);
    }

    public void init(ObjectTableModel<X> filteredModel) {
        init(new JTable(filteredModel),
             filteredModel);        
    }
    
    public abstract JPanel getInputsPanel();
    public abstract void filter();
    public abstract void edit(Collection<X> toEdit);
    
    /* returns the elements of the filtered list that should be added 
 		when the "edit" button is pressed */
    public Collection<X> getFilteredForEdit() {
        ArrayList<X> output = new ArrayList<X>();
        int[] inds = filteredList.getSelectedRows();
        for (int i = 0; i < inds.length; i++) {
            X o = filteredModel.getObject(inds[i]);        
            output.add(o);
        }
        return output;
    }
    
    //Double click handling
    public void filterToEdit(MouseEvent e) {
        if(e.getButton() == MouseEvent.BUTTON1 && e.getClickCount() == 2) {
        	ArrayList<X> output = new ArrayList<X>();
            int row = filteredList.rowAtPoint(e.getPoint());
            X x = filteredModel.getObject(row);
            output.add(x);
            edit(output);
        }
    }
    
    public void close() {}
    public boolean isClosed() {return true;}

}
