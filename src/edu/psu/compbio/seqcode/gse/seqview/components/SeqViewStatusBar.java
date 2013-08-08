package edu.psu.compbio.seqcode.gse.seqview.components;

import java.awt.Font;
import java.awt.Graphics;
import java.util.EventObject;

import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.BevelBorder;

import edu.psu.compbio.seqcode.gse.utils.Listener;

public class SeqViewStatusBar extends JPanel implements Listener<EventObject>{

	private static final long serialVersionUID = 1986253667638898752L;
	private JLabel statusLabel;
	private SeqViewStatus status;
	
	public SeqViewStatusBar(SeqViewStatus status){
		statusLabel =new JLabel("");
		statusLabel.setFont(new Font("SansSerif", Font.PLAIN, 10));
		this.status = status;
		status.addEventListener(this);
		statusLabel.setText(status.getMessage());
		statusLabel.setForeground(status.getColor());
		
		this.setBorder(new BevelBorder(BevelBorder.LOWERED));
		this.setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
		add(statusLabel);
	}
	
	public void paintComponent(Graphics g) {
		super.paintComponent(g);  
	}
	
	public synchronized void eventRegistered(EventObject e) {
		if (e.getSource() == status){
			statusLabel.setText(status.getMessage());
			statusLabel.setForeground(status.getColor());
			repaint();
		}
	}
}
