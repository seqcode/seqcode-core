package edu.psu.compbio.seqcode.gse.seqview.components;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;

import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.BevelBorder;

public class SeqViewStatusBar extends JPanel{

	private JLabel statusLabel;
	private Color statusColor = Color.green;
	
	public SeqViewStatusBar(){
		statusLabel =new JLabel("");
		statusLabel.setFont(new Font("SansSerif", Font.PLAIN, 10));
		this.setBorder(new BevelBorder(BevelBorder.LOWERED));
		this.setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
		
		add(statusLabel);
	}
	
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		g.setColor(statusColor);
	    //g.fillOval(3, 3, getHeight()-6, getHeight()-6);  
	}
	
	public void updateStatus(String text, Color indicator){
		statusLabel.setText(text);
		statusColor = indicator;
	}
}
