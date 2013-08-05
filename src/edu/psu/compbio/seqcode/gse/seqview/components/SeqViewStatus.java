package edu.psu.compbio.seqcode.gse.seqview.components;

import java.awt.Color;
import java.util.EventObject;
import java.util.HashSet;

import edu.psu.compbio.seqcode.gse.utils.EventSource;
import edu.psu.compbio.seqcode.gse.utils.Listener;

public class SeqViewStatus implements EventSource<EventObject>{
	private String statusMessage="";
	private Color statusColor = Color.black;
	private HashSet<Listener<EventObject>> listeners;
	
	public SeqViewStatus(){
		listeners = new HashSet<Listener<EventObject>>();
	}
	
	public void setStatus(String m, Color c){
		statusMessage=m;
		statusColor=c;
		notifyListeners();
	}
	public String getMessage(){return statusMessage;}
	public Color getColor(){return statusColor;}
	
	
	public synchronized void notifyListeners(EventObject obj) {
        for (Listener<EventObject> l : listeners) {
            l.eventRegistered(obj);
        }
    }
    
    public void notifyListeners() { 
        notifyListeners(new EventObject(this));
    }
    
    public void addEventListener(Listener<EventObject> l) {
        listeners.add(l);
    }
    public void removeEventListener(Listener<EventObject> l) {
        listeners.remove(l);
    }
    public boolean hasListeners() {return listeners.size() > 0;}
}
