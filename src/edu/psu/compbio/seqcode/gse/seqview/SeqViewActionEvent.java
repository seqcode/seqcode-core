/*
 * Created on Mar 21, 2006
 */
package edu.psu.compbio.seqcode.gse.seqview;

import java.util.EventObject;

/**
 * @author tdanford
 */
public class SeqViewActionEvent extends EventObject {
    
    public static final int MOVE_LEFT = 0;
    public static final int MOVE_RIGHT = 1;
    public static final int ZOOM_IN = 2;
    public static final int ZOOM_OUT = 3;
    
    public static final int NEXT_CHROM = 4;
    public static final int PREV_CHROM = 5;
    
    private Object data;
    private int type;

    /**
     * @param src The source of the event (Swing component, usually).
     * @param t The "type" of event.
     * @param d Any associated data, to the event.
     */
    public SeqViewActionEvent(Object src, int t, Object d) {
        super(src);
        type = t;
        data = d;
    }
    
    /**
     * @param src The source of the event (Swing component, usually).
     * @param t The "type" of event.
     */    
    public SeqViewActionEvent(Object src, int t) { 
        super(src);
        type = t;
        data = null;
    }
    
    public int getType() { return type; }
    public Object getData() { return data; }
}
