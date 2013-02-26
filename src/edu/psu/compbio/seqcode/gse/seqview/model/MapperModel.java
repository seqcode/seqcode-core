package edu.psu.compbio.seqcode.gse.seqview.model;

import java.util.*;

import edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper;
import edu.psu.compbio.seqcode.gse.utils.*;

/* see ExpanderModel for an overview */

public class MapperModel<IN,OUT> extends SeqViewModel implements Runnable {

    private boolean newinput;
    private Mapper<IN,OUT> mapper;
    private IN input;
    private OUT result;

    public MapperModel(Mapper<IN,OUT> mapper) {
        super();
        this.mapper = mapper;
        newinput = false;
    }

    public synchronized void run() {
        while(keepRunning()) {
            try {
                if (!newinput) {
                    wait();
                }
            } catch (InterruptedException ex) {

            }
            if (newinput) {
                try {
                    result = mapper.execute(input);
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
                newinput = false;
                notifyListeners();
            }
        }
    }
    
    public synchronized void setInput(IN i) {
        if (newinput == false) {
            if (!i.equals(input)) {
                input = i;
                newinput = true;
            } else {
                notifyListeners();
            }
        }
    }

    public OUT getResults() {
        if (newinput == false) {
            return result;
        } else {
            // or should we throw some exception?
            return null;
        }
    }

    public boolean isReady() {return !newinput;}

    

    

}
