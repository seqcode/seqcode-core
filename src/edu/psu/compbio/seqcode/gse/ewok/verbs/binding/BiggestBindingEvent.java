package edu.psu.compbio.seqcode.gse.ewok.verbs.binding;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Distiller;

public class BiggestBindingEvent<X extends BindingEvent> implements Distiller<X,X> {

    private X biggest;
    public BiggestBindingEvent() {
        biggest = null;
    }
    public X execute(X input) {
        if (biggest == null || 
            input.getSize() > biggest.getSize()) {
            biggest = input;
        }
        return biggest;
    }
    public X getCurrent() {
        return biggest;
    }
    public void reset() {
        biggest = null;
    }
}
