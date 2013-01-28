package edu.psu.compbio.seqcode.gse.ewok.types;

import javax.swing.*;

import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.types.*;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class ValueWrapper implements SelfDescribingConstant { 
	private EchoType vclass;
	private Object value;

	public ValueWrapper(Object v) { 
		value = v;
		vclass = new ClassType(v.getClass());
	}

	public Object getConstantValue() { return value; }
	public EchoType getConstantClass() { return vclass; }
	public void setConstantValue(Object v) { value = v; }
    
    public String toString() { return value.toString(); }

}
