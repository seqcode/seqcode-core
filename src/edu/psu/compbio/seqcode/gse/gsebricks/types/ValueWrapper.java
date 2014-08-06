package edu.psu.compbio.seqcode.gse.gsebricks.types;


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
