package org.seqcode.gse.gsebricks.types;

public class SequenceType extends StructuredType {
	
	private EchoType innerType;
	
	public SequenceType(EchoType inner) { 
		super("Sequence", "base", inner);
		innerType = inner;
	}
	
	public EchoType getInnerType() { return innerType; }
}
