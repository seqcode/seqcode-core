package edu.psu.compbio.seqcode.gse.ewok.types;

public class CollectionType extends StructuredType {
	
	private EchoType innerType;
	
	public CollectionType(EchoType inner) { 
		super("Collection", "base", inner);
		innerType = inner;
	}
	
	public EchoType getInnerType() { return innerType; }
}
