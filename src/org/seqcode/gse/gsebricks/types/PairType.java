/*
 * Created on Apr 17, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package org.seqcode.gse.gsebricks.types;


public class PairType extends StructuredType {
    
    private EchoType firstType, secondType;

    public PairType(EchoType f, EchoType s) { 
        super("Pair", "first", f, "second", s);
        firstType = f; secondType = s;
    }
    
    public EchoType getFirstType() { return firstType; }
    public EchoType getSecondType() { return secondType; }
}
