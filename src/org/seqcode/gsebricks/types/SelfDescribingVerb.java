/*
 * Created on Feb 16, 2007
 */
package org.seqcode.gsebricks.types;



public interface SelfDescribingVerb extends SelfDescribingParameterized {    
    public String[] getInputNames();
    public EchoType[] getInputClasses();
    public EchoType getOutputClass();
}
