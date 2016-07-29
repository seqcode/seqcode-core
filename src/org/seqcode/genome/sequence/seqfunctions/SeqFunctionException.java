package org.seqcode.genome.sequence.seqfunctions;

public class SeqFunctionException extends Exception{
	public SeqFunctionException() { super(); }
	public SeqFunctionException(String message) { super(message); }
	public SeqFunctionException(String message, Throwable cause) { super(message, cause); }
	public SeqFunctionException(Throwable cause) { super(cause); }
}
