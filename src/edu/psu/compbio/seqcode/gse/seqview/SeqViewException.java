package edu.psu.compbio.seqcode.gse.seqview;

/**
 * This class can be used to signal error cases in SeqView and 
 * as a wrapper for exceptions that may get generated within the program  
 * @author Bob
 *
 */
public class SeqViewException extends Exception {

	/**
	   * Constructs a <code>WarpDriveException</code> with no detail message.
	   */
	  public SeqViewException() {
	    super();
	  }

	  /**
	   * Constructs a <code>WarpDriveException</code> with specified message.
	   * @param   msg   the detail message.
	   */
	  public SeqViewException(String msg) {
	    super(msg);
	  }

	  /**
	   * Constructs a <code>WarpDriveException</code> with an Exception Object.
	   * @param   e   the Exception Object.
	   */
	  public SeqViewException(Exception e) {
	    super(e.getMessage());
	  }
}
