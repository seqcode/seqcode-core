/*
 * Author: tdanford
 * Date: Jun 12, 2008
 */
package edu.psu.compbio.seqcode.gse.utils.io.parsing.sexprs;

import java.util.*;

public interface SExpr {
	public boolean isCompound();
	public int length();
	public SExpr subExpr(int i);
}