package org.seqcode.data.seqdata.tools.editor;

import org.seqcode.data.seqdata.SeqExpt;
import org.seqcode.viz.components.ObjectTableModel;

public class SeqExptTableModel extends ObjectTableModel<SeqExpt> {

	private int findNewIndex(SeqExpt bs) {
		String n = bs.getName();
		for (int i = 0; i < getSize(); i++) {
			SeqExpt os = getObject(i);
			String on = os.getName();
			if (n.equals(on)) {
				return i;
			}
		}
		return getSize();
	}

	public int getColumnCount() {
		return 6;
	}

	public Class getColumnClass(int i) {
		if (i == 0) {
			return String.class;
		}
		if (i == 1) {
			return String.class;
		}
		if (i == 2) {
			return String.class;
		}
		if (i == 3) {
			return String.class;
		}
		if (i == 4) {
			return String.class;
		}
		if (i == 5) {
			return String.class;
		}
		return null;
	}

	public String getColumnName(int i) {
		if (i == 0) {
			return "ExptType";
		}
		if (i == 1) {
			return "Lab";
		}
		if (i == 2) {
			return "ExptCondition";
		}
		if (i == 3) {
			return "ExptTarget";
		}
		if (i == 4) {
			return "CellLine";
		}
		if (i == 5) {
			return "Replicate";
		}
		return null;
	}

	public Object getValueAt(int rowIndex, int c) {
		if (c == 0) {
			return getObject(rowIndex).getExptType().getName();
		}
		if (c == 1) {
			return getObject(rowIndex).getLab().getName();
		}
		if (c == 2) {
			return getObject(rowIndex).getExptCondition().getName();
		}
		if (c == 3) {
			return getObject(rowIndex).getExptTarget().getName();
		}
		if (c == 4) {
			return getObject(rowIndex).getCellLine().getName();
		}
		if (c == 5) {
			return getObject(rowIndex).getReplicate();
		}
		return null;
	}

}
