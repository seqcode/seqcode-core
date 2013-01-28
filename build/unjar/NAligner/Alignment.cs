#region NAligner Copyright
/*
 * NAligner
 * C# port of JAligner API, http://jaligner.sourceforge.net
 * 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
#endregion

using System;

namespace NAligner
{
	/// <summary> Holds the output of a pairwise sequences alignment. </summary>
	public class Alignment
	{
		/// <returns> Returns the extend. </returns>
		public virtual float Extend
		{
			get { return extend; }
			set { this.extend = value; }

		}

		/// <returns> Returns the matrix. </returns>
		public virtual Matrix Matrix
		{
			get { return matrix; }
			set { this.matrix = value; }

		}

		/// <returns> Returns the name1. </returns>
		public virtual System.String Name1
		{
			get { return name1 == null || name1.Length == 0 ? DEFAULT_SEQUENCE1_NAME : name1; }
			set { this.name1 = value; }

		}

		/// <returns> Returns the name2. </returns>
		public virtual System.String Name2
		{
			get { return name2 == null || name2.Length == 0 ? DEFAULT_SEQUENCE2_NAME : name2; }

			set { this.name2 = value; }

		}

		/// <returns> Returns the open. </returns>
		public virtual float Open
		{
			get { return open; }
			set { this.open = value; }

		}

		/// <returns> Returns the score. </returns>
		public virtual float Score
		{
			get { return score; }
			set { this.score = value; }

		}

		/// <returns> Returns the sequence1. </returns>
		public virtual char[] Sequence1
		{
			get { return sequence1; }
			set { this.sequence1 = value; }

		}

		/// <returns> Returns the sequence2. </returns>
		public virtual char[] Sequence2
		{
			get { return sequence2; }
			set { this.sequence2 = value; }

		}

		/// <returns> Returns the start1. </returns>
		public virtual int Start1
		{
			get { return start1; }
			set { this.start1 = value; }

		}

		/// <returns> Returns the start2. </returns>
		public virtual int Start2
		{
			get { return start2; }
			set { this.start2 = value; }

		}

		/// <returns> Returns the gaps. </returns>
		public virtual int Gaps
		{
			get { return gaps; }
			set { this.gaps = value; }

		}

		/// <returns> Returns the identity. </returns>
		public virtual int Identity
		{
			get { return identity; }
			set { this.identity = value; }

		}

		/// <returns> Returns the markupLine. </returns>
		public virtual char[] MarkupLine
		{
			get { return markupLine; }
			set { this.markupLine = value; }

		}

		/// <returns> Returns the similarity. </returns>
		public virtual int Similarity
		{
			get { return similarity; }
			set { this.similarity = value; }

		}

		/// <summary> Returns a summary for alignment</summary>
		public virtual string Summary
		{
			get
			{
				throw new NotImplementedException("Unconverted part of the package.");

//				System.Text.StringBuilder buffer = new System.Text.StringBuilder();
//				//UPGRADE_ISSUE: Class 'java.text.DecimalFormat' was not converted. 'ms-help://MS.VSCC.2003/commoner/redir/redirect.htm?keyword="jlca1000_javatextDecimalFormat_3"'
//				//UPGRADE_ISSUE: Constructor 'java.text.DecimalFormat.DecimalFormat' was not converted. 'ms-help://MS.VSCC.2003/commoner/redir/redirect.htm?keyword="jlca1000_javatextDecimalFormat_3"'
//				DecimalFormat f1 = new DecimalFormat("0.00");
//				//UPGRADE_ISSUE: Class 'java.text.DecimalFormat' was not converted. 'ms-help://MS.VSCC.2003/commoner/redir/redirect.htm?keyword="jlca1000_javatextDecimalFormat_3"'
//				//UPGRADE_ISSUE: Constructor 'java.text.DecimalFormat.DecimalFormat' was not converted. 'ms-help://MS.VSCC.2003/commoner/redir/redirect.htm?keyword="jlca1000_javatextDecimalFormat_3"'
//				DecimalFormat f2 = new DecimalFormat("0.00%");
//
//				int length = Sequence1.Length;
//
//				buffer.Append("Sequence #1: " + Name1);
//				buffer.Append(Commons.LineSeparator);
//				buffer.Append("Sequence #2: " + Name2);
//				buffer.Append(Commons.LineSeparator);
//				buffer.Append("Matrix: " + (matrix.Id == null ? "" : matrix.Id));
//				buffer.Append(Commons.LineSeparator);
//				buffer.Append("Gap open: " + open);
//				buffer.Append(Commons.LineSeparator);
//				buffer.Append("Gap extend: " + extend);
//				buffer.Append(Commons.LineSeparator);
//				buffer.Append("Length: " + length);
//				buffer.Append(Commons.LineSeparator);
//				//UPGRADE_WARNING: Data types in Visual C# might be different.  Verify the accuracy of narrowing conversions. 'ms-help://MS.VSCC.2003/commoner/redir/redirect.htm?keyword="jlca1042_3"'
//				buffer.Append("Identity: " + identity + "/" + length + " (" + f2.FormatDouble(identity/(float) length) + ")");
//				buffer.Append(Commons.LineSeparator);
//				//UPGRADE_WARNING: Data types in Visual C# might be different.  Verify the accuracy of narrowing conversions. 'ms-help://MS.VSCC.2003/commoner/redir/redirect.htm?keyword="jlca1042_3"'
//				buffer.Append("Similarity: " + similarity + "/" + length + " (" + f2.FormatDouble(similarity/(float) length) + ")");
//				buffer.Append(Commons.LineSeparator);
//				//UPGRADE_WARNING: Data types in Visual C# might be different.  Verify the accuracy of narrowing conversions. 'ms-help://MS.VSCC.2003/commoner/redir/redirect.htm?keyword="jlca1042_3"'
//				buffer.Append("Gaps: " + gaps + "/" + length + " (" + f2.FormatDouble(gaps/(float) length) + ")");
//				buffer.Append(Commons.LineSeparator);
//				buffer.Append("Score: " + f1.FormatDouble(score));
//				buffer.Append(Commons.LineSeparator);
//
//				return buffer.ToString();
			}

		}

		/// <summary> Gap character</summary>
		public const char GAP = '-';

		/// <summary> Default name for sequence #1</summary>
		private const System.String DEFAULT_SEQUENCE1_NAME = "jaligner_1";
		/// <summary> Default name for sequence #2</summary>
		private const System.String DEFAULT_SEQUENCE2_NAME = "jaligner_2";

		/// <summary> Scoring matrix</summary>
		private Matrix matrix;

		/// <summary> Gap open cost</summary>
		private float open;

		/// <summary> Gap extend cost</summary>
		private float extend;

		/// <summary> Alignment score</summary>
		private float score;

		/// <summary> Aligned sequence #1</summary>
		private char[] sequence1;

		/// <summary> Name of sequence #1</summary>
		private System.String name1;

		/// <summary> Alignment start location in sequence #1</summary>
		private int start1;

		/// <summary> Aligned sequence #2</summary>
		private char[] sequence2;

		/// <summary> Name of sequence #2</summary>
		private System.String name2;

		/// <summary> Alignment start location in sequence #2</summary>
		private int start2;

		/// <summary> Markup line</summary>
		private char[] markupLine;

		/// <summary> Count of identical locations</summary>
		private int identity;

		/// <summary> Count of similar locations</summary>
		private int similarity;

		/// <summary> Count of gap locations</summary>
		private int gaps;


		/// <summary> Constructor for Alignment</summary>
		public Alignment() : base()
		{
		}
	}
}