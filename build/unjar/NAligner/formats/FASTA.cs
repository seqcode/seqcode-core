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

using System.Text;
using NAligner;

namespace NAligner.formats
{
	/// <summary> <a href="http://www.ncbi.nlm.nih.gov/BLAST/fasta.html">FASTA</a> formatter</summary>
	public class FASTA : AlignmentFormatter
	{
		/// <summary> Number of characters per line</summary>
		private const int LINE_WIDTH = 60;

		/// <summary> Constructor for FASTA.</summary>
		public FASTA() : base()
		{
			Id = "FASTA";
		}

		/// <summary> Returns the name, description and sequence combined in one string.
		/// The length of each line in the sequence is FASTA.LINE_LENGTH</summary>
		public virtual string Format(Sequence sequence)
		{
			StringBuilder buffer = new StringBuilder(">");
			buffer.Append(sequence.Id == null ? "" : sequence.Id);
			buffer.Append("\n");
			for (int i = 0, n = sequence.Length; i*LINE_WIDTH < n; i++)
			{
				for (int j = i*LINE_WIDTH, m = (i + 1)*LINE_WIDTH < n ? (i + 1)*LINE_WIDTH : n; j < m; j++)
				{
					buffer.Append(sequence.Subsequence(j, 1));
				}
				buffer.Append("\n");
			}
			return buffer.ToString();
		}

		/// <returns> FASTA format of the input alignment</returns>
		public override string Format(Alignment alignment)
		{
			System.Text.StringBuilder buffer = new System.Text.StringBuilder();
			System.Text.StringBuilder s1 = new System.Text.StringBuilder();
			System.Text.StringBuilder s2 = new System.Text.StringBuilder();
			s1.Append(alignment.Sequence1);
			s2.Append(alignment.Sequence2);
			buffer.Append(Format(new Sequence(s1.ToString(), alignment.Name1, "", SequenceType.Protein)));
			buffer.Append(Format(new Sequence(s2.ToString(), alignment.Name2, "", SequenceType.Protein)));
			return buffer.ToString();
		}
	}
}