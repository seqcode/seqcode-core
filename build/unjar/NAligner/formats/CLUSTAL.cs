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
	/// <summary> <p>CLUSTAL format.</p>
	/// Example:
	/// <code>
	/// CLUSTAL_FORMAT W(1.60) multiple sequence alignment
	/// 
	/// 
	/// JC2395          NVSDVNLNK---YIWRTAEKMK---ICDAKKFARQHKIPESKIDEIEHNSPQDAAE----
	/// KPEL_DROME      MAIRLLPLPVRAQLCAHLDAL-----DVWQQLATAVKLYPDQVEQISSQKQRGRS-----
	/// FASA_MOUSE      NASNLSLSK---YIPRIAEDMT---IQEAKKFARENNIKEGKIDEIMHDSIQDTAE----
	/// 
	/// 
	/// JC2395          -------------------------QKIQLLQCWYQSHGKT--GACQALIQGLRKANRCD
	/// KPEL_DROME      -------------------------ASNEFLNIWGGQYN----HTVQTLFALFKKLKLHN
	/// FASA_MOUSE      -------------------------QKVQLLLCWYQSHGKS--DAYQDLIKGLKKAECRR
	/// 
	/// 
	/// JC2395          IAEEIQAM
	/// KPEL_DROME      AMRLIKDY
	/// FASA_MOUSE      TLDKFQDM
	/// </code>
	/// 
	/// </summary>
	/// <author>  Ahmed Moustafa (ahmed@users.sf.net) </author>
	public class CLUSTAL : AlignmentFormatter
	{
		/// <summary> Name width </summary>
		private const int NAME_WIDTH = 36;

		/// <summary> Sequence width</summary>
		private const int SEQUENCE_WIDTH = 50;

		/// <summary> CLUSTAL header</summary>
		private const string HEADER = "CLUSTAL_FORMAT W(1.60) multiple sequence alignment\n\n";

		/// <summary> Constructor</summary>
		public CLUSTAL() : base()
		{
			Id = "CLUSTAL";
		}

		/// <summary> Returns CLUSTAL format</summary>
		/// <param name="names">array of the names of the sequences.</param>
		/// <param name="sequences">array of the sequences</param>
		public virtual string Format(string[] names, string[] sequences)
		{
			StringBuilder buffer = new StringBuilder(HEADER);
			int maxSequenceLength = 0;
			for (int i = 0; i < sequences.Length; i++)
			{
				if (sequences[i].Length > maxSequenceLength)
				{
					maxSequenceLength = sequences[i].Length;
				}
			}

			for (int i = 0; i*SEQUENCE_WIDTH < maxSequenceLength; i++)
			{
				for (int j = 0; j < sequences.Length; j++)
				{
					buffer.Append(NAME_WIDTH <= names[j].Length ? names[j].Substring(0, (NAME_WIDTH - 1) - (0)) : names[j]);
					for (int k = names[j].Length; k < NAME_WIDTH; k++)
					{
						buffer.Append(" ");
					}
					if (names[j].Length >= NAME_WIDTH)
					{
						buffer.Append(" ");
					}
					buffer.Append(sequences[j].Substring(i*SEQUENCE_WIDTH, (((i + 1)*SEQUENCE_WIDTH) < sequences[j].Length ? (i + 1)*SEQUENCE_WIDTH : sequences[j].Length) - (i*SEQUENCE_WIDTH)));
					if (j < sequences.Length)
					{
						buffer.Append("\n");
					}
				}
				if ((i + 1)*SEQUENCE_WIDTH < maxSequenceLength)
				{
					buffer.Append("\n\n");
				}
			}
			return buffer.ToString();
		}

		/// <summary> Returns CLUSTAL format of the alignment</summary>
		public override string Format(Alignment alignment)
		{
			string[] sequences = new string[] {new string(alignment.Sequence1), new string(alignment.Sequence2)};
			string[] names = new string[] {alignment.Name1, alignment.Name2};
			return Format(names, sequences);
		}
	}
}