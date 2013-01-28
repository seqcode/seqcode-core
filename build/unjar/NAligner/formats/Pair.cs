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
using NAligner.util;

namespace NAligner.formats
{
	/// <summary> 
	/// <a href="http://www.hgmp.mrc.ac.uk/Software/EMBOSS/Themes/AlignExamples/align.pair">Pair</a> format.</summary>
	public class Pair : AlignmentFormatter
	{
		/// <summary> Name width </summary>
		private const int NAME_WIDTH = 13;

		/// <summary> Position width</summary>
		private const int POSITION_WIDTH = 6;

		/// <summary> Sequence width</summary>
		private const int SEQUENCE_WIDTH = 50;

		/// <summary> Space</summary>
		private const string BLANK = " ";

		/// <summary> Constructor </summary>
		public Pair() : base()
		{
			Id = "Pair";
		}

		/// <summary> Formats an alignment object to the Pair_FORMAT format</summary>
		/// <param name="alignment">alignment object to be formated </param>
		/// <returns> string of the alignment pair-formatted </returns>
		public override string Format(Alignment alignment)
		{
			char[] sequence1 = alignment.Sequence1;
			char[] sequence2 = alignment.Sequence2;
			char[] markup = alignment.MarkupLine;

			int length = sequence1.Length > sequence2.Length ? sequence2.Length : sequence1.Length;

			string name1 = AdjustName(alignment.Name1);
			string name2 = AdjustName(alignment.Name2);

			StringBuilder buffer = new StringBuilder();

			StringBuilder preMarkup = new StringBuilder();
			for (int j = 0; j < NAME_WIDTH + 1 + POSITION_WIDTH + 1; j++)
			{
				preMarkup.Append(BLANK);
			}

			int oldPosition1, position1 = 1 + alignment.Start1;
			int oldPosition2, position2 = 1 + alignment.Start2;

			char[] subsequence1;
			char[] subsequence2;
			char[] submarkup;
			int line;

			char c1, c2;

			for (int i = 0; i*SEQUENCE_WIDTH < length; i++)
			{
				oldPosition1 = position1;
				oldPosition2 = position2;

				line = ((i + 1)*SEQUENCE_WIDTH) < length ? (i + 1)*SEQUENCE_WIDTH : length;

				subsequence1 = new char[line - i*SEQUENCE_WIDTH];
				subsequence2 = new char[line - i*SEQUENCE_WIDTH];
				submarkup = new char[line - i*SEQUENCE_WIDTH];

				for (int j = i*SEQUENCE_WIDTH, k = 0; j < line; j++, k++)
				{
					subsequence1[k] = sequence1[j];
					subsequence2[k] = sequence2[j];
					submarkup[k] = markup[j];
					c1 = subsequence1[k];
					c2 = subsequence2[k];
					if (c1 == c2)
					{
						position1++;
						position2++;
					}
					else if (c1 == Alignment.GAP)
					{
						position2++;
					}
					else if (c2 == Alignment.GAP)
					{
						position1++;
					}
					else
					{
						position1++;
						position2++;
					}
				}

				buffer.Append(name1);
				buffer.Append(BLANK);
				buffer.Append(AdjustPosition((oldPosition1).ToString()));
				buffer.Append(BLANK);
				buffer.Append(subsequence1);
				buffer.Append(BLANK);
				buffer.Append(AdjustPosition(((position1 - 1)).ToString()));
				buffer.Append(Commons.LineSeparator);

				buffer.Append(preMarkup);
				buffer.Append(submarkup);
				buffer.Append(Commons.LineSeparator);

				buffer.Append(name2);
				buffer.Append(BLANK);
				buffer.Append(AdjustPosition((oldPosition2).ToString()));
				buffer.Append(BLANK);
				buffer.Append(subsequence2);
				buffer.Append(BLANK);
				buffer.Append(AdjustPosition(((position2 - 1)).ToString()));
				buffer.Append(Commons.LineSeparator);

				buffer.Append(Commons.LineSeparator);
			}
			return buffer.ToString();
		}

		private string AdjustName(string name)
		{
			System.Text.StringBuilder buffer = new System.Text.StringBuilder();

			if (name.Length > NAME_WIDTH)
			{
				buffer.Append(name.Substring(0, (NAME_WIDTH) - (0)));
			}
			else
			{
				buffer.Append(name);
				for (int j = buffer.Length; j < NAME_WIDTH; j++)
				{
					buffer.Append(BLANK);
				}
			}
			return buffer.ToString();
		}

		private string AdjustPosition(string position)
		{
			System.Text.StringBuilder buffer1 = new System.Text.StringBuilder();
			System.Text.StringBuilder buffer2 = new System.Text.StringBuilder();

			if (position.Length > POSITION_WIDTH)
			{
				buffer1.Append(position.Substring(position.Length - POSITION_WIDTH, (position.Length) - (position.Length - POSITION_WIDTH)));
			}
			else
			{
				buffer1.Append(position);
			}

			for (int j = 0, n = POSITION_WIDTH - buffer1.Length; j < n; j++)
			{
				buffer2.Append(BLANK);
			}

			buffer2.Append(buffer1.ToString());
			return buffer2.ToString();
		}
	}
}