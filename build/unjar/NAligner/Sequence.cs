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
using System.IO;
using System.Text;
using NAligner.util;

namespace NAligner
{
	public enum SequenceType
	{
		Protein,

		Nucleic
	}

	/// <summary> A basic (nucleic or protein) sequence. It's a 
	/// wrapper to <c>string</c>. </summary>
	/// <remarks>The <see cref="Sequence"/> can be parsed if encoded with
	/// <a href="http://www.ncbi.nlm.nih.gov/BLAST/fasta.html">FASTA</a>.</remarks>
	[Serializable]
	public class Sequence
	{
		/// <summary> Returns the sequence id</summary>
		public virtual string Id
		{
			get { return id; }
			set { this.id = value; }

		}

		/// <summary> Returns the sequence description</summary>
		public virtual string Description
		{
			get { return description; }
			set { this.description = value; }

		}

		/// <summary> Returns the sequence type (nucleic or protein) </summary>
		public virtual SequenceType Type
		{
			get { return type; }
			set { this.type = value; }

		}

//		/// <summary> Sequence type nucleic.</summary>
//		public const int NUCLEIC = 0;
//
//		/// <summary> Sequence type protein.</summary>
//		public const int PROTEIN = 1;

		/// <summary> List of amino acids</summary>
		private string aalist;

		/// <summary> Sequence id.</summary>
		private string id = null;

		/// <summary> Sequence description.</summary>
		private string description = null;

		/// <summary> Sequence type.</summary>
		private SequenceType type = SequenceType.Protein;

		/// <summary> Constructor</summary>
		public Sequence() : base()
		{
		}

		/// <summary> Constructor</summary>
		public Sequence(string sequence, string id, string description, SequenceType type) : base()
		{
			this.aalist = sequence;
			this.id = id;
			this.description = description;
			this.type = type;
		}

		/// <summary> Gets or sets the amino acid sequence</summary>
		public virtual string AAList
		{
			get { return aalist; }
			set { this.aalist = value; }
		}

		public virtual int Length
		{
			get { return this.aalist.Length; }
		}

		/// <summary> Returns a subsequence</summary>
		/// <param name="index">start index </param>
		/// <param name="length">length of subsequence </param>
		/// <returns> subsequence </returns>
		public virtual string Subsequence(int index, int length)
		{
			return this.aalist.Substring(index, (index + length) - (index));
		}

		/// <summary> Returns the acid at specific location in the sequence</summary>
		/// <param name="index">index </param>
		/// <returns> acid at index </returns>
		public virtual char AcidAt(int index)
		{
			return this.aalist[index];
		}

		/// <summary> Returns the sequence as an array of characters.</summary>
		public virtual char[] ToArray()
		{
			return this.aalist.ToCharArray();
		}

		#region I/O

		/// <summary> Returns a parsed Sequence from a FASTA string.</summary>
		/// <param name="stringToParse">FASTA string to parse</param>
		public static Sequence Parse(string stringToParse)
		{
			stringToParse = stringToParse.Replace("\r\n", "\n");

			string sequenceName = null;
			string sequenceDescription = null;

			if (stringToParse.StartsWith(">"))
			{
				// FASTA format
				int index = stringToParse.IndexOf("\n");

				if (index == - 1)
				{
					throw new System.Exception("Invalid sequence");
				}

				string first = stringToParse.Substring(1, (index) - (1));
				stringToParse = stringToParse.Substring(index);

				index = 0;
				for (int i = 0; i < first.Length && first[i] != ' ' && first[i] != '\t'; i++, index++)
				{
					// Skip white spaces
				}
				sequenceName = first.Substring(0, (index) - (0));
				sequenceDescription = index + 1 > first.Length ? "" : first.Substring(index + 1);
			}
			else
			{
				// Plain format ... nothing to do here
			}

			Sequence s = new Sequence(PrepareAndValidate(stringToParse), sequenceName, sequenceDescription, SequenceType.Protein);

			return s;
		}

		/// <summary> Returns a <see cref="Sequence"/> parsed and loaded from a file</summary>
		/// <param name="file">to parse </param>
		/// <returns> parsed sequence  </returns>
		public static Sequence Parse(FileInfo file)
		{
			string sequenceName = null;
			string sequenceDescription = null;

			StreamReader reader = new StreamReader(file.FullName);
			StringBuilder buffer = new StringBuilder();

			// Read & parse the first line
			string line = reader.ReadLine();

			if (line.StartsWith(">"))
			{
				// FASTA sequence

				line = line.Substring(1).Trim();
				int index = 0;
				for (int i = 0; i < line.Length && line[i] != ' ' && line[i] != '\t'; i++, index++)
				{
					// Skip white spaces
				}

				sequenceName = line.Substring(0, (index) - (0));
				sequenceDescription = index + 1 > line.Length ? "" : line.Substring(index + 1);
			}
			else
			{
				// Plain sequence
				buffer.Append(PrepareAndValidate(line));
			}

			// Read the remaining the file (the actual sequence)
			while ((line = reader.ReadLine()) != null)
			{
				buffer.Append(PrepareAndValidate(line));
			}
			reader.Close();

			Sequence s = new Sequence(buffer.ToString(), sequenceName, sequenceDescription, SequenceType.Protein);
			return s;
		}

		/// <summary> Removes whitespaces from a sequence and validates other characters.</summary>
		/// <param name="stringToParse">sequence to be prepared</param>
		/// <returns> prepared array of characters </returns>
		private static string PrepareAndValidate(string stringToParse)
		{
			StringBuilder buffer = new StringBuilder();
			string copy = stringToParse.Trim().ToUpper();

			for (int i = 0, n = copy.Length; i < n; i++)
			{
				switch (copy[i])
				{
						// skip whitespaces
					case (char) (9):
					case (char) (10):
					case (char) (13):
					case (char) (32):
						break;

						// add a valid character

					case 'A':
					case 'B':
					case 'C':
					case 'D':
					case 'E':
					case 'F':
					case 'G':
					case 'H':
					case 'I':
					case 'K':
					case 'L':
					case 'M':
					case 'N':
					case 'P':
					case 'Q':
					case 'R':
					case 'S':
					case 'T':
					case 'U':
					case 'V':
					case 'W':
					case 'Y':
					case 'Z':
					case 'X':
					case '-':
					case '*':
						buffer.Append(copy[i]);
						break;

						// throw an exception for anything else

					default:
						throw new NAlignerException("Invalid sequence character: " + copy[i]);

				}
			}
			return buffer.ToString();
		}

		#endregion
	}
}