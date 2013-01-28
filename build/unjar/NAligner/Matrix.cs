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
using NAligner.util;

namespace NAligner
{
	/// <summary> Scoring matrix. </summary>
	[Serializable]
	public class Matrix
	{
		/// <summary> The starter character of a comment line.</summary>
		private const char COMMENT_STARTER = '#';

		/// <summary> The size of the scoring matrix. It is the number of the 
		/// characters in the ASCII table. It is more than the 20 amino acids just 
		/// to save the processing time of the mapping.</summary>
		private const int SIZE = 127;

		/// <summary> Matrix id (or name)</summary>
		private string id = null;

		/// <summary> Scores</summary>
		private float[,] scores = null;

		public Matrix(string id, float[,] scores)
		{
			this.id = id;
			this.scores = scores;
		}

		/// <summary> Loads scoring matrix from a file.</summary>
		/// <returns> loaded matrix </returns>
		public static Matrix Load(string matrixPath)
		{
			char[] acids = new char[SIZE];

			// Initialize the acids array to null values (ascii = 0)
			for (int i = 0; i < SIZE; i++)
			{
				acids[i] = (char) (0);
			}

			float[,] scores = new float[SIZE, SIZE];

			StreamReader reader = new StreamReader(matrixPath);

			// Skip the comment lines
			string line;
			while ((line = reader.ReadLine()) != null && line.Trim()[0] == COMMENT_STARTER)
				;

			// Read the headers line (the letters of the acids)
			Tokenizer tokenizer = new Tokenizer(line.Trim());
			for (int j = 0; tokenizer.HasMoreTokens(); j++)
			{
				acids[j] = tokenizer.NextToken()[0];
			}

			// Read the scores
			while ((line = reader.ReadLine()) != null)
			{
				tokenizer = new Tokenizer(line.Trim());
				char acid = tokenizer.NextToken()[0];
				for (int i = 0; i < SIZE; i++)
				{
					if (acids[i] != 0)
					{
						scores[acid, acids[i]] = System.Single.Parse(tokenizer.NextToken());
					}
				}
			}
			return new Matrix(matrixPath, scores);
		}

		/// <returns> Returns the id. </returns>
		public virtual string Id
		{
			get { return this.id; }

		}

		/// <returns> Returns the scores. </returns>
		public virtual float[,] Scores
		{
			get { return this.scores; }

		}


		/// <returns> score </returns>
		public virtual float GetScore(char a, char b)
		{
			return this.scores[a,b];
		}
	}
}