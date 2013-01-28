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

namespace NAligner
{
	/// <summary> A cell in a similarity matrix, to hold row, column and score.</summary>
	public class Cell
	{
		/// <summary> Gets or sets the column. </summary>
		public virtual int Column
		{
			get { return this.column; }
			set { this.column = value; }

		}

		/// <summary>Gets or sets the row.</summary>
		public virtual int Row
		{
			get { return this.row; }
			set { this.row = value; }

		}

		/// <summary>Gets or sets the score.</summary>
		public virtual float Score
		{
			get { return this.score; }
			set { this.score = value; }
		}

		/// <summary> Row of the cell</summary>
		private int row;

		/// <summary> Column of the cell</summary>
		private int column;

		/// <summary> Alignment score at this cell</summary>
		private float score;

		/// <summary> Constructor</summary>
		public Cell() : base()
		{
			this.row = 0;
			this.column = 0;
			this.score = System.Single.NegativeInfinity;
		}

		/// <summary> Sets the row, column and score of the cell.</summary>
		public virtual void Set(int row, int column, float score)
		{
			this.row = row;
			this.column = column;
			this.score = score;
		}
	}
}