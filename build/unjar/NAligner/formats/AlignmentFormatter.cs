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

using NAligner;

namespace NAligner.formats
{
	/// <summary> Abstract format </summary>
	public abstract class AlignmentFormatter
	{
		/// <summary> Gets or sets format id</summary>
		public virtual string Id
		{
			get { return this.id == null ? this.GetType().FullName : this.id; }
			set { this.id = value; }

		}

		/// <summary> Formatter id </summary>
		private string id = null;

		/// <summary> Formats alignment</summary>
		/// <returns> formatted alignment </returns>
		public abstract string Format(Alignment alignment);
	}
}