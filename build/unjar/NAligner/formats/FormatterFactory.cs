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

using System.Collections;

namespace NAligner.formats
{
	/// <summary> Singleton formatter factory. </summary>
	public class FormatterFactory
	{
		/// <summary> Returns the singleton instance for <see cref="FormatterFactory"/>.</summary>
		public static FormatterFactory Default
		{
			get
			{
				if (singletonInstance == null)
				{
					singletonInstance = new FormatterFactory();
				}
				return singletonInstance;
			}

		}

		/// <summary> Returns a list of registered formatters</summary>
		public virtual ICollection Formatters
		{
			get { return formats.Keys; }
		}

		private static FormatterFactory singletonInstance = null;

		/// <summary>Table containing the <see cref="AlignmentFormatter"/>s.</summary>
		private Hashtable formats = new Hashtable();

		/// <summary> Hidden constructor </summary>
		private FormatterFactory() : base()
		{
		}

		/// <summary> Registers a new formatter.</summary>
		public virtual void RegisterFormatter(AlignmentFormatter alignmentFormatter)
		{
			formats[alignmentFormatter.Id] = alignmentFormatter;
		}

		/// <summary> Returns a <see cref="AlignmentFormatter"/>.</summary>
		/// <param name="id">formatter id </param>
		public virtual AlignmentFormatter GetFormatter(string id)
		{
			return (AlignmentFormatter) formats[id];
		}
	}
}