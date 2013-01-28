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

namespace NAligner.util
{
	/// <summary> Comparator to sort the scoring matrices by their names. </summary>
	public class AlphanumericComparator : IComparer
	{
		public virtual int Compare(object o1, object o2)
		{
			string s1 = (string) o1;
			string s2 = (string) o2;

			int index1 = FirstDigitIndex(s1);
			int index2 = FirstDigitIndex(s2);

			if (index1 == - 1 || index2 == - 1)
			{
				return string.Compare(s1, s2, true);
			}
			else
			{
				string s3 = s1.Substring(0, (index1) - (0));
				string s4 = s2.Substring(0, (index2) - (0));
				if (s3.ToUpper().Equals(s4.ToUpper()))
				{
					return System.Int32.Parse(s1.Substring(index1)).CompareTo(
							System.Int32.Parse(s2.Substring(index2)));
				}
				else
				{
					return string.Compare(s1, s2, true);
				}
			}
		}

		/// <summary> Returns the index of the first digit in a String.
		/// If there are no digits, returns -1.</summary>
		private int FirstDigitIndex(string s)
		{
			for (int i = 0; i < s.Length; i++)
			{
				if (s[i] >= '0' && s[i] <= '9')
				{
					return i;
				}
			}
			return - 1;
		}
	}
}