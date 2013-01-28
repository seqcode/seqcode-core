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

namespace NAligner.util
{
	/// <summary> Global constants and varilables  </summary>
	public abstract class Commons
	{
		/// <summary> Returns system file separator.</summary>
		/// <returns> file separator </returns>
		public static string FileSeparator
		{
			get { return fileSeparator; }

		}

		/// <summary> Returns system line separator.</summary>
		/// <returns> line separator </returns>
		public static string LineSeparator
		{
			get { return lineSeparator; }

		}

		/// <summary> Returns user's current directory.</summary>
		/// <returns> user's current directory </returns>
		public static string UserDirectory
		{
			get { return userDirectory; }

		}

		/// <summary> Returns the current release version of JAligner</summary>
		public static string CurrentRelease
		{
			get { return CURRENT_RELEASE; }

		}

		/// <summary> Current release version of JAligner</summary>
		private const string CURRENT_RELEASE = "0.9";

		/// <summary> Default home directory</summary>
		private const string DEFAULT_USER_DIRECTORY = ".";

		/// <summary> Default file separator</summary>
		private const string DEFAULT_FILE_SEPARATOR = "/";

		/// <summary> Default line separator</summary>
		private const string DEFAULT_LINE_SEPARATOR = "\r\n";

		/// <summary> User home directory</summary>
		private static string userDirectory = DEFAULT_USER_DIRECTORY;

		/// <summary> Line separator</summary>
		private static string fileSeparator = DEFAULT_FILE_SEPARATOR;

		/// <summary> Line separator</summary>
		private static string lineSeparator = DEFAULT_LINE_SEPARATOR;

		static Commons()
		{
			{
				try
				{
					userDirectory = System.Environment.CurrentDirectory;
				}
				catch (System.Exception e)
				{
					Console.Error.WriteLine("Failed getting user current directory: " + e.ToString());
				}
			}
			{
				try
				{
					fileSeparator = System.IO.Path.DirectorySeparatorChar.ToString();
				}
				catch (System.Exception e)
				{
					System.Console.Error.WriteLine("Failed getting system file separator: " + e.ToString());
				}
			}
			{
				try
				{
					lineSeparator = System.Environment.NewLine;
				}
				catch (System.Exception e)
				{
					System.Console.Error.WriteLine("Failed getting system line separator: " + e.ToString());
				}
			}
		}
	}
}