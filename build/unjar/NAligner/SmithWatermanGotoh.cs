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

#if TESTS
using System.IO;
using NUnit.Framework;
#endif

namespace NAligner
{
	/// <summary> An implementation of the Smith-Waterman algorithm with 
	/// Gotoh's improvement for biological local pairwise sequence alignment. </summary>
	public class SmithWatermanGotoh
	{
		/// <summary> Hidden constructor</summary>
		private SmithWatermanGotoh() : base()
		{
		}

		/// <summary> Aligns two sequences by Smith-Waterman algorithm</summary>
		/// <param name="s1">sequene #1 </param>
		/// <param name="s2">sequene #2 </param>
		/// <param name="matrix">scoring matrix </param>
		/// <param name="o">open gap penalty </param>
		/// <param name="e">extend gap penalty </param>
		/// <returns> alignment object contains the two aligned sequences, 
		/// the alignment score and alignment statistics</returns>
		/// <seealso cref="Sequence"/>
		/// <seealso cref="Matrix"/>
		public static Alignment Align(Sequence s1, Sequence s2, Matrix matrix, float o, float e)
		{		
			float[,] scores = matrix.Scores;

			SmithWatermanGotoh sw = new SmithWatermanGotoh();

			int m = s1.Length + 1;
			int n = s2.Length + 1;

			byte[] pointers = new byte[m * n];

			// Initializes the boundaries of the traceback matrix to STOP.
			for (int i = 0, k = 0; i < m; i++, k += n) 
			{
				pointers[k] = Directions.STOP;
			}
			for (int j = 1; j < n; j++) 
			{
				pointers[j] = Directions.STOP;
			}

			short[] sizesOfVerticalGaps = new short[m * n];
			short[] sizesOfHorizontalGaps = new short[m * n];
			for (int i = 0, k = 0; i < m; i++, k += n) 
			{
				for (int j = 0; j < n; j++) 
				{
					sizesOfVerticalGaps[k + j] = sizesOfHorizontalGaps[k + j] = 1;
				}
			}

			Cell cell = sw.Construct(s1, s2, scores, o, e, pointers,
				sizesOfVerticalGaps, sizesOfHorizontalGaps);

			Alignment alignment = sw.Traceback(s1, s2, matrix, pointers, cell,
				sizesOfVerticalGaps, sizesOfHorizontalGaps);
			
			alignment.Name1 = s1.Id;
			alignment.Name2 = s2.Id;
			alignment.Matrix = matrix;
			alignment.Open = o;
			alignment.Extend = e;

			return alignment;

		}

		/// <summary> Constructs directions matrix for the traceback </summary>
		/// <param name="s1">sequence #1 </param>
		/// <param name="s2">sequence #2 </param>
		/// <param name="matrix">scoring matrix </param>
		/// <param name="o">open gap penalty </param>
		/// <param name="e">extend gap penalty </param>
		/// <returns> The cell where the traceback starts. </returns>
		private Cell Construct(Sequence s1, Sequence s2, float[,] matrix, float o,
			float e, byte[] pointers, short[] sizesOfVerticalGaps, short[] sizesOfHorizontalGaps)
		{
		
			char[] a1 = s1.ToArray();
			char[] a2 = s2.ToArray();

			int m = s1.Length + 1;
			int n = s2.Length + 1;

			float f; // score of alignment x1...xi to y1...yi if xi aligns to yi
			float[] g = new float[n]; // score if xi aligns to a gap after yi
			float h; // score if yi aligns to a gap after xi
			float[] v = new float[n]; // best score of alignment x1...xi to y1...yi
			float vDiagonal;

			g[0] = float.NegativeInfinity;
			h = float.NegativeInfinity;
			v[0] = 0;

			for (int j = 1; j < n; j++) 
			{
				g[j] = float.NegativeInfinity;
				v[j] = 0;
			}

			float similarityScore, g1, g2, h1, h2;

			Cell cell = new Cell();

			for (int i = 1, k = n; i < m; i++, k += n) 
			{
				h = float.NegativeInfinity;
				vDiagonal = v[0];
				for (int j = 1, l = k + 1; j < n; j++, l++) 
				{
					similarityScore = matrix[a1[i - 1], a2[j - 1]];

					// Fill the matrices
					f = vDiagonal + similarityScore;

					g1 = g[j] - e;
					g2 = v[j] - o;
					if (g1 > g2) 
					{
						g[j] = g1;
						sizesOfVerticalGaps[l] = (short) (sizesOfVerticalGaps[l - n] + 1);
					} 
					else 
					{
						g[j] = g2;
					}

					h1 = h - e;
					h2 = v[j - 1] - o;
					if (h1 > h2) 
					{
						h = h1;
						sizesOfHorizontalGaps[l] = (short) (sizesOfHorizontalGaps[l - 1] + 1);
					} 
					else 
					{
						h = h2;
					}

					vDiagonal = v[j];
					v[j] = Max(f, g[j], h, 0);

					// Determine the traceback direction
					if (v[j] == 0) 
					{
						pointers[l] = Directions.STOP;
					} 
					else if (v[j] == f) 
					{
						pointers[l] = Directions.DIAGONAL;
					} 
					else if (v[j] == g[j]) 
					{
						pointers[l] = Directions.UP;
					} 
					else 
					{
						pointers[l] = Directions.LEFT;
					}

					// Set the traceback start at the current cell i, j and score
					if (v[j] > cell.Score) 
					{
						cell.Set(i, j, v[j]);
					}
				}
			}

			return cell;

		}

		/// <summary> Returns the alignment of two sequences based on the passed array of pointers</summary>
		/// <param name="s1">sequence #1 </param>
		/// <param name="s2">sequence #2 </param>
		/// <param name="m">scoring matrix </param>
		/// <param name="cell">The cell where the traceback starts. </param>
		/// <returns> <see cref="Alignment"/> with the two aligned sequences and alignment score. </returns>
		/// <seealso cref="Cell"/>
		/// <seealso cref="Alignment"/>
		private Alignment Traceback(Sequence s1, Sequence s2, Matrix m,
			byte[] pointers, Cell cell, short[] sizesOfVerticalGaps, short[] sizesOfHorizontalGaps)
		{
		
			char[] a1 = s1.ToArray();
			char[] a2 = s2.ToArray();
		
			float[,] scores = m.Scores;

			int n = s2.Length + 1;

			Alignment alignment = new Alignment();
			alignment.Score = cell.Score;

			int maxlen = s1.Length + s2.Length; // maximum length after the
			// aligned sequences

			char[] reversed1 = new char[maxlen]; // reversed sequence #1
			char[] reversed2 = new char[maxlen]; // reversed sequence #2
			char[] reversed3 = new char[maxlen]; // reversed markup

			int len1 = 0; // length of sequence #1 after alignment
			int len2 = 0; // length of sequence #2 after alignment
			int len3 = 0; // length of the markup line

			int identity = 0; // count of identitcal pairs
			int similarity = 0; // count of similar pairs
			int gaps = 0; // count of gaps

			char c1, c2;

			int i = cell.Row; // traceback start row
			int j = cell.Column; // traceback start col
			int k = i * n;

			bool stillGoing = true; // traceback flag: true -> continue & false
			// -> stop

			while (stillGoing) 
			{
				switch (pointers[k + j]) 
				{
					case Directions.UP:

						for (int l = 0, len = sizesOfVerticalGaps[k + j]; l < len; l++) 
						{
							reversed1[len1++] = a1[--i];
							reversed2[len2++] = Alignment.GAP;
							reversed3[len3++] = Markups.GAP;
							k -= n;
							gaps++;
						}
						break;

					case Directions.DIAGONAL:
						c1 = a1[--i];
						c2 = a2[--j];
						k -= n;
						reversed1[len1++] = c1;
						reversed2[len2++] = c2;
						if (c1 == c2) 
						{
							reversed3[len3++] = Markups.IDENTITY;
							identity++;
							similarity++;
						} 
						else if (scores[c1,c2] > 0) 
						{
							reversed3[len3++] = Markups.SIMILARITY;
							similarity++;
						} 
						else 
						{
							reversed3[len3++] = Markups.MISMATCH;
						}
						break;

					case Directions.LEFT:
						for (int l = 0, len = sizesOfHorizontalGaps[k + j]; l < len; l++) 
						{
							reversed1[len1++] = Alignment.GAP;
							reversed2[len2++] = a2[--j];
							reversed3[len3++] = Markups.GAP;
							gaps++;
						}
						break;

					case Directions.STOP:
						stillGoing = false;
						break;
				}
			}


			alignment.Sequence1 = Reverse(reversed1, len1);
			alignment.Start1 = i;
			alignment.Sequence2 = Reverse(reversed2, len2);
			alignment.Start2 = j;
			alignment.MarkupLine = Reverse(reversed3, len3);
			alignment.Identity = identity;
			alignment.Gaps = gaps;
			alignment.Similarity = similarity;

			return alignment;

		}

		/// <summary> Returns the maximum of 4 float numbers.</summary>
		private static float Max(float a, float b, float c, float d)
		{
			return Math.Max(Math.Max(a, b), Math.Max(c, d));
		}

		/// <summary> Reverses an array of chars</summary>
		private static char[] Reverse(char[] a, int len)
		{
			// TODO: replace this method by System.Array.Reverse
			char[] b = new char[len];
			for (int i = len - 1, j = 0; i >= 0; i--, j++)
			{
				b[j] = a[i];
			}
			return b;
		}

#if TESTS
		/// <summary>Testing the class <see cref="SmithWatermanGotoh"/>.</summary>
		[TestFixture]
			public class Tests
		{
			string NL = System.Environment.NewLine;

			[Test] public void Align()
			{
				string matrixPath = @"../../matrices/PAM250";
				if(!File.Exists(matrixPath)) matrixPath = @"PAM250";

				Matrix pam250 = Matrix.Load(matrixPath);

				float opengap = 15.0f;
				float extendgap = 3.0f;

				string str1 = @">100K_RAT  100 kDa protein (EC 6.3.2.-)." + NL +
					"MMSARGDFLN YALSLMRSHN DEHSDVLPVL DVCSLKHVAY VFQALIYWIK AMNQQTTLDT" +
					"PQLERKRTRE LLELGIDNED SEHENDDDTS QSATLNDKDD ESLPAETGQN HPFFRRSDSM" +
					"TFLGCIPPNP FEVPLAEAIP LADQPHLLQP NARKEDLFGR PSQGLYSSSA GSGKCLVEVT" +
					"MDRNCLEVLP TKMSYAANLK NVMNMQNRQK KAGEDQSMLA EEADSSKPGP SAHDVAAQLK" +
					"SSLLAEIGLT ESEGPPLTSF RPQCSFMGMV ISHDMLLGRW RLSLELFGRV FMEDVGAEPG" +
					"SILTELGGFE VKESKFRREM EKLRNQQSRD LSLEVDRDRD LLIQQTMRQL NNHFGRRCAT" +
					"TPMAVHRVKV TFKDEPGEGS GVARSFYTAI AQAFLSNEKL PNLDCIQNAN KGTHTSLMQR" +
					"LRNRGERDRE REREREMRRS SGLRAGSRRD RDRDFRRQLS IDTRPFRPAS EGNPSDDPDP" +
					"LPAHRQALGE RLYPRVQAMQ PAFASKITGM LLELSPAQLL LLLASEDSLR ARVEEAMELI" +
					"VAHGRENGAD SILDLGLLDS SEKVQENRKR HGSSRSVVDM DLDDTDDGDD NAPLFYQPGK" +
					"RGFYTPRPGK NTEARLNCFR NIGRILGLCL LQNELCPITL NRHVIKVLLG RKVNWHDFAF" +
					"FDPVMYESLR QLILASQSSD ADAVFSAMDL AFAVDLCKEE GGGQVELIPN GVNIPVTPQN" +
					"VYEYVRKYAE HRMLVVAEQP LHAMRKGLLD VLPKNSLEDL TAEDFRLLVN GCGEVNVQML" +
					"ISFTSFNDES GENAEKLLQF KRWFWSIVER MSMTERQDLV YFWTSSPSLP ASEEGFQPMP" +
					"SITIRPPDDQ HLPTANTCIS RLYVPLYSSK QILKQKLLLA IKTKNFGFV";

				string str2 = @">104K_THEPA  104 kDa microneme-rhoptry antigen." + NL +
					"MKFLILLFNI LCLFPVLAAD NHGVGPQGAS GVDPITFDIN SNQTGPAFLT AVEMAGVKYL" +
					"QVQHGSNVNI HRLVEGNVVI WENASTPLYT GAIVTNNDGP YMAYVEVLGD PNLQFFIKSG" +
					"DAWVTLSEHE YLAKLQEIRQ AVHIESVFSL NMAFQLENNK YEVETHAKNG ANMVTFIPRN" +
					"GHICKMVYHK NVRIYKATGN DTVTSVVGFF RGLRLLLINV FSIDDNGMMS NRYFQHVDDK" +
					"YVPISQKNYE TGIVKLKDYK HAYHPVDLDI KDIDYTMFHL ADATYHEPCF KIIPNTGFCI" +
					"TKLFDGDQVL YESFNPLIHC INEVHIYDRN NGSIICLHLN YSPPSYKAYL VLKDTGWEAT" +
					"THPLLEEKIE ELQDQRACEL DVNFISDKDL YVAALTNADL NYTMVTPRPH RDVIRVSDGS" +
					"EVLWYYEGLD NFLVCAWIYV SDGVASLVHL RIKDRIPANN DIYVLKGDLY WTRITKIQFT" +
					"QEIKRLVKKS KKKLAPITEE DSDKHDEPPE GPGASGLPPK APGDKEGSEG HKGPSKGSDS" +
					"SKEGKKPGSG KKPGPAREHK PSKIPTLSKK PSGPKDPKHP RDPKEPRKSK SPRTASPTRR" +
					"PSPKLPQLSK LPKSTSPRSP PPPTRPSSPE RPEGTKIIKT SKPPSPKPPF DPSFKEKFYD" +
					"DYSKAASRSK ETKTTVVLDE SFESILKETL PETPGTPFTT PRPVPPKRPR TPESPFEPPK" +
					"DPDSPSTSPS EFFTPPESKR TRFHETPADT PLPDVTAELF KEPDVTAETK SPDEAMKRPR" +
					"SPSEYEDTSP GDYPSLPMKR HRLERLRLTT TEMETDPGRM AKDASGKPVK LKRSKSFDDL" +
					"TTVELAPEPK ASRIVVDDEG TEADDEETHP PEERQKTEVR RRRPPKKPSK SPRPSKPKKP" +
					"KKPDSAYIPS ILAILVVSLI VGIL";


				Sequence seq1 = Sequence.Parse(str1);
				Sequence seq2 = Sequence.Parse(str2);

				Alignment alignment = SmithWatermanGotoh.Align(seq1, seq2, pam250, opengap, extendgap);

				Assert.AreEqual(98.0f, alignment.Score, "Unexpected score.");
			}

			[Test] public void Align2()
			{
				string matrixPath = @"../../matrices/TEST1";
				if(!File.Exists(matrixPath)) matrixPath = @"TEST1";

				Matrix customMatrix = Matrix.Load(matrixPath);

				Sequence seq1 = Sequence.Parse(@"AACCC"); 
				Sequence seq2 = Sequence.Parse(@"CCACC");

				float opengap = 2.0f;
				float extendgap = 1.0f;

				Alignment alignment = SmithWatermanGotoh.Align(seq1, seq2, customMatrix, opengap, extendgap);

				Assert.AreEqual(14.0f, alignment.Score, "Unexpected score.");
			}
		}
#endif
	}
}