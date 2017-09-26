package org.genvisis.gwas;

import java.io.File;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Set;

import org.genvisis.bioinformatics.Sequence;
import org.genvisis.common.Aliases;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.StrandOps;
import org.genvisis.seq.manage.StrandOps.CONFIG;

public class AlleleVerification {
	public static void verifyAlleles(String filename, String refFile, String freqFile, String posFile,
																	 boolean reorder,
																	 Logger log) {
		try {
			String[] header = Files.getHeaderOfFile(filename, log);
			String[] dataHeader = header;
			int[] indices = ext.indexFactors(new String[][] {Aliases.MARKER_NAMES, Aliases.ALLELES[0],
																											 Aliases.ALLELES[1], Aliases.CHRS},
																			 header, false, true, false);
			String[][] data;
			if (indices[3] == -1) {
				log.report("Couldn't find chromosome in " + filename + ".");
				data = HashVec.loadFileToStringMatrix(filename, false,
																							null);
				data = addPosition(data, posFile, log);
				indices[3] = data[0].length - 1;
				dataHeader = ArrayUtils.addStrToArray("Chr", dataHeader);
				data[0] = dataHeader;
			} else {
				data = HashVec.loadFileToStringMatrix(filename, false, null);
			}

			header = Files.getHeaderOfFile(freqFile, log);
			int[] cols = ext.indexFactors(new String[][] {Aliases.MARKER_NAMES, Aliases.ALLELE_FREQS},
																		header, false, true, false);
			Hashtable<String, String> frequencies = HashVec.loadFileToHashString(freqFile, cols[0],
																																					 new int[] {cols[1]},
																																					 null, true);

			boolean[] rowsToKeep = new boolean[data.length];
			rowsToKeep[0] = true;
			int opflip = 0;
			int opp = 0;
			int flip = 0;
			int same = 0;
			int err = 0;
			int ambg = 0;
			Hashtable<String, String> ref = null;
			String currentChr = null;
			for (int i = 1; i < data.length; i++) {
				String[] line = data[i];
				String[] alleles = new String[] {line[indices[1]].toUpperCase(),
																				 line[indices[2]].toUpperCase()};
				String chr = line[indices[3]] == null ? "" : line[indices[3]];

				String r = refFile;

				// only read in a new reference file if we have to
				if ((refFile.contains("#") && !chr.equals(currentChr)) || ref == (null)) {
					currentChr = chr;
					r = ext.replaceAllWith(r, "#", chr);
					if (!new File(r).exists())
						continue;
					header = Files.getHeaderOfFile(r, new Logger());
					cols = ext.indexFactors(new String[][] {Aliases.MARKER_NAMES, Aliases.ALLELES[0],
																									Aliases.ALLELES[1]},
																	header, false, true, false);

					ref = HashVec.loadFileToHashString(r, cols[0], new int[] {cols[1], cols[2]}, "\t", true);
				}

				// 1000G uses the format chr:position:a1:a2 OR rsid:position:a1:a2, so we need to convert to
				// chr:position or rsid
				Set<String> keys = ref.keySet();
				Hashtable<String, String> tempref = new Hashtable<String, String>();
				for (String k : keys) {
					String newKey = k;
					if (k.contains(":")) {

						String[] sp = k.split(":");
						if (k.startsWith("rs")) {
							newKey = sp[0];
						} else if (k.startsWith(chr + ":")) {
							newKey = sp[0] + ":" + sp[1];
						} else {
							newKey = chr + ":" + sp[1];
						}
					}
					String value = ref.get(k);
					tempref.put(newKey, value);
				}
				ref = tempref;


				String s = ref.get(line[indices[0]]);

				String[] ref_alleles;
				double freq;
				if (s != null) {
					String[] temp = s.split("\t");
					ref_alleles = new String[] {temp[0], temp[1]};
					freq = frequencies.get(line[indices[0]]) != null ? Double.parseDouble(frequencies.get(line[indices[0]]))
																													 : 1.0;
				} else {
					// keep this marker as is bc we can't find a ref for it
					rowsToKeep[i] = true;
					continue;
				}

				// write a method to consider freq that returns the same thing basically
				CONFIG config = StrandOps.determineStrandConfig(alleles, ref_alleles);

				switch (config) {
					case STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND:
						opflip++;
						if (reorder) {
							data[i][indices[1]] = Sequence.flip(alleles[1]);
							data[i][indices[2]] = Sequence.flip(alleles[0]);
							rowsToKeep[i] = true;
							break;
						}
					case STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:
						flip++;
						// flip strand, keep allele order
						data[i][indices[1]] = Sequence.flip(alleles[0]);
						data[i][indices[2]] = Sequence.flip(alleles[1]);
						rowsToKeep[i] = true;
						break;
					case STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND:
						opp++;
						if (reorder) {
							data[i][indices[1]] = alleles[1];
							data[i][indices[2]] = alleles[0];
							rowsToKeep[i] = true;
							break;
						}
					case STRAND_CONFIG_SAME_ORDER_SAME_STRAND:
						same++;
						// keep alleles as is
						rowsToKeep[i] = true;
						break;
					case STRAND_CONFIG_AMBIGUOUS:
						ambg++;
						if (freq < 0.3 || (freq > 0.5 && freq - 0.5 < 0.3)) {
							rowsToKeep[i] = true;
							if (alleles[0].equals(ref_alleles[1]) && reorder) {
								data[i][indices[1]] = alleles[1];
								data[i][indices[2]] = alleles[0];
							}
						} else {
							log.report("MAF is too high for marker " + line[0]);
						}
						break;
					case STRAND_CONFIG_DIFFERENT_ALLELES:
					case STRAND_CONFIG_BOTH_NULL:
					case STRAND_CONFIG_SPECIAL_CASE:
					default:
						err++;
						// drop this allele pair bc something's wrong
						log.report("Unable to map marker " + line[indices[0]]);
						break;
				}
			}

			data = Matrix.prune(data, rowsToKeep, null, log);
			Files.writeMatrix(data, ext.rootOf(filename, false) + "_allele_verified.txt", "\t");

			log.report("Opposite and flipped: " + opflip + "\tFlipped: " + (flip - opflip)
								 + "\tOpposite: " + opp + "\tSame: " + (same - opp) + "\tErr: " + err
								 + "\tAmbiguous: " + ambg);
		} catch (Exception e) {
			log.report("Something went wrong.");
			e.printStackTrace();
		}
	}

	private static String[][] addPosition(String[][] data, String posFile,
																				Logger log) throws Elision {

		if (posFile == null) {
			log.reportError("Unable to load chromosomes. Verification may be slower.");
			data = ArrayUtils.append(data, new String[data.length][data[0].length + 1]);
			return data;
		}
		log.report("Attempting to get chr mapping from posfile.");

		String[] keys = Matrix.extractColumn(data, 0);
		String[] header = Files.getHeaderOfFile(posFile, log);
		int[] indices = ext.indexFactors(new String[][] {Aliases.MARKER_NAMES, Aliases.CHRS},
																		 header, false, false, false);
		String[][] positions = HashVec.loadFileToStringMatrix(posFile, true, indices);

		String[][] results = Files.combineInMemory(keys, positions, "NA", true, true, log);
		data = ArrayUtils.append(data, results);

		// sort our data so that we can go through the ref faster
		Arrays.sort(data, new Comparator<String[]>() {
			@Override
			public int compare(final String[] entry1, final String[] entry2) {
				final String chr1 = entry1[entry1.length - 1];
				final String chr2 = entry2[entry2.length - 1];
				// put missing chrs at the end
				if (chr1.equals("NA") && !chr2.equals("NA")) {
					return -1;
				} else if (chr2.equals("NA")) {
					return 1;
				}

				// if chromosomes are not equal, sort by chrom
				return chr1.compareTo(chr2);
			}
		});

		return data;
	}

	public static void main(String[] args) {
		String filename = "C:/Users/MR/bin/parkinsons/replication_InvVar1_mapped.txt";
		String refFile = "C:/Users/MR/Google Drive/statgen/1000G_work/Frequencies/chr#_eu_unrel.frq.xln";
		String freqFile = "C:/Users/MR/Google Drive/statgen/ParkinsonDisease/ReplicationData/freqs_genotypes_freq.xln";
		String posFile = "C:/Users/MR/bin/parkinsons/data/1000G_PD.map";

		String usage = "(1) Name of the file with alleles to verify. (eg file=Metal.tbl)\n"
									 + "(2) Name of the reference file containing expected alleles (eg ref=chr#_eu_unrel.frq.xln (# may be used to indicate a chromosome number))\n"
									 + "(3) Name of the file containing MAF. May be the same as one of the other files (eg freq=freqs_genotypes_freq.xln)\n"
									 + "(Optional) Name of file containing chromosome mappings for the given markers (eg pos=1000G_PD.map)";

		for (String arg : args) {
			if (arg.startsWith("file="))
				filename = ext.parseStringArg(arg);
			else if (arg.startsWith("ref="))
				refFile = ext.parseStringArg(arg);
			else if (arg.startsWith("freq="))
				freqFile = ext.parseStringArg(arg);
			else if (arg.startsWith("pos="))
				posFile = ext.parseStringArg(arg);
			else {
				System.out.println(usage);
				System.exit(0);
			}
		}

		verifyAlleles(filename, refFile, freqFile, posFile, false, new Logger());
	}
}
