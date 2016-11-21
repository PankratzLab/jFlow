// -Xmx1024M
package org.genvisis.parse;

import java.util.List;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class LookupTable {
	public static void fromParameters(String filename, Logger log) {
		String hitsFile;
		String[] line, hits;
		int col;
		boolean ignoreCase, ignoreFirstLine, commaDelimited, tabDelimited;
		String outfile;
		List<String> params;
		String head, missingValue;
		boolean finalHeader, hideIndex, lessMemoryButSlower, keepIntermediateFiles;

		params = Files.parseControlFile(filename, "lookup",
																		new String[] {"topHits.txt 0 out=results.xln head=IID missingValue=. ignoreCase noFinalHeader hideIndex lessMemoryButSlower keepIntermediateFiles",
																									"F7_SingleSNP.csv , 1 3=MAF",
																									"plink.bim noHeader 1 0=Chr 3=loc_36.1",
																									"freq.frq 1 4", "plink.assoc 1 4 5 8 3",
																									"logistic.xls 0 4 6 7 9 10 12",
																									"#logistic.xls 0 4 6 7 9 10 12 13 15",
																									"#additive.assoc.logistic !4=ADD 1 6=AdditiveOR 8=AdditivePvalue",
																									"#dominant.assoc.logistic !4=DOM 1 6=DominantOR 8=DominantPvalue",
																									"#recessive.assoc.logistic !4=REC 1 6=Recessive 8=RecessivePvalue",
																									"#genotypic.assoc.logistic !4=GENO_2DF 1 8=GenotypicPvalue",
																									"#missingness.txt 0 3=MissingRate 4=HWE_pvalue 7=MissingByPhenotype_pvalue 6=MissingByHaplotype_minimum_pvalue",
																									"missing.lmiss 1 4=MissingRate",
																									"hardy.hwe !2=UNAFF 1 8=HWE_pvalue",
																									"test.missing.missing 1 4=MissingByPhenotype_pvalue",
																									"gender.missing 1 4=MissingByGender_pvalue",
																									"gender.assoc 1 5=Freq_Male 4=Freq_Female 8=AssocByGender_pvalue",
																									"mishap.missing.hap !1=HETERO 0 7=MissingByHaplotype_minimum_pvalue",
																									"#example.txt 0 5=Needed;0 fail",
																									"# annotation.csv 0 1=CommentWithQuotes doNotSimplifyQuotes"},
																		log);
		if (params != null) {
			line = params.remove(0).trim().split("[\\s]+");
			hitsFile = line[0];
			col = 0;
			ignoreCase = false;
			commaDelimited = line[0].endsWith(".csv");
			tabDelimited = false;
			ignoreFirstLine = false;
			finalHeader = true;
			hideIndex = false;
			lessMemoryButSlower = false;
			keepIntermediateFiles = false;
			head = "SNP";
			missingValue = ".";
			outfile = ext.rootOf(hitsFile) + "_described.xln";
			for (int j = 1; j < line.length; j++) {
				if (line[j].equalsIgnoreCase("ignoreCase")) {
					ignoreCase = true;
				} else if (line[j].startsWith("ignore") || line[j].equalsIgnoreCase("header")) {
					ignoreFirstLine = true;
				} else if (line[j].equalsIgnoreCase("nofinalheader")) {
					finalHeader = false;
				} else if (line[j].equalsIgnoreCase("hideIndex")) {
					hideIndex = true;
				} else if (line[j].startsWith("head=")) {
					head = line[j].split("=")[1];
				} else if (line[j].startsWith("missingValue=")) {
					missingValue = line[j].split("=")[1];
				} else if (line[j].startsWith("out=")) {
					outfile = line[j].split("=")[1];
				} else if (line[j].equals(",")) {
					commaDelimited = true;
				} else if (line[j].equals("tab")) {
					tabDelimited = true;
					commaDelimited = false;
				} else if (line[j].equals("lessMemoryButSlower")) {
					lessMemoryButSlower = true;
				} else if (line[j].equals("keepIntermediateFiles")) {
					keepIntermediateFiles = true;
				} else {
					col = Integer.parseInt(line[j]);
				}
			}

			log.report("Memory available: " + ext.prettyUpSize(Runtime.getRuntime().maxMemory(), 1) + "");
			log.report("Loading keys from '" + hitsFile + "'");
			hits = HashVec.loadFileToStringArray(	hitsFile, false, ignoreFirstLine, new int[] {col}, true,
																						false,
																						commaDelimited	? ","
																														: (tabDelimited ? "\t" : "[\\s]+"));
			if (lessMemoryButSlower) {
				Files.combineWithLessMemory(hits, Array.toStringArray(params), new String[params.size()][],
																		head, missingValue, outfile, log, ignoreCase, finalHeader,
																		hideIndex, keepIntermediateFiles);
			} else {
				Files.combine(hits, Array.toStringArray(params), new String[params.size()][], head,
											missingValue, outfile, log, ignoreCase, finalHeader, hideIndex);
			}
		}
	}
}
