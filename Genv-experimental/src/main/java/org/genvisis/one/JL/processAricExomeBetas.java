package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import org.genvisis.cnv.analysis.pca.BetaOptimizer;
import org.genvisis.cnv.analysis.pca.BetaOptimizer.MarkerRsFormat;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.StrandOps;
import org.genvisis.seq.manage.StrandOps.CONFIG;

public class processAricExomeBetas {

	public static void main(String[] args) {

		String flipFile = "/home/pankrat2/shared/MitoPipeLineResources/betas/ExomeBetas/ExomeChipFlip.txt";

		Hashtable<String, String> flipLook = HashVec.loadFileToHashString(flipFile, 0, new int[] { 1, 2 }, "\t", false);
		Project proj = new Project("/home/pankrat2/lanej/projects/aric_exome.properties", false);
		String out = proj.PROJECT_DIRECTORY.getValue() + "betaOpti/";
		new File(out).mkdirs();
		proj.AB_LOOKUP_FILENAME.setValue(out + "AB_LookupBeta.dat");
		ABLookup abLookup = new ABLookup();
		if (!Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
			abLookup.parseFromAnnotationVCF(proj);
			abLookup.writeToFile(proj.AB_LOOKUP_FILENAME.getValue(), proj.getLog());
		}
		MarkerSet markerSet = proj.getMarkerSet();
		abLookup = new ABLookup(markerSet.getMarkerNames(), proj.AB_LOOKUP_FILENAME.getValue(), true, true, proj.getLog());
		String mapSer = out + "rsIDMap.ser";
		if (!Files.exists(mapSer)) {
			BetaOptimizer.mapToRsIds(proj, abLookup, Resources.genome(GENOME_BUILD.HG19, proj.getLog()).getDBSNP().get(), markerSet.getMarkerNames(), mapSer, proj.getLog());
		}
		ArrayList<MarkerRsFormat> markerRsFormats = MarkerRsFormat.readSerial(mapSer, proj.getLog());
		ArrayList<String> rsOut = new ArrayList<String>();
		rsOut.add("MarkerName\tPos\trsID\trsRef\trsAlt\tmarkerA\tmarkerB\tconfig\tsiteType");
		for (int i = 0; i < markerRsFormats.size(); i++) {
			MarkerRsFormat m = markerRsFormats.get(i);
			rsOut.add(m.getMarkerName() + "\t" + m.getPosMarker() + "\t" + m.getRs() + "\t" + Array.toStr(m.getDbSnpAlleles()) + "\t" + Array.toStr(m.getMarkerAlleles()) + "\t" + m.getConfig() + "\t" + m.getType());
		}
		Files.writeArrayList(rsOut, out + "rsmatch.txt");

		Hashtable<String, Integer> match = new Hashtable<String, Integer>();
		for (int i = 0; i < markerRsFormats.size(); i++) {
			match.put(markerRsFormats.get(i).getMarkerName(), i);
		}

		ArrayList<String> betas = new ArrayList<String>();
		betas.add("/home/pankrat2/shared/skatMeta/exome_chip_hematology/WBC_TOTAL/Whites/SingleSNP/Whites_WBC_TOTAL_SingleSNP.csv");
		betas.add("/home/pankrat2/shared/skatMeta/exome_chip_hematology/WBC_TOTAL/Blacks/SingleSNP/Blacks_WBC_TOTAL_SingleSNP.csv");

		String[] others = Files.list("/home/pankrat2/shared/MitoPipeLineResources/betas/ExomeBetas/", "", ".csv", true, false, true);
		for (int i = 0; i < others.length; i++) {
			betas.add(others[i]);
		}
		// betas.add("/home/pankrat2/shared/skatMeta/exome_chip_hematology/WBC_TOTAL/Asians/SingleSNP/Asians_WBC_TOTAL_SingleSNP.csv");

		for (int i = 0; i < betas.size(); i++) {
			String beta = betas.get(i);
			System.out.println(beta);

			Logger log = proj.getLog();
			String[] betaHeader = Files.getHeaderOfFile(beta, ",", new Logger());
			for (int j = 0; j < betaHeader.length; j++) {
				betaHeader[j] = betaHeader[j].replaceAll("\"", "");
			}
			int markerIndex = ext.indexOfStr("Name", betaHeader);
			int betaIndex = ext.indexOfStr("beta", betaHeader);
			int pvalIndex = ext.indexOfStr("p", betaHeader);

			System.out.println(markerIndex + "\t" + Array.toStr(betaHeader));
			String outDir = "/home/pankrat2/shared/MitoPipeLineResources/betasExome2flip/";
			new File(outDir).mkdirs();
			String outBeta = outDir + ext.rootOf(beta) + "matched.txt";
			String outBetaFinal = outDir + ext.rootOf(beta) + "matched.final.beta";

			try {
				// PrintWriter writer = new PrintWriter(new FileWriter(outBeta));
				PrintWriter writerFinal = new PrintWriter(new FileWriter(outBetaFinal));

				BufferedReader reader = Files.getAppropriateReader(beta);
				// writer.println(Array.toStr(betaHeader) + "\tMarkerName\tPos\trsID\trsRef\trsAlt\tmarkerA\tmarkerB\tconfig\tsiteType");
				// writerFinal.println("rsID\tref\talt\ta2\tbeta\tp\tA\tB\tF1\tF2\tconfig\toriginalbeta");
				writerFinal.println("rsID\tref\talt\tbeta\tp");

				reader.readLine();
				while (reader.ready()) {

					String[] line = ext.splitCommasIntelligently(reader.readLine().trim(), true, log);

					if (!match.containsKey(line[markerIndex])) {
						System.err.println(line[markerIndex]);
						System.exit(1);
					}
					MarkerRsFormat m = markerRsFormats.get(match.get(line[markerIndex]));
					String[] aNo = m.getMarkerAlleles();
					String effAllele = flipLook.get(line[markerIndex]);
					double betaVal = Double.NaN;

					if (!effAllele.contains("0")) {
						CONFIG config = StrandOps.determineStrandConfig(m.getMarkerAlleles(), effAllele.split("\t"));
						try {
							betaVal = Double.parseDouble(line[betaIndex]);

						} catch (NumberFormatException nfe) {

						}
						switch (config) {
						case STRAND_CONFIG_BOTH_NULL:
							System.out.println(config + "\t" + line[markerIndex] + "\t" + effAllele + "\t" + Array.toStr(m.getDbSnpAlleles()) + "\t" + Array.toStr(m.getMarkerAlleles()) + "\t" + m.getConfig() + "\t" + m.getRs());

							break;
						case STRAND_CONFIG_DIFFERENT_ALLELES:
							System.out.println(config + "\t" + line[markerIndex] + "\t" + effAllele + "\t" + Array.toStr(m.getDbSnpAlleles()) + "\t" + Array.toStr(m.getMarkerAlleles()) + "\t" + m.getConfig() + "\t" + m.getRs());

							break;
						case STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND:// opposit order of the effect allele (forward Strand)
						case STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND:
							betaVal = -1 * betaVal;
							break;
						case STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:
						case STRAND_CONFIG_SAME_ORDER_SAME_STRAND:
							break;
						case STRAND_CONFIG_SPECIAL_CASE:
							System.out.println(config + "\t" + line[markerIndex] + "\t" + effAllele + "\t" + Array.toStr(m.getDbSnpAlleles()) + "\t" + Array.toStr(m.getMarkerAlleles()) + "\t" + m.getConfig() + "\t" + m.getRs());

							break;
						case STRAND_CONFIG_UNKNOWN:
							break;
						default:
							break;

						}
					}

					// line[betaIndex] = betaVal + "";
					// writer.println(Array.toStr(line) + "\t" + m.getMarkerName() + "\t" + m.getPosMarker() + "\t" + m.getRs() + "\t" + Array.toStr(m.getDbSnpAlleles()) + "\t" + Array.toStr(m.getMarkerAlleles()) + "\t" + m.getConfig() + "\t" + m.getType());
					if (m.isValidMatch() && Double.isFinite(betaVal)) {
						// if (m.flipBetas()) {
						// betaVal = -1 * betaVal;
						// }
						// writerFinal.println(m.getRs() + "\t" + Array.toStr(effAllele.split("\t")) + "\t" + betaVal + "\t" + line[pvalIndex] + "\t" + Array.toStr(m.getMarkerAlleles()) + "\t" + effAllele + "\t" + m.getConfig() + "\t" + line[betaIndex]);
						writerFinal.println(m.getRs() + "\t" + Array.toStr(effAllele.split("\t")) + "\t" + line[betaIndex] + "\t" + line[pvalIndex]);

					}

				}
				writerFinal.close();
				reader.close();
				// writer.close();

			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + beta + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + beta + "\"");
				return;
			}
		}
	}
}
