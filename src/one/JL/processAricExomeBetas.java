package one.JL;

import htsjdk.tribble.annotation.Strand;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import one.JL.BetaOptimizer.MarkerRsFormat;
import seq.manage.StrandOps;
import seq.manage.StrandOps.CONFIG;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;
import cnv.filesys.ABLookup;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.manage.Resources.GENOME_BUILD;
import cnv.manage.Resources.GENOME_RESOURCE_TYPE;

public class processAricExomeBetas {

	public static void main(String[] args) {

		String flipFile = "/panfs/roc/groups/5/pankrat2/shared/MitoPipeLineResources/ExomeCHip/alleleExome.txt";

		Hashtable<String, String> flipLook = HashVec.loadFileToHashString(flipFile, 0, new int[] { 1, 2 }, "\t", true);
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
			BetaOptimizer.mapToRsIds(markerSet, abLookup, GENOME_RESOURCE_TYPE.DB_SNP147.getResource(GENOME_BUILD.HG19).getResource(proj.getLog()), markerSet.getMarkerNames(), mapSer, proj.getLog());
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

		String beta = "/home/pankrat2/shared/skatMeta/exome_chip_hematology/WBC_TOTAL/Whites/SingleSNP/Whites_WBC_TOTAL_SingleSNP.csv";
		Logger log = proj.getLog();
		String[] betaHeader = Files.getHeaderOfFile(beta, ",", new Logger());
		for (int i = 0; i < betaHeader.length; i++) {
			betaHeader[i] = betaHeader[i].replaceAll("\"", "");
		}
		int markerIndex = ext.indexOfStr("Name", betaHeader);
		int betaIndex = ext.indexOfStr("beta", betaHeader);
		int pvalIndex = ext.indexOfStr("p", betaHeader);

		System.out.println(markerIndex + "\t" + Array.toStr(betaHeader));
		String outDir = "/home/pankrat2/shared/MitoPipeLineResources/betas/";
		new File(outDir).mkdirs();
		String outBeta = outDir + ext.rootOf(beta) + "matched.txt";
		String outBetaFinal = outDir + ext.rootOf(beta) + "matched.final.txt";

		try {
			// PrintWriter writer = new PrintWriter(new FileWriter(outBeta));
			PrintWriter writerFinal = new PrintWriter(new FileWriter(outBetaFinal));

			BufferedReader reader = Files.getAppropriateReader(beta);
			// writer.println(Array.toStr(betaHeader) + "\tMarkerName\tPos\trsID\trsRef\trsAlt\tmarkerA\tmarkerB\tconfig\tsiteType");
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
					if (m.flipBetas()) {
						betaVal = -1 * betaVal;
					}
					writerFinal.println(m.getRs() + "\t" + Array.toStr(m.getDbSnpAlleles()) + "\t" + betaVal + "\t" + line[pvalIndex] + "\t" + Array.toStr(m.getMarkerAlleles()) + "\t" + m.getConfig() + "\t" + Array.toStr(m.getDbSnpAlleles()) + "\t" + effAllele);
				}

			}
			writerFinal.close();
			reader.close();
			// writer.close();
			System.exit(1);

		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + beta + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + beta + "\"");
			return;
		}

		Project proj2 = new Project("/home/pankrat2/lanej/projects/Aric_gw6.properties", false);
		String out2 = proj2.PROJECT_DIRECTORY.getValue() + "betaOpti/";
		new File(out2).mkdirs();
		proj2.AB_LOOKUP_FILENAME.setValue(out2 + "AB_LookupBeta.dat");
		ABLookup abLookup2 = new ABLookup();
		if (!Files.exists(proj2.AB_LOOKUP_FILENAME.getValue())) {
			abLookup2.parseFromAnnotationVCF(proj2);
			abLookup2.writeToFile(proj2.AB_LOOKUP_FILENAME.getValue(), proj2.getLog());
		}
		MarkerSet markerSet2 = proj2.getMarkerSet();
		abLookup2 = new ABLookup(markerSet2.getMarkerNames(), proj2.AB_LOOKUP_FILENAME.getValue(), true, true, proj2.getLog());
		String mapSer2 = out2 + "rsIDMap.ser";
		BetaOptimizer.run(markerSet2, abLookup2, GENOME_RESOURCE_TYPE.DB_SNP147.getResource(GENOME_BUILD.HG19).getResource(proj2.getLog()), proj2.getNonCNMarkers(), mapSer2, outBeta, proj2.getLog());

	}
}
