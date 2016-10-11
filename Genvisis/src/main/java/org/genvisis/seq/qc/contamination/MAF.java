package org.genvisis.seq.qc.contamination;

import java.io.FileWriter;
import java.io.PrintWriter;

import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MDL;
import org.genvisis.common.Logger;

public class MAF {

	public static void computeMAF(Project proj) {
		Logger log;

		log = proj.getLog();
		MDL mdl = new MDL(proj, proj.getMarkerSet(), proj.getMarkerNames(), 4, 4);
		String output = proj.PROJECT_DIRECTORY.getValue() + "mafs.txt";
		// String sampOutliers = proj.getDir(Project.SAMPLE_DIRECTORY)+"outliers.ser";
		// Hashtable<String, Float> outSamps = (Hashtable<String, Float>)
		// Files.readSerial(sampOutliers);
		// System.out.println(outSamps.get("482175\tF10139D\tlrr"));
		// TransposeData.r
		// System.exit(1);
		// for(String sampKey:outSamps.keySet()){
		// System.out.println(sampKey);
		// }

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.println("NAME\tMAF\tCALLRATE\tA_AF\tB_AF");
			mdl.setDebugMode(true);
			mdl.setReportEvery(10000);
			while (mdl.hasNext()) {
				MarkerData markerData = mdl.next();
				try {
					double maf = markerData.getMAF(null, null, null, 0, log);
					// int[] counts = markerData.getGenotypeCounts(null, null, null, 0, log);
					// AlleleFreq.calcMAF(counts[0], counts[1], counts[2])
					byte[] genotypes = markerData.getAbGenotypes();
					int noCall = 0;
					for (byte genotype : genotypes) {
						if (genotype < 0) {
							noCall++;
						}
					}
					int called = genotypes.length - noCall;
					double callRate = (double) called / genotypes.length;
					writer.println(markerData.getMarkerName() + "\t" + maf + "\t" + callRate);
				} catch (ArrayIndexOutOfBoundsException aoe) {
					byte[] genos = markerData.getAbGenotypes();
					for (int i = 0; i < genos.length; i++) {
						if (genos[i] == 3) {
							System.out.println(markerData.getMarkerName() + "\t" + proj.getSamples()[i]);
						}
					}
				}

			}
			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + output);
			proj.getLog().reportException(e);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;

		String usage = "\n"	+ "one.JL.MAF requires 0-1 arguments\n" + "   (1) filename (i.e. proj="
										+ filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			computeMAF(new Project(filename, false));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
