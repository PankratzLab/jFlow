package org.genvisis.one.ben;

import java.io.File;
import java.io.PrintWriter;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;

public class FCSBooleanComparator {

	String srcDir;
	String srcLbl;

	List<String> compDirs = new ArrayList<>();
	List<String> compLbls = new ArrayList<>();

	String out;

	Logger log = new Logger();

	public void run() {
		String[] dirs1 = new File(srcDir).list((d, f) -> {
			return new File(srcDir + f).isDirectory();
		});
		List<String[]> compSubDirs = new ArrayList<String[]>();
		for (int i = 0; i < compDirs.size(); i++) {
			String compDir = compDirs.get(i);
			compSubDirs.add(new File(compDir).list((d, f) -> {
				return new File(compDir + f).isDirectory();
			}));
		}

		HashSet<String> sharedDirs = new HashSet<>();
		for (String s : dirs1) {
			sharedDirs.add(s);
		}
		HashSet<String> dirSet2 = new HashSet<>();
		for (String[] subDirs : compSubDirs) {
			for (String s : subDirs) {
				dirSet2.add(s);
			}
		}
		sharedDirs.retainAll(dirSet2);

		PrintWriter diff = Files.getAppropriateWriter(out + "diff.txt");
		for (String s : dirs1) {
			if (!sharedDirs.contains(s)) {
				diff.println(s);
			}
		}
		diff.println("--");
		for (String[] sub : compSubDirs) {
			for (String s : sub) {
				if (!sharedDirs.contains(s)) {
					diff.println(s);
				}
			}
		}
		diff.close();

		log.reportTime("Found " + sharedDirs.size() + " shared samples across " + (compDirs.size() + 1)
									 + " revisions.");

		String outFile = out + "out_" + srcLbl;
		for (String s : compLbls) {
			outFile += "_" + s;
		}

		PrintWriter writer1 = Files.getAppropriateWriter(outFile + ".indiv.xln");
		String hdr = "Gate\tSample\tTP\tFP\tTN\tFN";
		writer1.println(hdr);

		// revision dir -> sample dir -> gate -> conf matrix
		Map<String, Map<String, Map<String, int[]>>> sampConfMap = new HashMap<>();
		for (int i = 0; i < compDirs.size(); i++) {
			String compDir = compDirs.get(i);
			Map<String, Map<String, int[]>> compMap = new HashMap<>();
			for (String s : sharedDirs) {
				compMap.put(s, processDir(writer1, compDir, compLbls.get(i), s + "/"));
			}
			sampConfMap.put(compDir, compMap);
		}
		writer1.close();

		// shared dirs is shared samples
		HashSet<String> allGates = new HashSet<String>();
		for (Map<String, Map<String, int[]>> sampGateMatrs : sampConfMap.values()) {
			if (sampGateMatrs == null)
				continue;
			for (Map<String, int[]> gateMatrs : sampGateMatrs.values()) {
				if (gateMatrs == null) {
					continue;
				}
				allGates.addAll(gateMatrs.keySet());
			}
		}
		List<String> orderedGates = new ArrayList<>(allGates);
		allGates = null;

		PrintWriter writer = Files.getAppropriateWriter(outFile + ".xln");

		String[] hdrLbls = {
												"TP_mn",
												"FP_mn",
												"TN_mn",
												"FN_mn",
												"P_mn",
												"R_mn",
												"F_mn",
												"A_mn"
		};

		hdr = "Gate";
		for (int c = 0; c < compDirs.size(); c++) {
			for (String h : hdrLbls) {
				hdr += "\t" + srcLbl + "->" + compLbls.get(c) + "_" + h;
			}
		}
		writer.println(hdr);

		double N = 0;
		long[] v = {0, 0, 0, 0};
		double[] prf = {0, 0, 0};
		double[] prf1 = {0, 0, 0};
		DecimalFormat aveForm = new DecimalFormat("#.######");
		aveForm.setRoundingMode(RoundingMode.HALF_EVEN);
		for (String gate : orderedGates) {
			StringBuilder ln = new StringBuilder(gate);
			for (int c = 0; c < compDirs.size(); c++) {
				N = 0;
				v = new long[] {0, 0, 0, 0};
				Map<String, Map<String, int[]>> sampGateMtrs = sampConfMap.get(compDirs.get(c));
				for (Map<String, int[]> gateConfs : sampGateMtrs.values()) {
					if (gateConfs.containsKey(gate)) {
						N++;
						// tp, fp, tn, fn
						v[0] += gateConfs.get(gate)[0];
						v[1] += gateConfs.get(gate)[1];
						v[2] += gateConfs.get(gate)[2];
						v[3] += gateConfs.get(gate)[3];
						prf1 = getPRF(gateConfs.get(gate));
						prf[0] += prf1[0];
						prf[1] += prf1[1];
						prf[2] += prf1[2];
					}
				}

				ln.append("\t")
					.append(aveForm.format(v[0] / N))
					.append(aveForm.format(v[1] / N))
					.append(aveForm.format(v[2] / N))
					.append(aveForm.format(v[3] / N))
					.append(aveForm.format(prf[0] / N))
					.append(aveForm.format(prf[1] / N))
					.append(aveForm.format(prf[2] / N))
					.append(aveForm.format((v[0] + v[2]) / (v[0] + v[1] + v[2] + v[3])));
			}
			ln.replace(0, gate.length(), gate + "_" + N);

			writer.println(ln.toString());
		}

		writer.close();
	}

	private Map<String, int[]> processDir(PrintWriter writer, String compDir, String compLbl,
																				String subDir) {
		String f1 = srcDir + subDir
								+ (subDir.endsWith("/") ? subDir.substring(0, subDir.length() - 1) : subDir)
								+ ".incl.xln";
		String f2 = compDir + subDir
								+ (subDir.endsWith("/") ? subDir.substring(0, subDir.length() - 1) : subDir)
								+ ".incl.xln";

		int cnt1 = Files.countLines(f1, 0);
		int cnt2 = Files.countLines(f2, 0);

		if (cnt1 != cnt2) {
			log.reportError("Files did not match line counts: " + f1 + " {" + cnt1 + "}; " + f2 + " {"
											+ cnt2 + "}");
			return null;
		}

		String[] hdr1 = Files.getHeaderOfFile(f1, log);
		String[] hdr2 = Files.getHeaderOfFile(f2, log);

		int ind;
		String h1, h2;
		HashMap<String, int[]> columnInds = new HashMap<>();
		for (int i1 = 0; i1 < hdr1.length; i1++) {
			h1 = (ind = hdr1[i1].indexOf('(')) > 0 ? hdr1[i1].substring(0, ind).trim() : hdr1[i1].trim();
			for (int i2 = 0; i2 < hdr2.length; i2++) {
				h2 = (ind = hdr2[i2].indexOf('(')) > 0 ? hdr2[i2].substring(0, ind).trim()
																							 : hdr2[i2].trim();

				if (h1.equals(h2)) {
					columnInds.put(h2, new int[] {i1, i2});
					break;
				}
			}
		}

		log.reportTime("Found " + columnInds.size() + " shared columns for " + subDir);

		HashMap<String, int[]> confMap = new HashMap<>();

		ArrayList<String> cols = new ArrayList<>(columnInds.keySet());

		String[] data1;
		String[] data2;
		for (String col : cols) {
			if (col.equalsIgnoreCase("index")) {
				continue;
			}
			int[] inds = columnInds.get(col);
			data1 = HashVec.loadFileToStringArray(f1, true, new int[] {inds[0]}, false);
			data2 = HashVec.loadFileToStringArray(f2, true, new int[] {inds[1]}, false);
			int tp = 0;
			int fp = 0;
			int tn = 0;
			int fn = 0;

			for (int l = 0; l < data1.length; l++) {
				boolean t = data1[l].equalsIgnoreCase("true");
				if (data1[l].equalsIgnoreCase(data2[l])) {
					if (t) {
						tp += 1;
					} else {
						tn += 1;
					}
				} else {
					if (t) {
						fn += 1;
					} else {
						fp += 1;
					}
				}
			}

			confMap.put(col, new int[] {tp, fp, tn, fn});
			writer.println(col + "\t" + subDir + "\t" + tp + "\t" + fp + "\t" + tn + "\t" + fn);
		}
		writer.flush();

		return confMap;
	}

	private double[] getPRF(int[] conf) {
		double tp = conf[0];
		double fp = conf[1];
		double fn = conf[3];
		double p;
		double r;
		double f;

		p = (((double) tp) / (double) (tp + fp));
		r = (((double) tp) / (double) (tp + fn));
		f = 2d * ((double) (p * r)) / (double) (p + r);

		return new double[] {p, r, f};
	}

	private double[] getPRF(double[] conf) {
		double tp = conf[0];
		double fp = conf[1];
		double fn = conf[3];
		double p;
		double r;
		double f;

		p = (((double) tp) / (double) (tp + fp));
		r = (((double) tp) / (double) (tp + fn));
		f = 2d * ((double) (p * r)) / (double) (p + r);

		return new double[] {p, r, f};
	}

	public static void main(String[] args) {
		int numArgs = args.length;

		String f = null;
		FCSBooleanComparator comp = new FCSBooleanComparator();

		String usage = "\n" +
									 "org.genvisis.one.ben.FCSBooleanComparator requires 0-1 arguments\n" +
									 "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
					|| args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("in=")) {
				f = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				comp.out = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (args.length == 0 || numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		String[][] d = HashVec.loadFileToStringMatrix(f, false, new int[] {0, 1});
		comp.srcLbl = d[0][0];
		comp.srcDir = d[0][1];

		for (int i = 1; i < d.length; i++) {
			comp.compLbls.add(d[i][0]);
			comp.compDirs.add(d[i][1]);
		}

		comp.run();

	}



}
