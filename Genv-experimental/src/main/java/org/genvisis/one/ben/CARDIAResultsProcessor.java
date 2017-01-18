package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class CARDIAResultsProcessor {
	static final String REGRESSION_TEMPLATE = "regression_chr#.out_add.out.txt";
	static final String POS_TEMPLATE = "chr#_positions.proc.xln";
	static final String OUTFILE_TEMPLATE = "processed_chr#.out";
	static final String HEADER = "Markername,Chr,Pos,N,Effect_allele,Other_allele,EAF,Imp_info,Beta,SE,Pvalue";

	public static void run(String regDir, String posDir) {

		for (int i = 1; i < 23; i++) {
			BufferedReader readerR;
			try {
				readerR = Files.getAppropriateReader(regDir + REGRESSION_TEMPLATE.replace("#", "" + i));
				BufferedReader readerP = Files.getAppropriateReader(posDir
																														+ POS_TEMPLATE.replace("#", "" + i));
				PrintWriter writer = Files.getAppropriateWriter(regDir
																												+ OUTFILE_TEMPLATE.replace("#", "" + i));

				String lineR = readerR.readLine();
				String lineP = readerP.readLine();
				String[] hdr = lineR.split("[\\s]+");
				hdr = Array.insertStringAt("Chr", hdr, 1);
				hdr = Array.insertStringAt("Pos", hdr, 2);
				writer.println(Array.toStr(hdr));
				while ((lineR = readerR.readLine()) != null && (lineP = readerP.readLine()) != null) {
					String[] partsR = lineR.split("[\\s]+");
					String[] partsP = lineP.split("[\\s]+");

					if (!partsR[0].equals(partsP[0])) {
						// TODO error on mismatched RS
						System.err.println("ERROR - MISMATCHED RS/SNP NAMES!");
						break;
					}

					StringBuilder sb = new StringBuilder();
					sb.append(partsP[0]).append(",").append(partsP[1]).append(",").append(partsP[2])
						.append(",").append(partsR[7]).append(",").append(partsR[1]).append(",")
						.append(partsR[2]).append(",").append(partsR[4]).append(",").append(partsR[6])
						.append(",").append(partsR[9]).append(",").append(partsR[10]).append(",")
						.append(partsR[11]);

					writer.println(sb.toString());

				}

				writer.flush();
				writer.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dirR = "";
		String dirP = "";

		String usage = "\n"	+ "one.ben.CARDIAResultsProcessor requires 0-1 arguments\n"
										+ "   (1) regression results directory (i.e. dirR=" + dirR + " (default))\n"
										+ "   (2) SNP positions directory (i.e. dirP=" + dirP + " (default))\n" + "";


		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("dirR=")) {
				dirR = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("dirP=")) {
				dirP = arg.split("=")[1];
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
			run(ext.verifyDirFormat(dirR), ext.verifyDirFormat(dirP));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
