package org.genvisis.widgets;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Date;
import java.util.Vector;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;


public class diff {
	public static void addFileToQueue(String filename) {
		String[] files;
		Calendar midnight;
		long fingerprint;
		long[] prints;

		midnight = Calendar.getInstance();
		midnight.set(Calendar.HOUR_OF_DAY, 0);
		midnight.set(Calendar.MINUTE, 0);
		midnight.set(Calendar.SECOND, 0);
		midnight.set(Calendar.MILLISECOND, 0);

		filename = ext.removeDirectoryInfo(filename);

		fingerprint = new Date().getTime() - midnight.getTimeInMillis();
		fingerprint = fingerprint * 1000 + (int) (Math.random() * 1000);
		// fingerprint = (10000+Long.parseLong(ext.right(fingerprint+"",
		// 5)))*1000+(int)(Math.random()*1000);
		// System.out.println(fingerprint);
		Files.write(filename, fingerprint + ".difftemp");
		for (int i = 0; i < 15; i++) {
			try {
				Thread.sleep(200);
			} catch (InterruptedException ie) {
			}
			files = Files.list("./", ".difftemp", false);
			if (files.length > 1) {
				if (files.length > 2) {
					System.err.println("Error - Select only 2 files at a time");
					System.err.println("        ...proceeding with first two files and then exiting");
				}
				prints = new long[files.length];
				for (int j = 0; j < files.length; j++) {
					prints[j] = Long.parseLong(files[j].substring(0, files[j].indexOf(".")));
				}
				Arrays.sort(prints);
				if (fingerprint == prints[0]) {
					runDiff(filename, HashVec.loadFileToStringArray(prints[1]	+ ".difftemp", false,
																													null, false)[0]);
					deleteAll();
					return;
				} else {
					try {
						Thread.sleep(300);
					} catch (InterruptedException ie) {
					}
					deleteAll();
					return;
				}
			} else {
				if (i == 3) {
					System.err.println("diff requires 2 files... still waiting for another file to open...");
				}
				if (i == 14) {
					System.err.println("Error - second file never opened, highlight 2 (ctrl+left click) and try again");
				}
			}
		}
	}

	public static void deleteAll() {
		String[] files;

		files = Files.list("./", ".difftemp", false);
		for (String file : files) {
			new File(file).delete();
		}
	}

	public static void runDiff(String file1, String file2) {
		PrintWriter writer;
		String trav1, trav2;
		BufferedReader reader1, reader2;
		Vector<String> buffer1, buffer2;
		boolean nothingChanged;
		int count1, count2;

		try {
			nothingChanged = true;
			reader1 = new BufferedReader(new FileReader(file1));
			reader2 = new BufferedReader(new FileReader(file2));
			writer = null;
			trav1 = trav2 = "";
			count1 = count2 = 0;
			buffer1 = new Vector<String>();
			buffer2 = new Vector<String>();
			while (trav1 != null || trav2 != null) {
				if (reader1.ready()) {
					trav1 = reader1.readLine();
					count1++;
				} else {
					trav1 = null;
				}

				if (reader2.ready()) {
					trav2 = reader2.readLine();
					count2++;
				} else {
					trav2 = null;
				}

				if (trav1 != null && trav2 != null && !trav1.equals(trav2)	|| trav1 == null && trav2 != null
						|| trav1 != null && trav2 == null || buffer1.size() + buffer2.size() > 0) {
					if (nothingChanged) {
						writer = new PrintWriter(new FileWriter("diff bw '"	+ file1 + "' & '" + file2
																										+ "'.out"));
						nothingChanged = false;
					}
					if (trav1 != null) {
						buffer1.addElement(trav1);
					}
					if (trav2 != null) {
						buffer2.addElement(trav2);
					}
					int pos1 = 0;
					while (pos1 < buffer1.size()) {
						for (int pos2 = 0; pos2 < buffer2.size(); pos2++) {
							if (!buffer1.elementAt(pos1).equals(buffer2.elementAt(pos2))
									&& (trav1 != null && trav2 != null	|| pos1 != buffer1.size() - 1
											|| pos2 != buffer2.size() - 1)) {
								continue;
							}
							if (pos1 == pos2) {
								for (int i = 0; i < (trav1 != null || trav2 != null ? pos1 : pos1 + 1); i++) {
									int offsets[] = whatsTheDiff(buffer1.elementAt(i), buffer2.elementAt(i));
									writer.println(count1 + " <<del<< " + buffer1.elementAt(i));
									writer.println(count2 + " >>ins>> " + buffer2.elementAt(i));
									if (offsets[0] == 1) {
										writer.println("  Swapped '"
																			+ buffer1.elementAt(i).substring(offsets[1], offsets[2])
																		+ "' for '"
																		+ buffer2.elementAt(i).substring(offsets[1], offsets[3]) + "'");
									}
								}

								for (int i = 0; i <= pos1; i++) {
									buffer1.removeElementAt(0);
									buffer2.removeElementAt(0);
								}

								pos2 = buffer2.size();
								continue;
							}
							for (int i = 0; i < (trav1 != null || trav2 != null ? pos1 : pos1 + 1); i++) {
								writer.println(count1 + " <<del<< " + buffer1.elementAt(i));
							}

							for (int i = 0; i <= pos1; i++) {
								buffer1.removeElementAt(0);
							}

							for (int i = 0; i < (trav1 != null || trav2 != null ? pos2 : pos2 + 1); i++) {
								writer.println(count2 + " >>ins>> " + buffer2.elementAt(i));
							}

							for (int i = 0; i <= pos2; i++) {
								buffer2.removeElementAt(0);
							}

							pos2 = buffer2.size();
						}

						pos1++;
					}
				}
			}
			for (int i = 0; i < buffer1.size(); i++) {
				writer.println("<<del<< " + buffer1.elementAt(i));
			}

			for (int i = 0; i < buffer2.size(); i++) {
				writer.println(">>ins>> " + buffer2.elementAt(i));
			}

			reader1.close();
			reader2.close();
			if (nothingChanged) {
				System.err.println("These files (" + file1 + " & " + file2 + ") are identical");
				ext.waitForResponse();
			} else {
				writer.close();
			}
		} catch (Exception e) {
			System.err.println("Error diff'ing " + file1 + " and " + file2);
			e.printStackTrace();
		}
	}

	public static int[] whatsTheDiff(String str1, String str2) {
		int offs[];
		for (offs = new int[4]; offs[1] < str1.length()	&& offs[1] < str2.length()
														&& str1.charAt(offs[1]) == str2.charAt(offs[1]); offs[1]++) {
			;
		}
		offs[2] = str1.length();
		for (offs[3] = str2.length(); offs[2] > offs[1]	&& offs[3] > offs[1]
																	&& str1.charAt(offs[2] - 1) == str2.charAt(offs[3]
																																							- 1); offs[3]--) {
			offs[2]--;
		}

		if (offs[2] == offs[1] && offs[3] == offs[1]) {
			offs[0] = 2;
		} else if ((double) (offs[2] - offs[1]) / (double) str1.length() < 0.5D
								&& (double) (offs[3] - offs[1]) / (double) str2.length() < 0.5D) {
			offs[0] = 1;
		} else {
			offs[0] = 0;
		}
		return offs;
	}

	public static void main(String args[]) {
		int numArgs = args.length;
		String file1 = "file1.txt";
		String file2 = "file2.txt";

		String usage = "\n"	+ "consol.diff compares 2 files and reports back where they are different\n"
										+ "     requires 2 arguments, the 2 filenames (defaults: '" + file1 + "' and '"
										+ file2 + "')\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				ext.waitForResponse();
				return;
			}
		}

		try {
			if (numArgs == 1) {
				addFileToQueue(args[0]);
			} else if (numArgs == 2) {
				runDiff(args[0], args[1]);
			} else {
				System.err.println(usage);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		ext.waitForResponse();
	}
}
