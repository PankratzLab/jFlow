package link.bat;

import java.io.*;
import java.util.*;

import common.*;

public class pdt {
	public static final int maxFams = 1000;

	public static final String DEFAULT_DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\David Ginsburg\\pdt_equalFreq\\";

	public static void createFiles(int chr, int markerIncrement, String templateFile) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		Vector<String> v = new Vector<String>();
		Vector<String[]> pedinfo = new Vector<String[]>();
		String[] line;
		String chrome = ext.chrome(chr);
		int numMarkers = -1;
		int numFams;

		try {
			reader = new BufferedReader(new FileReader(templateFile));
			while (reader.ready()) {
				pedinfo.add(Array.subArray(reader.readLine().split("[\\s]+"), 0, 10));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+templateFile+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+templateFile+"\" (possible there are not 10 columns for every row?)");
			ioe.printStackTrace();
			System.exit(2);
		}

		try {
			reader = new BufferedReader(new FileReader("map"+chrome+".dat"));
			numMarkers = Integer.parseInt(reader.readLine().trim().split("[\\s]+")[0])-1;
			for (int i = 0; i<6; i++) {
				reader.readLine();
			}
			for (int i = 1; i<=numMarkers; i++) {
				writer = new PrintWriter(new FileWriter(i+".dat"));
				writer.println("2 0 0 5");
				writer.println("0 0.0 0.0 0");
				writer.println("1 2");
				writer.println("1 2 # TRAIT");
				writer.println("0.99500 0.00500");
				writer.println("1");
				writer.println("0.03 0.80 0.80");
				line = reader.readLine().split("[\\s]+");
				writer.println(line[0]+" "+line[1]+" "+line[2]+line[3]);
				writer.println(reader.readLine());
				writer.println("0 0");
				writer.println("10.0");
				writer.println("1 0.1 0.45");
				writer.close();
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+"map"+chrome+".dat"+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+"map"+chrome+".dat"+"\"");
			System.exit(2);
		}

		writer = new PrintWriter(new FileWriter("run_pdt"+chrome));
		writer.println("#/bin/sh");
		writer.println();
		for (int i = 1; i<=numMarkers; i++) {
			writer.println("pdt "+i+".dat "+i+".ped "+i+".out > /dev/null");
		}
		writer.close();
		Files.chmod("run_pdt"+chrome);

		try {
			reader = new BufferedReader(new FileReader("re_chrom"+chrome+".pre"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				HashVec.addIfAbsent(line[0], v);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+"re_chrom"+chrome+".pre"+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+"re_chrom"+chrome+".pre"+"\"");
			System.exit(2);
		}

		numFams = v.size();
		numFams = numFams<=maxFams?numFams:maxFams;
		if (numMarkers<markerIncrement) {
			splitThis("re_chrom"+chrome+".pre", 0, numMarkers, numFams, pedinfo);
		} else {
			for (int count = 0; count<numMarkers; count += markerIncrement) {
				splitThis("re_chrom"+chrome+".pre", count, count+markerIncrement<=numMarkers?count+markerIncrement:numMarkers, numFams, pedinfo);
			}
		}
	}

	public static void splitThis(String prefile, int startAt, int stopAt, int numFams, Vector<String[]> pedinfo) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		Vector<String[]> data = new Vector<String[]>();
		String[] dataline;
		String line[];

		try {
			reader = new BufferedReader(new FileReader(prefile));

			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				dataline = new String[(stopAt-startAt)*2];
				for (int i = 0; i<stopAt-startAt; i++) {
					dataline[i*2+0] = line[(startAt+i)*2+6+0];
					dataline[i*2+1] = line[(startAt+i)*2+6+1];
				}
				data.add(dataline);
			}
			reader.close();
			if (data.size()!=pedinfo.size()) {
				System.err.println("Error - pedfile template does not have the same number of lines as the prefile");
				System.exit(1);
			}

			for (int i = 0; i<stopAt-startAt; i++) {
				writer = new PrintWriter(new FileWriter((startAt+i+1)+".ped"));
//				int count = 0;
				for (int j = 0; j<data.size(); j++) {
					dataline = data.elementAt(j);
					writer.println(Array.toStr(pedinfo.elementAt(j))+"\t"+dataline[i*2+0]+"\t"+dataline[i*2+1]);
//					count++;
				}
				writer.close();
			}
		} catch (IOException ioe) {
			System.err.println("Error - IOException in splitting file: "+prefile);
		}
	}

	public static void sumPDT(String dir, int start, int stop) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String temp;
		Vector<String> markerNames = new Vector<String>();
		String chrome;
		int numMarkers = -1;
		String[] dists = null;
		double posi;
		String sum, avg;

		try {
			writer = new PrintWriter(new FileWriter(dir+"pdt_summary.xls"));
			writer.println("marker\tchr\tposition\tsumPDT\tavePDT");

			for (int i = start; i<=stop; i++) {
				chrome = ext.chrome(i);
				try {
					reader = new BufferedReader(new FileReader("map"+chrome+".dat"));
					numMarkers = Integer.parseInt(reader.readLine().trim().split("[\\s]+")[0])-1;
					for (int j = 0; j<6; j++) {
						reader.readLine();
					}
					for (int j = 1; j<=numMarkers; j++) {
						markerNames.add(reader.readLine().split("[\\s]+")[3]);
						reader.readLine();
					}
					reader.readLine();
					dists = reader.readLine().split("[\\s]+");
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+"map"+chrome+".dat"+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+"map"+chrome+".dat"+"\"");
					System.exit(2);
				}

				posi = 0;
				for (int j = 1; j<=numMarkers; j++) {
					avg = sum = ".";
					try {
						reader = Files.getReader("chrom"+chrome+"/"+j+".out", dir);
						while (reader.ready()) {
							temp = reader.readLine();
							if (temp.indexOf("GLOBAL_SUM_PDT")>0) {
								sum = temp.split("[\\s]+")[3];
							} else if (temp.indexOf("GLOBAL_AVE_PDT")>0) {
								avg = temp.split("[\\s]+")[3];
							}
						}
						if (avg.equals(".")||sum.equals(".")) {
							System.err.println("Missing data for marker "+markerNames.elementAt(j-1)+" ("+j+".out)");
						}
						reader.close();
						posi += Double.parseDouble(j!=1?dists[j-1]:"0");
						writer.println(markerNames.elementAt(j-1)+"\t"+i+"\t"+posi+"\t"+sum+"\t"+avg);

						reader.close();
					} catch (Exception e) {
						System.err.println(e.getMessage());
						e.printStackTrace();
					}
				}
			}

			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to file");
			e.printStackTrace();
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		int chr = 1, start = 1, stop = 1;
		int markerIncrement = 49;
		String template = "template.ped";
		String dir = DEFAULT_DIR;
		boolean summarize = true;

		String usage = "\n"+"pdt requires 0-1 arguments\n"+"   (1) chromosome (i.e. chr="+chr+" (default is all chromosomes))\n"+"   (2) marker increment (i.e. inc="+markerIncrement+" (default))\n"+"   (3) template for additional .ped info (i.e. ped="+template+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("chr=")) {
				chr = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("inc=")) {
				markerIncrement = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("ped=")) {
				template = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (args.length==0) {
			System.err.println("No chromosome specified; doing chromosome "+chr);
		}
		try {
			if (summarize) {
				sumPDT(dir, start, stop);
			} else {
				createFiles(chr, markerIncrement, template);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
