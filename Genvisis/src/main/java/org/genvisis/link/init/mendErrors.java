package link.init;

import java.io.*;
import java.util.*;

import common.*;

public class mendErrors {
	public static boolean COMPLETE = false;

	public mendErrors(int start, int stop) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer;
		StringTokenizer st;
		String temp, chrome, fam, data, marker, id, win; // what I need
		Vector<String> listOfErrors = new Vector<String>();

		for (int chromosome = start; chromosome<=stop; chromosome++) {
			try {
				System.out.print("Mending errors in chromosome "+chromosome);
				temp = "";
				listOfErrors.removeAllElements();
				chrome = (chromosome<10)?"0"+chromosome:""+chromosome;
				writer = new PrintWriter(new FileWriter("logfile of errors.out", true));

				if ((new File("marker.database")).exists()) {
					if (!new File("log.prn").exists()) {
						System.err.println("Error - marker.database was in the directory, so I assume we're in a Windows environment and can't run pedcheck.");
						System.err.println("        But could not find "+"log.prn"+" in current directory, which is what is usually used instead.");
						System.exit(2);
					}
					reader = new BufferedReader(new FileReader("log.prn"));
				} else {
					if ((new File("error.dat")).exists()) {
						(new File("error.dat")).delete();
					}
					Process process = null;
					Runtime runtime = Runtime.getRuntime();
					// process = runtime.exec("pedcheck -p chrom"+chrome+".pre
					// "+((chromosome == 23)?" -n names -x":"-d
					// map"+chrome+".dat ")+" -4");
					process = runtime.exec("pedcheck "+((chromosome==23)?" -x":"")+" -p chrom"+chrome+".pre -d map"+chrome+".dat "+" -4");

					reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
				}

				while (!reader.readLine().startsWith(" CHECKING LEVEL 3  ERRORS. "));
				reader.readLine();
				reader.readLine();
				reader.readLine();
				reader.readLine();

				boolean done = false;
				temp = reader.readLine();
				if (temp.equals("")) {
					System.out.println("<<no errors>>done");
					continue;
				}

				while (!done) {
					st = new StringTokenizer(temp, " ."); // new Fam
					st.nextToken();
					st.nextToken();
					fam = st.nextToken();
					reader.readLine();
					temp = reader.readLine();

					while (!temp.startsWith(" ##### Pedigree")) {
						st = new StringTokenizer(temp); // new marker
						if (st.countTokens()==0) {
							done = true;
							break;
						}
						st.nextToken();
						marker = st.nextToken();
						st.nextToken();
						temp = st.nextToken();
						st = new StringTokenizer(temp, ".");
						temp = st.nextToken();
						win = fam+"\t"+temp;
						marker += " ("+temp+")";
						if (marker.length()<10) {
							marker += "\t";
						}

						do {
							temp = reader.readLine();
						} while (!temp.startsWith(" Person")&&!temp.startsWith(" LEVEL 3 FAILED"));
						if (temp.startsWith(" LEVEL 3 FAILED")) {
							writer.println(fam+"\t"+marker+"\t"+"failed");
							System.out.print("X");
							reader.readLine(); // always a blank
							temp = reader.readLine();
							if (temp.startsWith(" ##### Pedigree")) { // if
								// new
								// pedigree
								break; // then cycle
							}
							temp = reader.readLine();
							if (temp.startsWith(" Name:")) { // if new
								// marker, same
								// family
								continue; // then continue
							}
							if (temp.equals("")) { // if end of list
								done = true;
								break; // then quit
							}
							System.err.println("Error: something failed, but didn't want to unfail for some reason.");
						}

						boolean chosenOne = false;
						boolean found = false;
						while (!temp.equals("")&&!temp.startsWith(" ##### Pedigree")) {
							st = new StringTokenizer(temp, " :\t/");
							st.nextToken();
							temp = st.nextToken();
							id = temp;
							data = temp+"\t"+st.nextToken()+"/"+st.nextToken();

							temp = reader.readLine();
							do {
								st = new StringTokenizer(temp, " /\t");
								st.nextToken();
								st.nextToken();
								if (st.nextToken().equals("1.000")) {
									if (!found) {
										chosenOne = true;
									}
								}
								temp = reader.readLine();
							} while (!temp.startsWith(" Person")&&!temp.equals("")&&!temp.startsWith(" ##### Pedigree"));

							if (chosenOne&&!found) {
								System.out.print(".");
								win += "\t"+id;
								listOfErrors.add(win);
								writer.println(fam+"\t"+marker+"\t"+data+"\t\tBecause Pedcheck said so");
								found = true;
							}
						}
						if (!temp.startsWith(" ##### Pedigree")) {
							temp = reader.readLine();
						}
					}
				}

				reader.close();
				writer.close();

				if (!new File("pedcheck.err").exists()) {
					System.err.println("Error - could not find "+"pedcheck.err"+" in current directory");
					System.exit(2);
				}
				reader = new BufferedReader(new FileReader("pedcheck.err"));
				writer = new PrintWriter(new FileWriter("review chrom"+chrome+" errors.out", true), true);

				while (reader.ready()&&!reader.readLine().startsWith(" *********** LEVEL 1 ERRORS "));

				if (!reader.ready()) {
					writer.println("There was probably a level 3 error; can't say for sure");
				} else {
					reader.readLine();
					temp = reader.readLine();
					while (!temp.equals("")) {
						writer.println(temp);
						reader.readLine();
						writer.println(reader.readLine());
						reader.readLine();
						temp = reader.readLine();
						while (!temp.equals("")) {
							writer.println(temp);
							temp = reader.readLine();
						}
						writer.println();
						temp = reader.readLine();
					}
				}
				writer.println();
				writer.println();
				writer.println();

				reader.close();
				writer.close();
				String bakFilename = Files.getBakFilename("chromosome"+chromosome+".dat", super.getClass().getName());
				(new File("chromosome"+chromosome+".dat")).renameTo((new File(bakFilename)));

				reader = new BufferedReader(new FileReader(bakFilename));
				writer = new PrintWriter(new FileWriter("chromosome"+chromosome+".dat"));

				temp = reader.readLine();
				if (temp.equals("")) {
					temp = "placeholder line"; // for older file versions
				}
				String crap, lastFam = "";
				fam = "";
				for (int i = 0; i<listOfErrors.size(); i++) {
					st = new StringTokenizer(listOfErrors.elementAt(i));
					lastFam = fam;
					fam = st.nextToken();
					marker = st.nextToken();
					data = st.nextToken();

					if (fam.equals(lastFam)) {
						writer.println(temp);
						while (reader.ready()) {
							writer.println(reader.readLine());
						}
						reader.close();
						writer.close();
						(new File("chromosome"+chromosome+".dat")).renameTo((new File("temp")));
						reader = new BufferedReader(new FileReader("temp"));
						writer = new PrintWriter(new FileWriter("chromosome"+chromosome+".dat"));
						temp = reader.readLine();
					}

					st = new StringTokenizer(temp);
					crap = st.nextToken();
					while (!st.nextToken().equals(fam)) {
						writer.println(temp);
						temp = reader.readLine();
						st = new StringTokenizer(temp);
						crap = st.nextToken();
					}
					while (!st.nextToken().equals(data)) {
						writer.println(temp);
						temp = reader.readLine();
						st = new StringTokenizer(temp);
						crap = st.nextToken();
						st.nextToken();
					}
					writer.print(crap+"\t"+fam+"\t"+data);
					for (int j = 0; j<2*(Integer.valueOf(marker).intValue()-1); j++) {
						writer.print("\t"+st.nextToken());
					}
					writer.print("\t0\t0");
					st.nextToken();
					st.nextToken();
					while (st.hasMoreTokens()) {
						writer.print("\t"+st.nextToken());
					}
					writer.println();
					if (reader.ready()) {
						temp = reader.readLine();
					} else {
						temp = null;
					}
				}
				if (temp!=null) {
					writer.println(temp);
				}
				while (reader.ready()) {
					temp = reader.readLine();
					if (!temp.equals("null")) {
						writer.println(temp);
					}
				}

				reader.close();
				writer.close();
				System.out.println("done");
				(new File("pedcheck.err")).delete();
				(new File("temp")).delete();
			} catch (Exception e) {
				System.err.println();
				if (chromosome==23) {
					System.err.println("Got an error, probably need to check for males heterozygous for an X-linked marker.");
				} else {
					System.err.println("Got an error running chromosome "+chromosome+", probably couldn't run pedcheck. Have the map and pre files been made?");
					System.err.println("If that's not it, try deleting the files \'pedcheck.err\' and \'temp\'");
					System.err.println("If that's not it, make sure the first line of the .dat file is either blank or has more than one word.");
				}
				System.err.println();
			}
		}
	}

	public static void main(String[] args) throws IOException {
		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-?")||args[i].equals("-help")||args[i].equals("-h")) {
				System.out.println("Expecting 0-2 arguments (chromosome number to start [and stop], default all).");
				System.exit(1);
			}
		}

		// new mendErrors(23, 23);
		// new mendErrors(5, 5);
		// System.exit(1);

		try {
			if (args.length==0) {
				new mendErrors(1, 23);
			}
			if (args.length==1) {
				new mendErrors(Integer.valueOf(args[0]).intValue(), Integer.valueOf(args[0]).intValue());
			}
			if (args.length==2) {
				new mendErrors(Integer.valueOf(args[0]).intValue(), Integer.valueOf(args[1]).intValue());
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
