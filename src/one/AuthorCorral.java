package one;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import common.Files;
import common.ext;

public class AuthorCorral {
	
	/*
		Tab-Delimited Authorship file:
			[0] First
			[1] Middle
			[2] Last
			[3] Dept/Div/Inst
			[4] Institution
			[5] City, State, Country, Zip
			[6] E-mail
			[7+] contribution columns - title taken verbatim from column title
	*/
	public static void run(String inFile, String outFile) throws IOException {
		BufferedReader inReader = Files.getReader(inFile, false, true, true);
		String[] header = inReader.readLine().split("\t");
		
		ArrayList<String> authorNamesOrder = new ArrayList<String>();
		ArrayList<String> contribOrder = new ArrayList<String>();
		HashMap<String, String[]> fullnameToParts = new HashMap<String, String[]>();
		HashMap<String, ArrayList<String[]>> namesToDepts = new HashMap<String, ArrayList<String[]>>();
		HashMap<String, ArrayList<String>> contribsToInitials = new HashMap<String, ArrayList<String>>();
		
		ArrayList<String> errors = new ArrayList<String>();
		
		for (int i = 7; i < header.length; i++) {
			contribOrder.add(header[i].trim());
			contribsToInitials.put(header[i].trim(), new ArrayList<String>());
		}
		
		String temp = null;
		while ((temp = inReader.readLine()) != null) {
			String[] line = temp.split("\t");
			if (line.length < 3) {
				String err = "Insufficient data for <" + temp.trim() + ">";
				System.err.println("Error - " + err);
				errors.add(err);
			} else {
				String fullname = line[0].trim() + " " + (!"".equals(line[1].trim()) ? (line[1].trim() + " ") : "") + line[2].trim();
				
				int ord = authorNamesOrder.indexOf(fullname);
				if (ord == -1) {
					ord = authorNamesOrder.size();
					authorNamesOrder.add(fullname);
					fullnameToParts.put(fullname, new String[]{line[0].trim(), line[1].trim(), line[2].trim()});
				}
				
				ArrayList<String[]> myDepts = namesToDepts.get(fullname); 
				if (myDepts == null) {
					myDepts = new ArrayList<String[]>();
					namesToDepts.put(fullname, myDepts);
				}
				
				if (line.length < 7) {
					String err = "Insufficient data for <" + temp.trim() + ">";
					System.err.println("Error - " + err);
					errors.add(err);
				} else {
					if ("".equals(line[3].trim())) {
						String err = "[Dept/Div/Inst] is blank for <" + line[0].trim() + " " + line[1].trim() + " " + line[2].trim() + ">";
						System.err.println("Error - " + err);
						errors.add(err);
					}
					if ("".equals(line[4].trim())) {
						String err = "[Institution] is blank for <" + line[0].trim() + " " + line[1].trim() + " " + line[2].trim() + ">";
						System.err.println("Error - " + err);
						errors.add(err);
					}
					if ("".equals(line[5].trim())) {
						String err = "[City, State, Country, Zip] is blank for <" + line[0].trim() + " " + line[1].trim() + " " + line[2].trim() + ">";
						System.err.println("Error - " + err);
						errors.add(err);
					}
					if ("".equals(line[6].trim())) {
						String err = "[E-mail] is blank for <" + line[0].trim() + " " + line[1].trim() + " " + line[2].trim() + ">";
						System.err.println("Error - " + err);
						errors.add(err);
					}
					myDepts.add(new String[]{ext.removeQuotes(line[3].trim()).trim(), ext.removeQuotes(line[4].trim()).trim(), ext.removeQuotes(line[5].trim()).trim(), ext.removeQuotes(line[6].trim()).trim()});
				}
				
				boolean foundAtLeastOne = false;
				for (int i = 7; i < line.length; i++) {
					if (!"".equals(line[i].trim())) {
						contribsToInitials.get(header[i]).add(getInitials(line[0].trim(), line[1].trim(), line[2].trim()));
						foundAtLeastOne = true;
					}
				}
				if (!foundAtLeastOne) {
					String err = "[contributions] is blank for <" + line[0].trim() + " " + line[1].trim() + " " + line[2].trim() + ">";
					System.err.println("Error - " + err);
					errors.add(err);
				}
			}
		}
		inReader.close();
		
		PrintWriter outWriter = Files.getWriter(outFile);
		
		StringBuilder authorString = new StringBuilder();
		ArrayList<String> deptOrder = new ArrayList<String>(); 
		
		for (int i = 0; i < authorNamesOrder.size(); i++) {
			ArrayList<String[]> myDepts = namesToDepts.get(authorNamesOrder.get(i));
			if (myDepts.size() > 0 && !(myDepts.size() == 1 && myDepts.get(0)[0].equals("") && myDepts.get(0)[1].equals("") && myDepts.get(0)[2].equals(""))) {
				if (i == authorNamesOrder.size() - 1) {
					authorString.append("and ");
				}
				authorString.append(authorNamesOrder.get(i));
				for (String[] dept : myDepts) {
					StringBuilder deptStr = new StringBuilder();
					int nonEmptyCnt = 0;
					if (!"".equals(dept[0])) {
						deptStr.append(dept[0]);
						nonEmptyCnt++;
					}
					if (!"".equals(dept[1])) {
						if (nonEmptyCnt > 0 && !"".equals(dept[1])) {
							deptStr.append(", ");
						}
						deptStr.append(dept[1]);
						nonEmptyCnt++;
					}
					if (!"".equals(dept[2])) {
						if (nonEmptyCnt > 0 && !"".equals(dept[2])) {
							deptStr.append(", ");
						}
						deptStr.append(dept[2]);
						nonEmptyCnt++;
					}
					int index = deptOrder.indexOf(deptStr.toString());
					if (index == -1) {
						index = deptOrder.size();
						deptOrder.add(deptStr.toString()); 
					}
					authorString.append("'").append(index + 1);
				}
				if (i < authorNamesOrder.size() - 1) {
					authorString.append(", ");
				}
			} else {
				String err = "No affiliation data for <" + authorNamesOrder.get(i).trim() + ">";
				System.err.println("Error - " + err);
				errors.add(err);
			}
		}
		
		outWriter.println(authorString.toString());
		outWriter.println();
		for (int i = 0; i < deptOrder.size(); i++) {
			outWriter.println((i + 1) + " " + deptOrder.get(i));
		}
		outWriter.println();
		outWriter.println();
		
		outWriter.println("Author Contributions");
		StringBuilder contribString = new StringBuilder();
		for (String contrib : contribOrder) {
			ArrayList<String> contribInitials = contribsToInitials.get(contrib);
			for (int i = 0; i < contribInitials.size(); i++) {
				if (i == contribInitials.size() - 1) {
					contribString.append("and ");
				}
				contribString.append(contribInitials.get(i));
				if (i < contribInitials.size() - 1) {
					contribString.append(", ");
				} else {
					contribString.append(" ");
				}
			}
			contribString.append(contrib).append(". ");
		}
		contribString.append("All authors were given the opportunity to comment and provide revisions to the manuscript text.");
		outWriter.println(contribString.toString());

		outWriter.println();
		outWriter.println();
		
		outWriter.println("Email list");
		StringBuilder emailString = new StringBuilder();
		HashSet<String> noEmailSet = new HashSet<String>();
		int outerCount = 0;
		for (int i = 0; i < authorNamesOrder.size(); i++) {
			ArrayList<String[]> myDepts = namesToDepts.get(authorNamesOrder.get(i));
			int innerCount = 0;
			for (int j = 0; j < myDepts.size(); j++) {
				if (myDepts.get(j)[3].equals("")) {
					StringBuilder noEmailString = new StringBuilder();
					noEmailString.append(fullnameToParts.get(authorNamesOrder.get(i))[0])
									.append(" ")
									.append(fullnameToParts.get(authorNamesOrder.get(i))[2]);
					noEmailSet.add(noEmailString.toString());
					continue;
				}
				if (outerCount > 0 || innerCount > 0) {
					emailString.append("; ");
				}
				emailString.append(fullnameToParts.get(authorNamesOrder.get(i))[0])
							.append(" ")
							.append(fullnameToParts.get(authorNamesOrder.get(i))[2])
							.append(" <");
				emailString.append(myDepts.get(j)[3]);
				emailString.append(">");
				innerCount++;
			}
			if (innerCount > 0) {
				outerCount++;
			} else {
				StringBuilder noEmailString = new StringBuilder();
				noEmailString.append(fullnameToParts.get(authorNamesOrder.get(i))[0])
								.append(" ")
								.append(fullnameToParts.get(authorNamesOrder.get(i))[2]);
				noEmailSet.add(noEmailString.toString());
			}
		}
		outWriter.println(emailString.toString());
		outWriter.println();
		
		
		outWriter.println("No Email Found");
		StringBuilder noEmailString = new StringBuilder();
		int cnt = 0;
		for (String name : noEmailSet) {
			if (cnt > 0) {
				noEmailString.append("; ");
			}
			noEmailString.append(name);
			cnt++;
		}
		outWriter.println(noEmailString.toString());
		
		outWriter.println();
		outWriter.println();
		
		if (errors.size() > 0) {
			outWriter.println("Transcription Errors");
		}
		for (String error : errors) {
			outWriter.println(error);
		}
		
		outWriter.flush();
		outWriter.close();
	}
	
	private static String getInitials(String first, String midTemp, String last) {
		StringBuilder init = new StringBuilder();
		for (String s : first.split("[\\s]+")) {
			init.append(s.charAt(0));
		}
		boolean allUpper = true;
		String mid = midTemp.replaceAll("\\.", "");
		for (int i = 0; i < mid.length(); i++) {
			if (!Character.isUpperCase(mid.charAt(i))) {
				allUpper = false;
				break;
			}
		}
		if (allUpper) {
			init.append(mid);
		} else {
			for (String s : mid.split("[\\s]+")) {
				init.append(s.charAt(0));
			}
		}
		for (String s : last.split("[\\s]+")) {
			init.append(s.charAt(0));
		}
		return init.toString().trim();
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String inFile = "N:/statgen/authors/input2.xln";
		String outFile = "N:/statgen/authors/authorship.out";

		String usage =  "\n" + 
						"one.AuthorCorral requires 2 arguments\n" + 
						"   (1) input filename (i.e. file=" + inFile + " (default))\n" + "" + 
						"   (1) output filename (i.e. out=" + outFile + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				inFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outFile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			run(inFile, outFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
