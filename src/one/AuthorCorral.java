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
	
	private static String rtfFormatted(String str, boolean authStr) {
		StringBuilder sb = new StringBuilder();
		String PREFIX = "{\\rtlch\\fcs1 \\af31507 \\ltrch\\fcs0 ";
		sb.append(PREFIX);
		for (int i = 0; i < str.length(); i++) {
			sb.append(str.charAt(i));
		}
		sb.append("}");
		String ret = sb.toString();
		if (authStr) {
			ret = sb.toString().replaceAll("#", "}{\\\\rtlch\\\\fcs1 \\\\af0\\\\afs24 \\\\ltrch\\\\fcs0 \\\\f0\\\\fs24\\\\super\\\\insrsid6508612 \\\\hich\\\\af0\\\\dbch\\\\af31505\\\\loch\\\\f0 ");
			ret = ret.replaceAll("_", "}{\\\\rtlch\\\\fcs1 \\\\af31507 \\\\ltrch\\\\fcs0 ");
		}
		return ret;
	}
	
	static final String PARAGRAH = "{\\rtlch\\fcs1 \\af0\\afs24 \\ltrch\\fcs0 \\line }";
	
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
	public static void run(String inFile, String outFile, boolean rtfOutput, boolean printErr) throws IOException {
		BufferedReader inReader = Files.getReader(inFile, false, true, true);
		String[] header = inReader.readLine().split("\t");
		
		ArrayList<String> authorNamesOrder = new ArrayList<String>();
		ArrayList<String> contribOrder = new ArrayList<String>();
		HashMap<String, String[]> fullnameToParts = new HashMap<String, String[]>();
		HashMap<String, ArrayList<String[]>> namesToDepts = new HashMap<String, ArrayList<String[]>>();
		HashMap<String, ArrayList<String>> contribsToNames = new HashMap<String, ArrayList<String>>();
		
		ArrayList<String> errors = new ArrayList<String>();
		
		for (int i = 7; i < header.length; i++) {
			contribOrder.add(header[i].trim());
			contribsToNames.put(header[i].trim(), new ArrayList<String>());
		}
		
		String temp = null;
		while ((temp = inReader.readLine()) != null) {
			String[] line = temp.split("\t");
			if (line.length < 3) {
				String err = "Insufficient data for <" + temp.trim() + ">";
				System.err.println("Error - " + err);
				errors.add(err);
			} else {
				String fullname = line[0].trim() + " " + (!"".equals(line[1].trim()) ? (line[1].trim().charAt(0) + ". ") : "") + line[2].trim();
				
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
						if (contribsToNames.get(header[i]).indexOf(fullname) == -1) {
							contribsToNames.get(header[i]).add(fullname);//getInitials(line[0].trim(), line[1].trim(), line[2].trim()));
						}
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
				for (int d = 0; d < myDepts.size(); d++) {
					if (d == 0 && rtfOutput) {
						authorString.append("#");
					}
					StringBuilder deptStr = new StringBuilder();
					int nonEmptyCnt = 0;
					if (!"".equals(myDepts.get(d)[0])) {
						deptStr.append(myDepts.get(d)[0]);
						nonEmptyCnt++;
					}
					if (!"".equals(myDepts.get(d)[1])) {
						if (nonEmptyCnt > 0 && !"".equals(myDepts.get(d)[1])) {
							deptStr.append(", ");
						}
						deptStr.append(myDepts.get(d)[1]);
						nonEmptyCnt++;
					}
					if (!"".equals(myDepts.get(d)[2])) {
						if (nonEmptyCnt > 0 && !"".equals(myDepts.get(d)[2])) {
							deptStr.append(", ");
						}
						deptStr.append(myDepts.get(d)[2]);
						nonEmptyCnt++;
					}
					int index = deptOrder.indexOf(deptStr.toString());
					if (index == -1) {
						index = deptOrder.size();
						deptOrder.add(deptStr.toString()); 
					}
					authorString.append(index + 1);
					if (d == myDepts.size() - 1 && rtfOutput) {
						authorString.append("_");
					} else if (d >= 0 && d < myDepts.size() - 1) {
						authorString.append(",");
					}
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
		
		if (rtfOutput) {
			outWriter.println("{\\rtf1 ");
		}
		outWriter.println(rtfOutput ? rtfFormatted(authorString.toString(), true) : authorString.toString());
		outWriter.println(rtfOutput ? PARAGRAH : "");
		outWriter.println(rtfOutput ? PARAGRAH : "");
		for (int i = 0; i < deptOrder.size(); i++) {
			outWriter.println(rtfOutput ? rtfFormatted((i + 1) + " " + deptOrder.get(i), false) : (i + 1) + " " + deptOrder.get(i));
			if (rtfOutput) {
				outWriter.println(rtfOutput ? PARAGRAH : "");
			}
		}
		outWriter.println(rtfOutput ? PARAGRAH : "");
		outWriter.println(rtfOutput ? PARAGRAH : "");
		
		outWriter.println(rtfOutput ? rtfFormatted("Author Contributions", false) : "Author Contributions");
		outWriter.print(rtfOutput ? PARAGRAH : "");
		StringBuilder contribString = new StringBuilder();
		for (String contrib : contribOrder) {
			ArrayList<String> contribInitials = contribsToNames.get(contrib);
			for (int i = 0; i < contribInitials.size(); i++) {
				if (i == contribInitials.size() - 1) {
					contribString.append("and ");
				}
				String[] parts = fullnameToParts.get(contribInitials.get(i));
				contribString.append(getInitials(parts[0], parts[1], parts[2]));
				if (i < contribInitials.size() - 1) {
					contribString.append(", ");
				} else {
					contribString.append(" ");
				}
			}
			contribString.append(contrib).append(". ");
		}
		contribString.append("All authors were given the opportunity to comment and provide revisions to the manuscript text.");
		outWriter.println(rtfOutput ? rtfFormatted(contribString.toString(), false) : contribString.toString());

		outWriter.println(rtfOutput ? PARAGRAH : "");
		outWriter.println(rtfOutput ? PARAGRAH : "");
		
		outWriter.println(rtfOutput ? rtfFormatted("Email list", false) : "Email list");
		outWriter.print(rtfOutput ? PARAGRAH : "");
		StringBuilder emailString = new StringBuilder();
		HashSet<String> emailSet = new HashSet<String>();
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
				String email = fullnameToParts.get(authorNamesOrder.get(i))[0] +" "+ fullnameToParts.get(authorNamesOrder.get(i))[2] +" <"+	myDepts.get(j)[3] +">";
				if (!emailSet.contains(email)) {
					if (outerCount > 0 || innerCount > 0) {
						emailString.append("; ");
					}
				
					emailString.append(email);
					innerCount++;
					emailSet.add(email);
				}
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
		outWriter.println(rtfOutput ? rtfFormatted(emailString.toString(), false) : emailString.toString());
		outWriter.println(rtfOutput ? PARAGRAH : "");
		outWriter.println(rtfOutput ? PARAGRAH : "");
		
		outWriter.println(rtfOutput ? rtfFormatted("No Email Found", false) : "No Email Found");
		outWriter.print(rtfOutput ? PARAGRAH : "");
		StringBuilder noEmailString = new StringBuilder();
		int cnt = 0;
		for (String name : noEmailSet) {
			if (cnt > 0) {
				noEmailString.append("; ");
			}
			noEmailString.append(name);
			cnt++;
		}
		outWriter.println(rtfOutput ? rtfFormatted(noEmailString.toString(), false) : noEmailString.toString());

		outWriter.println(rtfOutput ? PARAGRAH : "");
		outWriter.println(rtfOutput ? PARAGRAH : "");
		
		if (errors.size() > 0) {
			outWriter.println(rtfOutput ? rtfFormatted("Transcription Errors", false) : "Transcription Errors");
			outWriter.print(rtfOutput ? PARAGRAH : "");
		}
		for (String error : errors) {
			outWriter.println(rtfOutput ? rtfFormatted(error, false) : error);
			outWriter.print(rtfOutput ? PARAGRAH : "");
		}
		if (rtfOutput) {
			outWriter.println("}");
		}
		
		outWriter.flush();
		outWriter.close();
	}
	
	private static String getInitials(String first, String midTemp, String last) {
		StringBuilder init = new StringBuilder();
		try {
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
		} catch (Exception e) {
			System.err.println("Error - problem getting initials from [" + first + "; " + midTemp + "; " + last + "]");
		}
		return init.toString().trim();
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String inFile = "N:/statgen/authors/input2h.xln";
		String outFile = "N:/statgen/authors/2nd_submission_authorship_final.rtf";
		boolean rtf = true;
		boolean err = true;
		
		String usage =  "\n" + 
						"one.AuthorCorral requires 2+ arguments\n" + 
						"   (1) input filename (i.e. file=" + inFile + " (default))\n" + "" + 
						"   (2) output filename (i.e. out=" + outFile + " (default))\n" + "" + 
						"   (3) OPTIONAL output in RTF format (used for superscripts) (i.e. rtf=" + rtf + " (default))\n" + 
						"   (4) OPTIONAL include errors in output (i.e. err=" + err + " (default))\n" + "";
		

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
			} else if (args[i].startsWith("rtf=")) {
				rtf = Boolean.valueOf(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("err=")) {
				err = Boolean.valueOf(args[i].split("=")[1]);
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
			run(inFile, outFile, rtf, err);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
