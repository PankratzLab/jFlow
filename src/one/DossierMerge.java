package one;

import java.io.*;
import java.util.*;

import common.*;

public class DossierMerge {
	public static void run(String dir, boolean rtfOutput) throws IOException {
		RealTextFormatWriter writer;
		String[] keys, bits, line;
		Hashtable<String, String> citationHash, timesCitedHash, impactHash, rolesHash, abbreviationsHash;
		HashSet<String> authorsToBoldHash, jointFirstAuthorsHash;
		String pmid, citation, journal, timesCited, impactFactor, role;
		Logger log;
		boolean jointFirstAuthor;
		
		writer = new RealTextFormatWriter(dir+"formattedDossier"+(rtfOutput?".rtf":".out"), rtfOutput);
		log = new Logger(dir+"formattedDossier.log");
		
//		citationHash = HashVec.loadFileToHashString(dir+"pubmed_result.utf8.csv", new int[] {9}, new int[] {2, 0, 3}, true, "\t", true, false, false);
		citationHash = HashVec.loadFileToHashString(dir+"personal bibliography.txt", new int[] {0}, new int[] {1, 2, 3}, false, "\t", true, false, false);
		abbreviationsHash = HashVec.loadFileToHashString(dir+"Abbreviations.txt", new int[] {0}, new int[] {1}, false, null, false, false, false);
		
		impactHash = HashVec.loadFileToHashString(dir+"impactFactors.dat", new int[] {0}, new int[] {2,1}, false, " in ", false, false, false);
		timesCitedHash = HashVec.loadFileToHashString(dir+"citations.dat", new int[] {0}, new int[] {1,2}, false, "\t", false, false, false);
		rolesHash = HashVec.loadFileToHashString(dir+"roles.dat", new int[] {0}, new int[] {1}, false, null, false, false, false);
		

		authorsToBoldHash = HashVec.loadToHashSet(HashVec.loadFileToStringArray(dir+"wordsToBold.dat", false, null, false));
		jointFirstAuthorsHash = HashVec.loadFileToHashSet(dir+"jointFirstAuthors.dat", false);

		new File(dir+"journalsInDatabase.xln").delete();
		
		keys = HashVec.getKeys(citationHash, true, true);
		for (int i = 0; i < keys.length; i++) {
			pmid = keys[keys.length-i-1];
			citation = citationHash.remove(pmid);
			impactFactor = impactHash.remove(pmid);
			timesCited = timesCitedHash.remove(pmid);
			role = rolesHash.remove(pmid);
			jointFirstAuthor = jointFirstAuthorsHash.contains(pmid);
			
			if (role == null || !role.equals("None")) {
				bits = citation.split("\t", -1);
				
				line = bits[0].substring(0, bits[0].length()-1) .split(",");
				for (int j = 0; j < line.length; j++) {
					String author = line[j].trim();
					if (authorsToBoldHash.contains(author)) {
						line[j] = "<b>"+line[j]+"</>";
					}
				}
				if (jointFirstAuthor) {
					line[0] += "<super>†</>";
					line[1] += "<super>†</>";
				}
				bits[0] = Array.toStr(line, ",")+".";
				
				if (bits[2].contains(" doi:")) {
					bits[2] = bits[2].substring(0, bits[2].indexOf(" doi:"))+(bits[2].contains("[Epub ahead of print]")?" [Epub ahead of print]":"");
				}

				journal = bits[2].substring(0, bits[2].indexOf("."));
				if (abbreviationsHash.containsKey(journal)) {
					ext.appendToFile(journal+"\t"+abbreviationsHash.get(journal), dir+"journalsInDatabase.xln");
					journal = abbreviationsHash.get(journal);
				} else if (!abbreviationsHash.contains(journal)) {
					ext.appendToFile(journal+"\t.", dir+"journalsInDaabase.xln");
					log.reportError("Error journal '"+journal+"' was not found as an abbreviation or as a full journal name");
				}
				bits[2] = "<u>"+journal+"</>"+bits[2].substring(bits[2].indexOf("."));
				
				writer.println((i+1)+"."+(rtfOutput?"\\tab ":"\t")+Array.toStr(bits, " "));
			
				if (jointFirstAuthor) {
					writer.println((rtfOutput?"\\bullet \\tab ":"\t")+"<super>†</>Joint first author");
				}
				
				if (impactFactor != null) {
					writer.println((rtfOutput?"\\bullet \\tab ":"\t")+"Journal Impact Factor: "+impactFactor);
				} else {
					writer.println((rtfOutput?"\\bullet \\tab ":"\t")+"Journal Impact Factor: Not available");
					log.reportError("Error - no impact factor available for "+pmid+" ("+citation+")");
				}
				
				if (timesCited != null) {
					line = timesCited.split("[\\s]+");
					writer.println((rtfOutput?"\\bullet \\tab ":"\t")+"Times Cited: "+line[0]+" (Scopus) / "+line[1]+" (Google Scholar)");
				} else {
					writer.println((rtfOutput?"\\bullet \\tab ":"\t")+"Times Cited: Not available");
					log.reportError("Error - no times cited available for "+pmid+" ("+citation+")");
				}
				
				if (role != null) {
					writer.println((rtfOutput?"\\bullet \\tab ":"\t")+"Role: "+role);
				} else {
					writer.println((rtfOutput?"\\bullet \\tab ":"\t")+"Role: ");
					log.reportError("Error - no role available for "+pmid+" ("+citation+")");
				}
				
				writer.newParagraph();
			}
		}

		writer.close();
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "D:/umn/promotion/citations/";
		boolean rtf = true;
		
		String usage =  "\n" + 
						"one.AuthorCorral requires 2+ arguments\n" + 
						"   (1) working directory (i.e. dir=" + dir + " (default))\n" + "" + 
						"   (2) output in RTF format (used for superscripts) (i.e. rtf=" + rtf + " (default))\n";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("rtf=")) {
				rtf = Boolean.valueOf(args[i].split("=")[1]);
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
			run(dir, rtf);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
