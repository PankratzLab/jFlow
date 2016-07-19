package one;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

public class ARIC_Manuscript {

	public static void main(String[] args) {
		Scanner reader;
		String line;
		String[] items, journalDetails;
		String pmid, title, journal, year;
		byte counter;
		BufferedWriter writer;
		
		try {
			counter = 0;
//			reader = new Scanner(new File("N:/statgen/ARIC_manuscripts/source/ansi.txt"));
			reader = new Scanner(new File("N:/statgen/ARIC_manuscripts/source/page0.txt"));
			writer = new BufferedWriter(new FileWriter("N:/statgen/ARIC_manuscripts/source/paperList_1.txt"));
			while(reader.hasNext()) {
				line = reader.nextLine();
//				if (line.startsWith("</div></form><div class=\"biblio-export\">Export 972 results")) {
					System.out.println(line);
//				}
				if (line.contains("class=\"biblio-title\"")) {
					counter ++;
					items = line.split("class=\"biblio-title\"");
					for (int i = 1; i < items.length; i++) {
						title = items[i].substring(items[i].indexOf(">") + 1, items[i].indexOf("</span>")).trim();
						if (title.endsWith(".")) {
							title = title.substring(0, title.length() - 1);
						}
						journal = items[i].split("</span></a>")[1].trim();
						journal = journal.substring(0, journal.indexOf("<span"));
						if (journal.startsWith(".")) {
							journal = journal.substring(1);
						}
						journalDetails = journal.split("[.]");
						journal = journalDetails[0].trim();
						year = journalDetails[1].split(";")[0].trim();
						pmid = items[i].substring(items[i].indexOf("href=\"http://www.ncbi.nlm.nih.gov/pubmed/") + 41, items[i].indexOf("?dopt=Abstract\""));
						writer.write(counter + "\t" + pmid + "\t" + journal + "\t" + year + "\t" + title + "\n");
					}
				}
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
