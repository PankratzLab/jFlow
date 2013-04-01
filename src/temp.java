import common.Internat;

import java.io.*;
//import java.util.*;
import common.*;

public class temp {
	public static void parseAll(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String temp;
		String root;
		String annotation;
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename)+"_parsed.out"));
			root = "";
			annotation = "nada";
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.startsWith("rootLink=")) {
					root = temp.substring(9);
				} else if (temp.startsWith("<tr>")) {
					annotation = "";
					reader.readLine();
					reader.readLine();
					temp = reader.readLine().trim();
					temp = temp.substring(4, temp.length()-5);
					annotation += temp;
					reader.readLine();
					temp = reader.readLine().trim();
					temp = temp.substring(4, temp.length()-5);
					if (temp.length() > 0) {
						annotation += " ("+temp+")";
					}
				} else if (temp.indexOf(".mp4'>") > 0) {
					temp = temp.substring(0, temp.indexOf(".mp4'>")+4);
					temp = temp.substring(temp.lastIndexOf("href='")+6);
					writer.println(ext.link(root, temp)+"\t"+ext.replaceWithLinuxSafeCharacters(annotation+" "+ext.removeDirectoryInfo(temp), false));
				}
				
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}

	public static void downloadAll(String filename, String dir) {
		String[][] files;
		
		files = HashVec.loadFileToStringMatrix(filename, false, new int[] {0,1}, "\t", false, 1000, false);
		
		for (int i = 0; i < files.length; i++) {
			Internat.downloadFile(files[i][0], dir+files[i][1]);
		}

	}
	
	public static void parseMerc(String dir) {
		BufferedReader reader;
		PrintWriter writer;
		String[] files, data;
		String temp;

		files = Files.list(dir, ".htm", false);
		try {
			writer = new PrintWriter(new FileWriter(dir+"all.xln"));
			for (int i = 0; i < files.length; i++) {
				try {
					reader = new BufferedReader(new FileReader(dir+files[i]));
					data = new String[7];
					data[0] = ext.rootOf(files[i]);
					while (reader.ready()) {
						temp = reader.readLine();
						if (temp.contains("Valor:")) {
							data[1] = temp.trim().split("[\\s]+")[1];
						}
						if (temp.contains("Spirit:")) {
							data[2] = temp.trim().split("[\\s]+")[1];
						}
						if (temp.contains("Element:")) {
							data[3] = temp.trim().split("[\\s]+")[1];
						}
						if (temp.contains("Troop:")) {
							data[4] = temp.trim().split("\\>")[2];
							data[4] = data[4].substring(0, data[4].indexOf("<")).trim();
						}
						if (temp.contains("B.M.:")) {
							data[5] = temp.trim().split("\\>")[2];
							data[5] = data[5].substring(0, data[5].indexOf("<")).trim();
						}
						if (temp.contains("Area:")) {
							data[6] = temp.trim().split("\\>")[2];
							data[6] = data[6].substring(0, data[6].indexOf("<")).trim();
						}
					}
					writer.println(Array.toStr(data));
					writer.flush();
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + dir+files[i] + "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + dir+files[i] + "\"");
					System.exit(2);
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + dir+"all.xln");
			e.printStackTrace();
		}
	}
	
	public static void downTV() {
		int inarow = 0;
		int count = 3202;
		
		while (inarow < 100) {
			Internat.downloadFile("http://www.tv.com/buffy-the-vampire-slayer/show/"+count+"/summary.html", "D:/tv/"+count+".html");
			if (new File("D:/tv/"+count+".html").length() == 0) {
				inarow++;
			} else {
				inarow = 0;
			}
			count++;
		}
	}

	public static void parseTV() {
		BufferedReader reader;
		PrintWriter writer;
		String temp, filename, rating, votes, title, start, end;		
		
		try {
			writer = new PrintWriter(new FileWriter("d:/tv.xln"));
			writer.println("Index\tTitle\tRating\tVotes\tStartDate\tEndDate\tYearsOn\tLink");
			for (int i = 10; i < 1000; i++) {
				filename = "D:/tv/"+i+".html";
				title = "??";
				rating = votes = start = end = ".";
				if (new File(filename).exists() && new File(filename).length() != 0) {
					try {
						reader = new BufferedReader(new FileReader(filename));
						while (reader.ready()) {
							temp = reader.readLine();
							if (temp.contains("ratingValue")) {
								try {
									temp = temp.substring(temp.indexOf(">")+1);
									rating = temp.substring(0, temp.indexOf("<"));
								} catch (Exception e) {
									rating = "X";
								}
							} else if (temp.contains("ratingCount")) {
								try {
									temp = temp.substring(temp.indexOf(">")+1);
									votes = temp.substring(0, temp.indexOf(" votes"));
								} catch (Exception e) {
									votes = "X";
								}
							} else if (temp.contains("og:title")) {
								try {
									temp = temp.substring(temp.indexOf("content=\"")+9);
									title = temp.substring(0, temp.indexOf("\""));
								} catch (Exception e) {
									title = "X";
								}
							} else if (temp.contains("startDate")) {
								try {
									temp = temp.substring(temp.indexOf("datetime=\"")+10);
									start = temp.substring(0, temp.indexOf("\""));
								} catch (Exception e) {
									start = "X";
								}
							} else if (temp.contains("endDate")) {
								try {
									temp = temp.substring(temp.indexOf("datetime=\"")+10);
									end = temp.substring(0, temp.indexOf("\""));
								} catch (Exception e) {
									end = "X";
								}
							}
						}
						reader.close();
						writer.println(i+"\t"+title+"\t"+rating+"\t"+votes+"\t"+start+"\t"+end+"\t"+(start.equals(".")?".":(end.equals(".")?2011:Integer.parseInt(end.substring(0,end.indexOf("-")))) - Integer.parseInt(start.substring(0,start.indexOf("-"))))+"\t=HYPERLINK(\"D:/tv/"+i+".html\", \""+i+".html\")");
					} catch (FileNotFoundException fnfe) {
						System.err.println("Error: file \"" + filename + "\" not found in current directory");
						System.exit(1);
					} catch (IOException ioe) {
						System.err.println("Error reading file \"" + filename + "\"");
						System.exit(2);
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + "d:/tv.xln");
			e.printStackTrace();
		}
		
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "source.txt";

		String usage = "\n"+
		".temp requires 0-1 arguments\n"+
		"   (1) filename (i.e. file=" + filename + " (default))\n"+
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
//			Internat.downloadFile("http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2708794/pdf/JPATH175000054.pdf", "D:/SH3GL2.pdf");
//			parseAll(filename);
//			downloadAll(ext.rootOf(filename)+"_parsed.out", "C:/Ezgi/");
//			downloadAll("catchup.txt", "C:/Ezgi/");
//			parseMerc("C:\\Users\\npankrat\\Downloads\\Batheo\\mercs\\");

//			downTV();
//			parseTV();
			
//			String temp = "<td>Troop: <a href=\"/wiki/Category:Peltast\" title=\"Category:Peltast\">Peltast</a>;";
//			String[] line = temp.trim().split("\\>");
//			for (int i = 0; i < line.length; i++) {
//				System.out.println("'"+line[i]+"'");
//			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
