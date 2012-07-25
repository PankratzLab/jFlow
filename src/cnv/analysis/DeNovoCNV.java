package cnv.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Vector;

import cnv.filesys.Project;

import common.StringVector;
import common.CmdLine;
import common.ext;

public class DeNovoCNV {

	private static Project proj;
	
	public static void loadData(String outfile) {
		BufferedReader reader;
		String[] line;
		PrintWriter writer;
//		HashVec trav;
//		String[] newID;
		StringVector fId, iId, faId, moId, dna;
		
//		trav = new HashVec();
		try {
//	        reader = Files.getReader(new FileReader(filename), proj.getJarStatus(), true, false);
			reader = new BufferedReader(new FileReader(outfile));
//			reader.readLine();
			//TODO detect the column number of FID IID FAID MOID
			fId = new StringVector();
			iId = new StringVector();
			faId = new StringVector();
			moId = new StringVector();
			dna = new StringVector();
            while (reader.ready()) {
            	line = reader.readLine().trim().split("\t",-1);
            	if (line[6]!=null && !line[6].equals(".")) {
            		fId.add(line[0]);
            		iId.add(line[0]+"\t"+line[1]);
            		faId.add(line[2]);
            		moId.add(line[3]);
            		dna.add(line[6]);
            	}
            }
            reader.close();

            writer = new PrintWriter(new FileWriter("D:/PedigreeTempOutput.txt"));
			writer.println("fId\tiId\tfaId\tmoId\tiDna\tfaDna\tmoDna");
			for (int i=0; i<iId.size(); i++) {
				if (iId.contains(fId.elementAt(i)+"\t"+ faId.elementAt(i)) && iId.contains(fId.elementAt(i)+"\t"+moId.elementAt(i))) {
					writer.println(iId.elementAt(i)
							+"\t"+faId.elementAt(i)
							+"\t"+moId.elementAt(i)
							+"\t"+dna.elementAt(i)
							+"\t"+dna.elementAt(iId.indexOf(fId.elementAt(i)+"\t"+faId.elementAt(i)))
							+"\t"+dna.elementAt(iId.indexOf(fId.elementAt(i)+"\t"+moId.elementAt(i))));
				}
			}
			writer.flush();
			writer.close();
        } catch (FileNotFoundException fnfe) {
        	System.err.println("Error: file \""+outfile+"\" not found in current directory");
        } catch (IOException ioe) {
            System.err.println("Error reading file \""+outfile+"\"");
        }
	}
	
	public static void runPennCnv(String infile) {
		BufferedReader reader;
		String[] line;
		
		if (!(new File("penn_output")).exists()) {
			new File("penn_output").mkdir();
		}
		PennCNV.populationBAF(proj);
		PennCNV.gcModel(proj, "D:/PennCNV_Related/GcModel/gc5Base_hg18.txt", "penn_output/gcmodel", 0);
		
		try {
			reader = new BufferedReader(new FileReader(infile));
			line = reader.readLine().trim().split("[\\s]+");
//			ext.checkHeader(line, expected, kill)
			while (reader.ready()) {
				line = reader.readLine().trim().split("\t",-1);
				for (int i=4; i<=6; i++) {
					if (!(new File("penn_data/"+line[i])).exists()) {
						cnv.analysis.AnalysisFormats.penncnv(proj, new String[] {line[i]}, null);
					}
				}
				if (System.getProperty("os.name").startsWith("Windows")) {
					//TODO How to get hmm, pfb, and gcmodel files?
					CmdLine.run("perl C:/penncnv/detect_cnv.pl -joint -hmm C:/penncnv/lib/hh550.hmm -pfb custom.pfb penn_data/"+line[4]+" penn_data/"+line[5]+" penn_data/"+line[6]+" -out "+line[4]+".jointcnv", "penn_output");
					(new File("penn_output/resultAll.jointcnv")).renameTo(new File("penn_output/resultAll_tmp.jointcnv"));
					common.Files.cat(new String[] {"penn_output/"+line[4]+".jointcnv", "penn_output/resultAll_tmp.jointcnv"}, "penn_output/resultAll.jointcnv", null, null);
					(new File("penn_output/"+line[4]+".jointcnv")).delete();
					if (!(new File("penn_output/resultAll_tmp.jointcnv")).delete()) {
						System.err.println("Error deleting the temporaryfile 'resultAll_tmp.jointcnv'.");
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+infile+"\" not found in current directory");
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+infile+"\"");
		}
		PennCNV.parseResults(proj, "penn_output/resultAll.jointcnv");
	}

	/**
	 * A script to select relevant sample data and move them to another folder
	 */
	public static void selectFilesToTestInPennCNV() {
		BufferedReader reader;
		String[] line;
		try {
			reader = new BufferedReader(new FileReader("D:/BOSS_outputFiles/TriosForDenovoCnv.txt"));
//			line = reader.readLine().trim().split("[\\s]+");
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().trim().split("\t",-1);
				for (int i=4; i<=6; i++) {
					if (!(new File("D:/BOSS/penn_data/TriosSamples/"+line[i])).exists()) {
//						Files.copy("penn_data/"+line[i], "penn_data/TriosSamples/"+line[i], REPLACE_EXISTING);
						Files.copy(Paths.get("D:/BOSS/penn_data/"+line[i]), Paths.get("D:/BOSS/penn_data/TriosSamples/"+line[i]));
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file D:/BOSS_outputFiles/TriosForDenovoCnv.txt not found in current directory.");
		} catch (IOException ioe) {
			System.err.println("Error reading file D:/BOSS_outputFiles/TriosForDenovoCnv.txt.");
		}
	}
	
	/**
	 * A script to select relevant sample data and move them to another folder
	 */
	public static void detectDenovoCnv() {

	    File folder = new File("D:/BOSS_outputFiles/joint_results/");
	    File[] listOfFiles = folder.listFiles();
		BufferedReader reader;
		String[] line;
		Vector<String[]> cnv;
		int j, k;
		boolean found;

	    for (int i = 0; i < listOfFiles.length; i++) {
	      if (listOfFiles[i].isFile() && listOfFiles[i].toString().contains(".jointcnv")) {
	        cnv = new Vector<String[]>();
			try {
				reader = new BufferedReader(new FileReader(listOfFiles[i].getPath()));
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					if(!line[0].contains("NOTICE")) {
						cnv.add(line);
					} else {
						break;
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file D:/BOSS_outputFiles/TriosForDenovoCnv.txt not found in current directory.");
			} catch (IOException ioe) {
				System.err.println("Error reading file D:/BOSS_outputFiles/TriosForDenovoCnv.txt.");
			}
			j = 0;
			while(j<cnv.size()) {
				if(cnv.elementAt(j)[7].contains("offspring")) {
					break;
				}
				j++;
			}
			k=j;
			while(j<cnv.size()) {
				found = false;
				for (int l=0; l<k; l++) {
					if (cnv.elementAt(j)[0].split(":")[0].equals(cnv.elementAt(l)[0].split(":")[0])
							&&	(   (Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[0]) >= Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[0])
									 &&  Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[0]) <= Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[1]) )
								||  (Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[1]) >= Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[0])
								     &&  Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[1]) <= Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[1]) )
								||  (Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[0]) >= Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[0])
								     &&  Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[0]) <= Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[1]) )
								||  (Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[1]) >= Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[0])
								     &&  Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[1]) <= Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[1]) )
								)) {
						found = true;
						break;
					}
				}
				if (!found) {
					System.out.println(listOfFiles[i].getName()+"\t"+cnv.elementAt(j)[0]);
				}
				j++;
			}


	      } else if (listOfFiles[i].isDirectory()) {
	        System.out.println("Directory " + listOfFiles[i].getName());
	      }
	    }


	}
	
	public static void main(String[] args) {
		proj = new Project("C:/workspace/Genvisis/projects/boss.properties", false);
		PennCNV.populationBAF(proj);
		PennCNV.gcModel(proj, "D:/PennCNV_Related/GcModel/gc5Base_hg18.txt", "D:/PennCNV_Related/GcModel/BOSS_Genvisis.gcmodel", 0);

//		loadData("D:/BOSS_pedigree_wDNA.txt");
//		runPennCnv("D:/TriosForDenovoCnv.txt");

//		selectFilesToTestInPennCNV();

//		detectDenovoCnv();

	}
	
}
