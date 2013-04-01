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

import common.HashVec;
import common.StringVector;
import common.CmdLine;
import common.ext;

public class DeNovoCNV {

	private static Project proj;
	
	public static void generatePedigreeOfTrios(String pedigreeFileName) {
		BufferedReader reader;
		String[] line;
		PrintWriter writer;
//		HashVec trav;
//		String[] newID;
		StringVector fId, iId, faId, moId, dna;
		
//		trav = new HashVec();
		try {
//	        reader = Files.getReader(new FileReader(filename), proj.getJarStatus(), true, false);
			reader = new BufferedReader(new FileReader(pedigreeFileName));
//			reader.readLine();
			//TODO detect the column number of FID IID FAID MOID
			fId = new StringVector();
			iId = new StringVector();
			faId = new StringVector();
			moId = new StringVector();
			dna = new StringVector();
            while (reader.ready()) {
            	line = reader.readLine().trim().split("\t",-1);
//            	if (line[6]!=null && !line[6].equals(".")) {
            	if (line[4]!=null && !line[4].equals(".")) {
            		fId.add(line[0]);
            		iId.add(line[0]+"\t"+line[1]);
            		faId.add(line[2]);
            		moId.add(line[3]);
//            		dna.add(line[6]);
            		dna.add(line[4]);
            	}
            }
            reader.close();

            writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(pedigreeFileName) + "PedigreeOfTrios.txt"));
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
        	System.err.println("Error: file \""+pedigreeFileName+"\" not found in current directory");
        } catch (IOException ioe) {
            System.err.println("Error reading file \""+pedigreeFileName+"\"");
        }
		System.out.println("Pedigree of Trios is ready at " + ext.parseDirectoryOfFile(pedigreeFileName) + "PedigreeOfTrios.txt");
	}
	
	/**
	 * 
	 * @param pedigreeOfTrio
	 */
	public static void runPennCnv(String pedigreeOfTrio, String cnvOutputDir, String sampleDir, String rawCnvFileNameNoPath) {
		BufferedReader reader;
		String[] line;
		
		if (!(new File(cnvOutputDir)).exists()) {
			new File(cnvOutputDir).mkdir();
		}
//		PennCNV.populationBAF(proj);
//		PennCNV.gcModel(proj, "D:/PennCNV_Related/GcModel/gc5Base_hg18.txt", "penn_output/gcmodel", 0);
		
		try {
			reader = new BufferedReader(new FileReader(pedigreeOfTrio));
			line = reader.readLine().trim().split("[\\s]+");
//			ext.checkHeader(line, expected, kill)
			while (reader.ready()) {
				line = reader.readLine().trim().split("\t",-1);
				for (int i=4; i<=6; i++) {
					if (!(new File(sampleDir + line[i])).exists()) {
						cnv.analysis.AnalysisFormats.penncnv(proj, new String[] {line[i]}, null);
					}
				}
				if (System.getProperty("os.name").startsWith("Windows")) {
					//TODO How to get hmm, pfb, and gcmodel files?
					CmdLine.run("perl C:/penncnv/detect_cnv.pl -joint -hmm C:/penncnv/lib/hh550.hmm -pfb C:/penncnv/lib/hh550.hg18.pfb -gcmodel C:/penncnv/lib/hh550.hg18.gcmodel " + sampleDir + line[4] + " " + sampleDir + line[5] + " " + sampleDir + line[6] + " -out " + line[4] + ".jointcnv", cnvOutputDir);
					CmdLine.run("perl C:/penncnv/detect_cnv.pl -trio -hmm C:/penncnv/lib/hh550.hmm -pfb C:/penncnv/lib/hh550.hg18.pfb -cnv " + cnvOutputDir + rawCnvFileNameNoPath + " " + sampleDir + line[4] + " " + sampleDir + line[5] + " " + sampleDir + line[6] + " -out " + line[4] + ".triocnv", cnvOutputDir);
//					(new File(cnvOutputDir+"resultAll.jointcnv")).renameTo(new File(cnvOutputDir + "resultAll_tmp.jointcnv"));
//					common.Files.cat(new String[] {cnvOutputDir+line[4]+".jointcnv", cnvOutputDir+"resultAll_tmp.jointcnv"}, cnvOutputDir+"resultAll.jointcnv", null, null);
//					(new File(cnvOutputDir + line[4] + ".jointcnv")).delete();
//					if (!(new File("penn_output/resultAll_tmp.jointcnv")).delete()) {
//						System.err.println("Error deleting the temporaryfile 'resultAll_tmp.jointcnv'.");
//					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+pedigreeOfTrio+"\" not found in current directory");
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+pedigreeOfTrio+"\"");
		}
		PennCNV.parseResults(proj, "penn_output/resultAll.jointcnv");
	}
	
	public static void batch(String pedigreeOfTrio) {
		String command;
		String[][] iterations;
		
		command = "perl C:/penncnv/detect_cnv.pl -joint -hmm C:/penncnv/lib/hh550.hmm -pfb C:/penncnv/lib/hh550.hg18.pfb -gcmodel C:/penncnv/lib/hh550.hg18.gcmodel [%1] [%2] [%3] -out [%1].jointcnv";
		
		iterations = HashVec.loadFileToStringMatrix(pedigreeOfTrio, true, new int[] {4, 5, 6}, false);
		
		common.Files.qsub("denovo", 20, command, iterations);
	}

	/**
	 * A script to select relevant sample data and move them to another folder
	 */
	public static void selectFilesToTestInPennCNV(String trioCnvFileName) {
		BufferedReader reader;
		String[] line;
		try {
			reader = new BufferedReader(new FileReader(trioCnvFileName));
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
//		proj = new Project("C:/workspace/Genvisis/projects/boss.properties", false);
//		PennCNV.populationBAF(proj);
//		PennCNV.gcModel(proj, "D:/PennCNV_Related/GcModel/gc5Base_hg18.txt", "D:/PennCNV_Related/GcModel/BOSS_Genvisis.gcmodel", 0);

//		generatePedigreeOfTrios("D:/GEDI/GEDI_CNV_Pedigree.txt");
//		runPennCnv("D:/GEDI/PedigreeOfTrios.txt", "D:/GEDI/cnv_output/", "D:/GEDI/samples_all/", "gedi_all.rawcnv");
		batch("D:/GEDI/PedigreeOfTrios.txt");

//		selectFilesToTestInPennCNV("D:/BOSS_outputFiles/TriosForDenovoCnv.txt");

//		detectDenovoCnv();

	}
	
}
