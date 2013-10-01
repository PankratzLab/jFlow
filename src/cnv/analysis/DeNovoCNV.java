package cnv.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import common.Files;
import java.nio.file.Paths;
import java.util.Vector;

import cnv.filesys.Project;

import common.HashVec;
import common.Logger;
import common.StringVector;
import common.CmdLine;
import common.ext;

public class DeNovoCNV {

//	private static Project proj;

	public static void generatePedigreeOfTrios(String pedigreeFileName) {
		BufferedReader reader;
		String[] line;
		PrintWriter writer;
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
//            		faId.add(line[0]+"\t"+line[2]);
//            		moId.add(line[0]+"\t"+line[3]);
            		dna.add(line[6]);
//            		dna.add(line[4]);
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

	// TODO 2-pass cleanup
//	public static void generatePedigreeOfTrios(String pedigreeFileName) {
//		BufferedReader reader;
//		String[] line;
//		PrintWriter writer;
//		StringVector fId, iId, faId, moId, dna;
//		
////		trav = new HashVec();
//		try {
////	        reader = Files.getReader(new FileReader(filename), proj.getJarStatus(), true, false);
//			reader = new BufferedReader(new FileReader(pedigreeFileName));
////			reader.readLine();
//			//TODO detect the column number of FID IID FAID MOID
//			fId = new StringVector();
//			iId = new StringVector();
//			faId = new StringVector();
//			moId = new StringVector();
//			dna = new StringVector();
//            while (reader.ready()) {
//            	line = reader.readLine().trim().split("\t",-1);
////            	if (line[6]!=null && !line[6].equals(".")) {
//            	if (line[4]!=null && !line[4].equals(".")) {
//            		fId.add(line[0]);
//            		iId.add(line[0]+"\t"+line[1]);
//            		faId.add(line[2]);
//            		moId.add(line[3]);
////            		dna.add(line[6]);
//            		dna.add(line[4]);
//            	}
//            }
//            reader.close();
//
//            writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(pedigreeFileName) + "PedigreeOfTrios.txt"));
//			writer.println("fId\tiId\tfaId\tmoId\tiDna\tfaDna\tmoDna");
//			for (int i=0; i<iId.size(); i++) {
//				if (iId.contains(fId.elementAt(i)+"\t"+ faId.elementAt(i)) && iId.contains(fId.elementAt(i)+"\t"+moId.elementAt(i))) {
//					writer.println(iId.elementAt(i)
//							+"\t"+faId.elementAt(i)
//							+"\t"+moId.elementAt(i)
//							+"\t"+dna.elementAt(i)
//							+"\t"+dna.elementAt(iId.indexOf(fId.elementAt(i)+"\t"+faId.elementAt(i)))
//							+"\t"+dna.elementAt(iId.indexOf(fId.elementAt(i)+"\t"+moId.elementAt(i))));
//				}
//			}
//			writer.flush();
//			writer.close();
//        } catch (FileNotFoundException fnfe) {
//        	System.err.println("Error: file \""+pedigreeFileName+"\" not found in current directory");
//        } catch (IOException ioe) {
//            System.err.println("Error reading file \""+pedigreeFileName+"\"");
//        }
//		System.out.println("Pedigree of Trios is ready at " + ext.parseDirectoryOfFile(pedigreeFileName) + "PedigreeOfTrios.txt");
//	}
	

	/**
	 * Call PennCNV to run trios analysis
	 * @param pedigreeForTrioAnalysis
	 * @param pennDataDir
	 * @param pennOutDir
	 * @param pennOutFileNameRoot
	 */
	public static void runPennCnv(Project proj, String pedigreeForTrioAnalysis, String pennDataDir, String pennOutDir, String pennOutFileNameRoot) {
		BufferedReader reader;
		String[] line;
		String command = null;
		
		if (!(new File(pennOutDir)).exists()) {
			new File(pennOutDir).mkdir();
		}
//		PennCNV.populationBAF(proj);
//		PennCNV.gcModel(proj, "D:/PennCNV_Related/GcModel/gc5Base_hg18.txt", "penn_output/gcmodel", 0);
		
		try {
			reader = new BufferedReader(new FileReader(pedigreeForTrioAnalysis));
			reader.readLine();
//			ext.checkHeader(line, expected, kill)
			while (reader.ready()) {
				line = reader.readLine().trim().split("\t",-1);
				for (int i=4; i<=6; i++) {
					if (!(new File(pennDataDir + line[i])).exists()) {
						cnv.analysis.AnalysisFormats.penncnv(proj, new String[] {line[i]}, null);
					}
				}
				if (System.getProperty("os.name").startsWith("Windows")) {
					//TODO How to get hmm, pfb, and gcmodel files?
					CmdLine.run("perl C:/penncnv/detect_cnv.pl -joint -hmm C:/penncnv/lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + proj.getProjectDir() + "custom.gcmodel " + pennDataDir + line[4] + " " + pennDataDir + line[5] + " " + pennDataDir + line[6] + " -out " + line[4] + ".jointcnv", pennOutDir);
					CmdLine.run("perl C:/penncnv/detect_cnv.pl -trio -hmm C:/penncnv/lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -cnv " + pennOutDir + pennOutFileNameRoot + " " + pennDataDir + line[4] + " " + pennDataDir + line[5] + " " + pennDataDir + line[6] + " -out " + line[4] + ".triocnv", pennOutDir);
//					(new File(cnvOutputDir+"resultAll.jointcnv")).renameTo(new File(cnvOutputDir + "resultAll_tmp.jointcnv"));
//					common.Files.cat(new String[] {cnvOutputDir+line[4]+".jointcnv", cnvOutputDir+"resultAll_tmp.jointcnv"}, cnvOutputDir+"resultAll.jointcnv", null, null);
//					(new File(cnvOutputDir + line[4] + ".jointcnv")).delete();
//					if (!(new File("penn_output/resultAll_tmp.jointcnv")).delete()) {
//						System.err.println("Error deleting the temporaryfile 'resultAll_tmp.jointcnv'.");
//					}
				} else {
//					Runtime.getRuntime().exec("perl /home/pankrat2/shared/penncnv/detect_cnv.pl -joint -hmm /home/pankrat2/shared/penncnv/lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + proj.getProjectDir() + "custom.gcmodel " + pennDataDir + line[4] + " " + pennDataDir + line[5] + " " + pennDataDir + line[6] + " -out " + line[4] + ".jointcnv");
//					Runtime.getRuntime().exec("perl /home/pankrat2/shared/penncnv/detect_cnv.pl -trio -hmm /home/pankrat2/shared/penncnv/lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -cnv " + pennOutDir + pennOutFileNameRoot + " " + pennDataDir + line[4] + " " + pennDataDir + line[5] + " " + pennDataDir + line[6] + " -out " + line[4] + ".triocnv");
//					Files.qsub("runPennCNV.qsub", "perl /home/pankrat2/shared/penncnv/detect_cnv.pl -joint -hmm /home/pankrat2/shared/penncnv/lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + proj.getProjectDir() + "custom.gcmodel " + pennDataDir + line[4] + " " + pennDataDir + line[5] + " " + pennDataDir + line[6] + " -out " + line[4] + ".jointcnv", -1, -1);
//					Files.qsub("runPennCNV.qsub", "perl /home/pankrat2/shared/penncnv/detect_cnv.pl -trio -hmm /home/pankrat2/shared/penncnv/lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -cnv " + pennOutDir + pennOutFileNameRoot + " " + pennDataDir + line[4] + " " + pennDataDir + line[5] + " " + pennDataDir + line[6] + " -out " + line[4] + ".triocnv", -1, -1);
					command += "perl /home/pankrat2/shared/penncnv/detect_cnv.pl -joint -hmm /home/pankrat2/shared/penncnv/lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + proj.getProjectDir() + "custom.gcmodel " + pennDataDir + line[4] + " " + pennDataDir + line[5] + " " + pennDataDir + line[6] + " -out " + line[4] + ".jointcnv\n"
								+ "perl /home/pankrat2/shared/penncnv/detect_cnv.pl -trio -hmm /home/pankrat2/shared/penncnv/lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -cnv " + pennOutDir + pennOutFileNameRoot + " " + pennDataDir + line[4] + " " + pennDataDir + line[5] + " " + pennDataDir + line[6] + " -out " + line[4] + ".triocnv\n";
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+pedigreeForTrioAnalysis+"\" not found in current directory");
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+pedigreeForTrioAnalysis+"\"");
		}

// TODO		
//		command = "./my_perl_script.pl -arg1 [%0] [%1] [%2] -pfb constantfile -gcmodel same_file";
//		String[][] iterations = new String[][] {{"fa1", "mo1", "child1"}, {"fa2", "mo2", "child2"}};
//		Files.qsub(root_batch_name, dirToSwitchToBeforeRunning, 20, command, iterations, 1000);
		
		Files.qsub("runPennCNV.qsub", command, -1, -1);
		
		// this method will create master.[root_batch_name]
		// qsub batch1.qsub
		// qsub batch2.qsub
		// qsub batch3.qsub
		
		// ./master.[root]
		
//		Files.batchIt(root_batch_name, init, numBatches, commands, iterations);
//		PennCNV.parseResults(proj, "penn_output/resultAll.jointcnv", false, new Logger());
	}
	
	public static void batch(String pedigreeOfTrio) {
		String command;
		String[][] iterations;
		
//		command = "perl C:/penncnv/detect_cnv.pl -joint -hmm C:/penncnv/lib/hh550.hmm -pfb C:/penncnv/lib/hh550.hg18.pfb -gcmodel C:/penncnv/lib/hh550.hg18.gcmodel [%1] [%2] [%3] -out [%1].jointcnv";
		command = "perl ../bin/penncnv/detect_cnv.pl -joint -hmm ../bin/penncnv/lib/hhall.hmm -pfb ../bin/custom.pfb -gcmodel ../bin/custom.gcmodel [%1] [%2] [%0] -out ../results/[%0].jointcnv -log ../results/[%0].log";
		
		iterations = HashVec.loadFileToStringMatrix(pedigreeOfTrio, true, new int[] {4, 5, 6}, false);
		
		common.Files.qsub("denovo", "/share/bulk/gedi/pankr018/denovo/penn_data", 65, command, iterations, -1);
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
//						java.nio.file.Files.copy("penn_data/"+line[i], "penn_data/TriosSamples/"+line[i], REPLACE_EXISTING);
						java.nio.file.Files.copy(Paths.get("D:/BOSS/penn_data/"+line[i]), Paths.get("D:/BOSS/penn_data/TriosSamples/"+line[i]));
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
//		proj = new Project("C:/workspace/Genvisis/projects/OSv2.properties", false);
//		PennCNV.populationBAF(proj);
//		PennCNV.gcModel(proj, "D:/PennCNV_Related/GcModel/gc5Base_hg18.txt", "D:/PennCNV_Related/GcModel/BOSS_Genvisis.gcmodel", 0);

//		generatePedigreeOfTrios("D:/GEDI/GEDI_CNV_Pedigree.txt");
//		runPennCnv("D:/GEDI/PedigreeOfTrios.txt", "D:/GEDI/cnv_output/", "D:/GEDI/samples_all/", "gedi_all.rawcnv");
//		batch("D:/GEDI/PedigreeOfTrios.txt");

//		selectFilesToTestInPennCNV("D:/BOSS_outputFiles/TriosForDenovoCnv.txt");

//		detectDenovoCnv();

		String projPropertyFileFullPath;
		Project proj;
		String gcBaseFileFullPath;
		String gcModelFileFullPath;
		int numWindowUnits;
		String pedigreeFullPath;
		String pennDataDir;
		String pennOutputDir;
		String pennOutputFilenameRoot;

		projPropertyFileFullPath = "C:/workspace/Genvisis/projects/OSv2.properties";
		gcBaseFileFullPath = "D:/PennCNV_Related/GcModel/gc5Base_hg19.txt";
		numWindowUnits = 0;
		pennOutputFilenameRoot = "gedi_all.rawcnv";

		String usage = "\n"+
				"cnv.analysis.DeNovoCNV requires 4 arguments\n"+
				"   (1) The project's property file's name (i.e. proj=" + projPropertyFileFullPath + " (not the default))\n"+
//				"   (2) The output gcBase file name  (i.e. gcbase=" + gcBaseFileFullPath + " (not the default))\n"+
				"   (3) Number of window units for GC model (i.e. nwin=" + numWindowUnits + " (not the default))\n"+
//				"   (4) The output gc model file name (i.e. gcmodel=" + gcModelFileFullPath + " (not the default))\n"+
//				"   (5) The output pedigree file name (i.e. pedigree=" + pedigreeFullPath + " (not the default))\n"+
//				"   (6) The penn data directory (i.e. penndatadir=" + pennDataDir + " (not the default))\n"+
//				"   (7) The penn output directory (i.e. pennoutdir=" + pennOutputDir + " (not the default))\n"+
//				"   (8) The final output file name (i.e. outfile=" + pennOutputFilenameRoot + " (not the default))\n"+
				"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				projPropertyFileFullPath = args[i].split("=")[1];
			} else if (args[i].startsWith("gcbase=")) {
				gcBaseFileFullPath = args[i].split("=")[1];
			} else if (args[i].startsWith("nwin=")) {
				numWindowUnits = Integer.parseInt(args[i].split("=")[1]);
			} else if (args[i].startsWith("gcmodel=")) {
				gcModelFileFullPath = args[i].split("=")[1];
			} else if (args[i].startsWith("pedigree=")) {
				pedigreeFullPath = args[i].split("=")[1];
			} else if (args[i].startsWith("penndatadir=")) {
				pennDataDir = args[i].split("=")[1];
			} else if (args[i].startsWith("pennoutdir=")) {
				pennOutputDir = args[i].split("=")[1];
			} else if (args[i].startsWith("outfile=")) {
				pennOutputFilenameRoot = args[i].split("=")[1];
			} else {
				System.err.println("Error - invalid argument: "+args[i]);
			}
		}
		
		proj = new Project(projPropertyFileFullPath, false);
		gcModelFileFullPath = proj.getProjectDir() + "custom.gcmodel";
		pedigreeFullPath = proj.getDir(Project.DATA_DIRECTORY) + "pedigree.dat";
		pennDataDir = proj.getProjectDir() + "penn_data/";
		pennOutputDir = proj.getProjectDir() + "penn_output/";

//		PennCNV.populationBAF(proj);
//		PennCNV.gcModel(proj, gcBaseFileFullPath, gcModelOutputFullPath, numWindowUnits);

//		generatePedigreeOfTrios(pedigreeFullPath);
		runPennCnv(proj, ext.parseDirectoryOfFile(pedigreeFullPath) + "PedigreeOfTrios.txt", pennDataDir, pennOutputDir, pennOutputFilenameRoot);
//		batch(pedigreeFullPath);
//
//		selectFilesToTestInPennCNV("D:/BOSS_outputFiles/TriosForDenovoCnv.txt");
//
//		detectDenovoCnv();
	}
	
}
