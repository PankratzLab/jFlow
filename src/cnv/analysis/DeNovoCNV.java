package cnv.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import common.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import sun.util.logging.resources.logging;

import cnv.filesys.Project;
import cnv.var.CNVariant;
import cnv.var.SampleData;

import common.Array;
import common.HashVec;
import common.Logger;
import common.Positions;
import common.StringVector;
import common.CmdLine;
import common.ext;

public class DeNovoCNV {

//	private static Project proj;

	public static void generateTriosPedigree(String input_pedigreeFullPath, String output_trioPedigreeFullPath) {
		BufferedReader reader;
		String[] line;
		PrintWriter writer;
		StringVector fId, iId, faId, moId, dna;
		
//		trav = new HashVec();
		try {
//	        reader = Files.getReader(new FileReader(filename), proj.getJarStatus(), true, false);
			reader = new BufferedReader(new FileReader(input_pedigreeFullPath));
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

//          writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(input_pedigreeFullPath) + "PedigreeOfTrios.txt"));
            writer = new PrintWriter(new FileWriter(output_trioPedigreeFullPath));
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
        	System.err.println("Error: file \""+input_pedigreeFullPath+"\" not found in current directory");
        } catch (IOException ioe) {
            System.err.println("Error reading file \""+input_pedigreeFullPath+"\"");
        }
		System.out.println("Pedigree of Trios is ready at " + output_trioPedigreeFullPath);
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
	 * @param proj
	 * @param pedigreeForTrioAnalysis
	 * @param pennDataDir
	 * @param pennOutDir
	 * @param pennOutFileNameRoot
	 * @param pennBinDir
	 */
	public static void runPennCnv(Project proj, String pedigreeForTrioAnalysis, String pennDataDir, String gcModelFullPath, String pennOutDir, String pennOutFileNameRoot, String pennBinDir, int numQsubFiles, int numCommandsPerQsubFile) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String commands;
		Vector <String[]> iterationsVec;
		String[][] iterations;
		
		if (!(new File(pennOutDir)).exists()) {
			new File(pennOutDir).mkdir();
		}
//		PennCNV.populationBAF(proj);
//		PennCNV.gcModel(proj, "D:/PennCNV_Related/GcModel/gc5Base_hg18.txt", "penn_output/gcmodel", 0);

		iterationsVec = new Vector<String[]>();
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
				iterationsVec.add(new String[] {line[5], line[6], line[4]});
				writer = new PrintWriter(new FileOutputStream(pennDataDir + "lists/listSolo_" + line[4]));
				writer.println(line[5] + "\n" + line[6] + "\n" + line[4]);
				writer.close();
				writer = new PrintWriter(new FileOutputStream(pennDataDir + "lists/listTrio_" + line[4]));
				writer.println(line[5] + "\t" + line[6] + "\t" + line[4]);
				writer.close();
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+pedigreeForTrioAnalysis+"\" not found in current directory");
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+pedigreeForTrioAnalysis+"\"");
		}
		iterations = new String[iterationsVec.size()][3];
		for (int i = 0; i < iterations.length; i++) {
			iterations[i] = iterationsVec.elementAt(i);
		}

		if (numQsubFiles >= 0) {
			if (numCommandsPerQsubFile >= 0) {
				System.err.println("Warning: cannot specify numQsubFiles and numCommandsPerQsubFile at the same time. Program continues by ignoring numCommandsPerQsubFile.");
			}
		} else if (numCommandsPerQsubFile >= 0) {
			numQsubFiles = (int) Math.ceil((double) iterations.length * 3 / numCommandsPerQsubFile);
		} else {
			System.err.println("Error: must specify either numQsubFiles or numCommandsPerQsubFile. Program exits.");
			return;
		}

		if (System.getProperty("os.name").startsWith("Windows")) {
			for (int i = 0; i < iterations.length; i++) {
				CmdLine.run("perl " + pennBinDir + "detect_cnv.pl -test -hmm " + pennBinDir + "lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + gcModelFullPath + " " + pennDataDir + iterations[i][0] + " " + pennDataDir + iterations[i][1] + " " + pennDataDir + iterations[i][2] + " -out " + pennOutDir + iterations[i][2] + ".rawcnv", pennOutDir);
				CmdLine.run("perl " + pennBinDir + "detect_cnv.pl -trio -hmm " + pennBinDir + "lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + gcModelFullPath + " -cnv " + pennOutDir + iterations[i][0] + ".rawcnv " + pennDataDir + iterations[i][0] + " " + pennDataDir + iterations[i][1] + " " + pennDataDir + iterations[i][2] + " -out " + pennOutDir + iterations[i][2] + ".triocnv", pennOutDir);
				CmdLine.run("perl " + pennBinDir + "detect_cnv.pl -joint -hmm " + pennBinDir + "lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + gcModelFullPath + " " + pennDataDir + iterations[i][0] + " " + pennDataDir + iterations[i][1] + " " + pennDataDir + iterations[i][2] + " -out " + pennOutDir + iterations[i][2] + ".jointcnv", pennOutDir);
			}

//			(new File(cnvOutputDir+"resultAll.jointcnv")).renameTo(new File(cnvOutputDir + "resultAll_tmp.jointcnv"));
//			common.Files.cat(new String[] {cnvOutputDir+line[4]+".jointcnv", cnvOutputDir+"resultAll_tmp.jointcnv"}, cnvOutputDir+"resultAll.jointcnv", null, null);
//			(new File(cnvOutputDir + line[4] + ".jointcnv")).delete();
//			if (!(new File("penn_output/resultAll_tmp.jointcnv")).delete()) {
//				System.err.println("Error deleting the temporaryfile 'resultAll_tmp.jointcnv'.");
//			}

		} else {
//			Runtime.getRuntime().exec("/home/pankrat2/shared/penncnv/detect_cnv.pl -joint -hmm /home/pankrat2/shared/penncnv/lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + gcBaseFileFullPath + " " + pennDataDir + line[4] + " " + pennDataDir + line[5] + " " + pennDataDir + line[6] + " -out " + line[4] + ".jointcnv");
//			Runtime.getRuntime().exec("/home/pankrat2/shared/penncnv/detect_cnv.pl -trio -hmm /home/pankrat2/shared/penncnv/lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -cnv " + pennOutDir + pennOutFileNameRoot + " " + pennDataDir + line[4] + " " + pennDataDir + line[5] + " " + pennDataDir + line[6] + " -out " + line[4] + ".triocnv");

//			commands = "";
//			for (int i = 0; i < iterations.length; i++) {
//				commands += pennBinDir + "detect_cnv.pl -test -hmm " + pennBinDir + "lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + gcBaseFileFullPath + " " + pennDataDir + iterations[i][0] + " " + pennDataDir + iterations[i][1] + " " + pennDataDir + iterations[i][2] + " -out " + pennOutDir + iterations[i][2] + ".rawcnv\n"
//							+ pennBinDir + "detect_cnv.pl -trio -hmm " + pennBinDir + "lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + gcBaseFileFullPath + " -cnv " + pennOutDir + iterations[i][0] + ".rawcnv " + pennDataDir + iterations[i][0] + " " + pennDataDir + iterations[i][1] + " " + pennDataDir + iterations[i][2] + " -out " + pennOutDir + iterations[i][2] + ".triocnv\n"
//							+ pennBinDir + "detect_cnv.pl -joint -hmm " + pennBinDir + "lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + gcBaseFileFullPath + " " + pennDataDir + iterations[i][0] + " " + pennDataDir + iterations[i][1] + " " + pennDataDir + iterations[i][2] + " -out " + pennOutDir + iterations[i][2] + ".jointcnv\n";
//			}
			if (proj.getBoolean(Project.PENNCNV_GZIP_YESNO)) {
				commands = pennBinDir + "detect_cnv.pl -test -hmm " + pennBinDir + "lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + gcModelFullPath + " -list " + pennDataDir + "lists/listSolo_[%2].txt -out " + pennOutDir + "[%2].rawcnv\n"
						 + pennBinDir + "detect_cnv.pl -trio -hmm " + pennBinDir + "lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + gcModelFullPath + " -cnv " + pennDataDir +"[%2].rawcnv -list " + pennDataDir + "/lists/listTrio_[%2].txt -out " + pennOutDir + "[%2].triocnv\n"
						 + pennBinDir + "detect_cnv.pl -joint -hmm " + pennBinDir + "lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + gcModelFullPath + " -list " + pennDataDir + "/lists/listTrio_[%2].txt -out " + pennOutDir + "[%2].jointcnv\n";
				Files.qsub("runPennCNV", proj.getProjectDir() + "penn_scripts/", numQsubFiles, commands, iterations, 2000, 24);
				Files.qsubMultiple("chunkCNV", Array.stringArraySequence(numQsubFiles, proj.getProjectDir() + "penn_scripts/runPennCNV_", ".qsub"), 2000, 24);

			} else {
				commands = pennBinDir + "detect_cnv.pl -test -hmm " + pennBinDir + "lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + gcModelFullPath + " " + pennDataDir + "[%0] " + pennDataDir + "[%1] " + pennDataDir +"[%2] -out " + pennOutDir + "[%2].rawcnv\n"
						 + pennBinDir + "detect_cnv.pl -trio -hmm " + pennBinDir + "lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + gcModelFullPath + " -cnv " + pennDataDir +"[%2].rawcnv " + pennDataDir + "[%0] " + pennDataDir + "[%1] " + pennDataDir + "[%2] -out " + pennOutDir + "[%2].triocnv\n"
						 + pennBinDir + "detect_cnv.pl -joint -hmm " + pennBinDir + "lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb -gcmodel " + gcModelFullPath + " " + pennDataDir + "[%0] " + pennDataDir + "[%1] " + pennDataDir + "[%2] -out " + pennOutDir + "[%2].jointcnv\n";
				Files.qsub("runPennCNV", proj.getProjectDir() + "penn_scripts/", numQsubFiles, commands, iterations, 2000, 24);
				Files.qsubMultiple("chunkCNV", Array.stringArraySequence(numQsubFiles, proj.getProjectDir() + "penn_scripts/runPennCNV_", ".qsub"), 2000, 24);
	
//				this method will create master.[root_batch_name]
//				qsub batch1.qsub
//				qsub batch2.qsub
//				qsub batch3.qsub
	
//				./master.[root]
			}
		}
	}

	public static void parsePennCnvResult(Project proj, String pennCnvResultDir, String pennCnvResultFileNameExt, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp, currentOffSpring;
		Vector<String> warnings = new Vector<String>();
		int[] position;
		String score;
		SampleData sampleData;
		String famIndPair;
		Hashtable<String, String> offspringCnv;
		String[] ids;
		String[] filenames;
		boolean offspringExisted;

		sampleData = proj.getSampleData(2, false);
		if (pennCnvResultFileNameExt.startsWith(".")) {
			pennCnvResultFileNameExt = pennCnvResultFileNameExt.substring(1);
		}
		try {
			filenames = Files.list(pennCnvResultDir, pennCnvResultFileNameExt, false);
			writer = new PrintWriter(new FileWriter(proj.getDir(Project.DATA_DIRECTORY) + "denovo_" + pennCnvResultFileNameExt.replace("cnv", "") + ".cnv"));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			offspringCnv = new Hashtable<String, String>();
			for (int i = 0; i < filenames.length; i++) {
//				currentOffspring = ext.rootOf(filenames[i]);
				offspringExisted = false;
				reader = new BufferedReader(new FileReader(pennCnvResultDir + filenames[i]));
				while (reader.ready()) {
					temp = reader.readLine();
					if (! temp.startsWith("NOTICE:")) {
						temp = PennCNV.translateDerivedSamples(temp);
						line = temp.trim().split("[\\s]+");
						if (line[7].equalsIgnoreCase("offspring")) {
							if (line[8].split("=")[1].startsWith("33")) {
								position = Positions.parseUCSClocation(line[0]);
								currentOffSpring = line[4];
								currentOffSpring = currentOffSpring.substring(currentOffSpring.lastIndexOf("/")+1);
								ids = sampleData.lookup(currentOffSpring);
								if (ids == null) {
//									if (! offspringCnv.containsKey(trav)) {
//										if (log != null) {
//											log.reportError("Error - '" + trav + "' was not found in " + proj.getFilename(Project.SAMPLE_DATA_FILENAME));
//										} else {
//											System.out.println("Error - '" + trav + "' was not found in " + proj.getFilename(Project.SAMPLE_DATA_FILENAME));
//										}
//									}
									famIndPair = currentOffSpring + "\t" + currentOffSpring;
								} else {
									famIndPair = ids[1];
								}
								if (! offspringExisted) {
									offspringCnv.put(currentOffSpring, famIndPair);
									offspringExisted = true;
								}
				
								if (line.length<8 || !line[7].startsWith("conf=")) {
									score = "127";
									if (!warnings.contains(currentOffSpring) && warnings.size() < 10) {
										if (log != null) {
											log.reportError("Warning - no conf estimates for " + currentOffSpring);
										} else {
											System.out.println("Warning - no conf estimates for " + currentOffSpring);
										}
										warnings.add(currentOffSpring);
									}
								} else {
									score = ext.formDeci(Double.parseDouble(line[7].substring(5)), 4, true);
								}
								writer.println(famIndPair + "\t" + position[0] + "\t" + position[1] + "\t" + position[2] + "\t" + line[3].substring(line[3].indexOf("=")+1) + "\t" + score + "\t" + line[1].substring(7));
							}
						}
					}
				}
				reader.close();
			}
			writer.close();

//			FilterCalls.stdFilters(dir, ext.rootOf(filename)+".cnv", MAKE_UCSC_TRACKS);

		} catch (FileNotFoundException fnfe) {
			return;
		} catch (IOException ioe) {
			return;
		}
	}

	public static void parsePennCnvResult(Project proj, String pennCnvResultDir, String triosPedigreeFullPath, String pennCnvResultFileNameExt, Logger log) {
		BufferedReader reader;
		PrintWriter writer1, writer2;
		String[] line;
		String temp, currentOffspring, currentFather = null, currentMother = null;
//		Vector<String> warnings = new Vector<String>(); // TODO check to see if -conf is an option for trio/joint calling
		int[] position;
		String score;
		SampleData sampleData;
		String famIndPair = null;
		Hashtable<String, String[]> triosPedigree;
		String[] ids;
		String[] filenames;
		String[] currentParents;

		sampleData = proj.getSampleData(2, false);
		if (pennCnvResultFileNameExt.startsWith(".")) {
			pennCnvResultFileNameExt = pennCnvResultFileNameExt.substring(1);
		}
		try {
			triosPedigree = new Hashtable<String, String[]>();
			reader = new BufferedReader(new FileReader(triosPedigreeFullPath));
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				triosPedigree.put(line[4], new String[] {line[5], line[6]});
			}

			filenames = Files.list(pennCnvResultDir, pennCnvResultFileNameExt, false);
			writer1 = new PrintWriter(new FileWriter(proj.getDir(Project.DATA_DIRECTORY) + "denovo_" + pennCnvResultFileNameExt.replace("cnv", "") + ".cnv"));
			writer1.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			writer2 = new PrintWriter(new FileWriter(proj.getDir(Project.DATA_DIRECTORY) + "denovo_" + pennCnvResultFileNameExt.replace("cnv", "") + "_list.txt"));
			writer2.println("DNA\tPosition\tComments");
			for (int i = 0; i < filenames.length; i++) {
				currentOffspring = null;
				reader = new BufferedReader(new FileReader(pennCnvResultDir + filenames[i]));
				while (reader.ready()) {
					temp = reader.readLine();
					if (! temp.startsWith("NOTICE:")) {
						temp = PennCNV.translateDerivedSamples(temp);
						line = temp.trim().split("[\\s]+");
						if (line[7].equalsIgnoreCase("offspring")) {
							if (line[8].split("=")[1].startsWith("33")) {
								position = Positions.parseUCSClocation(line[0]);
								if (currentOffspring == null) {
									currentOffspring = line[4];
									currentOffspring = currentOffspring.substring(currentOffspring.lastIndexOf("/") + 1);
									currentParents = triosPedigree.get(currentOffspring);
									currentFather = currentParents[0];
									currentMother = currentParents[1];
									ids = sampleData.lookup(currentOffspring);
									if (ids == null) {
//										if (! offspringCnv.containsKey(trav)) {
//											if (log != null) {
//												log.reportError("Error - '" + trav + "' was not found in " + proj.getFilename(Project.SAMPLE_DATA_FILENAME));
//											} else {
//												System.out.println("Error - '" + trav + "' was not found in " + proj.getFilename(Project.SAMPLE_DATA_FILENAME));
//											}
//										}
										famIndPair = currentOffspring + "\t" + currentOffspring;
									} else {
										famIndPair = ids[1];
									}
					
								}
								if (line.length<8 || !line[7].startsWith("conf=")) {
									score = "127";
//									if (!warnings.contains(trav) && warnings.size() < 10) {
//										if (log != null) {
//											log.reportError("Warning - no conf estimates for " + trav);
//										} else {
//											System.out.println("Warning - no conf estimates for " + trav);
//										}
//										warnings.add(trav);
//									}
								} else {
									score = ext.formDeci(Double.parseDouble(line[7].substring(5)), 4, true);
								}
								writer1.println(famIndPair + "\t" + position[0] + "\t" + position[1] + "\t" + position[2] + "\t" + line[3].substring(line[3].indexOf("=")+1) + "\t" + score + "\t" + line[1].substring(7));
								writer2.println(currentOffspring + "\t" + line[0] + "\toffspring");
								writer2.println(currentFather + "\t" + line[0] + "\tfather");
								writer2.println(currentMother + "\t" + line[0] + "\tmother");
							}
						}
					}
				}
				reader.close();
			}
			writer1.close();
			writer2.close();

//			FilterCalls.stdFilters(dir, ext.rootOf(filename)+".cnv", MAKE_UCSC_TRACKS);

		} catch (FileNotFoundException fnfe) {
			return;
		} catch (IOException ioe) {
			return;
		}
	}

	public static void generateTriosFamFile(Project proj, String triosPedigreeFileFullPath, boolean pedigreeHasHeader, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		SampleData sampleData;

		sampleData = proj.getSampleData(2, false);
		try {
			reader = new BufferedReader(new FileReader(triosPedigreeFileFullPath));
			if (pedigreeHasHeader) {
				reader.readLine();
			}
			writer = new PrintWriter(new FileWriter(proj.getDir(Project.DATA_DIRECTORY) + "trios.fam"));
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				writer.println(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + sampleData.getSexForIndividual(line[4]) + "\t1");
			}
			writer.close();
			reader.close();

		} catch (FileNotFoundException fnfe) {
			return;
		} catch (IOException ioe) {
			return;
		}
	}

//	/**
//	 * Does not work. Still work in progress.
//	 * @param proj
//	 * @param cnvFileFullPath
//	 * @param cnvFileHasHeader
//	 * @param log
//	 */
//	public static void generateDenovoList(Project proj, String cnvFileFullPath, boolean cnvFileHasHeader, Logger log) {
//		BufferedReader reader;
//		PrintWriter writer;
//		String[] line;
//		String temp, currentOffSpring;
//		int[] position;
//		String score;
//		SampleData sampleData;
//		String famIndPair;
//		Hashtable<String, String> offspringCnv;
//		String[] ids;
//		String[] filenames;
//		boolean offspringExisted;
//
//		sampleData = proj.getSampleData(2, false);
//		try {
//			reader = new BufferedReader(new FileReader(cnvFileFullPath));
//			if (cnvFileHasHeader) {
//				reader.readLine();
//			}
//			writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(cnvFileFullPath) + ext.rootOf(cnvFileFullPath) + "_list.txt"));
//			writer.write("DNA\tPosition\tComments");
//			while (reader.ready()) {
//				line = reader.readLine().split("\t");
////				if(offspringCnv.containsKey(line[4])) {
////					writer.println(offspringCnv.get(line[4]) + "\t" + sampleData.lookup(line[5])[2]  + "\t" + sampleData.lookup(line[6])[2] + "\t" + sampleData.getSexForIndividual(line[4]) + "\t1");
////					writer.println(offspringCnv.get(line[4]) + "\t" + sampleData.lookup(line[5])[2]  + "\t" + sampleData.lookup(line[6])[2] + "\t" + sampleData.getSexForIndividual(line[4]) + "\t1");
////					writer.println(offspringCnv.get(line[4]) + "\t" + sampleData.lookup(line[5])[2]  + "\t" + sampleData.lookup(line[6])[2] + "\t" + sampleData.getSexForIndividual(line[4]) + "\t1");
////				}
//			}
//			writer.close();
//			reader.close();
//
//		} catch (FileNotFoundException fnfe) {
//			return;
//		} catch (IOException ioe) {
//			return;
//		}
//	}

	public static void batch(String pedigreeOfTrio) {
		String command;
		String[][] iterations;
		
//		command = "perl C:/penncnv/detect_cnv.pl -joint -hmm C:/penncnv/lib/hh550.hmm -pfb C:/penncnv/lib/hh550.hg18.pfb -gcmodel C:/penncnv/lib/hh550.hg18.gcmodel [%1] [%2] [%3] -out [%1].jointcnv";
		command = "perl ../bin/penncnv/detect_cnv.pl -joint -hmm ../bin/penncnv/lib/hhall.hmm -pfb ../bin/custom.pfb -gcmodel ../bin/custom.gcmodel [%1] [%2] [%0] -out ../results/[%0].jointcnv -log ../results/[%0].log";
		
		iterations = HashVec.loadFileToStringMatrix(pedigreeOfTrio, true, new int[] {4, 5, 6}, false);
		
		common.Files.qsub("denovo", "/share/bulk/gedi/pankr018/denovo/penn_data", 65, command, iterations, 2500, 2);
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
	
	// TODO this whole method needs to be formed correctly before release
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
		String gcModelFullPath;
		int numWindowUnits;
		String pedigreeFullPath;
		String pennDataDir;
		String pennOutputDir;
		String pennOutputFilenameRoot;
		String pennBinDir;
		String trioPedigreeFullPath;
		Logger log;

//		projPropertyFileFullPath = "C:/workspace/Genvisis/projects/OSv2.properties";
		projPropertyFileFullPath = "/home/pankrat2/zxu/projects/gedi_gwas.properties";
		proj = new Project(projPropertyFileFullPath, false);
		pedigreeFullPath = proj.getDir(Project.DATA_DIRECTORY) + "pedigree.dat";
		trioPedigreeFullPath = ext.parseDirectoryOfFile(pedigreeFullPath) + "pedigreeTrios.dat";
		gcBaseFileFullPath = "/home/pankrat2/zxu/PennCNV_Related/GcBase/gc5Base_hg19.txt";
		numWindowUnits = 0;
		pennOutputFilenameRoot = "gedi_all.rawcnv";
		pennBinDir = "C:/penncnv/";
//		gcBaseFileFullPath = "D:/PennCNV_Related/GcBase/gc5Base_hg19.txt";

		String usage = "\n"+
				"cnv.analysis.DeNovoCNV requires 3 arguments\n"+
				"   (1) The project's property file's name (i.e. proj=" + projPropertyFileFullPath + " (not the default))\n"+
				"   (2) Direction of PennCNV software (i.e. pennBinDir=" + pennBinDir + " (not the default))\n"+
				"   (2) The output gcBase file name  (i.e. gcbase=" + gcBaseFileFullPath + " (not the default))\n"+
				"   (3) Number of window units for GC model (i.e. nwin=" + numWindowUnits + " (not the default))\n"+
//				"   (4) The output gc model file name (i.e. gcmodel=" + gcModelFileFullPath + " (not the default))\n"+
				"   (5) The output pedigree file name (i.e. pedigree=" + pedigreeFullPath + " (not the default))\n"+
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
			} else if (args[i].startsWith("pennBinDir=")) {
				pennBinDir = args[i].split("=")[1];
			} else if (args[i].startsWith("nwin=")) {
				numWindowUnits = Integer.parseInt(args[i].split("=")[1]);
			} else if (args[i].startsWith("gcbase=")) {
				gcBaseFileFullPath = args[i].split("=")[1];
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
		log = new Logger(proj.getProjectDir() + "Genvisis_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
		pedigreeFullPath = proj.getDir(Project.DATA_DIRECTORY) + "pedigree.dat";
		trioPedigreeFullPath = ext.parseDirectoryOfFile(pedigreeFullPath) + "pedigreeTrios.dat";
		gcModelFullPath = proj.getProjectDir() + "custom.gcmodel";
//		pedigreeFullPath = "N:/statgen/GEDI_GWAS/pedigree.dat";
//		trioPedigreeFullPath = "N:/statgen/GEDI_GWAS/pedigreeTrios.dat";
		pennDataDir = proj.getProjectDir() + "penn_data/";
		pennOutputDir = proj.getProjectDir() + "penn_output/";


//		generateTriosPedigree(pedigreeFullPath, trioPedigreeFullPath);
//		PennCNV.populationBAF(proj);
//		PennCNV.gcModel(proj, gcBaseFileFullPath, gcModelFullPath, numWindowUnits);
		runPennCnv(proj, trioPedigreeFullPath, pennDataDir, gcModelFullPath, pennOutputDir, pennOutputFilenameRoot, pennBinDir, -1, 10);
//		batch(pedigreeFullPath);
//		parsePennCnvResult(proj, pennOutputDir, trioPedigreeFullPath, "jointcnv", log);
//		parsePennCnvResult(proj, pennOutputDir, trioPedigreeFullPath, "triocnv", log);
//		generateTriosFamFile(proj, trioPedigreeFullPath, true, log);
//		generateDenovoList(proj, proj.getDir(Project.DATA_DIRECTORY) + "denovo_joint.cnv", true, log);

//		selectFilesToTestInPennCNV("D:/BOSS_outputFiles/TriosForDenovoCnv.txt");
//
//		detectDenovoCnv();
	}
	
}
