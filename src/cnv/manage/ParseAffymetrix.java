package cnv.manage;

import java.io.*;
import java.util.*;

import cnv.filesys.*;
import common.*;

public class ParseAffymetrix implements Runnable {
	public static final String[][] SNP_HEADER_OPTIONS = {{"SNP Name", "rsID"}};
	public static final String[] FIELDS = {"Sample ID", "Sample Name"};
	public static final String[][] SNP_TABLE_FIELDS = {{"Name"}, {"Chr", "Chromosome"}, {"Position"}};

	private Project proj;
	private String[] files;
	private String[] markerNames;
	private int[] keysKeys;
	private long fingerprint;

	public ParseAffymetrix(Project proj, String[] files, String[] markerNames, int[] keysKeys, long fingerprint, Hashtable<String,String> fixes) {
		this.proj = proj;
		this.files = files;
		this.markerNames = markerNames;
		this.keysKeys = keysKeys;
		this.fingerprint = fingerprint;
	}

	public void run() {
		BufferedReader reader;
		String[] line;
		int count, version, snpIndex, sampIndex, key;
		String trav;
		int[] dataIndices, genotypeIndices;
		boolean parseAtAt;

		FullSample samp;
		String sampleName;
		float[][] data;
		byte[][] genotypes;
		String idHeader;
		String delimiter;
		
//		try {
//			PrintWriter writer = new PrintWriter(new FileWriter(files[0]+"_list.xln"));
//			for (int j = 0; j<files.length; j++) {
//				writer.println(files[j]);
//            }
//	        writer.close();
//	        return;
//        } catch (Exception e) {
//	        System.err.println("Error writing to "+files[0]+"_list.xln");
//	        e.printStackTrace();
//        }

		idHeader = proj.getProperty(Project.ID_HEADER);
		delimiter = proj.getSourceFileDelimiter();
        try {
			for (int i = 0; i<files.length; i++) {
				try {
					System.out.println(ext.getTime()+"\t"+(i+1)+" of "+files.length);
					reader = new BufferedReader(new FileReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[i]));
					do {
						line = reader.readLine().trim().split(delimiter, -1);
					} while (reader.ready()&&(ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || ext.indexOfStr(idHeader, line)==-1));

					dataIndices = ext.indexFactors(FullSample.DATA_FIELDS, line, false, true, false, false);
					genotypeIndices = ext.indexFactors(FullSample.GENOTYPE_FIELDS, line, false, true, false, false);
					sampIndex = ext.indexFactors(new String[] {idHeader}, line, false, true)[0];
					snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, true)[0];
					
					if (dataIndices[3] == -1 || dataIndices[4] == -1) {
						System.err.println("Error - File format not consistent! At the very least the files need to contain "+Array.toStr(FullSample.DATA_FIELDS[3], "/")+" and "+Array.toStr(FullSample.DATA_FIELDS[4], "/"));
						return;
					}
					if (genotypeIndices[0] == -1 || genotypeIndices[1] == -1) {
						System.err.println("Error - File format not consistent! The files need to contain "+Array.toStr(FullSample.GENOTYPE_FIELDS[0], "/")+" and "+Array.toStr(FullSample.GENOTYPE_FIELDS[1], "/"));
						return;
					}

					sampleName = "";
					data = new float[FullSample.DATA_FIELDS.length][];
					for (int j = 0; j<data.length; j++) {
						if (dataIndices[j] != -1) {
							data[j] = new float[markerNames.length];
						}
                    }
					genotypes = new byte[2][];
					genotypes[0] = Array.byteArray(markerNames.length, Byte.MIN_VALUE);
					genotypes[1] = Array.byteArray(markerNames.length, Byte.MIN_VALUE);
					
					count = 0;
					parseAtAt = Boolean.parseBoolean(proj.getProperty(Project.PARSE_AT_AT_SYMBOL));
					while (reader.ready()) {
						line = reader.readLine().split(delimiter);
						if (parseAtAt&&line[sampIndex].indexOf("@")==-1) {
							System.err.println("Error - "+idHeader+" '"+line[sampIndex]+"' did not contain an @ sample");
							parseAtAt = false;
						}
						trav = parseAtAt?line[sampIndex].substring(0, line[sampIndex].indexOf("@")):line[sampIndex];

						if (count==0) {
							sampleName = trav;
						} else if (!trav.equals(sampleName)) {
							System.err.println("Found more than one ID in file "+files[i]+"(found "+trav+", expecting "+sampleName+")");
							return;
						} else if (count>markerNames.length) {
							System.err.println("Error - expecting only "+markerNames.length+" markers and found more than that in file "+files[i]);
							return;
						} else if (!markerNames[count].equals(line[snpIndex]) && !markerNames[count].startsWith("Blank")) {
							System.err.println("Found "+line[snpIndex]+" at marker #"+(count+1)+" in file "+files[i]+"; expecting "+markerNames[count]);
							return;
						}

						key = keysKeys[count];
						// System.out.println(count+"\t"+markerNames[count]+"\t"+key);
						for (int j = 0; j<FullSample.DATA_FIELDS.length; j++) {
							try {
								if (dataIndices[j] != -1) {
									data[j][key] = Float.parseFloat(line[dataIndices[j]]);
								}
							} catch (NumberFormatException nfe) {
								System.err.println("Error - failed to parse"+line[dataIndices[j]]+" into a valid "+FullSample.DATA_FIELDS[j]);
								return;
							}
						}

						genotypes[0][key] = (byte)ext.indexOfStr(line[genotypeIndices[0]]+line[genotypeIndices[1]], FullSample.ALLELE_PAIRS);
						if (genotypes[0][key] == -1) {
							if (ext.indexOfStr(line[genotypeIndices[0]]+line[genotypeIndices[1]], FullSample.ALT_NULLS) == -1) {
								System.err.println("Error - failed to lookup "+line[genotypeIndices[0]]+line[genotypeIndices[1]]+" for marker "+markerNames[count]+" of sample "+files[i]);
							} else {
								genotypes[0][key] = 0;
							}								
						}
						genotypes[1][key] = (byte)ext.indexOfStr(line[genotypeIndices[2]]+line[genotypeIndices[3]], FullSample.AB_PAIRS);

						count++;
					}
					reader.close();
					if (count!=markerNames.length) {
						System.err.println("Error - expecting "+markerNames.length+" markers and only found "+count+" in file "+files[i]);
						return;
					}

					version = 0;
					do {
						trav = sampleName+(version==0?"":"."+version);
						version++;
					} while (new File(proj.getDir(Project.SAMPLE_DIRECTORY, true)+trav+".fsamp").exists());

					samp = new FullSample(sampleName, fingerprint, data, genotypes);
					samp.serialize(proj.getDir(Project.SAMPLE_DIRECTORY, true)+trav+".fsamp");
					samp.convertToSample().serialize(proj.getDir(Project.IND_DIRECTORY, true)+trav+".samp");
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+files[i]+"\" not found in current directory");
					return;
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+files[i]+"\"");
					return;
				}
			}

			SampleList.generateSampleList(proj).writeToTextFile(proj.getProjectDir()+"ListOfSamples.txt");

			System.out.println(ext.getTime()+"\tfinished");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void createFiles(Project proj, int numThreads) {
		BufferedReader reader;
		String[] line, markerNames, files;
		ArrayList<String> alNames;
		int snpIndex;
		int[] keys, keysKeys, indices;
		long fingerprint;
		Thread[] threads;
		Vector<Vector<String>> fileCabinet;
		String trav;
		int count;
        Hashtable<String,String> fixes;
        String idHeader, delimiter;

		if (!proj.getDir(Project.SOURCE_DIRECTORY).equals("")&&!new File(proj.getDir(Project.SOURCE_DIRECTORY)).exists()) {
			System.err.println("Error - the Project source location is invalid: "+proj.getDir(Project.SOURCE_DIRECTORY));
			return;
		}
        
		if (!new File(proj.getFilename(Project.MARKER_POSITION_FILENAME, false, false)).exists()) {
			System.err.println("Error - missing markerPositions: "+proj.getFilename(Project.MARKER_POSITION_FILENAME, false, false));
			return;
		}

		delimiter = proj.getSourceFileDelimiter();
		idHeader = proj.getProperty(Project.ID_HEADER);
		System.out.println(ext.getTime()+"\tSearching for "+proj.getProperty(Project.SOURCE_FILENAME_EXTENSION)+" files in: "+proj.getDir(Project.SOURCE_DIRECTORY));
		files = Files.list(proj.getDir(Project.SOURCE_DIRECTORY), proj.getProperty(Project.SOURCE_FILENAME_EXTENSION), false);
		
		if (files.length == 0) {
			System.err.println("Error - no files to parse");
			return;
		}

		System.out.println("\t\tFound "+files.length+" file"+(files.length==1?"":"s")+" to parse");
		fixes = new Hashtable<String,String>();
		if (new File(proj.getProjectDir()+"fixes.dat").exists()) {
			System.out.println("Also found a 'fixes.dat' file in the project directory, which will be used to rename samples");
			fixes = HashVec.loadFileToHashString(proj.getProjectDir()+"fixes.dat", false);
		} else {
			System.out.println("Did not find a 'fixes.dat' file; assuming you don't want to rename any IDs");
		}
		alNames = new ArrayList<String>();
		try {
			reader = new BufferedReader(new FileReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[0]));
			do {
				line = reader.readLine().trim().split(delimiter, -1);
			} while (reader.ready()&&(ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || ext.indexOfStr(idHeader, line)==-1));

			// check immediately to make sure these fields are valid
			indices = ext.indexFactors(FullSample.DATA_FIELDS, line, false, true, true, false); // dataIndices
			if (indices[3] == -1 || indices[4] == -1) {
				System.err.println("Error - at the very least the files need to contain "+Array.toStr(FullSample.DATA_FIELDS[3], "/")+" and "+Array.toStr(FullSample.DATA_FIELDS[4], "/"));
				System.err.println("      - failed to see that in "+files[0]);
				System.err.println(Array.toStr(line));
				
				return;
			}
			indices = ext.indexFactors(FullSample.GENOTYPE_FIELDS, line, false, true, true, false); // genotypeIndices
			if (indices[0] == -1 || indices[1] == -1) {
				System.err.println("Error - the files need to contain "+Array.toStr(FullSample.GENOTYPE_FIELDS[0], "/")+" and "+Array.toStr(FullSample.GENOTYPE_FIELDS[1], "/"));
				return;
			}
			
			ext.indexFactors(new String[] {idHeader}, line, false, true); // sampIndex
			snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, true, true)[0];

			if (Boolean.parseBoolean(proj.getProperty(Project.LONG_FORMAT))) {
//				createFilesFromLongFormat(proj, files, idHeader, fixes, delimiter);
//				return;
			}

			count = 1;
			while (reader.ready()) {
				trav = reader.readLine().trim().split(delimiter)[snpIndex];
				if (trav.equals("") || trav.equals("0")) {
					trav = "Blank"+count;
					count++;
				}
				alNames.add(trav);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+proj.getDir(Project.SOURCE_DIRECTORY)+files[0]+"\" not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+proj.getDir(Project.SOURCE_DIRECTORY)+files[0]+"\"");
			return;
		}

//		new File(proj.getDir(Project.SAMPLE_DIRECTORY, true)).mkdirs();
//		new File(proj.getDir(Project.IND_DIRECTORY)).mkdirs();
//		new File(proj.getDir(Project.DATA_DIRECTORY)).mkdirs();

		markerNames = Array.toStringArray(alNames);
		keys = Markers.orderMarkers(markerNames, proj.getFilename(Project.MARKER_POSITION_FILENAME), proj.getFilename(Project.MARKERSET_FILENAME, true, true));
		if (keys == null) {
			return;
		}
		keysKeys = Sort.quicksort(keys); // very important
		fingerprint = proj.getMarkerSet().getFingerprint();
		System.out.println("There are "+markerNames.length+" markers being processed (fingerprint: "+fingerprint+")");
		
		
		fileCabinet = new Vector<Vector<String>>();
		for (int i = 0; i<numThreads; i++) {
			fileCabinet.add(new Vector<String>());
		}
		for (int i = 0; i<files.length; i++) {
			fileCabinet.elementAt(i%numThreads).add(files[i]);
		}
		threads = new Thread[numThreads];
		for (int i = 0; i<numThreads; i++) {
			threads[i] = new Thread(new ParseAffymetrix(proj, fileCabinet.elementAt(i).toArray(new String[fileCabinet.elementAt(i).size()]), markerNames, keysKeys, fingerprint, fixes));
//			threads[i] = new Thread(new ParseIllumina(proj, fileCabinet.elementAt(i).toArray(new String[fileCabinet.elementAt(i).size()]), null, null, null, 0));
			threads[i].start();
			try {
				Thread.sleep(100L);
			} catch (InterruptedException ex) {}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj;
		String filename = Project.DEFAULT_PROJECT;
		int numThreads = 1;

		String usage = "\n"+
		"cnv.manage.ParseAffymetrix requires 0-2 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) number of threads to use (i.e. threads="+numThreads+" (default))\n"+
		"";
		
		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				return;
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("threads=")) {
				numThreads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			return;
		}

		proj = new Project(filename, false);
		
		try {
			createFiles(proj, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
