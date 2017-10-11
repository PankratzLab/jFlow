package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.genvisis.cnv.manage.PlinkMergePrep;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.DosageData;
import org.genvisis.filesys.SnpMarkerSet;

public class MergeExtractPipeline {

	public static class DataSource {
		public String label;
		public String dataFile;
		public String mapFile;
		public String idFile;
		public boolean relD, relM, relI;

		public DataSource(String label, String dir, String d, String m, String i) {
			relD = Files.isRelativePath(d);
			relM = Files.isRelativePath(m);
			relI = Files.isRelativePath(i);
			dataFile = (relD ? dir : "") + d;
			mapFile = null == m || ".".equals(m) ? null : (relM ? dir : "") + m;
			idFile = null == i || ".".equals(i) ? null : (relI ? dir : "") + i;
			this.label = label;
		}
	}

	ArrayList<DataSource> dataSources;
	private String outFileD = null;
	private String outFileM = null;
	int outFormat = -1;
	String markersFile = null;
	String[] markers = null;
	int[][] markerLocations;
	String runDir = null;
	Logger log = null;
	String regionsFile = null;
	int[][] regions;
	String[] regionLabels;

	boolean permissiveLookup = true;
	boolean splitOutput = false;
	boolean overwrite = false;
	boolean renameMarkers = true;
	boolean bestGuessOutput = false;
	float bestGuessThreshold = 0;

	public MergeExtractPipeline() {
		dataSources = new ArrayList<MergeExtractPipeline.DataSource>();
	}

	private void initLog() {
		if (log == null) {
			log = new Logger();
		}
	}

	public MergeExtractPipeline setLogger(Logger log) {
		this.log = log;
		return this;
	}

	/**
	 * Either absolute path to file, or relative path to file in either previously-specified
	 * run-directory or current directory (depending on where the file exists)
	 *
	 * @param regionFile
	 * @return
	 */
	public MergeExtractPipeline setRegions(String regionFile) {
		BufferedReader reader;
		String line;
		String[] temp;
		boolean hasLabels = false;
		ArrayList<int[]> reg;
		ArrayList<String> lbl;

		initLog();

		if (markersFile != null) {
			String msg = "Error - cannot specify both markers-to-extract and regions-to-extract!";
			log.reportError(msg);
			throw new IllegalArgumentException(msg);
		}

		lbl = new ArrayList<String>();
		reg = new ArrayList<int[]>();
		regionsFile = regionFile;
		if (Files.isRelativePath(regionFile)) {
			if (runDir != null && !"".equals(runDir)) {
				if (Files.exists(runDir + regionFile)) {
					regionsFile = runDir + regionFile;
				}
			} else {
				if (Files.exists("./" + regionFile)) {
					regionsFile = "./" + regionFile;
				}
			}
		}
		try {
			reader = Files.getAppropriateReader(regionsFile);
			line = reader.readLine();
			if (line != null) {
				temp = line.split("\t");
				if (temp.length > 1) {
					hasLabels = true;
				}
				reg.add(Positions.parseUCSClocation(temp[0]));
				if (hasLabels) {
					lbl.add(temp[1]);
				} else {
					lbl.add(ext.replaceWithLinuxSafeCharacters(temp[0], true));
				}
			}
			while ((line = reader.readLine()) != null) {
				temp = line.split("\t");
				reg.add(Positions.parseUCSClocation(temp[0]));
				if (hasLabels) {
					lbl.add(temp[1]);
				} else {
					lbl.add(ext.replaceWithLinuxSafeCharacters(temp[0], true));
				}
			}
			reader.close();
		} catch (IOException e) {
			throw new IllegalArgumentException(e);
		}

		regions = reg.toArray(new int[reg.size()][]);
		regionLabels = lbl.toArray(new String[lbl.size()]);

		log.report("Set regions-to-extract file: " + regionsFile);

		return this;
	}

	public MergeExtractPipeline setMarkers(String mkrsFile) {
		initLog();

		if (regionsFile != null) {
			String msg = "Error - cannot specify both markers-to-extract and regions-to-extract!";
			log.reportError(msg);
			throw new IllegalArgumentException(msg);
		}

		markersFile = mkrsFile;

		if (Files.isRelativePath(mkrsFile)) {
			if (runDir != null && !"".equals(runDir)) {
				if (Files.exists(runDir + mkrsFile)) {
					markersFile = runDir + mkrsFile;
				}
			} else {
				if (Files.exists("./" + mkrsFile)) {
					markersFile = "./" + mkrsFile;
				}
			}
		}
		markers = HashVec.loadFileToStringArray(markersFile, false, new int[] {0}, true, false, "\t");
		SnpMarkerSet markerSet = new SnpMarkerSet(markers, false, log);
		markerSet.parseSNPlocations(log);
		markerLocations = markerSet.getChrAndPositionsAsInts();
		for (int[] markerLocation : markerLocations) {
			if (markerLocation[0] == 0 || markerLocation[1] == -1) {
				log.reportError("Error - one or more markers could not be located.  Merging/Extracting may fail for this reason.");
				break;
			}
		}
		log.report("Set markers-to-extract file: " + markersFile);
		return this;
	}

	public MergeExtractPipeline setRenameMarkers(boolean rename) {
		renameMarkers = rename;
		return this;
	}

	public MergeExtractPipeline setPermissiveMarkerFiltering(boolean permissiveLookup) {
		this.permissiveLookup = permissiveLookup;
		return this;
	}

	public MergeExtractPipeline setBestGuessOutput(boolean bestGuess) {
		bestGuessOutput = bestGuess;
		return this;
	}

	public MergeExtractPipeline setBestGuessThreshold(float bestThresh) {
		bestGuessThreshold = bestThresh;
		return this;
	}

	public MergeExtractPipeline setRunDirectory(String runDir, boolean create) {
		initLog();
		this.runDir = runDir;
		if (!Files.exists(runDir)) {
			if (create) {
				(new File(runDir)).mkdirs();
			} else {
				throw new IllegalArgumentException(
																					 "Error - specified run directory \""
																							 + runDir
																							 + "\" doesn't exist, and create flag wasn't set.  Please fix or create directory and try again.");
			}
		}
		return this;
	}

	public MergeExtractPipeline setOutputFiles(String dataFile, String mapFile) {
		if (dataFile == null && mapFile == null) {
			if (markersFile != null && !"".equals(markersFile) && Files.exists(markersFile)) {
				outFileD = ext.rootOf(markersFile, false) + ".db.xln.gz";
				outFileM = ext.rootOf(markersFile, false) + ".snp";
			} else if (regionsFile != null && !"".equals(regionsFile) && Files.exists(regionsFile)) {
				outFileD = ext.rootOf(regionsFile, false) + ".db.xln.gz";
				outFileM = ext.rootOf(regionsFile, false) + ".snp";
			} else {
				outFileD = (this.runDir == null ? "./" : this.runDir)
									 + Files.getNextAvailableFilename("mergeExtract.db.xln.gz");
				outFileM = ext.rootOf(outFileD, false) + ".snp";
			}
		} else {
			if (dataFile == null && mapFile != null && !"".equals(mapFile)) {
				outFileD = ext.rootOf(mapFile, false) + ".db.xln.gz";
				outFileM = mapFile;
			} else if (mapFile == null && dataFile != null && !"".equals(dataFile)) {
				outFileD = dataFile;
				outFileM = ext.rootOf(dataFile, false) + ".snps";
			} else {
				outFileD = dataFile;
				outFileM = mapFile;
			}
		}
		return this;
	}

	public MergeExtractPipeline setSplitOutputByRegions(boolean split) {
		splitOutput = split;
		return this;
	}

	/**
	 * Refers to DosageData formats
	 *
	 * @param format
	 * @return
	 */
	public MergeExtractPipeline setOutputFormat(int format) {
		outFormat = format;
		return this;
	}

	public MergeExtractPipeline setOverwrite() {
		overwrite = true;
		return this;
	}

	public MergeExtractPipeline addDataSource(DataSource ds) {
		dataSources.add(ds);
		return this;
	}

	public MergeExtractPipeline addDataSource(String lbl, String dir, String dataFile,
																						String mapFile,
																						String idFile) {
		dataSources.add(new DataSource(lbl, dir, dataFile, mapFile, idFile));
		return this;
	}

	private static FilenameFilter getFilter(final int[][] markerLocations, final int[][] regions,
																					final String dataFileExt, final int bpWindow,
																					final Logger log) {
		return new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				boolean keep = false;
				if (name.endsWith(dataFileExt)) {
					if (markerLocations == null && regions == null) {
						keep = true;
					} else {
						Matcher m = Pattern.compile(DosageData.CHR_REGEX).matcher(name);
						if (m.matches()) {
							byte chr = Positions.chromosomeNumber(m.group(1), false, log);
							int bpStart = -1, bpEnd = -1;
							String[] testPos = name.split("\\.");
							for (int i = 0; i < testPos.length; i++) {
								if (testPos[i].startsWith("chr")) {
									boolean logMsg = false;
									try {
										if (testPos.length > (i + 2)) {
											bpStart = Integer.parseInt(testPos[i + 1]);
											bpEnd = Integer.parseInt(testPos[i + 2]);
										} else {
											logMsg = true;
										}
									} catch (NumberFormatException e) {
										logMsg = true;
									}
									if (logMsg) {
										log.reportError("not enough information in file name to determine bpStart/bpEnd from filename \""
																		+ name
																		+ "\" - whole chromosomes will be included in data merge.");
									}
									break;
								}
							}
							if (markerLocations != null) {
								for (int[] markerLocation : markerLocations) {
									if (markerLocation[0] == chr) {
										if (bpStart != -1 && bpEnd != -1) {
											if (bpStart <= markerLocation[1] - bpWindow
													&& bpEnd >= markerLocation[1] + bpWindow) {
												keep = true;
												break;
											}
										} else {
											keep = true;
											break;
										}
									}
								}
							}
							if (!keep && regions != null) {
								for (int[] region : regions) {
									if (region[0] == chr) {
										if (bpStart != -1 && bpEnd != -1) {
											if (bpStart <= region[1] - bpWindow && bpEnd >= region[1] + bpWindow) {
												keep = true;
												break;
											}
										} else {
											keep = true;
											break;
										}
									}
								}
							}
						} else {
							keep = true;
						}
					}
				}
				return keep;
			}
		};
	}

	public MergeExtractPipeline addDataSources(final String lbl, final String dir,
																						 final String dataFileExt, final String mapFileExt,
																						 final String idFile, final int bpWindow) {
		initLog();
		String[] filesToAdd = (new File(dir)).list(getFilter(markerLocations, regions, dataFileExt,
																												 bpWindow, log));
		for (String file : filesToAdd) {
			addDataSource(lbl, dir, file,
										file.substring(0, file.length() - dataFileExt.length()) + mapFileExt, idFile);
		}
		return this;
	}

	private String getOutputDataFile() {
		return (Files.isRelativePath(outFileD) ? (runDir == null ? "./" : runDir) : "") + outFileD;
	}

	private String getOutputDataFile(String subDir) {
		return ext.verifyDirFormat((Files.isRelativePath(outFileD) ? (runDir == null ? "./" : runDir)
																															: "")
															 + subDir)
					 + outFileD;
	}

	private String getOutputMapFile() {
		return (Files.isRelativePath(outFileM) ? (runDir == null ? "./" : runDir) : "") + outFileM;
	}

	private String getOutputMapFile(String subDir) {
		return ext.verifyDirFormat((Files.isRelativePath(outFileM) ? (runDir == null ? "./" : runDir)
																															: "")
															 + subDir)
					 + outFileM;
	}

	public void run() {
		initLog();
		if (dataSources.size() == 0) {
			log.reportError("no data sources specified!");
			return;
		}
		if (getOutputDataFile() == null || "".equals(getOutputDataFile())) {
			log.reportError("output data file name not specified!");
			return;
		}
		if (Files.exists(getOutputDataFile()) && !overwrite) {
			log.reportError("output data file already exists!  Please remove file or set -overwrite flag, and re-run.");
			return;
		}
		if (getOutputMapFile() == null || "".equals(getOutputMapFile())) {
			log.reportError("output map file name not specified!");
			return;
		}
		if (Files.exists(getOutputMapFile()) && !overwrite) {
			log.reportError("output map file already exists!  Please remove file or set -overwrite flag, and re-run.");
			return;
		}
		if (regionsFile == null) {
			log.reportTimeWarning("no regions file specified; all data will be exported.");
		} else if (!Files.exists(regionsFile)) {
			log.reportError("specified regions file not found; please check file path and try again.");
			return;
		}
		if (markersFile == null) {
			log.reportTimeWarning("no markers file specified; all markers will be exported.");
		} else if (!Files.exists(markersFile)) {
			log.reportError("specified markers file not found; please check file path and try again.");
			return;
		}

		if (!checkMemory()) {
			return;
		}

		String dir = "./";
		if (runDir == null) {
			log.reportTimeWarning("no run directory specified, merging will take place in the current directory.");
		} else {
			dir = runDir;
		}
		dir = ext.verifyDirFormat((new File(dir)).getAbsolutePath());

		// // all relative paths are already absolute
		// for (int i = 0; i < dataSources.size(); i++) {
		// if (Files.isRelativePath(dataSources.get(i).dataFile)) {
		// dataSources.get(i).dataFile = dir + dataSources.get(i).dataFile;
		// }
		// if (Files.isRelativePath(dataSources.get(i).mapFile)) {
		// dataSources.get(i).mapFile = dir + dataSources.get(i).mapFile;
		// }
		// if (Files.isRelativePath(dataSources.get(i).idFile)) {
		// dataSources.get(i).idFile = dir + dataSources.get(i).idFile;
		// }
		// }

		ArrayList<String> tempMarkerFiles = new ArrayList<>();
		if (permissiveLookup && markers != null) {
			HashSet<String> allMarkersLookup = new HashSet<>();
			for (int i = 0; i < markers.length; i++) {
				allMarkersLookup.add(markers[i]);
				if (!markers[i].contains(":")) {
					String temp = markerLocations[i][0] + ":" + markerLocations[i][1];
					allMarkersLookup.add("chr" + temp);
					allMarkersLookup.add(temp);
				}
			}
			if (allMarkersLookup.size() > markers.length) {
				String tempFile = "./tempMarkers." + ext.replaceWithLinuxSafeCharacters(ext.getDate())
													+ ".snp";

				log.reportTime("Permissive Lookup flagged, creating a temporary marker lookup file at "
											 + tempFile + ".  This will be removed after completion.");

				Files.writeIterable(allMarkersLookup, tempFile);
				tempMarkerFiles.add(tempFile);
				setMarkers(tempFile);
			}
		}

		ArrayList<String> plinkRoots = new ArrayList<String>();
		ArrayList<String> plinkLabels = new ArrayList<String>();
		// discover and merge all plink files
		for (int i = dataSources.size() - 1; i >= 0; i--) {
			if (dataSources.get(i).dataFile.endsWith(".bed")
					|| dataSources.get(i).dataFile.endsWith(".ped")) {
				plinkRoots.add(ext.rootOf(dataSources.get(i).dataFile, false));
				plinkLabels.add(dataSources.get(i).label);
				dataSources.remove(i);
			}
		}
		if (plinkRoots.size() > 0) {
			String[] roots = ArrayUtils.toStringArray(plinkRoots);
			String[] lbls = ArrayUtils.toStringArray(plinkLabels);
			plinkRoots = null;
			String outRoot = dir + "plink_merged";
			String mergeCommand = PlinkMergePrep.merge(PlinkMergePrep.BEDBIMFAM, outRoot, overwrite,
																								 renameMarkers, regionsFile, markersFile, roots,
																								 lbls);

			log.report("Running PLINK merge command: " + mergeCommand);
			boolean result = CmdLine.runDefaults(mergeCommand, dir);
			if (!result) {
				log.reportError("PLINK merge command failed.  Please check output for errors and try again.");
				return;
			}
			dataSources.add(new DataSource(null, dir, PSF.Plink.getBED(outRoot),
																		 PSF.Plink.getBIM(outRoot), PSF.Plink.getFAM(outRoot))); // no
																																														 // prepend
																																														 // here,
																																														 // as
																																														 // we've
																																														 // already
																																														 // renamed
																																														 // the
																																														 // markers
																																														 // using
																																														 // PlinkMergePrep
			if (markersFile != null && !"".equals(markersFile)) {
				String newMkrFile = (new File(PlinkMergePrep.TEMP_MKR_FILE)).getAbsolutePath();
				tempMarkerFiles.add(newMkrFile);
				log.report("Setting markers file to temporarily generated plink-renamed markers file: "
									 + newMkrFile);
				markersFile = newMkrFile;
			}
		}

		System.gc();
		boolean verbose = false;

		log.reportTime("Merging data...");
		log.reportTime("Starting from data file: " + dataSources.get(0).dataFile);

		DosageData dd1 = new DosageData(dataSources.get(0).dataFile, dataSources.get(0).idFile,
																		dataSources.get(0).mapFile, regionsFile, markersFile,
																		renameMarkers ? dataSources.get(0).label : "", verbose, log);
		HashMap<String, HashMap<String, Annotation>> annotations = getAnnotations(dd1,
																																							getAnnotationLabels(dataSources.get(0).mapFile,
																																																	dataSources.get(0).label));
		for (int i = 1; i < dataSources.size(); i++) {
			System.gc();
			log.reportTime("... merging with data file: " + dataSources.get(i).dataFile);
			DosageData dd2 = new DosageData(dataSources.get(i).dataFile, dataSources.get(i).idFile,
																			dataSources.get(i).mapFile, regionsFile, markersFile,
																			renameMarkers ? dataSources.get(i).label : "", verbose, log);
			if (!dd2.isEmpty()) {
				String[] dd2Annots = getAnnotationLabels(dataSources.get(i).mapFile,
																								 dataSources.get(i).label);
				if (dd2Annots != null && dd2Annots.length > 0
						&& dd2.getMarkerSet().getAnnotation() != null) {
					combineAnnotations(annotations, dd2Annots, dd2.getMarkerSet().getAnnotation(),
														 dd2.getMarkerSet().getMarkerNames());
				}
				dd1 = DosageData.combine(dd1, dd2, DosageData.COMBINE_OP.EITHER_IF_OTHER_MISSING,
																 bestGuessOutput, bestGuessThreshold, log);
			}
			dd2 = null;
			System.gc();
		}

		log.reportTime("Writing final merged data to output file: " + getOutputDataFile());

		if (outFormat == -1) {
			outFormat = DosageData.determineType(getOutputDataFile());
		}

		String[] allMarkers = markers;
		if (renameMarkers && markers != null) {
			allMarkers = dd1.getMarkerSet().getMarkerNames();
		}

		if (splitOutput) {
			String outD, outM;
			for (int i = 0; i < regions.length; i++) {
				outD = getOutputDataFile(regionLabels[i]);
				outM = getOutputMapFile(regionLabels[i]);
				(new File(ext.parseDirectoryOfFile(outD))).mkdirs();
				(new File(ext.parseDirectoryOfFile(outM))).mkdirs();
				dd1.writeToFile(outD, outM, allMarkers, new int[][] {regions[i]}, true, true,
												DosageData.PARAMETERS[outFormat], bestGuessOutput, bestGuessThreshold, log);
			}
		} else {
			dd1.writeToFile(getOutputDataFile(), getOutputMapFile(), allMarkers, regions, true, true,
											DosageData.PARAMETERS[outFormat], bestGuessOutput, bestGuessThreshold, log);
		}

		writeAnnotations(allMarkers, annotations);

		for (String f : tempMarkerFiles) {
			new File(f).delete();
		}

		System.gc();
	}

	private void combineAnnotations(HashMap<String, HashMap<String, Annotation>> mkrAnnotations,
																	String[] annotLabels, String[][] annotation,
																	String[] markerNames) {
		for (int i = 0; i < markerNames.length; i++) {
			String snp = markerNames[i];
			String[] annArr = annotation[i];
			HashMap<String, Annotation> annMap = mkrAnnotations.get(snp);
			if (annMap == null) {
				annMap = new HashMap<String, MergeExtractPipeline.Annotation>();
				mkrAnnotations.put(snp, annMap);
			}
			System.out.println("Combining annotations for " + snp + "; existing: " + annMap.keySet()
												 + "; adding: " + ArrayUtils.toStr(annotLabels, ","));
			for (int j = 0; j < annotLabels.length; j++) {
				Annotation ann = annMap.get(annotLabels[j]);
				if (ann == null) {
					ann = new Annotation();
					ann.snp = snp;
					ann.annotationName = annotLabels[j];
					ann.annotation = annArr[j];
					annMap.put(annotLabels[j], ann);
				} else {
					if (ann.annotation == null || ext.isMissingValue(ann.annotation)) {
						ann.annotation = annArr[j];
					} else if (!ann.annotation.equals(annArr[j])) {
						System.err.println("Error - DUPLICATE ANNOTATION: " + ann.annotation + " | "
															 + annArr[j]);
					}
				}
			}
		}


	}

	private void writeAnnotations(String[] allMarkers,
																HashMap<String, HashMap<String, Annotation>> annotations) {
		String file = getOutputMapFile() + ".annot";
		PrintWriter writer = Files.getAppropriateWriter(file);

		HashSet<String> labelSet = new HashSet<String>();
		for (HashMap<String, Annotation> m : annotations.values()) {
			labelSet.addAll(m.keySet());
		}
		ArrayList<String> allLabels = new ArrayList<String>(labelSet);

		StringBuilder sb = new StringBuilder("SNP");
		for (String s : allLabels) {
			sb.append("\t").append(s);
		}

		writer.println(sb.toString());
		for (String s : allMarkers) {
			if (s == null || "".equals(s))
				continue;
			String[] annots = new String[allLabels.size()];
			HashMap<String, Annotation> ann = annotations.get(s);
			if (ann == null) {
				annots = ArrayUtils.stringArray(allLabels.size(), ".");
			} else {
				for (int i = 0; i < allLabels.size(); i++) {
					if (ann.get(allLabels.get(i)) == null) {
						annots[i] = ".";
					} else {
						annots[i] = ann.get(allLabels.get(i)).annotation;
					}
				}
			}
			writer.println(s + "\t" + ArrayUtils.toStr(annots, "\t"));
		}
		writer.flush();
		writer.close();

	}

	private String[] getAnnotationLabels(String mapFile, String prepend) {
		int type = SnpMarkerSet.determineType(mapFile);
		if (type == -1) {
			return new String[0];
		}
		String[] hdr = SnpMarkerSet.HEADERS[type];
		if (hdr == null || hdr.length <= 6) {
			return new String[0];
		}
		int[] indices = SnpMarkerSet.INDICES[type];
		String[] newHdr = new String[indices.length - 6];
		for (int i = 6; i < indices.length; i++) {
			newHdr[i - 6] = (prepend != null && !"".equals(prepend) ? prepend + "_" : "")
											+ hdr[indices[i]];
		}
		return newHdr;
	}

	static class Annotation {
		String snp;
		String annotationName;
		String annotation;

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((annotationName == null) ? 0 : annotationName.hashCode());
			result = prime * result + ((snp == null) ? 0 : snp.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Annotation other = (Annotation) obj;
			if (annotationName == null) {
				if (other.annotationName != null)
					return false;
			} else if (!annotationName.equals(other.annotationName))
				return false;
			if (snp == null) {
				if (other.snp != null)
					return false;
			} else if (!snp.equals(other.snp))
				return false;
			return true;
		}

	}

	private HashMap<String, HashMap<String, Annotation>> getAnnotations(DosageData dd,
																																			String[] annotationLabels) {
		HashMap<String, HashMap<String, Annotation>> annotations = new HashMap<String, HashMap<String, Annotation>>();
		String[][] ann = dd.getMarkerSet().getAnnotation();
		if (ann != null) {
			for (int i = 0; i < ann.length; i++) {
				String[] annots = dd.getMarkerSet().getAnnotation()[i];
				for (int a = 0; a < annots.length; a++) {
					Annotation ann1 = new Annotation();
					ann1.snp = dd.getMarkerSet().getMarkerNames()[i];
					ann1.annotationName = annotationLabels[a];
					ann1.annotation = annots[a];
					HashMap<String, Annotation> annMap = annotations.get(ann1.snp);
					if (annMap == null) {
						annMap = new HashMap<String, MergeExtractPipeline.Annotation>();
						annotations.put(ann1.snp, annMap);
					}
					annMap.put(ann1.annotationName, ann1);
				}
			}
		}
		return annotations;
	}

	private boolean checkMemory() {
		initLog();
		log.reportError("Warning - memory check not implemented yet - if not enough memory is provided, an Out Of Memory exception may occur.");
		return true;
	}

	public static ArrayList<DataSource> parseDataFile(String runDir, int[][] markerLocations,
																										int[][] regions, String data, int bpWindow,
																										Logger log) {
		BufferedReader reader;
		String line, file;
		String[] temp;

		if (null == data || "".equals(data)) {
			throw new IllegalArgumentException("Error - provided data file \"" + data
																				 + "\" doesn't exist.");
		}

		file = Files.isRelativePath(data) ? (Files.exists(runDir + data) ? runDir + data : "./" + data)
																		 : data;
		if (!Files.exists(file)) {
			throw new IllegalArgumentException("Error - provided data file \"" + file
																				 + "\" doesn't exist.");
		}

		ArrayList<DataSource> sources = new ArrayList<MergeExtractPipeline.DataSource>();
		try {
			reader = Files.getAppropriateReader(file);
			while ((line = reader.readLine()) != null) {
				if ("".equals(line)) {
					continue;
				}
				temp = line.split(ext.determineDelimiter(line));
				if (temp.length == 4) {
					log.report("Added data source: " + temp[1]);
					sources.add(new DataSource(temp[0], null, temp[1], temp[2], temp[3]));
				} else if (temp.length == 5) {
					String lbl = temp[0];
					String dir = temp[1];
					String dataFileExt = temp[2];
					String mapFileExt = temp[3];
					String idFile = temp[4];
					FilenameFilter ff = getFilter(markerLocations, regions, dataFileExt, bpWindow, log);
					// the program will die with only a "null" in the error log if the user does not have
					// permissions to access dir
					String[] filesToAdd = (new File(dir)).list(ff);
					log.report("Found " + filesToAdd.length + " files to add from " + dir);
					for (String fileToAdd : filesToAdd) {
						sources.add(new DataSource(lbl, dir, fileToAdd,
																			 fileToAdd.substring(0,
																													 fileToAdd.length()
																															 - dataFileExt.length())
																					 + mapFileExt,
																			 idFile));
						log.report("Added data source: " + fileToAdd);
					}
				} else {
					log.reportError("Error - skipping invalid entry in data file: " + line);
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e);
		}

		return sources;
	}

	// MergeExtractPipeline pipeline = new MergeExtractPipeline();
	// pipeline.setMarkers(markersFile);
	// pipeline.setRunDirectory("/scratch.global/cole0482/merge/", true);
	// pipeline.setOutputFormat(DosageData.DATABASE_DOSE_FORMAT);
	// pipeline.setOutputFiles(outFile, mapOutFile);
	// pipeline.setRenamePlinkMarkers(true);
	// pipeline.addDataSource("/scratch.global/cole0482/merge/blacks/", "gwas.bed", "gwas.bim",
	// "gwas.fam");
	// pipeline.addDataSource("/scratch.global/cole0482/merge/blacks/", "exome.bed", "exome.bim",
	// "exome.fam");
	// pipeline.addDataSource("/scratch.global/cole0482/merge/blacks/", "metab.bed", "metab.bim",
	// "metab.fam");
	// add more data;
	// pipeline.run();
	public static void main(String[] args) {
		int numArgs = args.length;
		String outfileD = null;
		String outfileM = null;
		String data = "./dataSources.txt";
		String markers = null;
		String regions = null;
		String rundir = "./";
		String logFile = null;
		boolean split = false;
		boolean overwrite = false;
		boolean rename = true;
		float bestThresh = 0;
		boolean bestGuess = false;
		boolean permissive = true;

		String usage = "\n"
									 + "filesys.MergeExtractPipeline requires 4+ arguments\n"
									 + "   (1) Run directory (output files and temporary files will be created here) (i.e. runDir="
									 + rundir
									 + " (default))\n"
									 + "   (2) File listing data sources (i.e. data="
									 + data
									 + " (default))\n"
									 + "          Example:\n"
									 + "          dataLabel1\tfullPathDataFile1\tFullPathMapFile1\tFullPathIdFile1\n"
									 + "          dataLabel2\tfullPathDataFile2\tFullPathMapFile2\tFullPathIdFile2\n"
									 + "          dataLabel3\tdir1\tdataFileExt1\tmapFileExt1\tidFile3\n"
									 + "          dataLabel4\tdir2\tdataFileExt2\tmapFileExt2\tidFile4\n"
									 + "   (3a) Regions-to-extract filename (i.e. regions="
									 + regions
									 + " (default))\n"
									 + "   (3b) Markers-to-extract filename (i.e. markers="
									 + markers
									 + " (default))\n"
									 + "          (Note: only one is allowed, either regions or markers, not both)\n"
									 + "   (4) Optional: output Data filename (i.e. outData="
									 + outfileD
									 + " (defaults to root of marker or region file + \".db.xln.gz\"))\n"
									 + "   (5) Optional: output Map filename (i.e. outMap="
									 + outfileM
									 + " (defaults to root of marker or region file + \".snp\"))\n"
									 + "   (6) Optional: Log file name (i.e. log="
									 + logFile
									 + " (default))\n"
									 + "   (7) Optional: Split output files (if region file is specified) (i.e. split="
									 + split
									 + " (default))\n"
									 + "   (8) Optional: enable/disable permissive marker lookups (creates a temporary file with the chr:pos identifier for each marker added) (i.e. permissive="
									 + permissive
									 + " (default))\n"
									 + "   (9) Optional: Overwrite files if they already exist (i.e. overwrite="
									 + overwrite
									 + " (default))\n"
									 + "   (10) Optional: Rename markers in any data files to dataLabel+MarkerName (i.e. rename="
									 + rename
									 + " (default))\n"
									 + "\n"
									 + "   (11) Optional: if source data contains genotype probabilities and output format is in dosage values, use \"Best Guess\" dosage rather than computed dosage value. The value specified is used such that if all genotype probability values are below the given value, the dosage value is set to missing. (i.e. best="
									 + bestThresh + " (not the default))\n"
									 + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("outData=")) {
				outfileD = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("outMap=")) {
				outfileM = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("runDir=")) {
				rundir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("markers=")) {
				markers = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("data=")) {
				data = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("regions=")) {
				regions = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("split=")) {
				split = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("overwrite=")) {
				overwrite = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("rename=")) {
				rename = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("permissive=")) {
				permissive = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("best=")) {
				bestGuess = true;
				bestThresh = ext.parseFloatArg(arg);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0 || args.length == 0) {
			System.err.println(usage);
			System.exit(1);
		}

		MergeExtractPipeline mep = new MergeExtractPipeline();
		mep.setLogger(logFile == null
																 ? new Logger(ext.rootOf(data, false)
																							+ ext.replaceWithLinuxSafeCharacters(ext.getDate()
																																									 + "_"
																																									 + ext.getTime())
																							+ ".log")
																 : new Logger(logFile));
		if (rundir != null) {
			mep.setRunDirectory(rundir, true);
		}
		if (regions != null) {
			mep.setRegions(regions);
		}
		if (markers != null) {
			mep.setMarkers(markers);
		}
		if (overwrite) {
			mep.setOverwrite();
		}
		mep.setSplitOutputByRegions(split);
		mep.setRenameMarkers(rename);
		mep.setOutputFiles(outfileD, outfileM);
		if (bestGuess) {
			mep.setBestGuessOutput(bestGuess);
			mep.setBestGuessThreshold(bestThresh);
		}
		mep.setPermissiveMarkerFiltering(permissive);
		mep.initLog();
		ArrayList<DataSource> dss = parseDataFile(mep.runDir, mep.markerLocations, mep.regions, data,
																							0,
																							mep.log);
		for (DataSource ds : dss) {
			mep.addDataSource(ds);
		}

		mep.run();
	}
}
