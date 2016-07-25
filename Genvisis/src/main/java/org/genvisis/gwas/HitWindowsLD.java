package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;

//TODO, this is a bit sprawling for a simple thing
public class HitWindowsLD {
	private static final String[] MARK_REPORT_HEADER = { "Target", "Target.chr", "Target.pos", "LD.Hits", "LD.refHits", "LD.minPos", "LD.minMarker", "LD.maxPos", "LD.maxMarker", "Target.inRef", "UCSC", "UCSC.ExcelHyperLink", "Window", "Source File" };
	private static final String[] REGION_REPORT_HEADER = { "Region", "Target.chr", "LD.Hits", "Region.minPos", "Region.minMarker", "Region.maxPos", "Region.maxMarker", "Region.markersWithLDHits", "UCSC", "UCSC.ExcelHyperLink", "Window", "Source File" };
	private static String[] FILE_TYPES = { ".ld", ".flt", ".other" };

	private String target;
	private int minHitPos;
	private int maxHitPos;
	private String minHitMarker;
	private String maxHitMarker;
	private ArrayList<String> markers;
	private int numHits;
	private int numRefHits;
	private int chr;
	private int window;
	private boolean inRef;
	private int pos;
	private boolean region;
	// private Logger log;
	private double filter;

	public HitWindowsLD(String target, int window, Logger log) {
		this.target = target;
		this.minHitMarker = "NA";
		this.maxHitMarker = "NA";
		this.numHits = 0;
		this.numRefHits = 0;
		this.window = window;
		this.inRef = false;
		// this.log = log;
	}

	public void setLoc(ReferenceMap referenceMap) {
		this.chr = referenceMap.getChr(target);
		this.pos = referenceMap.getPos(target);
		this.minHitPos = pos;
		this.maxHitPos = pos;
	}

	public double getFilter() {
		return filter;
	}

	public void setFilter(double filter) {
		this.filter = filter;
	}

	public int getWindow() {
		return window;
	}

	public int getChr() {
		return chr;
	}

	public boolean isInRef() {
		return inRef;
	}

	public void setInRef(boolean inRef) {
		this.inRef = inRef;
	}

	public String getTarget() {
		return target;
	}

	public boolean isRegion() {
		return region;
	}

	public int getMinHitPos() {
		return minHitPos;
	}

	public int getMaxHitPos() {
		return maxHitPos;
	}

	public String getMinHitMarker() {
		return minHitMarker;
	}

	public String getMaxHitMarker() {
		return maxHitMarker;
	}

	public int getNumHits() {
		return numHits;
	}

	public void setRegion(boolean region) {
		this.region = region;
		this.markers = new ArrayList<String>();
	}

	public ArrayList<String> getMarkers() {
		return markers;
	}

	public void checkHit(String hit1, int position1, String hit2, int position2, int checkchr, ReferenceMap referenceMap) {
		if (checkchr == chr && hit1.equals(target)) {
			checkMin(hit2, position2, referenceMap);
			checkMax(hit2, position2, referenceMap);
		} else if (checkchr == chr && hit2.equals(target)) {
			checkMin(hit1, position1, referenceMap);
			checkMax(hit1, position1, referenceMap);
		}
	}

	public void checkMin(String hit, int position, ReferenceMap referenceMap) {
		if (window == -1 || Math.abs((position - pos)) < window) {
			if (position < pos) {
				numHits++;
				if (referenceMap.inReference(hit)) {
					numRefHits++;
					if (region) {
						markers.add(hit);
					}
				}
			}
			if (position < minHitPos) {
				minHitPos = position;
				minHitMarker = hit;
			}
		}
	}

	public void checkMax(String hit, int position, ReferenceMap referenceMap) {
		if (window == -1 || Math.abs((position - pos)) < window) {
			if (position > pos) {
				numHits++;
				if (referenceMap.inReference(hit)) {
					numRefHits++;
					if (region) {
						markers.add(hit);
					}
				}
			}
			if (position > maxHitPos) {
				maxHitPos = position;
				maxHitMarker = hit;
			}
		}
	}

	public void setLiberal(ReferenceMap referenceMap) {
		if (inRef) {
			int previousAvailableIndex = referenceMap.getPreviousAvailableIndex(minHitMarker);
			int nextAvailableIndex = referenceMap.getNextAvailableIndex(maxHitMarker);
			int previousChr = referenceMap.getChrByIndex(previousAvailableIndex);
			int nextChr = referenceMap.getChrByIndex(nextAvailableIndex);
			int previousPos = referenceMap.getPosByIndex(previousAvailableIndex);
			int nextPos = referenceMap.getPosByIndex(nextAvailableIndex);
			if (previousChr == nextChr && nextChr == chr) {
				if (Math.abs(previousPos - pos) < window) {
					minHitPos = previousPos;
					minHitMarker = referenceMap.getNameByIndex(previousAvailableIndex);
				}
				if (Math.abs(nextPos - pos) < window) {
					maxHitPos = nextPos;
					maxHitMarker = referenceMap.getNameByIndex(nextAvailableIndex);
				}
			}
		}
	}

	public String getSummary(ReferenceMap referenceMap, boolean liberal) {
		int[] posUCSC = { chr, minHitPos, maxHitPos };
		String UCSC = Positions.getUCSCformat(posUCSC);
		String UCSCLink = Positions.getUCSClinkInExcel(posUCSC);
		String summary = "";
		if (!inRef) {
			String[] summ = new String[12];
			Arrays.fill(summ, "NA");
			summ[8] = "false";
			summ[11] = window + "";
			summary += target + "\t" + Array.toStr(summ);
		} else {
			if (minHitPos == pos) {
				minHitMarker = target;
			}
			if (maxHitPos == pos) {
				maxHitMarker = target;
			}
			if (liberal) {
				setLiberal(referenceMap);
			}
			summary += target + "\t" + chr + "\t" + pos + "\t" + numHits + "\t" + numRefHits + "\t" + minHitPos + "\t" + minHitMarker + "\t" + maxHitPos + "\t" + maxHitMarker + "\t" + inRef + "\t" + UCSC + "\t" + UCSCLink + "\t" + window;
		}
		return summary;
	}

	public static void analyzeLD(String dir, String[] targetFile, String refMap, int window, boolean liberal, String fileExt, boolean filtered, boolean region, String output, double filter, int filterOnColumn, Logger log) {
		analyzeLD(dir, targetFile, ReferenceMap.getMap(refMap, log), window, liberal, fileExt, filtered, region, output, filter, filterOnColumn, log);
	}

	public static void analyzeLD(String dir, String[] targetFiles, ReferenceMap referenceMap, int window, boolean liberal, String fileExt, boolean filtered, boolean region, String output, double filter, int filterOnColumn, Logger log) {
		int fileType = ext.indexOfStr(fileExt, FILE_TYPES);
		if (fileType == -1) {
			log.reportError("Warning - assuming haploview format");
			fileType = 2;
		}
		String[] files = Files.list(dir, fileExt, false);
		log.report("Info - found " + files.length + " files with extension " + fileExt + " in directory " + dir);
		for (int i = 0; i < files.length; i++) {
			log.report(ext.getTime() + " Info - checking file " + files[i]);
			for (int j = 0; j < targetFiles.length; j++) {
				analyzeTargetLDCombo(targetFiles[j], referenceMap, window, liberal, region, output, log, fileType, filter, filterOnColumn, dir + files[i]);
			}
		}
	}

	public static void determineLD(String targetFile, String mapFile, String output, String ldFile, int window, double filter, int column, Logger log) {
		int fileType = determineFileType(ldFile, log);
		analyzeTargetLDCombo(targetFile, ReferenceMap.getMap(mapFile, log), window, true, true, output, log, fileType, filter, column, ldFile);
	}

	public static void analyzeTargetLDCombo(String targetFile, ReferenceMap referenceMap, int window, boolean liberal, boolean region, String output, Logger log, int fileType, double filter, int filterOnColumn, String ldFile) {
		HitWindowsLD[] hitWindowsLD = prepWindows(targetFile, referenceMap, window, region, filter, log);
		harvestLD(ldFile, hitWindowsLD, referenceMap, fileType, filter, filterOnColumn, log);
		summarize(targetFile, ldFile, referenceMap, liberal, hitWindowsLD, region, output, log);
	}

	private static HitWindowsLD[] prepWindows(String targetFile, ReferenceMap referenceMap, int window, boolean region, double filter, Logger log) {
		return prepWindows(getTargetWindows(targetFile, window, log), referenceMap, window, region, filter, log);
	}

	public static void determineLD(String[] targets, String regionName, String mapFile, String output, String ldFile, int window, double filter, int column, boolean region, Logger log) {
		int fileType = determineFileType(ldFile, log);
		analyzeTargetLDCombo(targets, regionName, ReferenceMap.getMap(mapFile, log), window, true, region, output, log, fileType, filter, column, ldFile);
	}

	public static void analyzeTargetLDCombo(String[] targets, String regionName, ReferenceMap referenceMap, int window, boolean liberal, boolean region, String output, Logger log, int fileType, double filter, int filterOnColumn, String ldFile) {
		HitWindowsLD[] hitWindowsLD = prepWindows(getTargetWindows(targets, window, log), referenceMap, window, region, filter, log);
		harvestLD(ldFile, hitWindowsLD, referenceMap, fileType, filter, filterOnColumn, log);
		summarize(regionName, ldFile, referenceMap, liberal, hitWindowsLD, region, output, log);
	}

	public static HitWindowsLD[] prepWindows(HitWindowsLD[] hitWindowsLD, ReferenceMap referenceMap, int window, boolean region, double filter, Logger log) {
		for (int i = 0; i < hitWindowsLD.length; i++) {
			if (!referenceMap.inReference(hitWindowsLD[i].getTarget())) {
				log.reportError("Warning - marker " + hitWindowsLD[i].getTarget() + " was not present in the reference LD dataset/map file");
			} else {
				hitWindowsLD[i].setInRef(true);
				hitWindowsLD[i].setLoc(referenceMap);
				hitWindowsLD[i].setRegion(region);
				hitWindowsLD[i].setFilter(filter);
			}
		}
		return hitWindowsLD;
	}

	public static HitWindowsLD[] getTargetWindows(String targetFile, int window, Logger log) {
		return getTargetWindows(getTargets(targetFile, log), window, log);
	}

	public static HitWindowsLD[] getTargetWindows(String[] targets, int window, Logger log) {
		HitWindowsLD[] hitWindowsLD = new HitWindowsLD[targets.length];
		for (int i = 0; i < targets.length; i++) {
			hitWindowsLD[i] = new HitWindowsLD(targets[i], window, log);
		}
		return hitWindowsLD;
	}

	public static String[] getTargets(String targetFile, Logger log) {
		ArrayList<String> targets = new ArrayList<String>();
		BufferedReader reader;
		try {
			reader = Files.getAppropriateReader(targetFile);
			while (reader.ready()) {
				String[] line = reader.readLine().split("\t");
				targets.add(line[0]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + targetFile + "\" not found in current directory");
			return null;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + targetFile + "\"");
			return null;
		}
		return targets.toArray(new String[targets.size()]);
	}

	private static int determineFileType(String ldFile, Logger log) {
		for (int i = 0; i < FILE_TYPES.length; i++) {
			if (ldFile.endsWith(FILE_TYPES[i])) {
				return i;
			}
		}
		// defualt to other, which should be haploview
		return (FILE_TYPES.length - 1);
	}

	public static class ReferenceMap implements Serializable {
		private static final long serialVersionUID = 1L;
		private Hashtable<String, Integer> trackpos;
		private Hashtable<String, Integer> trackchr;
		private Hashtable<String, Integer> trackindex;
		private Hashtable<Integer, String> lookup;
		private int loadingChr;
		private int loadingPos;
		private Logger log;
		private int index;

		public ReferenceMap(Logger log) {
			this.trackpos = new Hashtable<String, Integer>();
			this.trackchr = new Hashtable<String, Integer>();
			this.trackindex = new Hashtable<String, Integer>();
			this.lookup = new Hashtable<Integer, String>();
			this.loadingChr = -1;
			this.loadingPos = -1;
			this.log = log;
			this.index = 0;
		}

		public void addMarker(String marker, int chr, int pos) {
			boolean add = true;
			if (trackpos.containsKey(marker)) {
				add = false;
				log.reportError("Error - a duplicate marker " + marker + " was found in the reference map, only retaining map info for the first seen");
			}
			if (chr < loadingChr) {
				add = false;
				log.reportError("Error - the reference map must be sorted in ascending order by chromosome");
				log.reportError("Error - the marker " + marker + " was not added to the reference hash");

			}
			if (chr == loadingChr && pos < loadingPos) {
				add = false;
				log.reportError("Error - the reference map must be sorted in ascending order by chromosome, and then by position");
				log.reportError("Error - the marker " + marker + " was not added to the reference hash");
			}
			if (loadingChr != chr) {
				loadingChr = chr;
			}
			if (add) {
				loadingPos = pos;
				trackpos.put(marker, pos);
				trackchr.put(marker, chr);
				trackindex.put(marker, index);
				lookup.put(index, marker);
				index++;
			}
		}

		public void serialize(String filename) {
			SerializedFiles.writeSerial(this, filename);
		}

		public static ReferenceMap load(String filename, boolean jar) {
			return (ReferenceMap) SerializedFiles.readSerial(filename, jar, true);
		}

		public boolean inReference(String marker) {
			return (trackpos.containsKey(marker) && trackchr.containsKey(marker));
		}

		public int getIndex(String marker) {
			return trackindex.get(marker);
		}

		public int getPreviousAvailableIndex(String marker) {
			int currentIndex = getIndex(marker);
			if (currentIndex >= 1 && trackchr.get(lookup.get(currentIndex)) == trackchr.get(lookup.get(currentIndex - 1))) {
				return getIndex(lookup.get(currentIndex - 1));
			} else {
				return getIndex(lookup.get(currentIndex));
			}
		}

		public int getNextAvailableIndex(String marker) {
			int currentIndex = getIndex(marker);
			if (currentIndex < (index - 1) && trackchr.get(lookup.get(currentIndex)) == trackchr.get(lookup.get(currentIndex + 1))) {
				return getIndex(lookup.get(currentIndex + 1));
			} else {
				return getIndex(lookup.get(currentIndex));
			}
		}

		public int getPosByIndex(int index) {
			if (lookup.containsKey(index)) {
				return trackpos.get(lookup.get(index));
			} else {
				log.reportError("Error -invalid index");
				return -1;
			}
		}

		public int getChrByIndex(int index) {
			if (lookup.containsKey(index)) {
				return trackchr.get(lookup.get(index));
			} else {
				log.reportError("Error -invalid index");
				return -1;
			}
		}

		public String getNameByIndex(int index) {
			if (lookup.containsKey(index)) {
				return lookup.get(index);
			} else {
				log.reportError("Error -invalid index");
				return "NA";
			}
		}

		public int getChr(String marker) {
			return trackchr.get(marker);
		}

		public int getPos(String marker) {
			return trackpos.get(marker);
		}

		public static ReferenceMap getMap(String mapFile, Logger log) {
			ReferenceMap referenceMap = new ReferenceMap(log);
			BufferedReader reader;
			String refMapSer = ext.rootOf(mapFile, false) + ".ser";
			if (Files.exists(refMapSer)) {
				log.report(ext.getTime() + " Info - Loading serialized map file " + refMapSer);
				referenceMap = ReferenceMap.load(refMapSer, false);
				log.report(ext.getTime() + " Info - Finished loading serialized map file " + refMapSer);
			} else {
				try {
					reader = Files.getAppropriateReader(mapFile);
					while (reader.ready()) {
						String[] line = reader.readLine().split("\t");
						referenceMap.addMarker(line[1], Integer.parseInt(line[0]), Integer.parseInt(line[3]));
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					log.reportError("Error: file \"" + mapFile + "\" not found in current directory");
					return null;
				} catch (IOException ioe) {
					log.reportError("Error reading file \"" + mapFile + "\"");
					return null;
				}
				referenceMap.serialize(refMapSer);
			}
			return referenceMap;
		}
	}

	public static void harvestLD(String ldfile, HitWindowsLD[] hitWindowsLD, ReferenceMap referenceMap, int fileType, double filter, int filterOnColumn, Logger log) {
		BufferedReader reader;
		try {
			reader = Files.getAppropriateReader(ldfile);
			while (reader.ready()) {
				String[] line = reader.readLine().split("[\\s]+");
				for (int i = 0; i < hitWindowsLD.length; i++) {
					if (fileType == 0) {
						try {
							double R2 = Double.parseDouble(line[7]);
							if (line[4].equals(line[1]) && R2 > hitWindowsLD[i].getFilter()) {
								hitWindowsLD[i].checkHit(line[3], Integer.parseInt(line[2]), line[6], Integer.parseInt(line[5]), Integer.parseInt(line[1]), referenceMap);
							}
						} catch (NumberFormatException e) {
							if (!line[7].equals("R2")) {
								log.report("Warning - skipping R2 value " + line[7] + " in file " + ldfile + ", not a number");
							}
						}
					} else if (fileType == 1) {
						hitWindowsLD[i].checkHit(line[2], Integer.parseInt(line[1]), line[4], Integer.parseInt(line[3]), Integer.parseInt(line[0]), referenceMap);
					}

					else if (fileType == 2) {
						try {
							double metric = Double.parseDouble(line[filterOnColumn - 1]);
							if (referenceMap.getChr(line[0]) == referenceMap.getChr(line[1]) && metric >= hitWindowsLD[i].getFilter()) {
								try {
									hitWindowsLD[i].checkHit(line[0], referenceMap.getPos(line[0]), line[1], referenceMap.getPos(line[1]), referenceMap.getChr(line[0]), referenceMap);
								} catch (NullPointerException npe) {
									log.reportError("Error - could not find both markers (" + line[0] + " and " + line[1] + " in reference map");
									log.reportException(npe);
								}
							}
						} catch (NumberFormatException e) {
							log.report("Warning - skipping value " + line[filterOnColumn - 1] + " in file " + ldfile + ", not a number");
						}
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + ldfile + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + ldfile + "\"");
			return;
		}
	}

	//TODO this is a crappy summary format...
	private static HitWindowsLD[] summarize(String targetFile, String ldFile, ReferenceMap referenceMap, boolean liberal, HitWindowsLD[] hitWindowsLD, boolean region, String output, Logger log) {
		PrintWriter writer;
		String source = ext.rootOf(ldFile);
		try {
			if (Files.exists(output + "." + source)) {
				writer = new PrintWriter(new FileWriter(output + "." + source, true));
			} else {
				writer = new PrintWriter(new FileWriter(output + "." + source, false));
			}
			writer.println("\n" + Array.toStr(MARK_REPORT_HEADER));
			for (int i = 0; i < hitWindowsLD.length; i++) {
				writer.println(hitWindowsLD[i].getSummary(referenceMap, liberal) + "\t" + source);
			}
			if (region) {
				writer.println("\n" + Array.toStr(REGION_REPORT_HEADER));
				writer.println(ext.rootOf(targetFile, true) + "\t" + getRegionSummary(hitWindowsLD, log) + "\t" + source);
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + output + "." + source + "\" could not be written to (it's probably open)");
			log.reportException(fnfe);
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + output + "." + source + "\"");
			log.reportException(ioe);
			System.exit(2);
		}
		return hitWindowsLD;
	}

	private static String getRegionSummary(HitWindowsLD[] hitWindowsLD, Logger log) {
		Hashtable<String, String> uniq = new Hashtable<String, String>();
		int uniqHits = 0;
		int regionMinPos = 0;
		int regionMaxPos = 0;
		String regionMinMark = "NA";
		String regionMaxMark = "NA";
		int numinLD = 0;
		int use = 0;
		int chr = 0;
		int window = 0;
		for (int i = 0; i < hitWindowsLD.length; i++) {
			if (!hitWindowsLD[i].inRef) {
				use++;
			} else {
				if (hitWindowsLD[i].getNumHits() > 0) {
					numinLD++;
				}
				if (i == use) {
					window = hitWindowsLD[i].getWindow();
					chr = hitWindowsLD[i].getChr();
					regionMinPos = hitWindowsLD[i].getMinHitPos();
					regionMaxPos = hitWindowsLD[i].getMaxHitPos();
					regionMinMark = hitWindowsLD[i].getMinHitMarker();
					regionMaxMark = hitWindowsLD[i].getMaxHitMarker();
				} else if (hitWindowsLD[i].getChr() != chr) {
					log.reportError("Error - mismatched chormosomes for this region, skipping " + hitWindowsLD[i].getTarget());
				} else {
					if (hitWindowsLD[i].getMinHitPos() < regionMinPos) {
						regionMinPos = hitWindowsLD[i].getMinHitPos();
						regionMinMark = hitWindowsLD[i].getMinHitMarker();
					}
					if (hitWindowsLD[i].getMaxHitPos() > regionMaxPos) {
						regionMaxPos = hitWindowsLD[i].getMaxHitPos();
						regionMaxMark = hitWindowsLD[i].getMaxHitMarker();
					}
				}
				ArrayList<String> marks = hitWindowsLD[i].getMarkers();
				for (int j = 0; j < marks.size(); j++) {
					if (!uniq.containsKey(marks.get(j))) {
						uniqHits++;
						uniq.put(marks.get(j), marks.get(j));
					}
				}
			}
		}
		int[] posUCSC = { chr, regionMinPos, regionMaxPos };
		String UCSC = Positions.getUCSCformat(posUCSC);
		String UCSCLink = Positions.getUCSClinkInExcel(posUCSC);
		String summary = chr + "\t" + uniqHits + "\t" + regionMinPos + "\t" + regionMinMark + "\t" + regionMaxPos + "\t" + regionMaxMark + "\t" + numinLD + " of " + hitWindowsLD.length + "\t" + UCSC + "\t" + UCSCLink + "\t" + window;
		return summary;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String refMap = "C:/data/gedi_exomechip/LD_windows/reference.bim";
		String logFile = "C:/data/gedi_exomechip/LD_windows/LD.log";
		String[] targetFiles = { "C:/data/gedi_exomechip/LD_windows/region1.txt", "C:/data/gedi_exomechip/LD_windows/region2.txt", "C:/data/gedi_exomechip/LD_windows/region3.txt" };
		String dir = "C:/data/gedi_exomechip/LD_windows/";
		String fileExt = ".hap";
		boolean filtered = true;
		int window = 100000;
		boolean liberal = true;
		boolean region = true;
		int column = 4;
		double filter = 0.3;
		String output = "C:/data/gedi_exomechip/LD_windows/testHap";
		String usage = "\n" + "gwas.HitWindowsLD requires 3 arguments\n" + "   (1) files (comma-separated) with a list of rs# targets (i.e. targets=" + Array.toStr(targetFiles) + ")\n" + "   (2) plink format map file for lookup (i.e. refMap=" + refMap + " (default))\n" + "   (3) directory containing LD input files (i.e. dir=" + dir + " (default))\n" + "   OPTIONAL: " + "   (4) file extension of LD input files (options are " + Array.toStr(FILE_TYPES) + ")(i.e. fileExt=" + fileExt + " (default))\n" + "   (5) window around hit to extend (i.e. window=" + window + " (default))\n" + "   (6) for \".ld\" and haploview format files, the cutoff for R2,D',or LOD  (i.e. filter=" + filter + " (default))\n" + "   (7) for haploview format files, the column to filter on (i.e. column=" + column + " (default))\n" + "   (8) output file base name (i.e. out=" + output + " (default))\n" + "   (9) do not treat each input target file as a region (i.e. -noRegion(not the default))\n";

		// String usage = "\n" +
		// "gwas.HitWindowsLD requires 3 arguments\n" +
		// "   (1) files (comma-separated) with a list of rs# targets (i.e. targets=" + Array.toStr(targetFiles)+")\n" +
		// "   (2) plink format map file for lookup (i.e. refMap=" + refMap + " (default))\n" +
		// "   (3) directory containing LD input files (i.e. dir=" + dir + " (default))\n" +
		// "   OPTIONAL: "+
		// "   (4) file extension of LD input files (options are " +Array.toStr(FILE_TYPES) +")(i.e. fileExt=" + fileExt + " (default))\n"+
		// "   (5) window around hit to extend (i.e. window=" + window + " (default))\n" +
		// "   (6) for \".ld\" and haploview format files, the cutoff for R2,D',or LOD  (i.e. filter=" + filter + " (default))\n" +
		// "   (7) for haploview format files, the column to filter on (i.e. column=" + column + " (default))\n"+
		// "   (8) output file base name (i.e. out=" + output + " (default))\n"+
		// "   (9) do not treat each input target file as a region (i.e. -noRegion(not the default))\n"
		// ;
		//
		for (int i = 0; i < args.length; i++) {
			if (args[i].startsWith("targets=")) {
				targetFiles = args[i].split("=")[1].split(",");
				numArgs--;
			} else if (args[i].startsWith("refMap=")) {
				refMap = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("dir=")) {
				dir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("fileExt=")) {
				fileExt = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("window=")) {
				window = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				output = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-raw")) {
				filtered = false;
				numArgs--;
			} else if (args[i].startsWith("-restrictive")) {
				liberal = false;
				numArgs--;
			} else if (args[i].startsWith("-noRegion")) {
				region = false;
				numArgs--;
			} else {
				System.err.println("Error - invalid argument " + args[i] + "\n" + usage);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		Logger log = new Logger(logFile, false);
		analyzeLD(dir, targetFiles, refMap, window, liberal, fileExt, filtered, region, output, filter, column, log);
	}
}

// public static Hashtable<String, String> getFilteredLDMap(String filteredLdFile, Hashtable<String, String> ldMap, Logger log) {
// if (ldMap == null) {
// ldMap = new Hashtable<String, String>();
// }
// BufferedReader reader;
// try {
// reader = Files.getAppropriateReader(filteredLdFile);
//
// while (reader.ready()) {
// String[] line = reader.readLine().split("\t");
// ldMap.put(line[2], line[4]);
// ldMap.put(line[4], line[2]);
// }
// reader.close();
// } catch (FileNotFoundException fnfe) {
// log.reportError("Error: file \"" + filteredLdFile + "\" not found in current directory");
// return null;
// } catch (IOException ioe) {
// log.reportError("Error reading file \"" + filteredLdFile + "\"");
// return null;
// }
// return ldMap;
// }
