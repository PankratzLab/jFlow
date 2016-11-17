package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.poi.util.IOUtils;
import org.genvisis.CLI;
import org.genvisis.cnv.Launch;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.ProgressMonitor;
import org.genvisis.common.ext;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class ImputationPrep {

	private static final String[] HRC_COLS = new String[] {"#CHROM", "POS", "ID", "REF", "ALT", "AF"};
	private static final Set<Character> VALID_ALLELES = ImmutableSet.of('A', 'T', 'G', 'C');
	private static final Map<Character, Character> PALINDROMIC_PAIRS = ImmutableMap.of(	'A', 'T', 'T',
																																											'A', 'G', 'C',
																																											'C', 'G');
	private static final double PALINDROME_MAF_LIMIT = 0.4;

	public static final String ARGS_PROJ = "proj";
	public static final String ARGS_TARGET_DIR = "dir";
	public static final String ARGS_REFERENCE_FILE = "refFile";

	private static final String DEFAULT_TARGET_DIR = "imputationPrep/";
	private static final String DEFAULT_REFERENCE_FILE = "HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz";


	private Project proj;
	private String referenceFile;
	private String targetDir;
	private Logger log;

	private Map<Byte, Map<Integer, Set<ReferencePosition>>> referencePositions;
	private String imputationRefAlleles;

	private static class ReferencePosition {
		private byte chr;
		private int position;
		private String id;
		private char ref;
		private char alt;
		private double altFreq;

		/**
		 * @param chr
		 * @param position
		 * @param id
		 * @param ref
		 * @param alt
		 * @param altFreq
		 */
		public ReferencePosition(	byte chr, int position, String id, char ref, char alt,
															double altFreq) {
			super();
			this.chr = chr;
			this.position = position;
			this.id = id;
			this.ref = ref;
			this.alt = alt;
			this.altFreq = altFreq;
		}

		public byte getChr() {
			return chr;
		}

		public int getPosition() {
			return position;
		}

		public String getId() {
			return id;
		}

		public char getRef() {
			return ref;
		}

		public char getAlt() {
			return alt;
		}

		public double getAltFreq() {
			return altFreq;
		}



	}


	/**
	 * @param proj
	 * @param referenceFile
	 * @param targetDir
	 * @param log
	 */
	public ImputationPrep(Project proj, String referenceFile, String targetDir, Logger log) {
		super();
		this.proj = proj;
		this.referenceFile = referenceFile;
		this.targetDir = targetDir;
		this.log = log;
		readRefFile();
	}

	private boolean readRefFile() {
		log.report("Parsing reference panel file");
		long time = System.currentTimeMillis();
		referencePositions = Maps.newHashMap();
		String taskName = "parseRefFile";
		proj.getProgressMonitor().beginIndeterminateTask(	taskName, "Parsing reference panel file",
																											ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
		BufferedReader reader = null;
		try {
			reader = Files.getAppropriateReader(referenceFile);
			String header = reader.readLine();
			if (header == null) {
				log.reportError("Reference file is empty");
				return false;
			}
			String delim = ext.determineDelimiter(header);
			int[] cols = ext.indexFactors(HRC_COLS, header.split(delim), false, log, true, false);
			while (reader.ready()) {
				String[] refLine = Array.subArray(reader.readLine().split(delim), cols);
				proj.getProgressMonitor().updateTask(taskName);
				byte chr = Positions.chromosomeNumber(refLine[0], log);
				int position;
				try {
					position = Integer.parseInt(refLine[1]);
				} catch (NumberFormatException e) {
					log.reportError("Imputation reference file ("	+ referenceFile
															+ ") contains a non-integer position: " + refLine[1]);
					return false;
				}
				String id = refLine[2];
				char ref = refLine[3].length() == 1	? refLine[3].toUpperCase().charAt(0)
																						: ABLookup.MISSING_ALLELE;
				char alt = refLine[4].length() == 1	? refLine[4].toUpperCase().charAt(0)
																						: ABLookup.MISSING_ALLELE;
				double altFreq;
				try {
					altFreq = Double.parseDouble(refLine[5]);
				} catch (NumberFormatException e) {
					log.reportError("Imputation reference file ("	+ referenceFile
															+ ") contains a non-numeric alternate allele frequency: "
															+ refLine[5]);
					return false;
				}
				Map<Integer, Set<ReferencePosition>> posMap = referencePositions.get(chr);
				if (posMap == null) {
					posMap = Maps.newHashMap();
					referencePositions.put(chr, posMap);
				}
				Set<ReferencePosition> refPosSet = posMap.get(position);
				if (refPosSet == null) {
					refPosSet = Sets.newHashSet();
					posMap.put(position, refPosSet);
				}
				refPosSet.add(new ReferencePosition(chr, position, id, ref, alt, altFreq));
			}
		} catch (IOException ioe) {
			log.reportIOException(referenceFile);
		} finally {
			IOUtils.closeQuietly(reader);
		}
		proj.getProgressMonitor().endTask(taskName);
		log.report("Finished parsing Reference File in " + ext.getTimeElapsed(time));
		return true;
	}

	private String generateFilteredPlinkset() {
		MarkerSet markerSet = proj.getMarkerSet();
		String[] markerNames = markerSet.getMarkerNames();
		ABLookup projAB =
										new ABLookup(markerNames, proj.AB_LOOKUP_FILENAME.getValue(), true, true, log);
		Set<String> keepMarkers = Sets.newHashSet();
		char[][] lookup = projAB.getLookup().clone();
		char[] refAlleles = new char[lookup.length];
		int mismatchPos = 0;
		int invalidAlleles = 0;
		int palindromes = 0;
		int alleleFlips = 0;
		int strandFlips = 0;
		int strandAlelleFlips = 0;
		int matches = 0;
		int mismatchAlleles = 0;

		log.report("Generating AB Lookup for Imputation");
		String taskName = "generateABLookup";
		proj.getProgressMonitor().beginDeterminateTask(	taskName, "Generating AB Lookup for Imputation",
																										markerNames.length,
																										ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
		for (int i = 0; i < markerNames.length; i++) {
			proj.getProgressMonitor().updateTask(taskName);
			Map<Integer, Set<ReferencePosition>> chrMap = referencePositions.get(markerSet.getChrs()[i]);
			if (chrMap == null) {
				log.reportError("Warning - Chr " + markerSet.getChrs()[i] + " missing from reference file");
				referencePositions.put(	markerSet.getChrs()[i],
																new HashMap<Integer, Set<ReferencePosition>>());
				chrMap = referencePositions.get(markerSet.getChrs()[i]);
			}
			Set<ReferencePosition> refMatches = chrMap.get(markerSet.getPositions()[i]);
			if (refMatches == null) {
				mismatchPos++;
				// TODO Maybe add checking by name, shouldn't be necessary with BLAST VCF though
			} else {
				refMatches = Sets.newHashSet(refMatches);
				char a = lookup[i][0];
				char b = lookup[i][1];
				refAlleles[i] = a;
				// Find the best matched ReferencePosition from all of the matches
				if (!validateAlleles(a, b, refMatches)) {
					invalidAlleles++;
				} else if (!validatePalindromes(refMatches)) {
					palindromes++;
				} else if (matches(a, b, refMatches)) {
					matches++;
					keepMarkers.add(markerNames[i]);
				} else if (matches(b, a, refMatches)) {
					alleleFlips++;
					keepMarkers.add(markerNames[i]);
					refAlleles[i] = b;
				} else {
					char aFlip = PALINDROMIC_PAIRS.get(a);
					char bFlip = PALINDROMIC_PAIRS.get(b);
					if (matches(aFlip, bFlip, refMatches)) {
						strandFlips++;
						keepMarkers.add(markerNames[i]);
						lookup[i] = new char[] {aFlip, bFlip};
						refAlleles[i] = aFlip;
					} else if (matches(bFlip, aFlip, refMatches)) {
						strandAlelleFlips++;
						keepMarkers.add(markerNames[i]);
						lookup[i] = new char[] {aFlip, bFlip};
						refAlleles[i] = bFlip;
					} else {
						mismatchAlleles++;
					}
				}
			}
		}
		proj.getProgressMonitor().endTask(taskName);

		log.report(mismatchPos + " positions did not match to a position in the reference set");
		log.report(invalidAlleles
								+ " positions had invalid allele codes in the project or reference set");
		log.report(mismatchAlleles + " positions had mismatched alleles and were excluded");
		log.report(palindromes
								+ " positions had A/T or G/C SNPs in the reference set and were excluded");
		log.report("A total of "	+ (mismatchPos + invalidAlleles + mismatchAlleles + palindromes)
								+ " positions were excluded");
		log.report("");
		log.report(matches + " positions matched the reference set");
		log.report(alleleFlips + " positions matched after flipping the alleles");
		log.report(strandFlips + " positions matched after flipping the strand");
		log.report(strandAlelleFlips + " positions matched after flipping the alleles and strand");
		log.report("A total of "	+ (matches + alleleFlips + strandFlips + strandAlelleFlips)
								+ " positions will be exported for imputation");

		String imputationABLookup = ext.addToRoot(proj.AB_LOOKUP_FILENAME.getDefaultValue(), "_imputation");
		new ABLookup(markerNames, lookup).writeToFile(proj.PROJECT_DIRECTORY.getValue()	+ targetDir
																									+ imputationABLookup, log);

		imputationRefAlleles = "refAlleles.txt";
		PrintWriter writer = Files.getAppropriateWriter(proj.PROJECT_DIRECTORY.getValue() + targetDir + imputationRefAlleles);
		for (int i = 0; i < refAlleles.length; i++) {
			writer.println(markerNames[i] + "\t" + refAlleles[i]);
		}
		writer.flush();
		writer.close();

		String imputationTargetMarkers = proj.PROJECT_DIRECTORY.getValue()	+ targetDir
																			+ "targetMarkers.txt";
		Files.writeIterable(keepMarkers, imputationTargetMarkers);

		log.report("Generating PLINK dataset for Imputation");
		String filteredPlinkroot = "plinkFiltered";
		if (!PlinkData.saveGenvisisToPlinkBedSet(	proj, targetDir + filteredPlinkroot,
																							proj.CLUSTER_FILTER_COLLECTION_FILENAME.getValue(),
																							imputationTargetMarkers, -1, true)) {
			return null;
		}

		return filteredPlinkroot;
	}

	private String dedupePlinkSNPs(String inputPlinkroot) {
		List<String> commandList;
		Set<String> inputs;
		Set<String> outputs;

		List<CmdLine.Command> commands = Lists.newLinkedList();



		String dupesFile = inputPlinkroot + ".dupvar";

		commandList = ImmutableList.of(	"plink2", "--bfile", inputPlinkroot, "--list-duplicate-vars",
																		"ids-only", "suppress-first", "--out", inputPlinkroot);
		inputs = ImmutableSet.copyOf(PSF.Plink.getPlinkBedBimFam(inputPlinkroot));
		outputs = ImmutableSet.of(dupesFile);
		commands.add(new CmdLine.Command(	commandList, inputs, outputs,
																			proj.PROJECT_DIRECTORY.getValue() + targetDir));

		String scrubbedPlinkroot = "scrubbed";
		commandList = ImmutableList.of(	"plink2", "--bfile", inputPlinkroot, "--extract", dupesFile,
																		"--make-bed", "--out", scrubbedPlinkroot);
		inputs = outputs;
		outputs = ImmutableSet.copyOf(PSF.Plink.getPlinkBedBimFam(scrubbedPlinkroot));
		commands.add(new CmdLine.Command(	commandList, inputs, outputs,
																			proj.PROJECT_DIRECTORY.getValue() + targetDir));

		// TODO Drop lower callrate duplicate

		for (CmdLine.Command command : commands) {
			if (!command.runCommand(true, true, false, false, log)) {
				return null;
			}
		}
		return scrubbedPlinkroot;
	}

	private List<String> exportPlinkToChrVCFs(String inputPlinkroot) {
		List<String> vcfs = Lists.newArrayList();
		Set<String> inputs = ImmutableSet.copyOf(PSF.Plink.getPlinkBedBimFam(inputPlinkroot));
		for (int i = 1; i <= 23; i++) {
			String outputRoot = proj.PROJECT_NAME.getValue() + "_imputation_chr" + i;
			String outputFile = outputRoot + ".vcf.gz";
			Set<String> outputs = ImmutableSet.of(outputFile);
			List<String> commandList = ImmutableList.of("plink2", "--bfile", inputPlinkroot,
																									"--a2-allele", imputationRefAlleles,
																									"--real-ref-alleles", "--chr",
																									Integer.toString(i), "--output-chr", "M",
																									"--recode", "vcf", "bgz", "--out", outputRoot);
			CmdLine.Command command = new CmdLine.Command(commandList, inputs, outputs,
																										proj.PROJECT_DIRECTORY.getValue() + targetDir);
			if (command.runCommand(true, true, false, true, log)) {
				vcfs.add(outputFile);
			} else {
				log.reportError("Could not generate VCF for chromosome " + i);
			}
		}
		return vcfs;
	}

	public boolean run() {
		String filteredPlinkroot = generateFilteredPlinkset();
		if (filteredPlinkroot == null
				|| !PSF.Plink.allFilesExist(proj.PROJECT_DIRECTORY.getValue()	+ targetDir
																		+ filteredPlinkroot, true)) {
			log.reportError("Could not generate a filtered plink dataset for imputation");
			return false;
		}

		String qcPlinkroot = Qc.fullGamut(proj.PROJECT_DIRECTORY.getValue()	+ targetDir,
																			filteredPlinkroot, true, log);
		if (qcPlinkroot == null || !PSF.Plink.allFilesExist(qcPlinkroot, true)) {
			log.reportError("GWAS QC failed");
			return false;
		}
		qcPlinkroot = qcPlinkroot.substring((proj.PROJECT_DIRECTORY.getValue() + targetDir).length());

		String dedupedPlinkroot = dedupePlinkSNPs(qcPlinkroot);
		if (dedupedPlinkroot == null
				|| !PSF.Plink.allFilesExist(proj.PROJECT_DIRECTORY.getValue()	+ targetDir
																		+ dedupedPlinkroot, true)) {
			log.reportError("Failed to remove duplicates from exported plink dataset");
			return false;
		}

		List<String> vcfs = exportPlinkToChrVCFs(dedupedPlinkroot);
		if (vcfs.isEmpty()) {
			log.reportError("No VCFs were exported");
			return false;
		}
		log.report("Succesfully exported "	+ vcfs.size() + " VCFs for imputation to "
								+ proj.PROJECT_DIRECTORY.getValue() + targetDir + ":\n\t"
								+ Joiner.on("\n\t").join(vcfs));
		return true;
	}

	private static boolean validateAlleles(char a, char b, Collection<ReferencePosition> refMatches) {
		if (!validAlleles(a, b)) {
			return false;
		}
		for (Iterator<ReferencePosition> i = refMatches.iterator(); i.hasNext();) {
			ReferencePosition refPos = i.next();
			char ref = refPos.getRef();
			char alt = refPos.getAlt();
			if (!validAlleles(ref, alt)) {
				i.remove();
			}
		}
		return !refMatches.isEmpty();
	}

	private static boolean validAlleles(char... alleles) {
		for (char allele : alleles) {
			if (!VALID_ALLELES.contains(allele)) {
				return false;
			}
		}
		return true;
	}

	private static boolean validatePalindromes(Collection<ReferencePosition> refMatches) {
		for (Iterator<ReferencePosition> i = refMatches.iterator(); i.hasNext();) {
			ReferencePosition refPos = i.next();
			char ref = refPos.getRef();
			char alt = refPos.getAlt();
			double maf = refPos.getAltFreq();
			if (maf > 0.5) {
				maf -= 0.5;
			}
			if (maf > PALINDROME_MAF_LIMIT && palindromic(ref, alt)) {
				i.remove();
			}
		}
		return !refMatches.isEmpty();
	}

	private static boolean palindromic(char a1, char a2) {
		Character a1Pair = PALINDROMIC_PAIRS.get(a1);
		return a1Pair != null && a1Pair.charValue() == a2;
	}

	private static boolean matches(char a, char b, Collection<ReferencePosition> refMatches) {
		for (ReferencePosition refPos : refMatches) {
			if (a == refPos.getRef() && b == refPos.getAlt()) {
				return true;
			}
		}
		return false;
	}

	public static void main(String... args) {
		Project proj;
		String projFile = Launch.getDefaultDebugProjectFile(false);
		String targetDir = DEFAULT_TARGET_DIR;
		String referenceFile = DEFAULT_REFERENCE_FILE;

		CLI c = new CLI(ImputationPrep.class);
		c.addArgWithDefault(ARGS_PROJ, "project properties filename", projFile);
		c.addArgWithDefault(ARGS_TARGET_DIR,
												"generate intermediates and outputs to this directory (relative to project directory)",
												targetDir);
		c.addArgWithDefault(ARGS_REFERENCE_FILE, "prepare for imputation based on this reference file",
												referenceFile);

		c.parseWithExit(args);

		proj = new Project(c.get(ARGS_PROJ), false);
		targetDir = c.get(ARGS_TARGET_DIR);
		referenceFile = c.get(ARGS_REFERENCE_FILE);
		if (Files.isRelativePath(referenceFile)) {
			referenceFile = Files.firstPathToFileThatExists(referenceFile,
																											proj.PROJECT_DIRECTORY.getValue(),
																											proj.PROJECT_DIRECTORY.getValue() + targetDir,
																											"");
		}
		if (referenceFile == null || !Files.exists(referenceFile)) {
			proj.getLog().reportError("Reference file could not be found");
			return;
		}

		new ImputationPrep(proj, referenceFile, targetDir, proj.getLog()).run();
	}
}
