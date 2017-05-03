package org.genvisis.cnv.analysis.pca;

import java.io.PrintWriter;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.stats.Maths;

import com.google.common.base.Joiner;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class PCImputeRace {

	public static class Sample {
		private final String fid;
		private final String iid;
		private double pc1;
		private double pc2;

		/**
		 * @param fid
		 * @param iid
		 * @param pc1
		 * @param pc2
		 */
		public Sample(String fid, String iid, double pc1, double pc2) {
			super();
			this.fid = fid;
			this.iid = iid;
			this.pc1 = pc1;
			this.pc2 = pc2;
		}

		public String getFid() {
			return fid;
		}

		public String getIid() {
			return iid;
		}

		public String getFidIid() {
			return fid + "\t" + iid;
		}

		public double getPc1() {
			return pc1;
		}

		public double getPc2() {
			return pc2;
		}

		public void setPc1(double pc1) {
			this.pc1 = pc1;
		}

		public void setPc2(double pc2) {
			this.pc2 = pc2;
		}



	}

	public static enum RACE {

		WHITE("White", 1), AA("African American", 2), HISPANIC("Hispanic", 3), ASIAN("Asian", 4);

		private final String description;
		private final int sampleDataClassNum;


		private RACE(String description, int sampleDataClassNum) {
			this.description = description;
			this.sampleDataClassNum = sampleDataClassNum;
		}

		public String getDescription() {
			return description;
		}

		public int getSampleDataClassNum() {
			return sampleDataClassNum;
		}

	}

	public static final String[] STEP_PCS_HEADER = {"FID", "IID", "PC1", "PC2"};
	public static final String[] CORRECTED_PCS_HEADER = {"FID", "IID", "PC1", "PC2", "%African",
																											 "%Asian", "%White"};
	public static final String[] IMPUTED_RACE_SAMPLE_DATA_HEADERS = new String[] {"Class=ImputedRace;1=White;2=African American;3=Hispanic;4=Asian",
																																								"% African",
																																								"% Asian",
																																								"% European"};

	private static final Comparator<Sample> PC1_COMPARATOR = new Comparator<Sample>() {
		@Override
		public int compare(Sample sample1, Sample sample2) {
			return Double.compare(sample1.getPc1(), sample2.getPc1());
		}
	};

	private static final Comparator<Sample> PC2_COMPARATOR = new Comparator<Sample>() {
		@Override
		public int compare(Sample sample1, Sample sample2) {
			return Double.compare(sample1.getPc2(), sample2.getPc2());
		}
	};
	private final Project proj;
	private final Collection<Sample> samples;
	private final Collection<Sample> eurSeeds;
	private final Collection<Sample> afrSeeds;
	private final Collection<Sample> asianSeeds;
	private final Logger log;


	/**
	 * @param proj
	 * @param samples
	 * @param eurSeeds
	 * @param afrSeeds
	 * @param asianSeeds
	 * @param log
	 */
	public PCImputeRace(Project proj, Collection<Sample> samples, Collection<Sample> eurSeeds,
											Collection<Sample> afrSeeds, Collection<Sample> asianSeeds, Logger log) {
		super();
		this.proj = proj;
		this.samples = samples;
		this.eurSeeds = eurSeeds;
		this.afrSeeds = afrSeeds;
		this.asianSeeds = asianSeeds;
		this.log = log;
	}

	public void correctPCsToRace(String outFile) {
		PrintWriter writer;

		log.report("Checking for PC1 and PC2 predictions of African and Asian");

		writeStepPCs(ext.rootOf(outFile, false) + "_start.mds");

		if (!maybeSwapPCs()) {
			return;
		}
		writeStepPCs(ext.rootOf(outFile, false) + "_step0.mds");

		shiftEuropeansOutOfFirstQuadrant();

		writeStepPCs(ext.rootOf(outFile, false) + "_step1.mds");

		forceRightAngleBetweenAsiansAndAfricans();

		writeStepPCs((ext.rootOf(outFile, false) + "_step2.mds"));

		scaleAfricanAndAsianTo1();

		log.report("Calculating Estimated Races");

		Map<Sample, Double> pctsAfrican = Maps.newHashMap();
		Map<Sample, Double> pctsAsian = Maps.newHashMap();
		Map<Sample, Double> pctsEuropean = Maps.newHashMap();
		Map<Sample, RACE> imputedRaces = Maps.newHashMap();

		for (Sample sample : samples) {
			double pctAfrican = Math.min(Math.max(sample.getPc1(), 0.0), 1.0);
			double pctAsian = Math.min(Math.max(sample.getPc2(), 0.0), 1.0);
			double pctEuropean = 1.0 - (pctAfrican + pctAsian);
			pctsAfrican.put(sample, pctAfrican);
			pctsAsian.put(sample, pctAsian);
			pctsEuropean.put(sample, pctEuropean);
			if (pctEuropean > 0.95) {
				imputedRaces.put(sample, RACE.WHITE);
			} else if (pctAsian > 0.999) {
				imputedRaces.put(sample, RACE.ASIAN);
			} else if (pctAfrican / Math.max(pctAsian, 0.015) > 10) {
				imputedRaces.put(sample, RACE.AA);
			} else {
				imputedRaces.put(sample, RACE.HISPANIC);
			}

		}

		log.report("Writing Results");

		SampleData sampleData = proj.getSampleData(false);

		Map<String, String> dataToAdd = Maps.newHashMap();
		for (Sample sample : samples) {
			dataToAdd.put(sampleData.lookup(sample.getFidIid())[0],
										Joiner.on('\t').join(imputedRaces.get(sample).getSampleDataClassNum(),
																				 pctsAfrican.get(sample), pctsAsian.get(sample),
																				 pctsEuropean.get(sample)));
		}


		sampleData.addData(dataToAdd, "DNA", IMPUTED_RACE_SAMPLE_DATA_HEADERS, ".", "\t", log);

		writer = Files.getAppropriateWriter(outFile);
		writer.println(ArrayUtils.toStr(CORRECTED_PCS_HEADER));
		Map<RACE, PrintWriter> raceWriters = generateRaceListWriters(outFile);

		for (Sample sample : samples) {
			writer.println(Joiner.on('\t').join(sample.getFid(), sample.getIid(), sample.getPc1(),
																					sample.getPc2(), pctsAfrican.get(sample),
																					pctsAsian.get(sample), pctsEuropean.get(sample)));
			raceWriters.get(imputedRaces.get(sample)).println(sample.getFidIid());
		}
		writer.close();
		for (PrintWriter raceWriter : raceWriters.values()) {
			raceWriter.close();
		}

	}

	private void writeStepPCs(String outFile) {
		PrintWriter writer = Files.getAppropriateWriter(outFile);
		writer.println(ArrayUtils.toStr(STEP_PCS_HEADER));

		for (Sample sample : samples) {
			writer.println(Joiner.on('\t').join(sample.getFid(), sample.getIid(), sample.getPc1(),
																					sample.getPc2(),
																					(1.0 - sample.getPc1() + sample.getPc2())));
		}
		writer.close();
	}

	private void shiftEuropeansOutOfFirstQuadrant() {
		log.report("Setting European Seeds to below 0,0");


		double maxEuropeanPC1 = Collections.max(eurSeeds, PC1_COMPARATOR).getPc1();
		double maxEuropeanPC2 = Collections.max(eurSeeds, PC2_COMPARATOR).getPc2();

		for (Sample sample : samples) {
			sample.setPc1(sample.getPc1() - maxEuropeanPC1);
			sample.setPc2(sample.getPc2() - maxEuropeanPC2);
		}

	}

	private void forceRightAngleBetweenAsiansAndAfricans() {
		log.report("Converting to Polar Coordinates");

		final Map<Sample, Double> polarThetas = Maps.newHashMap();
		Map<Sample, Double> polarRs = Maps.newHashMap();
		for (Sample sample : samples) {
			polarThetas.put(sample, Maths.cartesianToPolarTheta(sample.getPc1(), sample.getPc2()));
			polarRs.put(sample, Maths.cartesianToPolarR(sample.getPc1(), sample.getPc2()));
		}

		log.report("Forcing 90 degree difference between Asians and Africans");

		double maxAfricanTheta = polarThetas.get(Collections.max(afrSeeds, new Comparator<Sample>() {
			@Override
			public int compare(Sample sample1, Sample sample2) {
				return Double.compare(polarThetas.get(sample1), polarThetas.get(sample2));
			}
		}));

		double minAsianTheta = polarThetas.get(Collections.min(asianSeeds, new Comparator<Sample>() {
			@Override
			public int compare(Sample sample1, Sample sample2) {
				return Double.compare(polarThetas.get(sample1), polarThetas.get(sample2));
			}
		}));

		for (Sample sample : samples) {
			double theta = polarThetas.get(sample);
			if (theta <= minAsianTheta && theta >= maxAfricanTheta) {
				// If in target range, scale such that minAsianTheta becomes pi/2 (90) and maxAfricanTheta
				// becomes 0 (0)
				theta = ((theta - maxAfricanTheta) / (minAsianTheta - maxAfricanTheta)) * (Math.PI / 2.0);
			} else {
				// If out of target range, scale opposite such that minAsianTheta becomes pi/2 (90) and
				// maxAfricanTheta becomes 2*pi (360)
				if (theta < maxAfricanTheta) {
					// Values from pi (180) to 2*pi (360) are expressed from -pi (-180) to 0 (0), express all
					// values under maxAfricanTheta as >pi (180) for scaling
					theta = 2.0 * Math.PI + theta;
				}
				theta = ((theta - minAsianTheta) / (2.0 * Math.PI + maxAfricanTheta - minAsianTheta))
								* (3.0 * Math.PI / 2.0) + (Math.PI / 2.0);
			}
			polarThetas.put(sample, theta);
		}

		log.report("Converting back to Cartesian");

		for (Sample sample : samples) {
			double theta = polarThetas.get(sample);
			double r = polarRs.get(sample);
			sample.setPc1(Maths.polarToCartesianX(theta, r));
			sample.setPc2(Maths.polarToCartesianY(theta, r));
		}
	}

	private void scaleAfricanAndAsianTo1() {
		log.report("Forcing scale where African and Asian are at 1.0");

		double minAfricanPC1 = Collections.min(afrSeeds, PC1_COMPARATOR).getPc1();
		double minAsianPC2 = Collections.min(asianSeeds, PC2_COMPARATOR).getPc2();

		for (Sample sample : samples) {
			sample.setPc1(sample.getPc1() / minAfricanPC1);
			sample.setPc2(sample.getPc2() / minAsianPC2);
		}
	}

	private boolean maybeSwapPCs() {
		double europeanMeanPC1 = meanPC1(eurSeeds);
		double europeanMeanPC2 = meanPC2(eurSeeds);
		double africanMeanPC1 = meanPC1(afrSeeds);
		double africanMeanPC2 = meanPC2(afrSeeds);
		double asianMeanPC1 = meanPC1(asianSeeds);
		double asianMeanPC2 = meanPC2(asianSeeds);

		if (Math.abs(africanMeanPC1) > Math.abs(asianMeanPC1)
				&& Math.abs(asianMeanPC2) > Math.abs(africanMeanPC2)
				&& Math.abs(africanMeanPC1) > Math.abs(europeanMeanPC1)
				&& Math.abs(asianMeanPC2) > Math.abs(europeanMeanPC2)) {
			// PC1 = African, PC2 = Asian
		} else if (Math.abs(asianMeanPC1) > Math.abs(africanMeanPC1)
							 && Math.abs(africanMeanPC2) > Math.abs(asianMeanPC2)
							 && Math.abs(asianMeanPC1) > Math.abs(europeanMeanPC1)
							 && Math.abs(africanMeanPC2) > Math.abs(europeanMeanPC2)) {
			// PC1 = Asian, PC2 = African
			for (Sample sample : samples) {
				double pc1 = sample.getPc1();
				double pc2 = sample.getPc2();
				sample.setPc1(pc2);
				sample.setPc2(pc1);
			}

			double tempMean;

			tempMean = europeanMeanPC1;
			europeanMeanPC1 = europeanMeanPC2;
			europeanMeanPC2 = tempMean;

			tempMean = africanMeanPC1;
			africanMeanPC1 = africanMeanPC2;
			africanMeanPC2 = tempMean;

			tempMean = asianMeanPC1;
			asianMeanPC1 = asianMeanPC2;
			asianMeanPC2 = tempMean;
		} else {
			log.reportError("PC1 and PC2 do not appear to predict African and Asian, race cannot be imputed");
			return false;
		}

		if (africanMeanPC1 < 0) {
			for (Sample sample : samples) {
				sample.setPc1(sample.getPc1() * -1.0);
			}
		}

		if (asianMeanPC2 < 0) {
			for (Sample sample : samples) {
				sample.setPc2(sample.getPc2() * -1.0);
			}
		}
		return true;
	}

	private static double meanPC1(Collection<Sample> samples) {
		double sum = 0.0;
		int count = 0;
		for (Sample sample : samples) {
			double pc1 = sample.getPc1();
			if (Double.isFinite(pc1)) {
				count++;
				sum += pc1;
			}
		}
		if (count == 0)
			return Double.NaN;
		return sum / count;
	}

	private static double meanPC2(Collection<Sample> samples) {
		double sum = 0.0;
		int count = 0;
		for (Sample sample : samples) {
			double pc2 = sample.getPc2();
			if (Double.isFinite(pc2)) {
				count++;
				sum += pc2;
			}
		}
		if (count == 0)
			return Double.NaN;
		return sum / count;
	}


	private static Map<RACE, PrintWriter> generateRaceListWriters(String resultFile) {
		Map<RACE, String> raceFilenames = raceListFilenames(resultFile);
		Map<RACE, PrintWriter> raceWriters = Maps.newEnumMap(RACE.class);
		for (Map.Entry<RACE, String> raceFilename : raceFilenames.entrySet()) {
			RACE race = raceFilename.getKey();
			String filename = raceFilename.getValue();
			raceWriters.put(race, Files.getAppropriateWriter(filename));
		}
		return raceWriters;
	}

	private static Map<RACE, String> raceListFilenames(String resultFile) {
		Map<RACE, String> raceFilenames = Maps.newEnumMap(RACE.class);
		for (RACE race : RACE.values()) {
			String filename = formRaceListFilename(race, resultFile);
			raceFilenames.put(race, filename);
		}
		return raceFilenames;
	}

	public static String formRaceListFilename(RACE race, String resultFile) {
		String raceString = ext.replaceWithLinuxSafeCharacters(race.description);
		return ext.rootOf(resultFile, false) + "_" + raceString + "s.dat";
	}

	private static int countFounders(String plinkroot, String keepFile) {
		int founders = 0;
		Hashtable<String, String> plinkFam = HashVec.loadFileToHashString(plinkroot + ".fam",
																																			new int[] {0, 1},
																																			new int[] {2, 3}, false, "\t",
																																			false, false, false);
		Set<String> keeps = keepFile == null ? plinkFam.keySet()
																				 : HashVec.loadFileToHashSet(keepFile, new int[] {0, 1},
																																		 "\t", false);
		for (Entry<String, String> famEntry : plinkFam.entrySet()) {
			if (famEntry.getValue().equals("0\t0") && keeps.contains(famEntry.getKey())) {
				founders++;
			}
		}
		return founders;
	}

	public static void freqsByRace(String resultFile, String plinkroot, String outFile, Logger log) {
		String dir = ext.parseDirectoryOfFile(plinkroot, true);
		plinkroot = ext.rootOf(plinkroot);

		String overallFrqFile = ext.rootOf(resultFile, false) + "_all.frq";
		CmdLine.runDefaults("plink2 --noweb --bfile " + plinkroot + " --freq" + " --out "
												+ ext.rootOf(overallFrqFile, false), dir);
		String[] header = Files.getHeaderOfFile(overallFrqFile, log);
		int key = ext.indexOfStr("SNP", header);
		int[] targets = new int[] {ext.indexOfStr("A1", header), ext.indexOfStr("A2", header),
															 ext.indexOfStr("MAF", header)};
		Hashtable<String, String> overallFreq = HashVec.loadFileToHashString(overallFrqFile, key,
																																				 targets, "\t", true);

		Map<RACE, String> raceListFiles = raceListFilenames(resultFile);
		@SuppressWarnings("unchecked")
		Map<RACE, Map<String, String>> raceFreqs = Maps.newHashMap();

		PrintWriter writer = Files.getAppropriateWriter(outFile);
		writer.print("SNP\tA1\tA2\tOverall A1F (n=" + countFounders(dir + plinkroot, null) + ")");

		for (Map.Entry<RACE, String> raceListFileEntry : raceListFiles.entrySet()) {
			RACE race = raceListFileEntry.getKey();
			String raceListFile = raceListFileEntry.getValue();
			String raceFrqFile = ext.rootOf(raceListFile, false) + ".frq";
			CmdLine.runDefaults("plink2 --noweb --bfile " + plinkroot + " --keep " + raceListFile
													+ " --freq" + " --out " + ext.rootOf(raceFrqFile, false), dir);
			header = Files.getHeaderOfFile(raceFrqFile, log);
			key = ext.indexOfStr("SNP", header);
			targets = new int[] {ext.indexOfStr("A1", header), ext.indexOfStr("A2", header),
													 ext.indexOfStr("MAF", header)};
			raceFreqs.put(race, HashVec.loadFileToHashString(raceFrqFile, key, targets, "\t", true));

			writer.print("\t" + race + " A1F (n=" + countFounders(dir + plinkroot, raceListFile) + ")");
		}
		writer.println();

		for (Entry<String, String> overallEntry : overallFreq.entrySet()) {
			String marker = overallEntry.getKey();
			String[] data = overallEntry.getValue().split("\t");
			String A1 = data[0];
			String A2 = data[1];;
			String overallMAF;
			try {
				overallMAF = Double.toString(ext.roundToSignificantFigures(Double.parseDouble(data[2]), 4));
			} catch (NumberFormatException nfe) {
				log.reportError("Invalid MAF (" + data[2] + ") for SNP '" + marker + "'");
				overallMAF = ".";
			}

			writer.print(marker + "\t" + A1 + "\t" + A2 + "\t" + overallMAF);

			for (Map.Entry<RACE, Map<String, String>> raceFreqEntry : raceFreqs.entrySet()) {
				RACE race = raceFreqEntry.getKey();
				Map<String, String> raceFreq = raceFreqEntry.getValue();
				String a1f;
				if (marker == null) {
					log.reportError("SNP '" + marker + "' not found for " + raceListFiles.get(race));
					a1f = ".";
				} else {
					String[] raceData = raceFreq.get(marker).split("\t");
					String raceA1 = raceData[0];
					String raceA2 = raceData[1];
					try {
						double raceMaf = Double.parseDouble(raceData[2]);
						if (A1.equals(raceA1) && A2.equals(raceA2)) {
							a1f = Double.toString(ext.roundToSignificantFigures(raceMaf, 4));
						} else if (A1.equals(raceA2) && A2.equals(raceA1)) {
							a1f = Double.toString(ext.roundToSignificantFigures(1.0 - raceMaf, 4));
						} else {
							log.reportError("Alleles for SNP '" + marker + "' and " + raceListFiles.get(race)
															+ " (" + raceA1 + ", " + raceA2 + ") do not match overall alleles ("
															+ A1 + ", " + A2 + " )");
							a1f = ".";
						}
					} catch (NumberFormatException nfe) {
						log.reportError("Invalid MAF (" + raceData[2] + ") for SNP '" + marker + "' and "
														+ raceListFiles.get(race));
						a1f = ".";
					}
				}
				writer.print("\t" + a1f);
			}
			writer.println();
		}
		writer.close();
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj = null;
		String inFile = null;

		String usage = "\n" + "cnv.analysis.pca.PCImputeRace requires 2 arguments\n";
		usage += "   (1) Project Filename (i.e. proj=" + inFile + " (default))\n" + "";
		usage += "   (2) Input Filename (i.e. inFile=" + inFile + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				proj = new Project(arg.split("=")[1], false);
				numArgs--;
			} else if (arg.startsWith("inFile=")) {
				inFile = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			String[] input = HashVec.loadFileToStringArray(inFile, true, null, false);
			Set<Sample> samples = Sets.newHashSet();
			Set<Sample> europeans = Sets.newHashSet();
			Set<Sample> africans = Sets.newHashSet();
			Set<Sample> asians = Sets.newHashSet();

			for (int i = 0; i < input.length; i++) {
				String[] line = input[i].split("\t");
				String fid = line[0];
				String iid = line[1];
				double pc1 = Double.parseDouble(line[2]);
				double pc2 = Double.parseDouble(line[3]);
				Sample sample = new Sample(fid, iid, pc1, pc2);
				samples.add(sample);
				try {
					int race = Integer.parseInt(line[4]);
					switch (race) {
						case 1:
							europeans.add(sample);
							break;
						case 2:
							africans.add(sample);
							break;
						case 3:
						case 4:
							asians.add(sample);
							break;
						default:
							break;
					}
				} catch (NumberFormatException nfe) {
					// If not hapmap, don't add to a hapmap set
				}
			}

			PCImputeRace raceChecker = new PCImputeRace(proj, samples, europeans, africans, asians,
																									new Logger());
			raceChecker.correctPCsToRace(ext.rootOf(inFile, false) + "_Corrected_PCS.mds");

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
