package org.genvisis.cnv.analysis.pca;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;

import com.google.common.primitives.Doubles;

/**
 * @author lane0212 Class to compute percent race with genotyping pcs, hijacking the
 *         {@link PrincipalComponentsResiduals} methods
 */
public class PCPopulation {
	private final Project proj;
	private final PrincipalComponentsResiduals pResiduals;
	private final int clusterComponents;
	private final VcfPopulation vpop;
	private Population[] populations;
	// private boolean verbose;
	private final Logger log;

	public PCPopulation(Project proj, PrincipalComponentsResiduals pResiduals, VcfPopulation vpop,
											int clusterComponents, boolean verbose) {
		super();
		this.proj = proj;
		this.pResiduals = pResiduals;
		this.vpop = vpop;
		log = proj.getLog();
		this.clusterComponents = clusterComponents;
	}

	public TestSampleDistances[] computeDistance() {
		populations = getPopulations(pResiduals, vpop);
		determinePCClusters();
		Hashtable<String, TestSampleDistances> dists = new Hashtable<String, TestSampleDistances>();
		Hashtable<String, TestPopulationDistances> pdists =
																											new Hashtable<String, TestPopulationDistances>();
		for (int i = 0; i < populations.length; i++) {
			if (populations[i].isTest()) {
				for (int j = 0; j < populations.length; j++) {
					if (j != i) {
						populations[j].computeDistance(populations[i], clusterComponents, pResiduals, log);
						ArrayList<TestSampleDistance> popDists = populations[j].getDistances();
						for (TestSampleDistance testSampleDistance : popDists) {
							String sample = testSampleDistance.getSample();
							if (!dists.containsKey(sample)) {
								dists.put(testSampleDistance.getSample(), new TestSampleDistances(sample));
							}
							dists.get(sample).add(sample, testSampleDistance.getDistance(),
																		populations[j].getName());
						}
					}
				}
			} else {
				pdists.put(populations[i].getName(), new TestPopulationDistances(populations[i]));
				for (int j = 0; j < populations.length; j++) {
					if (!populations[j].isTest()) {
						pdists.get(populations[i].getName()).add(populations[j]);
					}
				}
			}
		}
		ArrayList<TestSampleDistances> finalizedDists =
																									new ArrayList<PCPopulation.TestSampleDistances>(dists.size());
		for (String ind : dists.keySet()) {
			TestSampleDistances curDist = dists.get(ind);
			curDist.computeNormDist(pdists, log);
			finalizedDists.add(curDist);

		}
		return finalizedDists.toArray(new TestSampleDistances[finalizedDists.size()]);
	}

	private void determinePCClusters() {
		for (int i = 0; i < populations.length; i++) {
			if (!populations[i].isTest()) {
				populations[i].determineClusters(clusterComponents, pResiduals, proj.getLog());
			}
		}
	}

	private static Population[] getPopulations(	PrincipalComponentsResiduals pResiduals,
																							VcfPopulation vpop) {
		ArrayList<Population> pops = new ArrayList<Population>(vpop.getSuperPop().size());
		if (vpop.valid()) {
			Hashtable<String, Integer> sampPCs = pResiduals.getSamplesInPc();
			for (String sPop : vpop.getSuperPop().keySet()) {
				boolean[] curMask = new boolean[sampPCs.size()];
				Arrays.fill(curMask, false);
				HashSet<String> vpopPCs = new HashSet<String>();
				for (String ind : vpop.getSuperPop().get(sPop)) {
					if (sampPCs.containsKey(ind)) {
						vpopPCs.add(ind);
						curMask[sampPCs.get(ind)] = true;
					}
				}
				if (Array.booleanArraySum(curMask) > 0) {
					pops.add(new Population(sPop, curMask, vpopPCs,
																	sPop.equals(VcfPopulation.DETERMINE_ANCESTRY)));
				} else {
					vpop.getLog().reportTimeWarning("Skipping population " + sPop + ", no PCs were found");
				}
			}
		}
		return pops.toArray(new Population[pops.size()]);

	}

	private static class TestPopulationDistances {
		private final Population population;
		private final ArrayList<Double> otherDistances;
		private final ArrayList<String> otherPopulations;

		public TestPopulationDistances(Population population) {
			super();
			this.population = population;
			otherDistances = new ArrayList<Double>();
			otherPopulations = new ArrayList<String>();
		}

		// public Population getPopulation() {
		// return population;
		// }

		public void add(Population otherPopulation) {
			if (otherPopulation.isTest()) {
				throw new IllegalStateException("This method must not be called on the test case");
			} else {
				double tmpDist = 0;
				for (int i = 0; i < population.getClusterCenters().length; i++) {
					tmpDist += Math.pow(population.getClusterCenters()[i]
															- otherPopulation.getClusterCenters()[i], 2);
				}
				tmpDist = Math.sqrt(tmpDist);
				otherDistances.add(tmpDist);
				otherPopulations.add(otherPopulation.getName());
			}
		}

		// public ArrayList<Double> getOtherDistances() {
		// return otherDistances;
		// }
		//
		// public void setOtherDistances(ArrayList<Double> otherDistances) {
		// this.otherDistances = otherDistances;
		// }
		//
		// public ArrayList<String> getOtherPopulations() {
		// return otherPopulations;
		// }
		//
		// public void setOtherPopulations(ArrayList<String> otherPopulations) {
		// this.otherPopulations = otherPopulations;
		// }

	}

	private static class TestSampleDistances {
		private final String sample;
		private ArrayList<Double> distances;
		private final ArrayList<Double> normDistances;

		private ArrayList<String> populations;

		public TestSampleDistances(String sample) {
			super();
			this.sample = sample;
			distances = new ArrayList<Double>();
			populations = new ArrayList<String>();
			normDistances = new ArrayList<Double>();

		}

		public String getSample() {
			return sample;
		}

		public void computeNormDist(Hashtable<String, TestPopulationDistances> pdists, Logger log) {
			double minDist = Double.MAX_VALUE;
			// String minDistPop = null;
			double[] allDists = new double[populations.size()];
			for (int i = 0; i < populations.size(); i++) {
				// String curPop = populations.get(i);

				// TestPopulationDistances curPopulationDistance = pdists.get(curPop);
				double distToCur = distances.get(i);
				allDists[i] = distToCur;

				if (distToCur < minDist) {
					minDist = distToCur;
					// minDistPop = curPop;
				}

				// for (int j = 0; j < populations.size(); j++) {
				// if (i != j) {
				// if (!curPopulationDistance.getOtherPopulations().get(j).equals(populations.get(j))) {
				// log.reportTimeError("Mismatched populations");
				// } else {
				//
				// double diff = curPopulationDistance.getOtherDistances().get(j)-distances.get(j);
				//
				// distToCur+=diff;
				// // System.out.println(sample + "\t" + curPop + "\t" + distToCur + "\t" +
				// curPopulationDistance.getOtherPopulations().get(j) + "\t" +
				// curPopulationDistance.getOtherDistances().get(j));
				// if (sample.startsWith("F10607D")) {
				// // System.out.println("SDF\t" + sample + "\tcurPop " + curPop + " > " +
				// populations.get(j) + "\t" + distances.get(j) + "\t" +
				// curPopulationDistance.getOtherDistances().get(j) + "\t" + diff);
				// //System.out.println(curPopulationDistance.getPopulation().getName() + "\t" +
				// curPopulationDistance.getOtherDistances().get(j)+"\t"+populations.get(j) + "\t" +
				// distances.get(j));
				// System.out.println("The difference from the
				// "+curPopulationDistance.getPopulation().getName() + " cluster to the
				// "+curPopulationDistance.getOtherPopulations().get(j)+ " cluster is "+
				// curPopulationDistance.getOtherDistances().get(j));
				// System.out.println("The distance from this sample to the "+populations.get(j)+" cluster
				// is"+distances.get(j));
				// System.out.println("The difference in distances is "+diff);
				// System.out.println();
				//
				// }
				//
				// // distToCur+=
				// //
				// // curPopulationDistance.getOtherPopulations().get(i);
				// }
				// }
				// }
				//

			}
			// System.out.println(Array.toStr(allDists));
			// int[] order = Sort.quicksort(allDists);
			// System.out.println(Array.toStr(order));
			// System.out.println(Array.toStr(allDists));

			// if (sample.startsWith("F10626D")) {
			// System.out.println(sample + "\t" + minDist + "\t" + minDistPop);
			// }
			// TestPopulationDistances clusterPop = pdists.get(minDistPop);
			for (int i = 0; i < populations.size(); i++) {
				double normDist = Double.NaN;

				// String curPop = populations.get(i);
				// if (curPop.equals(minDistPop)) {

				// if (sample.startsWith("F10626D")) {
				// System.out.println("The distance from this sample to the " + populations.get(i) + "
				// cluster is" + distances.get(i));
				// System.out.println(clusterPop.getOtherPopulations().get(i) + "\t" +
				// clusterPop.getOtherDistances().get(i));
				// System.out.println(clusterPop.getOtherDistances().get(i) - distances.get(i));
				// }
				// TestPopulationDistances clusterPop = pdists.get(curPop);
				boolean[] mask = new boolean[populations.size()];
				Arrays.fill(mask, true);
				mask[i] = false;
				double sumTotalDists = Array.sum(Array.subArray(getDistances(), mask));
				normDist = distances.get(i) / sumTotalDists;
				// System.out.println("The difference from the " + clusterPop.getPopulation().getName() + "
				// cluster to the " + clusterPop.getOtherPopulations().get(i) + " cluster is " +
				// clusterPop.getOtherDistances().get(i));
				// System.out.println("The distance from this sample to the " + populations.get(j) + "
				// cluster is" + distances.get(j));
				// System.out.println("The difference in distances is " + diff + " with percent diff ");
				// System.out.println("The current percentage " + curPop + " is " + (1 - normDist));


				// for (int j = 0; j < populations.size(); j++) {
				// if (i != j) {
				// double diff = distances.get(j) - clusterPop.getOtherDistances().get(j);
				// normDist += (diff / Array.sum(getDistances()));
				//
				// // double percentCurrent = diff
				// <0?(diff/distances.get(j)):(diff/clusterPop.getOtherDistances().get(j));
				// // normDist+=percentCurrent;
				//
				// // if(percentCurrent<1){
				// // normDist+=Math.abs(percentCurrent);
				// // if(percentCurrent>0){
				// // normDist += percentCurrent;
				// // }else{
				// // normDist +=1;
				// // }
				// // }
				//
				// if (sample.startsWith("F10626M")) {
				// System.out.println("The difference from the " + clusterPop.getPopulation().getName() + "
				// cluster to the " + clusterPop.getOtherPopulations().get(j) + " cluster is " +
				// clusterPop.getOtherDistances().get(j));
				//
				// System.out.println("The distance from this sample to the " + populations.get(j) + "
				// cluster is" + distances.get(j));
				// System.out.println("The difference in distances is " + diff + " with percent diff ");
				// // System.out.println(clusterPop.getOtherPopulations().get(j) + "\t" +
				// clusterPop.getOtherDistances().get(j));
				// System.out.println("The current percentage " + curPop + " is " + (1 - normDist));
				// }
				// }
				// }

				// if (sample.startsWith("D180")) {
				// System.out.println("The distance from this sample to the " + populations.get(i) + "
				// cluster is" + distances.get(i));
				// System.out.println(clusterPop.getOtherPopulations().get(i) + "\t" +
				// clusterPop.getOtherDistances().get(i));
				// System.out.println("The current percentage is" + normDist);
				// }
				// }
				normDistances.add(Math.max(Math.min(1 - normDist, 1), 0));
			}

		}

		public double[] getDistances() {
			return Doubles.toArray(distances);
		}

		public double[] getNormDistances() {
			return Doubles.toArray(normDistances);
		}

		public String[] getPopulations() {
			return populations.toArray(new String[populations.size()]);
		}

		public void add(String thisSample, double distance, String population) {
			if (!sample.equals(thisSample)) {
				distances = null;
				populations = null;
				throw new IllegalStateException("Un-matched samples " + sample + " and " + thisSample);
			} else {
				distances.add(distance);
				populations.add(population);
			}
		}

	}

	private static class TestSampleDistance {
		private final String sample;
		private final double distance;

		public TestSampleDistance(String sample, double distance) {
			super();
			this.sample = sample;
			this.distance = distance;
		}

		public String getSample() {
			return sample;
		}

		public double getDistance() {
			return distance;
		}

	}

	private static class Population {
		private final String name;
		private final boolean[] pcMatchedMask;
		private double[] clusterCenters;
		private ArrayList<TestSampleDistance> distances;
		private final Set<String> samples;
		private final boolean isTest;

		public Population(String name, boolean[] pcMatchedMask, Set<String> samples, boolean isTest) {
			super();
			this.name = name;
			this.pcMatchedMask = pcMatchedMask;
			this.isTest = isTest;
			this.samples = samples;
		}

		public String getName() {
			return name;
		}

		public boolean isTest() {
			return isTest;
		}

		// public boolean[] getPcMatchedMask() {
		// return pcMatchedMask;
		// }

		public double[] getClusterCenters() {
			return clusterCenters;
		}

		public Set<String> getSamples() {
			return samples;
		}

		public ArrayList<TestSampleDistance> getDistances() {
			return distances;
		}

		public void computeDistance(Population testPop, int clusterComponents,
																PrincipalComponentsResiduals pResiduals, Logger log) {
			if (!testPop.isTest()) {
				throw new IllegalStateException("This method must be called on the test case");
			} else {
				Set<String> testSamples = testPop.getSamples();
				distances = new ArrayList<TestSampleDistance>(testSamples.size());
				double[][] trimmedPcs = PrincipalComponentsResiduals.trimPcBasis(	clusterComponents,
																																					pResiduals.getPcBasis(),
																																					log);
				for (String testSamp : testSamples) {
					int index = pResiduals.getSamplesInPc().get(testSamp);
					double distTmp = 0;
					for (int i = 0; i < trimmedPcs.length; i++) {
						distTmp += Math.pow(trimmedPcs[i][index] - clusterCenters[i], 2);
					}
					distTmp = Math.sqrt(distTmp);
					distances.add(new TestSampleDistance(testSamp, distTmp));
				}
			}
		}

		public void determineClusters(int clusterComponents, PrincipalComponentsResiduals pResiduals,
																	Logger log) {
			if (isTest) {
				throw new IllegalStateException("This method cannot be called on the test case");
			} else {
				clusterCenters = new double[clusterComponents];
				double[][] trimmedPcs = PrincipalComponentsResiduals.trimPcBasis(	clusterComponents,
																																					pResiduals.getPcBasis(),
																																					log);
				for (int i = 0; i < clusterCenters.length; i++) {
					clusterCenters[i] = Array.median(Array.subArray(trimmedPcs[i], pcMatchedMask));
					log.reportTimeInfo(name + " Cluster " + i + " : " + clusterCenters[i]);
				}
			}
		}
	}

	public static void test(Project proj, String genoPCfile, int numComponents) {
		PrincipalComponentsResiduals pResiduals = new PrincipalComponentsResiduals(	genoPCfile,
																																								Files.getHeaderOfFile(genoPCfile,
																																																			proj.getLog()).length
																																														- 2,
																																								proj.getLog());
		VcfPopulation vpop = VcfPopulation.load(proj.SAMPLE_DATA_FILENAME.getValue(),
																						POPULATION_TYPE.PC_ANCESTRY, proj.getLog());
		vpop.report();
		if (vpop.valid()) {
			proj.getLog()
					.reportTimeInfo("Detected " + (vpop.getSuperPop().size() - 1) + " seed populations");
			if (numComponents < (vpop.getSuperPop().size() - 2)) {
				proj.getLog()
						.reportTimeWarning("Usually "	+ (vpop.getSuperPop().size() - 2)
																+ " populations should be clustered across at least "
																+ (vpop.getSuperPop().size() - 3) + " principal components");
			}
			PCPopulation pcPopulation = new PCPopulation(proj, pResiduals, vpop, numComponents, true);
			TestSampleDistances[] testSampleDistances = pcPopulation.computeDistance();
			String output = ext.rootOf(genoPCfile, false) + ".ancestry.txt";
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				writer.println("DNA\tSTUDY\t"	+ Array.toStr(testSampleDistances[0].getPopulations()) + "\t"
												+ Array.toStr(testSampleDistances[0].getPopulations()));
				for (TestSampleDistances testSampleDistance : testSampleDistances) {
					String sample = testSampleDistance.getSample();
					writer.println(sample	+ "\t" + vpop.getPopulationForInd(sample, RETRIEVE_TYPE.SUB)[0]
													+ "\t" + Array.toStr(testSampleDistance.getDistances()) + "\t"
													+ Array.toStr(testSampleDistance.getNormDistances()));
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + output);
				proj.getLog().reportException(e);
			}
		}

		// for (int i = 0; i < testSampleDistances.length; i++) {
		// System.out.println(testSampleDistances[i].getSample()+"\t"+Array.toStr(testSampleDistances[i].getPopulations()));
		// System.out.println(testSampleDistances[i].getSample()+"\t"+Array.toStr(testSampleDistances[i].getDistances()));
		// }
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String genoPCfile = null;
		int numComponents = 2;
		// String logfile = null;
		// Logger log;

		String usage = "\n" + "cnv.analysis.pca.PCPopulation requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) full path to a genotype pc file (i.e. file=" + filename + " (default))\n" + "";
		usage += "   (3) number of components (i.e. numComp=" + numComponents + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("pc=")) {
				genoPCfile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("numComp=")) {
				numComponents = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("log=")) {
				// logfile = args[i].split("=")[1];
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
			// log = new Logger(logfile);
			test(new Project(filename, false), genoPCfile, numComponents);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
