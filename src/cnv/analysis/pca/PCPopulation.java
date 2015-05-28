package cnv.analysis.pca;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import common.Array;
import common.Files;
import common.Logger;
import common.ext;
import cnv.filesys.Project;
import cnv.var.SampleData;

/**
 * @author lane0212 Class to compute percent race with genotyping pcs, hijacking the {@link PrincipalComponentsResiduals} methods
 */
public class PCPopulation {
	private Project proj;
	private PrincipalComponentsResiduals pResiduals;
	private int clusterComponents;
	private VcfPopulation vpop;
	private Population[] populations;
	private boolean verbose;
	private Logger log;

	public PCPopulation(Project proj, PrincipalComponentsResiduals pResiduals, VcfPopulation vpop, int clusterComponents, boolean verbose) {
		super();
		this.proj = proj;
		this.pResiduals = pResiduals;
		this.vpop = vpop;
		this.log = proj.getLog();
		this.clusterComponents = clusterComponents;
	}

	public TestSampleDistances[] computeDistance() {
		this.populations = getPopulations(pResiduals, vpop);
		determinePCClusters();
		Hashtable<String, TestSampleDistances> dists = new Hashtable<String, TestSampleDistances>();
		Hashtable<String, TestPopulationDistances> pdists = new Hashtable<String, TestPopulationDistances>();
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
							dists.get(sample).add(sample, testSampleDistance.getDistance(), populations[j].getName());
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
		ArrayList<TestSampleDistances> finalizedDists = new ArrayList<PCPopulation.TestSampleDistances>(dists.size());
		for (String ind : dists.keySet()) {
			TestSampleDistances curDist = dists.get(ind);

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

	private static Population[] getPopulations(PrincipalComponentsResiduals pResiduals, VcfPopulation vpop) {
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
					pops.add(new Population(sPop, curMask, vpopPCs, sPop.equals(VcfPopulation.DETERMINE_ANCESTRY)));
				} else {
					vpop.getLog().reportTimeWarning("Skipping population " + sPop + ", no PCs were found");
				}
			}
		}
		return pops.toArray(new Population[pops.size()]);

	}

	private static class TestPopulationDistances {
		private Population population;
		private ArrayList<Double> otherDistances;
		private ArrayList<String> otherPopulations;

		public TestPopulationDistances(Population population) {
			super();
			this.population = population;
			this.otherDistances = new ArrayList<Double>();
			this.otherPopulations = new ArrayList<String>();
		}

		public Population getPopulation() {
			return population;
		}

		public void add(Population otherPopulation) {
			if (otherPopulation.isTest()) {
				throw new IllegalStateException("This method must not be called on the test case");
			} else {
				double tmpDist = 0;
				for (int i = 0; i < population.getClusterCenters().length; i++) {
					tmpDist += Math.pow(population.getClusterCenters()[i] - otherPopulation.getClusterCenters()[i], 2);
				}
				tmpDist = Math.sqrt(tmpDist);
				otherDistances.add(tmpDist);
				otherPopulations.add(otherPopulation.getName());
			}
		}

		public ArrayList<Double> getOtherDistances() {
			return otherDistances;
		}

		public void setOtherDistances(ArrayList<Double> otherDistances) {
			this.otherDistances = otherDistances;
		}

		public ArrayList<String> getOtherPopulations() {
			return otherPopulations;
		}

		public void setOtherPopulations(ArrayList<String> otherPopulations) {
			this.otherPopulations = otherPopulations;
		}

	}

	private static class TestSampleDistances {
		private String sample;
		private ArrayList<Double> distances;
		private ArrayList<Double> normDistances;

		private ArrayList<String> populations;

		public TestSampleDistances(String sample) {
			super();
			this.sample = sample;
			this.distances = new ArrayList<Double>();
			this.populations = new ArrayList<String>();
			this.normDistances = new ArrayList<Double>();

		}

		public String getSample() {
			return sample;
		}

		public void computeNormDist(Hashtable<String, TestPopulationDistances> pdists) {
			for (int i = 0; i < populations.size(); i++) {
				String curPop = populations.get(i);
				TestPopulationDistances curPopulationDistance = pdists.get(curPop);
				double distToCur = distances.get(i);
				
			}

		}

		public double[] getDistances() {
			return Array.toDoubleArray(distances);
		}

		public String[] getPopulations() {
			return populations.toArray(new String[populations.size()]);
		}

		public void add(String thisSample, double distance, String population) {
			if (!sample.equals(thisSample)) {
				this.distances = null;
				this.populations = null;
				throw new IllegalStateException("Un-matched samples " + sample + " and " + thisSample);
			} else {
				distances.add(distance);
				populations.add(population);
			}
		}

	}

	private static class TestSampleDistance {
		private String sample;
		private double distance;

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
		private String name;
		private boolean[] pcMatchedMask;
		private double[] clusterCenters;
		private ArrayList<TestSampleDistance> distances;
		private Set<String> samples;
		private boolean isTest;

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

		public boolean[] getPcMatchedMask() {
			return pcMatchedMask;
		}

		public double[] getClusterCenters() {
			return clusterCenters;
		}

		public Set<String> getSamples() {
			return samples;
		}

		public ArrayList<TestSampleDistance> getDistances() {
			return distances;
		}

		public void computeDistance(Population testPop, int clusterComponents, PrincipalComponentsResiduals pResiduals, Logger log) {
			if (!testPop.isTest()) {
				throw new IllegalStateException("This method must be called on the test case");
			} else {
				Set<String> testSamples = testPop.getSamples();
				this.distances = new ArrayList<TestSampleDistance>(testSamples.size());
				double[][] trimmedPcs = PrincipalComponentsResiduals.trimPcBasis(clusterComponents, pResiduals.getPcBasis(), log);
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

		public void determineClusters(int clusterComponents, PrincipalComponentsResiduals pResiduals, Logger log) {
			if (isTest) {
				throw new IllegalStateException("This method cannot be called on the test case");
			} else {
				this.clusterCenters = new double[clusterComponents];
				double[][] trimmedPcs = PrincipalComponentsResiduals.trimPcBasis(clusterComponents, pResiduals.getPcBasis(), log);
				for (int i = 0; i < clusterCenters.length; i++) {
					clusterCenters[i] = Array.median(Array.subArray(trimmedPcs[i], pcMatchedMask));
					log.reportTimeInfo(name + " Cluster " + i + " : " + clusterCenters[i]);
				}
			}
		}
	}

	public static void test(Project proj, String genoPCfile) {
		PrincipalComponentsResiduals pResiduals = new PrincipalComponentsResiduals(genoPCfile, Files.getHeaderOfFile(genoPCfile, proj.getLog()).length - 2, proj.getLog());
		VcfPopulation vpop = VcfPopulation.load(proj.SAMPLE_DATA_FILENAME.getValue(), POPULATION_TYPE.PC_ANCESTRY, proj.getLog());
		vpop.report();
		PCPopulation pcPopulation = new PCPopulation(proj, pResiduals, vpop, 2, true);
		TestSampleDistances[] testSampleDistances = pcPopulation.computeDistance();
		String output = ext.rootOf(genoPCfile, false) + ".ancestry.txt";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.println("DNA\tSTUDY\t" + Array.toStr(testSampleDistances[0].getPopulations()));
			for (int i = 0; i < testSampleDistances.length; i++) {
				String sample = testSampleDistances[i].getSample();
				writer.println(sample + "\t" + vpop.getPopulationForInd(sample, RETRIEVE_TYPE.SUB)[0] + "\t" + Array.toStr(testSampleDistances[i].getDistances()));
			}

			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + output);
			proj.getLog().reportException(e);
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
		String logfile = null;
		Logger log;

		String usage = "\n" + "cnv.analysis.pca.PCPopulation requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) full path to a genotype pc file (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pc=")) {
				genoPCfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			log = new Logger(logfile);
			test(new Project(filename, false), genoPCfile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
