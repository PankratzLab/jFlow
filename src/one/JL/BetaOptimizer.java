package one.JL;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;

import seq.manage.StrandOps;
import seq.manage.VCOps;
import seq.manage.StrandOps.CONFIG;
import common.Array;
import common.Files;
import common.Logger;
import common.Positions;
import common.ext;
import cnv.filesys.ABLookup;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.manage.Resources.ARRAY_RESOURCE_TYPE;
import cnv.manage.Resources.GENOME_BUILD;
import cnv.manage.Resources.GENOME_RESOURCE_TYPE;
import cnv.manage.Resources.Resource;
import filesys.Segment;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Going to handle the beta-optimization via correlating effects to a previous meta-analysis
 *
 */
public class BetaOptimizer {
	// 

	private static final String RS_HEADER = "rsID";
	private static final String BETA_HEADER = "beta";

	public static void run(MarkerSet markerSet, ABLookup abLookup, String dbsnpVCF, String[] namesToQuery, String outpuDir, String betaFile, Logger log) {
		new File(outpuDir).mkdirs();
		String outSer = outpuDir + "rsIdLookup.ser";
		ArrayList<MarkerRsFormat> markerRsFormats = null;
		if (Files.exists(outSer)) {
			try {
				log.reportTimeInfo("Trying to load " + outSer);
				markerRsFormats = MarkerRsFormat.readSerial(outSer, log);
				log.reportTimeInfo("Loaded " + outSer);
			} catch (Exception e) {

			}
		}
		if (markerRsFormats == null) {
			markerRsFormats = mapToRsIds(markerSet, abLookup, dbsnpVCF, namesToQuery, outSer, log);
		}
		ArrayList<MetaBeta> metaBetas = loadBetas(markerRsFormats, betaFile, log);

		// for(Me)
	}

	private static ArrayList<MetaBeta> loadBetas(ArrayList<MarkerRsFormat> markerRsFormats, String betaFile, Logger log) {
		String[] header = Files.getHeaderOfFile(betaFile, log);
		int[] indices = ext.indexFactors(new String[] { RS_HEADER, BETA_HEADER }, header, false, false);
		if (Array.countIf(indices, -1) > 0) {
			log.reportTimeError("Did not detect proper header in " + betaFile + ", requires " + RS_HEADER + " AND " + BETA_HEADER);
			return null;
		} else {
			ArrayList<MetaBeta> metaBetas = new ArrayList<BetaOptimizer.MetaBeta>();
			Hashtable<String, Integer> index = new Hashtable<String, Integer>();
			for (int i = 0; i < markerRsFormats.size(); i++) {
				index.put(markerRsFormats.get(i).getRs(), i);
			}
			try {
				BufferedReader reader = Files.getAppropriateReader(betaFile);
				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("\t");
					String rsId = line[indices[0]];
					if (index.containsKey(rsId)) {
						try {
							double beta = Double.parseDouble(line[indices[1]]);
							if (Double.isFinite(beta)) {
								MarkerRsFormat current = markerRsFormats.get(index.get(rsId));
								if (current.flipBetas()) {
									beta *= -1;
								}
								MetaBeta me = new MetaBeta(current, beta);
								metaBetas.add(me);
							} else {
								log.reportTimeWarning("Invalid beta on line " + Array.toStr(line));
							}
						} catch (NumberFormatException nfe) {
							log.reportTimeWarning("Invalid beta on line " + Array.toStr(line));
						}
					}

				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + betaFile + "\" not found in current directory");
				return null;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + betaFile + "\"");
				return null;
			}
			return metaBetas;

		}

	}

	public enum SITE_TYPE {
		BIALLELIC, TRIALLELIC, UNKNOWN;
	}

	/**
	 * Map a projects markers to rsIds by chr,pos,and alleles
	 * 
	 */
	public static ArrayList<MarkerRsFormat> mapToRsIds(MarkerSet markerSet, ABLookup abLookup, String dbsnpVCF, String[] namesToQuery, String outSer, Logger log) {

		String[] markerNames = markerSet.getMarkerNames();
		int[] indices = ext.indexLargeFactors(namesToQuery, markerNames, true, log, true, false);
		Segment[] segs = new Segment[indices.length];
		for (int i = 0; i < indices.length; i++) {
			segs[i] = new Segment(markerSet.getChrs()[indices[i]], markerSet.getPositions()[indices[i]], markerSet.getPositions()[indices[i]] + 1);
		}

		VCFFileReader reader = new VCFFileReader(dbsnpVCF, true);
		ArrayList<MarkerRsFormat> markerRsFormats = new ArrayList<MarkerRsFormat>();
		log.reportTimeInfo("Attempting to look up " + namesToQuery.length + " markers as rsIds from " + dbsnpVCF);
		for (int i = 0; i < segs.length; i++) {
			Segment current = segs[i];
			String[] allelesMarker = new String[] { abLookup.getLookup()[indices[i]][0] + "".toUpperCase(), abLookup.getLookup()[indices[i]][1] + "".toUpperCase() };
			if (i % 10000 == 0 && i > 0) {
				log.reportTimeInfo("queried " + (i + 1) + " markers");
			}
			MarkerRsFormat markerRsFormat = new MarkerRsFormat(namesToQuery[i], current.getStart(), indices[i], "NA", new String[] { "N", "N" }, allelesMarker, CONFIG.STRAND_CONFIG_UNKNOWN, SITE_TYPE.UNKNOWN);

			if (Array.countIf(allelesMarker, "N") == 0) {
				CloseableIterator<VariantContext> vcIter = reader.query(Positions.getChromosomeUCSC(current.getChr(), false), current.getStart() - 2, current.getStop() + 2);
				boolean foundGood = false;
				while (!foundGood && vcIter.hasNext()) {
					VariantContext vc = vcIter.next();

					if (!vc.isPointEvent() || current.getStart() == vc.getStart() || Integer.parseInt(VCOps.getAnnotationsFor(new String[] { "RSPOS" }, vc, "-1")[0]) == current.getStart()) {
						String[][] allelesVC = new String[][] { { "N", "N" } };
						if (vc.isPointEvent()) {
							allelesVC = new String[vc.getAlternateAlleles().size()][];// tri allele possibility
							for (int j = 0; j < allelesVC.length; j++) {
								allelesVC[j] = new String[] { vc.getReference().getBaseString(), vc.getAlternateAllele(j).getBaseString() };
							}
						} else if (vc.isIndel()) {
							if (vc.getReference().getBaseString().length() > vc.getAlternateAllele(0).length()) {
								allelesVC = new String[][] { { "I", "D" } };
							} else {
								allelesVC = new String[][] { { "D", "I" } };
							}
						}
						// else {//MNP, CLUMPED, SV...should'nt be on arrays
						// throw new IllegalArgumentException("Unknown variant type " + vc.toStringWithoutGenotypes());
						//
						// }
						for (int j = 0; j < allelesVC.length; j++) {
							if (!foundGood) {
								CONFIG config = StrandOps.determineStrandConfig(allelesMarker, allelesVC[j]);
								markerRsFormat = new MarkerRsFormat(namesToQuery[i], current.getStart(), indices[i], vc.getID(), allelesVC[j], allelesMarker, config, vc.isBiallelic() ? SITE_TYPE.BIALLELIC : SITE_TYPE.TRIALLELIC);
								switch (config) {

								case STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND:
								case STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND:
								case STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:
								case STRAND_CONFIG_SAME_ORDER_SAME_STRAND:
									markerRsFormat.setValidMatch(true);
									foundGood = true;
									break;
								case STRAND_CONFIG_BOTH_NULL:
								case STRAND_CONFIG_DIFFERENT_ALLELES:
								case STRAND_CONFIG_SPECIAL_CASE:
								default:
									break;
								}
							}
						}
					}

				}
			}
			markerRsFormats.add(markerRsFormat);
		}
		reader.close();
		// log.reportTimeInfo("Detected " + markerRsFormats.size() + " valid rs ids from " + dbsnpVCF);
		MarkerRsFormat.writeSerial(markerRsFormats, outSer);
		return markerRsFormats;
	}

	public static class MarkerRsFormat implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private int projectIndex;
		private String markerName;
		private String rs;
		private String[] dbSnpAlleles;
		private String[] markerAlleles;
		private int posMarker;
		private CONFIG config;
		private boolean validMatch;
		private SITE_TYPE type;

		private MarkerRsFormat(String markerName, int posMarker, int projectIndex, String rs, String[] refs, String[] markers, CONFIG config, SITE_TYPE type) {
			super();
			this.markerName = markerName;
			this.posMarker = posMarker;
			this.projectIndex = projectIndex;
			this.rs = rs;
			this.dbSnpAlleles = refs;
			this.markerAlleles = markers;
			this.config = config;
			this.validMatch = false;
			this.type = type;
		}

		int getPosMarker() {
			return posMarker;
		}

		int getProjectIndex() {
			return projectIndex;
		}

		SITE_TYPE getType() {
			return type;
		}

		public String getRs() {
			return rs;
		}

		CONFIG getConfig() {
			return config;
		}

		private void setValidMatch(boolean validMatch) {
			this.validMatch = validMatch;
		}

		public String getMarkerName() {
			return markerName;
		}

		public String[] getDbSnpAlleles() {
			return dbSnpAlleles;
		}

		public String[] getMarkerAlleles() {
			return markerAlleles;
		}

		 boolean flipBetas() {
			if (!validMatch) {
				throw new IllegalArgumentException("Did not have valid rs id match, this method should not be used");
			}
			return config == CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND || config == CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND;
		}

		public static void writeSerial(ArrayList<MarkerRsFormat> markerRsFormats, String fileName) {
			Files.writeSerial(markerRsFormats, fileName, true);
		}

		@SuppressWarnings("unchecked")
		public static ArrayList<MarkerRsFormat> readSerial(String fileName, Logger log) {
			return (ArrayList<MarkerRsFormat>) Files.readSerial(fileName, false, log, false, true);
		}

	}

	/**
	 * Storing the betas of variants from a meta-analysis
	 *
	 */
	private static class MetaBeta {
		private MarkerRsFormat markerRsFormat;
		private double beta;

		public MetaBeta(MarkerRsFormat markerRsFormat, double beta) {
			super();
			this.markerRsFormat = markerRsFormat;
			this.beta = beta;
		}

	}

	public static void main(String[] args) {
		Project proj = new Project("C:/workspace/Genvisis/projects/OSv2_hg19.properties", false);
		String betaFile = null;
		if (!Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
			ABLookup abLookup = new ABLookup();
			abLookup.parseFromGenotypeClusterCenters(proj);
			abLookup.writeToFile(proj.AB_LOOKUP_FILENAME.getValue(), proj.getLog());
		}
		MarkerSet markerSet = proj.getMarkerSet();
		ABLookup abLookup = new ABLookup(markerSet.getMarkerNames(), proj.AB_LOOKUP_FILENAME.getValue(), true, true, proj.getLog());

		Resource dbsnp = GENOME_RESOURCE_TYPE.DB_SNP147.getResource(GENOME_BUILD.HG19);
		Resource test = ARRAY_RESOURCE_TYPE.AFFY_SNP6_MARKER_POSITIONS.getResource(GENOME_BUILD.HG18);
		test.getResource(proj.getLog());
		run(markerSet, abLookup, dbsnp.getResource(proj.getLog()), proj.getAutosomalMarkers(), proj.PROJECT_DIRECTORY.getValue() + "betaOpt/", betaFile, proj.getLog());
	}
}
