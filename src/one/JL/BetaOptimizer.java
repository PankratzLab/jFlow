package one.JL;

import java.util.ArrayList;

import seq.manage.StrandOps;
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
 * Going to handle the beta-optimization via correlating effects to a previous meta-analyis
 *
 */
public class BetaOptimizer {

	/**
	 * Map a projects markers to rsIds by chr,pos,and alleles
	 * 
	 */
	private static void mapToRsIds(MarkerSet markerSet, ABLookup abLookup, String dbsnpVCF, String[] namesToQuery, Logger log) {
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
			if (Array.countIf(allelesMarker, "N") == 0) {
				CloseableIterator<VariantContext> vcIter = reader.query(Positions.getChromosomeUCSC(current.getChr(), false), current.getStart(), current.getStop());

				while (vcIter.hasNext()) {
					VariantContext vc = vcIter.next();
					if (vc.isPointEvent()) {
						if (current.getStart() == vc.getStart()) {
							String[] allelesVC = new String[] { vc.getReference().getBaseString(), vc.getAlternateAllele(0).getBaseString() };
							CONFIG config = StrandOps.determineStrandConfig(allelesMarker, allelesVC);
							// System.out.println( vc.getID() + "\t" + i + "\t" + Array.toStr(allelesMarker) + "\t" + Array.toStr(allelesVC) + "\t" + config);
							switch (config) {

							case STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND:
							case STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND:
							case STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:
							case STRAND_CONFIG_SAME_ORDER_SAME_STRAND:
								// System.out.println(vc.getAttribute("RV")+ "\t" + Array.toStr(allelesMarker) + "\t" + Array.toStr(allelesVC) + "\t" + config);
								MarkerRsFormat markerRsFormat = new MarkerRsFormat(namesToQuery[i], indices[i], vc.getID(), allelesVC, allelesMarker, config);
								markerRsFormats.add(markerRsFormat);
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
		reader.close();
		log.reportTimeInfo("Detected " + markerRsFormats.size() + " valid rs ids from " + dbsnpVCF);
	}

	private static class MarkerRsFormat {
		private int projectIndex;
		private String markerName;
		private String rs;
		private String[] dbSnpAlleles;
		private String[] markerAlleles;
		private CONFIG config;

		public MarkerRsFormat(String markerName, int projectIndex, String rs, String[] refs, String[] markers, CONFIG config) {
			super();
			this.markerName = markerName;
			this.projectIndex = projectIndex;
			this.rs = rs;
			this.dbSnpAlleles = refs;
			this.markerAlleles = markers;
			this.config = config;
		}

		public boolean flipGenotypes() {
			return config == CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND || config == CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND;
		}
	}

	/**
	 * Storing the betas of variants from a meta-analysis
	 *
	 */
	private static class MetaBetas {

	}

	public static void main(String[] args) {
		Project proj = new Project("C:/workspace/Genvisis/projects/OSv2_hg19.properties", false);
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
		System.out.println(dbsnp.isAvailable(proj.getLog()));
		mapToRsIds(markerSet, abLookup, dbsnp.getResource(proj.getLog()), proj.getAutosomalMarkers(), proj.getLog());

	}

}
