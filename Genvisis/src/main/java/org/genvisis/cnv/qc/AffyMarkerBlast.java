package org.genvisis.cnv.qc;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.LinkedHashSet;

import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.PROBE_TAG;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.TOP_BOT;
import org.genvisis.cnv.annotation.markers.BlastParams;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.ExtProjectDataParser.ProjectDataParserBuilder;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;

import htsjdk.tribble.annotation.Strand;

public class AffyMarkerBlast extends MarkerBlast {

	private static final String PROBE_MARKER_NAME = "PROBESET_ID";
	private static final String PROBE_SEQUENCE = "PROBE_SEQUENCE";
	private static final String PROBE_STRAND = "TARGET_STRANDEDNESS";

	private static final String ANNOT_MARKER_NAME = "Probe Set ID";
	private static final String ANNOT_ALLELE_A = "Allele A";
	private static final String ANNOT_ALLELE_B = "Allele B";

	private final String probeFile;
	private final String annotFile;

	private AffyMarkerBlast(Project proj, int numThreads, String probeFile, String annotFile) {
		this(proj, getDefaultWordSize(proj), getDefaultWordSize(proj),
				 DEFAULT_MAX_ALIGNMENTS_REPORTED,
				 DEFAULT_REPORT_TO_TEMPORARY_FILE, DEFAULT_ANNOTATE_GC_CONTENT,
				 DEFAULT_DO_BLAST, numThreads, probeFile, annotFile);
	}

	/**
	 * @param proj
	 * @param blastWordSize
	 * @param reportWordSize
	 * @param maxAlignmentsReported
	 * @param reportToTmp
	 * @param annotateGCContent
	 * @param doBlast
	 * @param numThreads
	 * @param probeFile
	 * @param annotFile
	 */
	private AffyMarkerBlast(Project proj, int blastWordSize, int reportWordSize,
													int maxAlignmentsReported, boolean reportToTmp, boolean annotateGCContent,
													boolean doBlast, int numThreads, String probeFile, String annotFile) {
		super(proj, blastWordSize, reportWordSize, maxAlignmentsReported, reportToTmp,
					annotateGCContent, doBlast, numThreads);
		if (proj.getArrayType() != ARRAY.AFFY_GW6 && proj.getArrayType() != ARRAY.AFFY_GW6_CN) {
			proj.getLog()
					.reportError("Array type was set to " + proj.getArrayType() + " and this file is for "
											 + ARRAY.AFFY_GW6 + " or " + ARRAY.AFFY_GW6_CN);
			throw new IllegalArgumentException("Cannot generate " + AffyMarkerBlast.class.getName()
																				 + " for non-Affy array");
		}
		this.probeFile = probeFile;
		this.annotFile = annotFile;
	}

	@Override
	protected String getSourceString() {
		return probeFile + "," + annotFile;
	}

	@Override
	protected void generateNaiveMarkerSet() {
		if (!Files.exists(proj.MARKERSET_FILENAME.getValue())) {
			proj.MARKERSET_FILENAME.setValue(proj.MARKERSET_FILENAME.getDefaultValue());
			String[] markerNames = null;
			if (!Files.exists(proj.MARKER_POSITION_FILENAME.getValue())) {
				proj.getLog()
						.reportError("Marker Names not found at " + proj.MARKER_POSITION_FILENAME.getValue()
												 + ", must be provided for Affy arrays");
			}
			Markers.orderMarkers(markerNames, proj.MARKER_POSITION_FILENAME.getValue(true, false),
													 proj.MARKERSET_FILENAME.getValue(true, false), proj.getLog());

		}
	}

	@Override
	protected String getNameBase() {
		return ext.rootOf(probeFile) + "_" + ext.rootOf(annotFile);
	}

	@Override
	protected MarkerFastaEntry[] getMarkerFastaEntries(BlastParams params,
																										 boolean alleleLookup) {
		int seqLength = proj.ARRAY_TYPE.getValue().getProbeLength();
		ExtProjectDataParser probeFileParser;
		try {
			probeFileParser = probeFileParser().build(proj, probeFile);
		} catch (FileNotFoundException e) {
			proj.getLog().reportFileNotFound(probeFile);
			e.printStackTrace();
			return null;
		}
		ExtProjectDataParser annotFileParser;
		try {
			annotFileParser = annotFileParser().build(proj, annotFile);
		} catch (FileNotFoundException e) {
			proj.getLog().reportFileNotFound(annotFile);
			e.printStackTrace();
			return null;
		}

		probeFileParser.determineIndicesFromTitles();
		probeFileParser.loadData();
		annotFileParser.determineIndicesFromTitles();
		annotFileParser.loadData();
		ArrayList<MarkerFastaEntry> entries = new ArrayList<MarkerFastaEntry>(ArrayUtils.booleanArraySum(probeFileParser.getDataPresent()));
		MarkerSetInfo markerSet = proj.getMarkerSet();
		// SequenceLookup sequenceLookup = new SequenceLookup(proj.getLog());
		ReferenceGenome referenceGenome = Files.exists(proj.getReferenceGenomeFASTAFilename()) ? new ReferenceGenome(proj.getReferenceGenomeFASTAFilename(),
																																																								 proj.getLog())
																																													 : null;
		// ABLookup abLookup =
		// new ABLookup( markerSet.getMarkerNames(), proj.AB_LOOKUP_FILENAME.getValue(),
		// true, true, proj.getLog());
		//
		// Hashtable<String, Integer> indices = proj.getMarkerIndices();
		if (referenceGenome == null) {
			proj.getLog().reportTimeWarning("A reference genome was not found");
		}
		String refStrandTitle = getRefStrand(probeFileParser);
		if (refStrandTitle == null && params != null) {
			params.setNotes("Warning, the positive and negative strands of the probe design are actually forward and reverse designations, due to parsing IlmnIDs ");
		}
		for (int i = 0; i < probeFileParser.getDataPresent().length; i++) {
			if (probeFileParser.getDataPresent()[i]) {
				if (alleleLookup) {
					if ((i + 1) % 10000 == 0 && referenceGenome != null) {
						proj.getLog().reportTimeInfo("Loaded " + (i + 1) + " reference alleles from "
																				 + referenceGenome.getReferenceFasta());
					}
				}
				String markerName = probeFileParser.getDataToLoad()[i];
				Segment markerSegment = new Segment(markerSet.getChrs()[i], markerSet.getPositions()[i],
																						markerSet.getPositions()[i]);
				String[] allSeqs = ArrayUtils.unique(probeFileParser.getStringDataForTitle(PROBE_SEQUENCE)[i].split("\t"));
				LinkedHashSet<String> collapsed = new LinkedHashSet<String>();
				for (String seq : allSeqs) {
					collapsed.add(seq);
				}
				String[] tmpSeq = ArrayUtils.toStringArray(collapsed);

				if (tmpSeq.length != 2) {
					proj.getLog()
							.reportError("Marker " + markerName + " did not have 2 unique probe designs");
					proj.getLog().reportError("found the following " + markerName + "\t"
																		+ ArrayUtils.toStr(tmpSeq));
					return null;
				} else {
					if (tmpSeq[0].length() != seqLength || tmpSeq[1].length() != seqLength) {
						proj.getLog()
								.reportError("Sequence " + tmpSeq[0] + " or " + tmpSeq[1]
														 + "  did not have length "
														 + proj.ARRAY_TYPE.getValue().getProbeLength());
						return null;
					}
					int interrogationPosition = -1;
					for (int j = 0; j < tmpSeq[0].length(); j++) {
						if (tmpSeq[0].charAt(j) != tmpSeq[1].charAt(j)) {
							if (interrogationPosition != -1) {
								proj.getLog().reportError("Multiple interrogation position for " + markerName);
								return null;
							}
							interrogationPosition = j;
						}
					}
					String[] affyStrandtmp = probeFileParser.getStringDataForTitle(PROBE_STRAND)[i].split("\t");
					if (ArrayUtils.unique(affyStrandtmp).length != 1) {
						proj.getLog()
								.reportError("Multiple strands detected " + ArrayUtils.toStr(affyStrandtmp));
						return null;
					}
					String affyStrand = affyStrandtmp[0];
					Strand strand = null;
					if (affyStrand.equals("f")) {// PLUS seems to be for cnvi probes
						strand = Strand.POSITIVE;
					} else if (affyStrand.equals("r")) {// MINUS seems to be for cnvi probes
						strand = Strand.NEGATIVE;
					} else {
						proj.getLog().reportError("Invalid AffyStrand " + affyStrand);
						return null;
					}
					// probeFileParser and annotFileParser should both be indexed by marker indices
					String aS = annotFileParser.getStringDataForTitle(ANNOT_ALLELE_A)[i];
					String bS = annotFileParser.getStringDataForTitle(ANNOT_ALLELE_B)[i];
					// if (!Arrays.equals( new String[] {StrandOps.flipIfNeeded(aS, strand, false),
					// StrandOps.flipIfNeeded(bS, strand, false)},
					// new String[] {abLookup.getLookup()[indices.get(markerName)][0]
					// + "",
					// abLookup.getLookup()[indices.get(markerName)][1]
					// + ""})) {
					// System.err.println(markerName+ "\t" + ArrayUtils.toStr(tmpSeq) + "\t"
					// + ArrayUtils.toStr(new String[] {abLookup.getLookup()[indices.get(markerName)][0]
					// + "",
					// abLookup.getLookup()[indices.get(markerName)][1]
					// + "",
					// aS, bS, strand.toString()}));
					//
					// }
					AlleleParser alleleParser = new AlleleParser(markerName, markerSegment, aS, bS,
																											 tmpSeq[1], referenceGenome);
					if (alleleLookup) {
						alleleParser.parse(proj.getArrayType(), strand, null);
					}
					entries.add(new MarkerFastaEntry(markerName + PROBE_TAG.A.getTag(), tmpSeq[0],
																					 tmpSeq[1],
																					 strand, interrogationPosition, markerSegment,
																					 TOP_BOT.NA, TOP_BOT.NA, alleleParser.getA(),
																					 alleleParser.getB(), alleleParser.getRef(),
																					 alleleParser.getAlts()));

					entries.add(new MarkerFastaEntry(markerName + PROBE_TAG.B.getTag(), tmpSeq[1],
																					 tmpSeq[1],
																					 strand, interrogationPosition, markerSegment,
																					 TOP_BOT.NA, TOP_BOT.NA, alleleParser.getA(),
																					 alleleParser.getB(), alleleParser.getRef(),
																					 alleleParser.getAlts()));
				}
			}
		}

		proj.getLog().reportTimeInfo("Found " + entries.size() + " marker sequences");
		return entries.toArray(new MarkerFastaEntry[entries.size()]);
	}


	protected ProjectDataParserBuilder probeFileParser() {
		ExtProjectDataParser.ProjectDataParserBuilder builder = new ExtProjectDataParser.ProjectDataParserBuilder();
		builder.separator("\t");
		builder.dataKeyColumnName(PROBE_MARKER_NAME);
		builder.stringDataTitles(PROBE_SEQUENCE, PROBE_STRAND);
		builder.concatMultipleStringEntries(true);
		builder.sampleBased(false);
		builder.treatAllNumeric(false);
		builder.requireAll(false);
		builder.verbose(false);
		return builder;
	}


	protected ProjectDataParserBuilder annotFileParser() {
		ExtProjectDataParser.ProjectDataParserBuilder builder = new ExtProjectDataParser.ProjectDataParserBuilder();
		builder.separator(",");
		builder.dataKeyColumnName(ANNOT_MARKER_NAME);
		builder.headerFlags(ANNOT_MARKER_NAME, ANNOT_ALLELE_A, ANNOT_ALLELE_B);
		builder.stringDataTitles(ANNOT_ALLELE_A, ANNOT_ALLELE_B);
		builder.sampleBased(false);
		builder.treatAllNumeric(false);
		builder.requireAll(false);
		builder.verbose(false);
		return builder;
	}

}
