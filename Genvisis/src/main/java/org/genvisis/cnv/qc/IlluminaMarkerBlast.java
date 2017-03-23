package org.genvisis.cnv.qc;

import java.io.FileNotFoundException;
import java.util.ArrayList;

import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.PROBE_TAG;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.TOP_BOT;
import org.genvisis.cnv.annotation.markers.BlastParams;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;

import htsjdk.tribble.annotation.Strand;

public class IlluminaMarkerBlast extends MarkerBlast {

	private final String manifestFile;

	public IlluminaMarkerBlast(Project proj, int numThreads, String manifestFile) {
		this(proj, getDefaultWordSize(proj), getDefaultWordSize(proj),
				 DEFAULT_MAX_ALIGNMENTS_REPORTED,
				 DEFAULT_REPORT_TO_TEMPORARY_FILE, DEFAULT_ANNOTATE_GC_CONTENT,
				 DEFAULT_DO_BLAST, numThreads, manifestFile);

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
	 * @param manifestFile
	 */
	public IlluminaMarkerBlast(Project proj, int blastWordSize, int reportWordSize,
														 int maxAlignmentsReported, boolean reportToTmp,
														 boolean annotateGCContent, boolean doBlast, int numThreads,
														 String manifestFile) {
		super(proj, blastWordSize, reportWordSize, maxAlignmentsReported, reportToTmp,
					annotateGCContent, doBlast, numThreads);
		this.manifestFile = manifestFile;
	}

	@Override
	protected String getSourceString() {
		return manifestFile;
	}

	@Override
	protected String getNameBase() {
		return ext.rootOf(manifestFile, true);
	}

	@Override
	protected void generateNaiveMarkerSet() {
		if (!Files.exists(proj.MARKERSET_FILENAME.getValue())) {
			proj.MARKERSET_FILENAME.setValue(proj.MARKERSET_FILENAME.getDefaultValue());
			String[] markerNames = null;
			if (!Files.exists(proj.MARKER_POSITION_FILENAME.getValue())) {
				proj.MARKER_POSITION_FILENAME.setValue(proj.MARKER_POSITION_FILENAME.getDefaultValue());
				markerNames = extractMarkerPositionsFromManifest(manifestFile, proj.getArrayType(),
																												 FILE_SEQUENCE_TYPE.MANIFEST_FILE,
																												 proj.MARKER_POSITION_FILENAME.getValue(),
																												 ",", proj.getLog());
			}
			Markers.orderMarkers(markerNames, proj.MARKER_POSITION_FILENAME.getValue(true, false),
													 proj.MARKERSET_FILENAME.getValue(true, false), proj.getLog());

		}
	}

	@Override
	protected MarkerFastaEntry[] getMarkerFastaEntries(BlastParams params,
																										 boolean alleleLookup) {
		ExtProjectDataParser.ProjectDataParserBuilder builder = formatParser(proj,
																																				 FILE_SEQUENCE_TYPE.MANIFEST_FILE,
																																				 manifestFile);
		MarkerFastaEntry[] fastaEntries = null;
		try {
			int seqLength = proj.ARRAY_TYPE.getValue().getProbeLength();

			ExtProjectDataParser parser = builder.build(proj, manifestFile);
			parser.determineIndicesFromTitles();
			parser.loadData();
			ArrayList<MarkerFastaEntry> entries = new ArrayList<MarkerFastaEntry>(ArrayUtils.booleanArraySum(parser.getDataPresent()));
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
			String refStrandTitle = getRefStrand(parser);
			if (refStrandTitle == null && params != null) {
				params.setNotes("Warning, the positive and negative strands of the probe design are actually forward and reverse designations, due to parsing IlmnIDs ");
			}
			for (int i = 0; i < parser.getDataPresent().length; i++) {
				if (parser.getDataPresent()[i]) {
					if (alleleLookup) {
						if ((i + 1) % 10000 == 0 && referenceGenome != null) {
							proj.getLog().reportTimeInfo("Loaded " + (i + 1) + " reference alleles from "
																					 + referenceGenome.getReferenceFasta());
						}
					}
					String seqA = null;
					String seqB = null;
					String markerName = parser.getDataToLoad()[i];
					Segment markerSegment = new Segment(markerSet.getChrs()[i], markerSet.getPositions()[i],
																							markerSet.getPositions()[i]);

					// builder.stringDataTitles(new String[] { "AlleleA_ProbeSeq", "AlleleB_ProbeSeq",
					// "SNP", "RefStrand", "IlmnStrand", "SourceStrand", "SourceSeq" });
					seqA = parser.getStringDataForTitle("AlleleA_ProbeSeq")[i];
					seqB = parser.getStringDataForTitle("AlleleB_ProbeSeq")[i];

					String refStrand = null;
					TOP_BOT topBotProbe = TOP_BOT.valueOf(parser.getStringDataForTitle("IlmnStrand")[i].toUpperCase());
					TOP_BOT topBotRef = TOP_BOT.valueOf(parser.getStringDataForTitle("SourceStrand")[i].toUpperCase());
					if (refStrandTitle != null) {
						refStrand = parser.getStringDataForTitle(refStrandTitle)[i];
					} else {
						refStrand = parseStrandFromIlmnID(parser.getStringDataForTitle("IlmnID")[i].split("_"),
																							topBotProbe);
					}
					Strand strand = null;
					if (refStrand.equals("+")) {
						strand = Strand.POSITIVE;
					} else if (refStrand.equals("-")) {
						strand = Strand.NEGATIVE;
					} else {
						proj.getLog().reportError("Invalid RefStrand " + refStrand);
						return null;
					}
					if (seqA.length() != seqLength) {
						proj.getLog()
								.reportError("Sequence " + seqA + " did not have length "
														 + proj.ARRAY_TYPE.getValue().getProbeLength() + " " + markerName);
						return null;
					}

					String[] snp = parser.getStringDataForTitle("SNP")[i].replaceAll("\\[", "")
																															 .replaceAll("\\]", "").split("/");
					AlleleParser alleleParser = new AlleleParser(markerName, markerSegment, snp[0], snp[1],
																											 seqB, referenceGenome);
					if (alleleLookup) {
						alleleParser.parse(proj.getArrayType(), strand,
															 parser.getStringDataForTitle("SourceSeq")[i]);
					}
					if (seqB.equals("")) {// not ambiguous
						entries.add(new MarkerFastaEntry(markerName + PROBE_TAG.BOTH.getTag(), seqA, seqB,
																						 strand, seqLength, markerSegment, topBotProbe,
																						 topBotRef, alleleParser.getA(), alleleParser.getB(),
																						 alleleParser.getRef(), alleleParser.getAlts()));

					} else {// interrogationPosition is the last bp for Ambiguous SNPS
						entries.add(new MarkerFastaEntry(markerName + PROBE_TAG.A.getTag(), seqA, seqB,
																						 strand, seqLength - 1, markerSegment, topBotProbe,
																						 topBotRef, alleleParser.getA(), alleleParser.getB(),
																						 alleleParser.getRef(), alleleParser.getAlts()));
						entries.add(new MarkerFastaEntry(markerName + PROBE_TAG.B.getTag(), seqB, seqB,
																						 strand, seqLength - 1, markerSegment, topBotProbe,
																						 topBotRef, alleleParser.getA(), alleleParser.getB(),
																						 alleleParser.getRef(), alleleParser.getAlts()));
						if (seqB.length() != seqLength) {
							proj.getLog()
									.reportError("Sequence " + seqB + " did not have length "
															 + proj.ARRAY_TYPE.getValue().getProbeLength() + " " + markerName);
							return null;
						}
					}


				}
			}
			fastaEntries = entries.toArray(new MarkerFastaEntry[entries.size()]);
		} catch (FileNotFoundException e) {
			proj.getLog().reportFileNotFound(manifestFile);
			e.printStackTrace();
		}
		proj.getLog().reportTimeInfo("Found " + fastaEntries.length + " marker sequences");
		return fastaEntries;
	}

}
