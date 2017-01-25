package org.genvisis.cnv.qc;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.LinkedHashSet;
import java.util.List;

import org.genvisis.cnv.annotation.markers.AnalysisParams;
import org.genvisis.cnv.annotation.markers.AnnotationData;
import org.genvisis.cnv.annotation.markers.AnnotationFileWriter;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.PROBE_TAG;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.TOP_BOT;
import org.genvisis.cnv.annotation.markers.BlastAnnotationWriter;
import org.genvisis.cnv.annotation.markers.BlastParams;
import org.genvisis.cnv.annotation.markers.LocusAnnotation;
import org.genvisis.cnv.annotation.markers.LocusAnnotation.Builder;
import org.genvisis.cnv.annotation.markers.MarkerGCAnnotation;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.analysis.Blast;
import org.genvisis.seq.analysis.Blast.BlastWorker;
import org.genvisis.seq.analysis.Blast.FastaEntry;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.StrandOps;

import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;

/**
 * As starts handling more than just annotatating with Blast results, we should maybe re-name this
 * class
 *
 */
public class MarkerBlast {

	public static final int DEFAULT_MAX_ALIGNMENTS_REPORTED = 20;
	public static final boolean DEFAULT_REPORT_TO_TEMPORARY_FILE = true;
	public static final boolean DEFAULT_ANNOTATE_GC_CONTENT = true;
	public static final boolean DEFAULT_DO_BLAST = true;

	public enum FILE_SEQUENCE_TYPE {
																	// /**
																	// * like HumanOmni2-5-8-v1-2-A-Strand-Report-FDT.txt
																	// */
																	// STRAND_REPORT,
																	/**
																	 * Like HumanExome-12-v1-0-B.csv
																	 */
																	MANIFEST_FILE,
																	/**
																	 * Like GenomeWideSNP_6.probe_tab
																	 */
																	AFFY_ANNOT;
	}

	public static MarkerBlastResult blastEm(Project proj, String fileSeq, FILE_SEQUENCE_TYPE type,
																					int numThreads) {
		int wordSize = getDefaultWordSize(proj);
		return blastEm(	proj, fileSeq, type, wordSize, wordSize, DEFAULT_MAX_ALIGNMENTS_REPORTED,
										numThreads, DEFAULT_REPORT_TO_TEMPORARY_FILE, DEFAULT_ANNOTATE_GC_CONTENT,
										DEFAULT_DO_BLAST);
	}

	/**
	 * {@link MarkerBlast#blastEm(Project, String, FILE_SEQUENCE_TYPE, int, int, int, int, boolean, boolean, boolean)}
	 */
	public static MarkerBlastResult blastEm(Project proj, String fileSeq, FILE_SEQUENCE_TYPE type,
																					int blastWordSize, int reportWordSize,
																					int maxAlignmentsReported, int numThreads,
																					boolean reportToTmp, boolean annotateGCContent,
																					boolean doBlast) {
		String fastaDb = proj.getReferenceGenomeFASTAFilename();
		if (!Files.exists(fastaDb) && doBlast) {
			proj.getLog().reportError("Unable to find or obtain reference genome");
			return null;
		} else {
			double evalueCutoff = 10000;
			BlastParams blastParams = new BlastParams(fileSeq, fastaDb, maxAlignmentsReported,
																								reportWordSize, blastWordSize,
																								ext.getTimestampForFilename(), evalueCutoff,
																								proj.getMarkerSet().getFingerprint(), "",
																								proj.getLog());
			Blast blast = new Blast(fastaDb, blastWordSize, reportWordSize, proj.getLog(), true, true);
			blast.setEvalue(evalueCutoff);// we rely on the wordSize instead
			String dir = proj.PROJECT_DIRECTORY.getValue() + "Blasts/";
			new File(dir).mkdirs();

			String root = dir	+ ext.rootOf(fileSeq, true) + ".blasted.ws." + blastWordSize + ".rep."
										+ reportWordSize;
			String output = root + ".blasted";

			String[] tmps = new String[numThreads];
			for (int i = 0; i < numThreads; i++) {
				tmps[i] = root + ".tmp" + i;
				if (!Files.exists(tmps[i])) {
					tmps[i] = tmps[i] + ".gz";
				}
				// proj.getLog().reportTimeInfo("REmember gz");
			}

			if (doBlast && !Files.exists("", tmps)) {
				MarkerFastaEntry[] fastaEntries = getMarkerFastaEntries(proj, fileSeq, type, null, false);
				List<MarkerFastaEntry[]> splits =
																				Array.splitUpArray(fastaEntries, numThreads, proj.getLog());

				ArrayList<BlastWorker> workers = new ArrayList<Blast.BlastWorker>();
				if (fastaEntries != null && fastaEntries.length > 0) {
					for (int i = 0; i < splits.size(); i++) {
						if (!Files.exists(tmps[i])) {
							workers.add(new BlastWorker(blast, splits.get(i), reportToTmp ? tmps[i] : null));
						} else {
							proj.getLog().reportTimeWarning("Skipping index " + i + ", " + tmps[i] + " exists");
						}
					}
				}

				if (workers.size() > 0) {
					WorkerHive<Blast.BlastResultsSummary[]> hive = new WorkerHive<Blast.BlastResultsSummary[]>(	numThreads,
																																																			10,
																																																			proj.getLog());
					hive.addCallables(workers.toArray(new BlastWorker[workers.size()]));
					hive.setReportEvery(1);
					hive.execute(true);
					hive.getResults();
				}
			} else {
				for (String tmp : tmps) {
					if (doBlast) {
						proj.getLog().reportFileExists(tmp);
					}
				}
			}
			proj.getLog()
					.reportTimeWarning("Assuming that all sequences have length "
																+ proj.getArrayType().getProbeLength() + " based on array type "
															+ proj.getArrayType());
			MarkerBlastResult result = new MarkerBlastResult(	output, blastWordSize, reportWordSize,
																												proj.getArrayType().getProbeLength());
			result.setTmpFiles(tmps);

			// TODO, revert this to addition mode again
			if (Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue())) {
				new File(proj.BLAST_ANNOTATION_FILENAME.getValue()).delete();
			}
			MarkerFastaEntry[] entries = getMarkerFastaEntries(proj, fileSeq, type, blastParams, true);

			proj.getLog().reportTimeInfo("Summarizing blast results to "
																		+ proj.BLAST_ANNOTATION_FILENAME.getValue());
			BlastAnnotationWriter blastAnnotationWriter =
																									new BlastAnnotationWriter(proj,
																																						new AnalysisParams[] {blastParams},
																																						proj.BLAST_ANNOTATION_FILENAME.getValue(),
																																						doBlast	? tmps
																																										: new String[] {},
																																						reportWordSize,
																																						proj.getArrayType()
																																								.getProbeLength(),
																																						proj.getArrayType()
																																								.getProbeLength(),
																																						maxAlignmentsReported);

			if (proj.getArrayType() != ARRAY.ILLUMINA) {
				proj.getLog()
						.reportTimeWarning("Did not detect array type "	+ ARRAY.ILLUMINA
																+ " , probe sequence annotation may not reflect the true design since multiple designs may be reported");
			}
			blastAnnotationWriter.setMarkerFastaEntries(entries);
			blastAnnotationWriter.summarizeResultFiles(false);
			blastAnnotationWriter.close();
			if (annotateGCContent) {
				annotateGCContent(proj, fileSeq, type);
			}

			// } else {
			// proj.getLog().reportFileExists(proj.BLAST_ANNOTATION_FILENAME.getValue());
			// }
			return result;
		}
	}

	public static int getDefaultWordSize(Project proj) {
		int reportWordSize;
		switch (proj.getArrayType()) {
			case AFFY_GW6:
				reportWordSize = 15;
				break;
			case AFFY_GW6_CN:
				reportWordSize = 15;
				break;
			case ILLUMINA:
				reportWordSize = 25;
				break;
			default:
				reportWordSize = -1;
				proj.getLog().reportError("Invalid array type " + proj.getArrayType());
				break;
		}
		return reportWordSize;
	}

	/**
	 * Adds the gc content annotation to the summary file, does some repeat loading but oh well.
	 */
	private static void annotateGCContent(Project proj, String fileSeq, FILE_SEQUENCE_TYPE type) {
		MarkerFastaEntry[] fastaEntries = getMarkerFastaEntries(proj, fileSeq, type, null, false);
		String[] markerNames = proj.getMarkerNames();
		MarkerSet markerSet = proj.getMarkerSet();
		Hashtable<String, Integer> indices = proj.getMarkerIndices();
		// ReferenceGenome referenceGenome = new ReferenceGenome(fastaDb, proj.getLog());
		LocusAnnotation[] gcAnnotations = new LocusAnnotation[proj.getMarkerNames().length];
		for (MarkerFastaEntry fastaEntrie : fastaEntries) {
			String marker = fastaEntrie.getName();
			PROBE_TAG tag;
			switch (proj.getArrayType()) {
				// we have to remove the "_A" or "_B" from affy markers
				case AFFY_GW6:
					tag = PROBE_TAG.parseMarkerTag(marker, proj.getLog());
					marker = marker.substring(0, marker.length() - tag.getTag().length());
					break;
				case AFFY_GW6_CN:
					tag = PROBE_TAG.parseMarkerTag(marker, proj.getLog());
					marker = marker.substring(0, marker.length() - tag.getTag().length());
					break;
				case ILLUMINA:
					tag = PROBE_TAG.parseMarkerTag(marker, proj.getLog());
					marker = marker.substring(0, marker.length() - tag.getTag().length());
					break;
				default:
					proj.getLog().reportError("Invalid Array type " + proj.getArrayType());
					break;
			}
			int index = indices.get(marker);
			double gcContent = fastaEntrie.getGCMinusInterrogationPosition();
			Builder builder = new Builder();
			AnnotationData annotationData = MarkerGCAnnotation.getGCAnnotationDatas();
			annotationData.setData(gcContent + "");
			builder.annotations(new AnnotationData[] {annotationData});
			MarkerGCAnnotation markerGCAnnotation =
																						new MarkerGCAnnotation(	builder, marker,
																																		fastaEntrie.getMarkerSegment());
			gcAnnotations[index] = markerGCAnnotation;
		}
		AnnotationFileWriter writer = new AnnotationFileWriter(	proj, null,
																														new AnnotationData[] {MarkerGCAnnotation.getGCAnnotationDatas()},
																														proj.BLAST_ANNOTATION_FILENAME.getValue(),
																														false) {};
		for (int i = 0; i < gcAnnotations.length; i++) {
			if (gcAnnotations[i] != null) {// currently some markers may not be represented, such as affy
																			// CN markers
				writer.write(gcAnnotations[i], false, true);

			} else {
				double gcContent = Double.NaN;
				Builder builder = new Builder();
				AnnotationData annotationData = MarkerGCAnnotation.getGCAnnotationDatas();
				annotationData.setData(gcContent + "");
				builder.annotations(new AnnotationData[] {annotationData});
				MarkerGCAnnotation blank = new MarkerGCAnnotation(builder, markerNames[i],
																													new Segment(markerSet.getChrs()[i],
																																			markerSet.getPositions()[i],
																																			markerSet.getPositions()[i]));
				writer.write(blank, false, true);
			}
		}
		writer.close();
	}

	public static class MarkerBlastResult {
		private final String output;
		private String[] tmpFiles;
		private final int blastWordSize;
		private final int reportWordSize;
		private final int sequenceSize;

		public MarkerBlastResult(	String output, int blastWordSize, int reportWordSize,
															int sequenceSize) {
			super();
			this.output = output;
			this.blastWordSize = blastWordSize;
			this.reportWordSize = reportWordSize;
			this.sequenceSize = sequenceSize;
		}

		public int getSequenceSize() {
			return sequenceSize;
		}

		public void setTmpFiles(String[] tmpFiles) {
			this.tmpFiles = tmpFiles;
		}

		public String[] getTmpFiles() {
			return tmpFiles;
		}

		public int getBlastWordSize() {
			return blastWordSize;
		}

		public int getReportWordSize() {
			return reportWordSize;
		}

		public String getOutput() {
			return output;
		}
	}

	private static String getRefStrand(ExtProjectDataParser parser) {
		String strand = null;
		if (parser.hasStringDataForTitle("RefStrand")) {
			strand = "RefStrand";
		} else if (parser.hasStringDataForTitle("GenomicStrand")) {
			strand = "GenomicStrand";
		}
		return strand;
	}

	/**
	 *
	 * {@link MarkerBlast#getMarkerFastaEntries(Project, String, FILE_SEQUENCE_TYPE, BlastParams, boolean)}
	 */
	private static MarkerFastaEntry[] getMarkerFastaEntries(Project proj, String strandReportFile,
																													FILE_SEQUENCE_TYPE type,
																													BlastParams params,
																													boolean alleleLookup) {
		ExtProjectDataParser.ProjectDataParserBuilder builder = formatParser(	proj, type,
																																					strandReportFile);
		MarkerFastaEntry[] fastaEntries = null;
		try {
			int seqLength = proj.ARRAY_TYPE.getValue().getProbeLength();

			ExtProjectDataParser parser = builder.build(proj, strandReportFile);
			parser.determineIndicesFromTitles();
			parser.loadData();
			ArrayList<MarkerFastaEntry> entries =
																					new ArrayList<MarkerFastaEntry>(Array.booleanArraySum(parser.getDataPresent()));
			MarkerSet markerSet = proj.getMarkerSet();
			// SequenceLookup sequenceLookup = new SequenceLookup(proj.getLog());
			ReferenceGenome referenceGenome =
																			Files.exists(proj.getReferenceGenomeFASTAFilename())	? new ReferenceGenome(proj.getReferenceGenomeFASTAFilename(),
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
							proj.getLog().reportTimeInfo("Loaded "	+ (i + 1) + " reference alleles from "
																						+ referenceGenome.getReferenceFasta());
						}
					}
					String seqA = null;
					String seqB = null;
					String markerName = parser.getDataToLoad()[i];
					Segment markerSegment = new Segment(markerSet.getChrs()[i], markerSet.getPositions()[i],
																							markerSet.getPositions()[i]);
					if (type != FILE_SEQUENCE_TYPE.AFFY_ANNOT) {

						// builder.stringDataTitles(new String[] { "AlleleA_ProbeSeq", "AlleleB_ProbeSeq",
						// "SNP", "RefStrand", "IlmnStrand", "SourceStrand", "SourceSeq" });
						seqA = parser.getStringDataForTitle("AlleleA_ProbeSeq")[i];
						seqB = parser.getStringDataForTitle("AlleleB_ProbeSeq")[i];

						String refStrand = null;
						TOP_BOT topBotProbe =
																TOP_BOT.valueOf(parser.getStringDataForTitle("IlmnStrand")[i].toUpperCase());
						TOP_BOT topBotRef =
															TOP_BOT.valueOf(parser.getStringDataForTitle("SourceStrand")[i].toUpperCase());
						if (refStrandTitle != null) {
							refStrand = parser.getStringDataForTitle(refStrandTitle)[i];
						} else {
							refStrand =
												parseStrandFromIlmnID(parser.getStringDataForTitle("IlmnID")[i].split("_"),
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
									.reportError("Sequence "	+ seqA + " did not have length "
																+ proj.ARRAY_TYPE.getValue().getProbeLength() + " " + markerName);
							return null;
						}

						String[] snp = parser.getStringDataForTitle("SNP")[i]	.replaceAll("\\[", "")
																																	.replaceAll("\\]", "").split("/");
						AlleleParser alleleParser = new AlleleParser(	markerName, markerSegment, snp[0], snp[1],
																													seqB, referenceGenome);
						if (alleleLookup) {
							alleleParser.parse(	proj.getArrayType(), strand,
																	parser.getStringDataForTitle("SourceSeq")[i]);
						}
						if (seqB.equals("")) {// not ambiguous
							entries.add(new MarkerFastaEntry(markerName	+ PROBE_TAG.BOTH.getTag(), seqA, seqB,
																								strand, seqLength, markerSegment, topBotProbe,
																								topBotRef, alleleParser.getA(), alleleParser.getB(),
																								alleleParser.getRef(), alleleParser.getAlts()));

						} else {// interrogationPosition is the last bp for Ambiguous SNPS
							entries.add(new MarkerFastaEntry(markerName	+ PROBE_TAG.A.getTag(), seqA, seqB,
																								strand, seqLength - 1, markerSegment, topBotProbe,
																								topBotRef, alleleParser.getA(), alleleParser.getB(),
																								alleleParser.getRef(), alleleParser.getAlts()));
							entries.add(new MarkerFastaEntry(markerName	+ PROBE_TAG.B.getTag(), seqB, seqB,
																								strand, seqLength - 1, markerSegment, topBotProbe,
																								topBotRef, alleleParser.getA(), alleleParser.getB(),
																								alleleParser.getRef(), alleleParser.getAlts()));
							if (seqB.length() != seqLength) {
								proj.getLog()
										.reportError("Sequence "	+ seqB + " did not have length "
																	+ proj.ARRAY_TYPE.getValue().getProbeLength() + " " + markerName);
								return null;
							}
						}

					} else {

						String[] allSeqs = Array.unique(parser.getStringDataForTitle("PROBE_SEQUENCE")[i].split("\t"));
						LinkedHashSet<String> collapsed = new LinkedHashSet<String>();
						for (String seq : allSeqs) {
							collapsed.add(seq);
						}
						String[] tmpSeq = Array.toStringArray(collapsed);

						if (tmpSeq.length != 2) {
							proj.getLog()
									.reportError("Marker " + markerName + " did not have 2 unique probe designs");
							proj.getLog()
									.reportError("found the following " + markerName + "\t" + Array.toStr(tmpSeq));
							return null;
						} else {
							if (tmpSeq[0].length() != seqLength || tmpSeq[1].length() != seqLength) {
								proj.getLog()
										.reportError("Sequence "	+ tmpSeq[0] + " or " + tmpSeq[1]
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
							String[] affyStrandtmp = parser.getStringDataForTitle("TARGET_STRANDEDNESS")[i].split("\t");
							if (Array.unique(affyStrandtmp).length != 1) {
								proj.getLog()
										.reportError("Multiple strands detected " + Array.toStr(affyStrandtmp));
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
							String aS = tmpSeq[0].substring(interrogationPosition, interrogationPosition + 1);
							String bS = tmpSeq[1].substring(interrogationPosition, interrogationPosition + 1);
							// if (!Arrays.equals( new String[] {StrandOps.flipIfNeeded(aS, strand, false),
							// StrandOps.flipIfNeeded(bS, strand, false)},
							// new String[] {abLookup.getLookup()[indices.get(markerName)][0]
							// + "",
							// abLookup.getLookup()[indices.get(markerName)][1]
							// + ""})) {
							// System.err.println(markerName+ "\t" + Array.toStr(tmpSeq) + "\t"
							// + Array.toStr(new String[] {abLookup.getLookup()[indices.get(markerName)][0]
							// + "",
							// abLookup.getLookup()[indices.get(markerName)][1]
							// + "",
							// aS, bS, strand.toString()}));
							//
							// }
							AlleleParser alleleParser = new AlleleParser(	markerName, markerSegment, aS, bS,
																														tmpSeq[1], referenceGenome);
							if (alleleLookup) {
								alleleParser.parse(proj.getArrayType(), strand, null);
							}
							entries.add(new MarkerFastaEntry(markerName	+ PROBE_TAG.A.getTag(), tmpSeq[0],
																								tmpSeq[1], strand, interrogationPosition,
																								markerSegment, TOP_BOT.NA, TOP_BOT.NA,
																								alleleParser.getA(), alleleParser.getB(),
																								alleleParser.getRef(), alleleParser.getAlts()));
							entries.add(new MarkerFastaEntry(markerName	+ PROBE_TAG.B.getTag(), tmpSeq[1],
																								tmpSeq[1], strand, interrogationPosition,
																								markerSegment, TOP_BOT.NA, TOP_BOT.NA,
																								alleleParser.getA(), alleleParser.getB(),
																								alleleParser.getRef(), alleleParser.getAlts()));
							throw new IllegalArgumentException("Sorry, Affymetrix AB designation is not currently working, feel free to remove this exception if you just want blast results, but make sure to not rely on ABs - also shout at @jlanej");
						}
					}
				}
			}
			fastaEntries = entries.toArray(new MarkerFastaEntry[entries.size()]);
		} catch (FileNotFoundException e) {
			proj.getLog().reportFileNotFound(strandReportFile);
			e.printStackTrace();
		}
		proj.getLog().reportTimeInfo("Found " + fastaEntries.length + " marker sequences");
		return fastaEntries;
	}

	private static String parseStrandFromIlmnID(String[] fullId, TOP_BOT topBotRef) {
		String refStrand = null;
		String tb = fullId[fullId.length - 3];
		String fr = fullId[fullId.length - 2];
		if ((!tb.equals("B") && !tb.equals("T") && !tb.equals("M"))
				|| (!fr.equals("F") && !fr.equals("R"))) {
			throw new IllegalStateException("Invalid IlmnID parsing asumption, splitting with _ for "
																				+ Array.toStr(fullId, "_") + "\t" + tb + "\t" + fr
																			+ "\nPlease update TOP/BOT strand designation");
		}

		switch (topBotRef) {
			case BOT:
			case MINUS:
				if (((tb.equals("B") || tb.equals("M")) && fr.equals("F"))
						|| tb.equals("T") && fr.equals("R")) {
					refStrand = "+";
				} else if (((tb.equals("B") || tb.equals("M")) && fr.equals("R"))
										|| tb.equals("T") && fr.equals("F")) {
					refStrand = "-";
				}

				break;
			case PLUS:
			case TOP:
				if ((tb.equals("B") && fr.equals("R")) || tb.equals("T") && fr.equals("F")) {
					refStrand = "+";
				} else if ((tb.equals("B") && fr.equals("F")) || tb.equals("T") && fr.equals("R")) {
					refStrand = "-";
				}
				break;
			case NA:
				break;
			default:
				break;

		}
		if (refStrand == null) {
			throw new IllegalStateException("Could not parse ref TOP/BOT  "	+ topBotRef + " using IlmnID"
																			+ Array.toStr(fullId, "_") + " to a strand");
		}
		return refStrand;
	}

	/**
	 * @author lane0212 Parse alleles to ref/alt/A/B with strand flip for alts if needed
	 */
	static class AlleleParser {
		private final Segment loc;
		private final String aS;
		private final String bS;
		private final String pseqB;
		private final String markerName;
		private final ReferenceGenome referenceGenome;
		private Allele A;
		private Allele B;
		private Allele ref;
		private Allele[] alts;

		private AlleleParser(	String markerName, Segment loc, String aS, String bS, String pseqB,
													ReferenceGenome referenceGenome) {
			super();
			this.markerName = markerName;
			this.loc = loc;
			this.aS = aS;
			this.pseqB = pseqB;
			this.bS = bS;
			this.referenceGenome = referenceGenome;
		}

		public Allele getA() {
			return A;
		}

		public Allele getB() {
			return B;
		}

		public Allele getRef() {
			return ref;
		}

		public Allele[] getAlts() {
			return alts;
		}

		private boolean isIndel() {
			if (aS.equals("I") || bS.equals("I")) {
				return true;
			}
			if (aS.equals("D") || bS.equals("D")) {
				return true;
			}
			return false;
		}

		private void parseIndel(ARRAY array, Strand strand, String sourceSeq) {
			if (sourceSeq == null) {
				throw new IllegalArgumentException("Must have source seq to get indel");
			} else if (!pseqB.equals("")) {
				throw new IllegalArgumentException("Assumptions violated for single probe indel ");

			} else if (!(aS.equals("I") && bS.equals("D")) && !(aS.equals("D") && bS.equals("I"))) {
				throw new IllegalArgumentException("Assumptions violated for insertion deletion codes "	+ aS
																						+ " and " + bS);

			} else {

				int start = sourceSeq.indexOf("[") + 1;
				int stop = sourceSeq.indexOf("]");

				String indel = sourceSeq.substring(start, stop);
				if (!indel.startsWith("-/") && !indel.startsWith("+/")) {
					throw new IllegalArgumentException("Indel source seq must start with - or + for "
																							+ indel);
				}
				indel = indel.substring(2);
				int len = indel.length() - 1;
				if (loc.getChr() > 0) {

					Segment newQuery = new Segment(loc.getChr(), loc.getStart() - 1, loc.getStop() + len);// need
																																																// to
																																																// grab
																																																// preceeding
																																																// ref
																																																// bp
					String[] tmp = referenceGenome == null	? Array.stringArray(indel.length() + 1, "N")
																									: referenceGenome.getSequenceFor(newQuery);
					if (tmp.length - 1 != indel.length()) {
						throw new IllegalStateException("Invalid reference indel query; REF -> "
																						+ Array.toStr(tmp) + "\t indel -> " + indel);
					}
					boolean insertionisRef = true;
					for (int i = 0; i < indel.length(); i++) {
						if (!tmp[i + 1].equalsIgnoreCase(indel.charAt(i) + "")
								&& !tmp[i + 1].equalsIgnoreCase(StrandOps.flipIfNeeded(indel.charAt(i)	+ "", strand,
																																				false))) {
							insertionisRef = false;
							break;
						}
					}
					ref = Allele.create(Array.toStr(tmp, ""), true);
					if (insertionisRef) {
						if (aS.equals("I")) {

							String tmpA = Array.toStr(tmp, "");
							String tmpB = tmp[0];
							A = Allele.create(StrandOps.flipsIfNeeded(tmpA, strand, false), false);
							B = Allele.create(StrandOps.flipsIfNeeded(tmpB, strand, false), false);
							alts = new Allele[1];
							alts[0] = Allele.create(tmp[0], false);

						} else {
							String tmpA = tmp[0];
							String tmpB = Array.toStr(tmp, "");
							A = Allele.create(StrandOps.flipsIfNeeded(tmpA, strand, false), false);
							B = Allele.create(StrandOps.flipsIfNeeded(tmpB, strand, false), false);
							alts = new Allele[1];
							alts[0] = Allele.create(tmp[0], false);
						}
					} else {
						if (aS.equals("I")) {
							String tmpA = StrandOps.flipIfNeeded(tmp[0], strand, false) + indel;
							String tmpB = StrandOps.flipIfNeeded(tmp[0], strand, false);
							A = Allele.create(tmpA, false);
							B = Allele.create(tmpB, false);
							alts = new Allele[2];
							alts[0] = Allele.create(tmp[0], false);
							alts[1] =
											Allele.create(tmp[0] + StrandOps.flipsIfNeeded(indel, strand, false), false);
						} else {
							String tmpA = StrandOps.flipIfNeeded(tmp[0], strand, false);
							String tmpB = StrandOps.flipIfNeeded(tmp[0], strand, false) + indel;
							A = Allele.create(tmpA, false);
							B = Allele.create(tmpB, false);
							alts = new Allele[2];
							alts[0] = Allele.create(tmp[0], false);
							alts[1] =
											Allele.create(tmp[0] + StrandOps.flipsIfNeeded(indel, strand, false), false);
						}
					}
				} else {

					ref = Allele.create("NN", true);
					if (aS.equals("I")) {
						String tmpA = "N" + indel;
						String tmpB = "N";
						A = Allele.create(tmpA, false);
						B = Allele.create(tmpB, false);
						alts = new Allele[2];
						alts[0] = Allele.create("N" + StrandOps.flipsIfNeeded(indel, strand, false), false);
						alts[1] = Allele.create("N", false);
					} else {
						String tmpA = "N";
						String tmpB = "N" + indel;
						A = Allele.create(tmpA, false);
						B = Allele.create(tmpB, false);
						alts = new Allele[2];
						alts[0] = Allele.create("N" + StrandOps.flipsIfNeeded(indel, strand, false), false);
						alts[1] = Allele.create("N", false);
					}
				}
			}
		}

		private void parse(ARRAY array, Strand strand, String sourceSeq) {
			if (loc.getChr() > 0) {
				if (!isIndel()) {
					String[] tmp;
					if (loc.getStart() > referenceGenome.getContigLength(loc) || referenceGenome == null) {
						tmp = new String[] {"N"};
					} else {
						tmp = referenceGenome.getSequenceFor(loc);

					}
					if (tmp.length != 1) {// don't think we need multiple
						throw new IllegalArgumentException("base query must be length one ("
																									+ loc.getUCSClocation() + " returned "
																								+ Array.toStr(tmp));
					} else if (array.isCNOnly(markerName)) {

						ref = Allele.create(tmp[0], true);
						alts = new Allele[1];
						alts[0] = Allele.create("<" + array + "_CNP>", false);// Symbolic allele
						A = Allele.create(alts[0], false);
						B = Allele.create(alts[0], false);

					} else {
						String flipA = StrandOps.flipIfNeeded(aS, strand, false);
						String flipB = StrandOps.flipIfNeeded(bS, strand, false);

						ref = Allele.create(tmp[0], true);
						A = Allele.create(aS, false);
						B = Allele.create(bS, false);
						if (Allele.create(flipA, false).equals(ref, true)) {// A is ref
							alts = new Allele[1];
							alts[0] = Allele.create(flipB);
						} else if (Allele.create(flipB, false).equals(ref, true)) {// B is ref
							alts = new Allele[1];
							alts[0] = Allele.create(flipA);
						} else {// neither are ref
							alts = new Allele[2];
							alts[0] = Allele.create(flipA, false);
							alts[1] = Allele.create(flipB, false);
						}
					}
				} else {
					parseIndel(array, strand, sourceSeq);
				}

			} else {
				if (array.isCNOnly(markerName)) {

					ref = Allele.create("N", true);
					alts = new Allele[1];
					alts[0] = Allele.create("<" + array + "_CN>", false);// Symbolic allele
					A = Allele.create(alts[0], false);
					B = Allele.create(alts[0], false);

				} else if (isIndel()) {
					parseIndel(array, strand, sourceSeq);
				} else {
					ref = Allele.create("N", true);
					alts = new Allele[2];
					alts[0] = Allele.create(StrandOps.flipIfNeeded(aS, strand, false), false);
					alts[1] = Allele.create(StrandOps.flipIfNeeded(bS, strand, false), false);
					A = Allele.create(aS, false);
					B = Allele.create(bS, false);
				}
			}
		}
	}

	private static ExtProjectDataParser.ProjectDataParserBuilder formatParser(Project proj,
																																						FILE_SEQUENCE_TYPE type,
																																						String strandReportFile) {
		ExtProjectDataParser.ProjectDataParserBuilder builder =
																													new ExtProjectDataParser.ProjectDataParserBuilder();
		switch (type) {
			case MANIFEST_FILE:
				if (proj.getArrayType() != ARRAY.ILLUMINA) {
					proj.getLog().reportError("Array type was set to "	+ proj.getArrayType()
																		+ " and this file is for " + ARRAY.ILLUMINA);
					builder = null;
					break;
				}
				builder.separator(",");
				builder.dataKeyColumnName("Name");
				String[] headerFlags = new String[] {"Name", "AlleleA_ProbeSeq"};
				builder.headerFlags(headerFlags);
				String ref = null;
				String[] header =
												Files.getLineContaining(strandReportFile, ",", headerFlags, proj.getLog());
				if (header != null) {
					if (ext.indexOfStr("RefStrand", header) >= 0) {
						ref = "RefStrand";
					} else if (ext.indexOfStr("GenomicStrand", header) >= 0) {
						ref = "GenomicStrand";
					}

				} else {
					proj.getLog().reportError("Header of " + strandReportFile + " not found");
					return null;
				}
				String[] dataTitles = new String[] {"AlleleA_ProbeSeq", "AlleleB_ProbeSeq", "SNP",
																						"IlmnStrand", "SourceStrand", "SourceSeq", "IlmnID"};
				if (ref != null) {
					dataTitles = Array.concatAll(dataTitles, new String[] {ref});
				} else {
					proj.getLog()
							.reportTimeWarning(strandReportFile
																	+ " did not have a column for RefStrand, will determine strand from IlmnID but these will be inacurate, forward assigned to +");
					// proj.getLog().reportTimeError("Actully we are stopping, get the strands from
					// http://www.well.ox.ac.uk/~wrayner/strand/");
					// return null;
				}
				builder.stringDataTitles(dataTitles);

				break;

			case AFFY_ANNOT:
				if (proj.getArrayType() != ARRAY.AFFY_GW6 && proj.getArrayType() != ARRAY.AFFY_GW6_CN) {
					proj.getLog()
							.reportError("Array type was set to "	+ proj.getArrayType() + " and this file is for "
														+ ARRAY.AFFY_GW6 + " or " + ARRAY.AFFY_GW6_CN);
					builder = null;
					break;
				}
				builder.separator("\t");
				builder.dataKeyColumnName("PROBESET_ID");
				builder.stringDataTitles(new String[] {"PROBE_SEQUENCE", "TARGET_STRANDEDNESS"});
				builder.concatMultipleStringEntries(true);

			default:
				break;

		}
		builder.sampleBased(false);
		builder.treatAllNumeric(false);
		builder.requireAll(false);
		builder.verbose(false);
		return builder;
	}

	public static class MarkerFastaEntry extends FastaEntry {

		private final int interrogationPosition;// Illumina = sequence length, affy = somewhere in the
																						// middle
		private final Segment markerSegment;
		private final String seqB;
		private final Strand strand;
		private final TOP_BOT topBotProbe;
		private final TOP_BOT topBotRef;
		private final Allele A;
		private final Allele B;
		private final Allele ref;
		private final Allele[] alts;

		// IlmnStrand
		public MarkerFastaEntry(String name, String seqA, String seqB, Strand strand,
														int interrogationPosition, Segment markerSegment, TOP_BOT topBotProbe,
														TOP_BOT topBotRef, Allele A, Allele B, Allele ref, Allele[] alts) {
			super(name, seqA);
			this.seqB = seqB;
			this.strand = strand;
			this.interrogationPosition = interrogationPosition;
			this.markerSegment = markerSegment;
			this.topBotProbe = topBotProbe;
			this.topBotRef = topBotRef;
			this.A = A;
			this.B = B;
			this.ref = ref;
			this.alts = alts;
		}

		public String getSeqB() {
			return seqB;
		}

		public Allele getRef() {
			return ref;
		}

		public Allele[] getAlts() {
			return alts;
		}

		public Strand getStrand() {
			return strand;
		}

		public Allele getA() {
			return A;
		}

		public Allele getB() {
			return B;
		}

		public Segment getMarkerSegment() {

			return markerSegment;
		}

		public int getInterrogationPosition() {
			return interrogationPosition;
		}

		public TOP_BOT getTopBotProbe() {
			return topBotProbe;
		}

		public TOP_BOT getTopBotRef() {
			return topBotRef;
		}

		public double getGCMinusInterrogationPosition() {
			int GorC = 0;
			int total = 0;
			for (int i = 0; i < sequence.length(); i++) {
				if (i != interrogationPosition) {
					total++;
					if ((sequence.charAt(i) + "").equalsIgnoreCase("G")
							|| (sequence.charAt(i) + "").equalsIgnoreCase("C")) {
						GorC++;
					}
				}
			}
			return (double) GorC / total;
		}

	}

	private static Project prepareDummyProject(String csv, FILE_SEQUENCE_TYPE type) {
		if (Files.exists(csv)) {
			String dir = ext.parseDirectoryOfFile(csv);
			Project proj = new Project();
			proj.PROJECT_DIRECTORY.setValue(dir);
			Logger log = proj.getLog();
			switch (type) {
				case AFFY_ANNOT:
					System.err.println("Not implemented for AFFY, yet");
					return null;
				// proj.ARRAY_TYPE.setValue(ARRAY.AFFY_GW6);
				case MANIFEST_FILE:
					proj.ARRAY_TYPE.setValue(ARRAY.ILLUMINA);
					log.reportTimeWarning("Extracting marker positions from " + csv);
					String[] markerNames = Markers.extractMarkerPositionsFromManifest(csv,
																																						proj.getArrayType(),
																																						proj.GENOME_BUILD_VERSION.getValue(),
																																						type,
																																						proj.MARKER_POSITION_FILENAME.getValue(),
																																						log);
					Markers.orderMarkers(	markerNames, proj.MARKER_POSITION_FILENAME.getValue(),
																proj.MARKERSET_FILENAME.getValue(true, true), proj.getLog());

					break;
				default:
					break;
			}
			return proj;
		} else {
			System.err.println("could not find " + csv);
			return null;
		}
	}



	// public static void test() {
	// String file = "/home/pankrat2/shared/aric_exome_chip/Blast/HumanExome-12-v1-0-B.csv";
	// String fastaDb = "/home/pankrat2/public/bin/ref/hg19_canonical.fa";
	// blastEm(new Project("/home/pankrat2/lanej/projects/aric_exome.properties", false), file,
	// FILE_SEQUENCE_TYPE.MANIFEST_FILE, fastaDb, 2);
	// }

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		// String fastaDb = "hg19_canonical.fa";
		String fileSeq = "HumanExome-12-v1-0-B.csv";
		int numThreads = 4;
		int blastWordSize = -1;
		int reportWordSize = -1;
		int maxAlignmentsReported = DEFAULT_MAX_ALIGNMENTS_REPORTED;
		boolean annotateGCContent = DEFAULT_ANNOTATE_GC_CONTENT;
		boolean doBlast = DEFAULT_DO_BLAST;
		// boolean report = true;
		FILE_SEQUENCE_TYPE fSequence_TYPE = FILE_SEQUENCE_TYPE.MANIFEST_FILE;
		String usage = "\n" + "cnv.qc.MarkerBlast requires 3 arguments\n";
		usage += "   (1) Project file name (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) full path to an Illumina manifest file  (i.e. fileSeq="	+ fileSeq
							+ " (default))\n" + "";
		usage += "   (3) number of threads to use  (i.e. "	+ PSF.Ext.NUM_THREADS_COMMAND + numThreads
							+ " (default))\n" + "";
		usage += "   (4) word size for initial match in blast db  (i.e. blastWordSize="	+ blastWordSize
							+ " (defaults are array dependent))\n" + "";
		usage += "   (5) number of base pairs with an exact match to report (i.e. reportWordSize="
							+ reportWordSize + " (defaults are array dependent))\n" + "";
		// usage += " (6) report results to temporary files (i.e. -report (not the default))\n" + "";
		usage += "   (6) the maximum number of alignments to summarize (i.e. maxAlignmentsReported="
							+ maxAlignmentsReported + " (default))\n" + "";
		usage +=
					"   (7) annotate the summary file with GC content, both the probe's gc content and the regions gc content  (i.e. annoGC="
							+ annotateGCContent + " (default))\n";
		usage += "   (8) sequence file type  (i.e. seqType="	+ fSequence_TYPE + " (default))\n"
							+ ", Options are:\n ";
		usage +=
					"   (9) skip the blast analysis and simply add the manifest annotations (i.e. -noBlast (not the default))\n"
							+ ", Options are:\n ";
		for (int i = 0; i < FILE_SEQUENCE_TYPE.values().length; i++) {
			usage += FILE_SEQUENCE_TYPE.values()[i] + "\n";
		}

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("fileSeq=")) {
				fileSeq = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("blastWordSize=")) {
				blastWordSize = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("reportWordSize=")) {
				reportWordSize = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("seqType=")) {
				fSequence_TYPE = FILE_SEQUENCE_TYPE.valueOf(ext.parseStringArg(	arg,
																																				FILE_SEQUENCE_TYPE.MANIFEST_FILE.toString()));
				numArgs--;
			} else if (arg.startsWith("annoGC=")) {
				annotateGCContent = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("maxAlignmentsReported=")) {
				maxAlignmentsReported = ext.parseIntArg(arg);
				numArgs--;

			} else if (arg.startsWith("-noBlast")) {
				doBlast = false;
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
			Project proj = null;

			if (filename == null) {
				proj = prepareDummyProject(fileSeq, fSequence_TYPE);
			} else {
				proj = new Project(filename, false);
			}
			if (proj == null) {
				System.err.println("Invalid project");
			}
			if (reportWordSize == -1) {
				switch (proj.getArrayType()) {
					case AFFY_GW6:
						reportWordSize = 15;
						break;
					case AFFY_GW6_CN:
						reportWordSize = 15;
						break;
					case ILLUMINA:
						reportWordSize = 25;
						break;
					default:
						proj.getLog().reportError("Invalid array type " + proj.getArrayType());
						break;
				}
				proj.getLog().reportTimeInfo("report word size updated to default "	+ reportWordSize
																			+ " for array " + proj.getArrayType());
			}
			if (blastWordSize == -1) {
				blastWordSize = getDefaultWordSize(proj);
				if (blastWordSize != -1) {
					proj.getLog().reportTimeInfo("blast word size updated to default "	+ blastWordSize
																				+ " for array " + proj.getArrayType());
				} else {
					return;
				}

			}
			blastEm(proj, fileSeq, fSequence_TYPE, blastWordSize, reportWordSize, maxAlignmentsReported,
							numThreads, DEFAULT_REPORT_TO_TEMPORARY_FILE, annotateGCContent, doBlast);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}

// String alts = null;
// for (int j = 0; j < tmpSeq[0].length(); j++) {
// if (tmpSeq[0].charAt(j) != tmpSeq[1].charAt(j)) {
// if (alts != null) {
// proj.getLog().reportTimeError("found multiple probe differences for " + parser.getDataToLoad()[i]
// + "\t" + Array.toStr(tmpSeq));
// } else {
// alts = sequenceLookup.getNucForCombo("" + tmpSeq[0].charAt(j) + tmpSeq[1].charAt(j));
// seq = tmpSeq[0].substring(0, j) + alts + tmpSeq[0].substring(j + 1, tmpSeq[0].length());
// System.out.println(alts + "\t" + parser.getDataToLoad()[i] + "\t" + tmpSeq[0]);
// System.out.println(alts + "\t" + parser.getDataToLoad()[i] + "\t" + tmpSeq[1]);
// System.out.println(alts + "\t" + parser.getDataToLoad()[i] + "\t" + seq + "\t" +
// tmpSeq[0].length() + "\t" + seq.length());
// // }
// // }
// // }
// private static boolean isPerfectlyGoodMatch(Blast.BlastResultsSummary bSummary, byte chr) {
// if (bSummary.getNumPerfectMatches() != 1) {
// return false;
// }
//
// return true;
//
// }

//
// private static void summarizeOneHitWonders(Project proj, ArrayList<Blast.BlastResultsSummary[]>
// bSummaries, String outputFile, Logger log) {
// try {
// PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
// Hashtable<String, Integer> indices = proj.getMarkerIndices();
// byte[] chrs = proj.getMarkerSet().getChrs();
// int[] pos = proj.getMarkerSet().getPositions();
// for (int i = 0; i < bSummaries.size(); i++) {
// for (int j = 0; j < bSummaries.get(i).length; j++) {
// if (bSummaries.get(i)[j].getNumPerfectMatches() == 1) {
// if (Array.countIf(bSummaries.get(i)[j].getPercentIdentityHistogram().getCounts(), 0) ==
// bSummaries.get(i)[j].getPercentIdentityHistogram().getCounts().length - 1) {
// int currentMarkerIndex = indices.get(bSummaries.get(i)[j].getName());
// Segment markerSegment = new Segment(chrs[currentMarkerIndex], pos[currentMarkerIndex] - 1,
// pos[currentMarkerIndex] + 1);
// Segment blastPerfectSegment = bSummaries.get(i)[j].getPerfectMatchSegment();
// if (blastPerfectSegment != null && blastPerfectSegment.overlaps(markerSegment)) {
// writer.println(bSummaries.get(i)[j].getName() + "\t" + chrs[currentMarkerIndex] + "\t" +
// pos[currentMarkerIndex] + "\t" + blastPerfectSegment.getChr() + "\t" +
// blastPerfectSegment.getStart() + "\t" + blastPerfectSegment.getStop());
// }
// }
// }
// }
// }
//
// writer.close();
// } catch (Exception e) {
// log.reportError("Error writing to " + outputFile);
// log.reportException(e);
// }
// }
//
// private static void summarize(Project proj, ArrayList<Blast.BlastResultsSummary[]> bSummaries,
// String outputFile, Logger log) {
// try {
// PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
// Hashtable<String, Integer> indices = proj.getMarkerIndices();
// byte[] chrs = proj.getMarkerSet().getChrs();
// writer.print("Marker\tDesignContig\tNumTotalContigs");
// double[] bins = bSummaries.get(0)[0].getPercentIdentityHistogram().getBins();
// HashSet<String> allContigs = new HashSet<String>();
// ArrayList<String> uniqContigs = new ArrayList<String>();
// for (int i = 0; i < bSummaries.size(); i++) {
// for (int j = 0; j < bSummaries.get(i).length; j++) {
// allContigs.addAll(bSummaries.get(i)[j].getHitCounts().keySet());
// }
// }
// uniqContigs.addAll(allContigs);
// for (int i = 0; i < bins.length; i++) {
// writer.print("\t" + ext.formDeci(bins[i], 0) + "_PercentIdentity");
// }
// for (String contig : uniqContigs) {
// writer.print("\t" + contig);
// }
// writer.println();
// for (int i = 0; i < bSummaries.size(); i++) {
// for (int j = 0; j < bSummaries.get(i).length; j++) {
// int[] counts = bSummaries.get(i)[j].getPercentIdentityHistogram().getCounts();
// int contigHits = bSummaries.get(i)[j].getHitCounts().size();
// String chr = Positions.getChromosomeUCSC(chrs[indices.get(bSummaries.get(i)[j].getName())],
// true);
// writer.print(bSummaries.get(i)[j].getName() + "\t" + chr + "\t" + contigHits + "\t" +
// Array.toStr(counts));
// for (String contig : uniqContigs) {
// if (bSummaries.get(i)[j].getHitCounts().containsKey(contig)) {
// writer.print("\t" + bSummaries.get(i)[j].getHitCounts().get(contig));
// } else {
// writer.print("\t" + "0");
// }
// }
// writer.println();
// }
//
// }
//
// writer.close();
// } catch (Exception e) {
// log.reportError("Error writing to " + outputFile);
// log.reportException(e);
// }
// }

// case STRAND_REPORT:
// if (proj.getArrayType() != ARRAY.ILLUMINA) {
// proj.getLog().reportTimeError("Array type was set to " + proj.getArrayType() + " and this file is
// for " + ARRAY.ILLUMINA);
// return null;
//
// }
// builder.dataKeyColumnName("SNP_Name");
// builder.commentString("##");
// builder.stringDataTitles(new String[] { "Design_Seq" });
// break;
