package cnv.annotation;

import java.io.File;

import seq.manage.ReferenceGenome;
import seq.manage.VCOps;
import common.Array;
import common.Files;
import common.ext;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.WritingFilePrep;

/**
 * @author lane0212 Basically we hijack a vcf file for storing our annotations...seems like most of the things we want to annotate will have a chr, start and stop
 */
public abstract class AnnotationFileWriter extends AnnotationFile implements WritingFilePrep {

	private VariantContextWriter writer;
	private boolean overWriteExisting;
	private boolean additionMode;
	private String tmpFile;
	private CloseableIterator<VariantContext> additionReader;

	public AnnotationFileWriter(Project proj, Annotation[] annotations, String annotationFilename, boolean overWriteExisting) {
		super(proj, annotationFilename);
		setAnnotations(annotations);
		this.overWriteExisting = overWriteExisting;
		this.additionMode = false;
		this.tmpFile = null;
		this.additionReader = null;
		init();
	}

	@Override
	public boolean overwriteExisting() {
		return overWriteExisting;
	}

	/**
	 * @param locusAnnotation
	 *            we convert this to a minimal {@link VariantContext} and write to the underlying {@link VariantContextWriter}
	 */
	public void write(LocusAnnotation locusAnnotation) {
		if (writer != null) {
			VariantContext vcAnno = locusAnnotation.annotateLocus();
			if (additionReader != null) {
				if (!additionReader.hasNext()) {
					String error = "Mismatched number of entries in " + annotationFilename + " , cancelling addition...";
					proj.getLog().reportTimeError(error);
					throw new IllegalStateException(error);
				} else {
					VariantContext vcAdd = additionReader.next();
					if (vcAdd.hasSameAllelesAs(vcAnno) && VCOps.getSegment(vcAdd).equals(VCOps.getSegment(vcAnno)) && vcAdd.getID().equals(vcAnno.getID())) {
						VariantContextBuilder builder = new VariantContextBuilder(vcAdd);
						for (String att : vcAnno.getAttributes().keySet()) {
							builder.attribute(att, vcAnno.getAttributes().get(att));
						}
						vcAnno = builder.make();
					} else {
						String error = "Entries must be in the exact same order to be added together, ";
						error += vcAdd.toStringWithoutGenotypes() + " from " + annotationFilename;
						error += vcAnno.toStringWithoutGenotypes() + " being added";
						proj.getLog().reportTimeError(error);
						throw new IllegalStateException(error);
					}

				}
			}

			writer.add(vcAnno);
		} else {
			String error = "annotation writer has not been intialized";
			proj.getLog().reportTimeError(error);
			throw new IllegalStateException(error);
		}
	}

	@Override
	public boolean validate() {
		if (!Files.exists(annotationFilename) || overwriteExisting()) {
			additionMode = false;
		} else {
			additionMode = true;
			proj.getLog().reportTimeInfo("Attempting to initialize addition mode");
			String ts = ext.getTimestampForFilename();
			Files.copyFileUsingFileChannels(new File(annotationFilename), new File(annotationFilename + "." + ts + ".bak"), proj.getLog());
			Files.copyFileUsingFileChannels(new File(annotationFilename + ".tbi"), new File(annotationFilename + ".tbi." + ts + ".bak"), proj.getLog());
		}

		boolean valid = true;
		if (additionMode) {
			try {
				this.additionReader = new VCFFileReader(annotationFilename, false).iterator();
				this.tmpFile = getTmpFile(annotationFilename);
			} catch (Exception e) {
				proj.getLog().reportTimeError("Trying to initialize addition mode, but " + annotationFilename + " did not pass vcf file checks");
				proj.getLog().reportException(e);
				valid = false;
			}
		}

		if (valid) {
			valid = Files.exists(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue());
			if (valid) {
				valid = annotations != null;
				if (!valid) {
					proj.getLog().reportTimeError("Must provided annotation array");
				}
			} else {
				proj.getLog().reportTimeError("Could not find required file " + proj.REFERENCE_GENOME_FASTA_FILENAME.getValue());
			}
		}

		return valid;
	}

	@Override
	public void init() {
		if (validate()) {
			VCFHeader vcfHeader = additionMode ? new VCFFileReader(annotationFilename, true).getFileHeader() : new VCFHeader();// take previous header if needed
			for (int i = 0; i < annotations.length; i++) {
				if (!vcfHeader.hasInfoLine(annotations[i].getName())) {
					VCFInfoHeaderLine vHeaderLine = null;
					if (annotations[i].getCount() != null) {
						vHeaderLine = new VCFInfoHeaderLine(annotations[i].getName(), annotations[i].getCount(), annotations[i].getType(), annotations[i].getDescription());
					} else {
						vHeaderLine = new VCFInfoHeaderLine(annotations[i].getName(), annotations[i].getNumber(), annotations[i].getType(), annotations[i].getDescription());
					}

					vcfHeader.addMetaDataLine(vHeaderLine);
				} else {
					proj.getLog().reportTimeWarning("Detected that info line " + annotations[i].getName() + " is already present, any new data added will overwrite previous");
				}
			}

			VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(additionMode ? tmpFile : annotationFilename);// tmp file if needed
			builder.clearOptions();
			builder.setOption(Options.INDEX_ON_THE_FLY);
			builder.setOption(Options.DO_NOT_WRITE_GENOTYPES);

			String refGenome = proj.REFERENCE_GENOME_FASTA_FILENAME.getValue();
			proj.getLog().reportTimeInfo("Using reference genome" + refGenome);

			SAMSequenceDictionary samSequenceDictionary = new ReferenceGenome(refGenome, proj.getLog()).getIndexedFastaSequenceFile().getSequenceDictionary();
			MarkerSet markerSet = proj.getMarkerSet();

			if (samSequenceDictionary.getSequenceIndex("chrXY") == -1 && markerSet.getIndicesByChr()[25].length > 0) {// not always present in ref
				int[] xyLen = Array.subArray(markerSet.getPositions(), markerSet.getIndicesByChr()[25]);
				proj.getLog().reportTimeInfo("Since the project contained markers designated as pseudo-autosomal, a chrXY contig is being added");
				samSequenceDictionary.addSequence(new SAMSequenceRecord("chrXY", Array.max(xyLen)));
			}
			builder.setReferenceDictionary(samSequenceDictionary);
			vcfHeader.setSequenceDictionary(samSequenceDictionary);
			this.writer = builder.build();
			writer.writeHeader(vcfHeader);

		} else {
			proj.getLog().reportTimeError("Could not intialize annotation file " + annotationFilename);
		}
	}

	private static String getTmpFile(String annotationFilename) {
		return annotationFilename + ext.getTimestampForFilename() + "tmp.gz";
	}

	@Override
	public void close() {
		if (writer != null) {
			writer.close();
		}
		if (additionReader != null) {
			additionReader.close();
			proj.getLog().reportTimeInfo("Copying temporary file " + tmpFile + " to " + annotationFilename + ", " + tmpFile + " can be deleted on successful completion");
			Files.copyFileUsingFileChannels(new File(tmpFile), new File(annotationFilename), proj.getLog());
			Files.copyFileUsingFileChannels(new File(tmpFile + ".tbi"), new File(annotationFilename + ".tbi"), proj.getLog());
			new File(tmpFile).delete();
			new File(tmpFile + ".tbi").delete();
		}
	}
}
