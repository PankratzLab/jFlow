package cnv.annotation;

import seq.manage.ReferenceGenome;
import common.Array;
import common.Files;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
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

	public AnnotationFileWriter(Project proj, Annotation[] annotations, String annotationFilename, boolean overWriteExisting) {
		super(proj, annotationFilename);
		setAnnotations(annotations);
		this.overWriteExisting = overWriteExisting;
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
			writer.add(vcAnno);
		} else {
			String error = "annotation writer has not been intialized";
			proj.getLog().reportTimeError(error);
			throw new IllegalStateException(error);
		}
	}

	@Override
	public boolean validate() {
		boolean valid = !Files.exists(annotationFilename) || overwriteExisting();
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
		} else {
			proj.getLog().reportTimeError("Detected that " + annotationFilename + " exists and the overwrite option was not flagged");
		}
		return valid;
	}

	@Override
	public void init() {
		if (validate()) {
			VCFHeader vcfHeader = new VCFHeader();
			for (int i = 0; i < annotations.length; i++) {
				VCFInfoHeaderLine vHeaderLine = new VCFInfoHeaderLine(annotations[i].getName(), 1, annotations[i].getType(), annotations[i].getDescription());
				vcfHeader.addMetaDataLine(vHeaderLine);
			}
			// if(analysisInfo!=null){
			// for (int i = 0; i < analysisInfo.length; i++) {
			// VCFHeaderLine vcfHeaderLine = new VCFInfoHeaderLine(name, count, type, description)
			//
			// }
			// }

			VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(annotationFilename);
			builder.clearOptions();
			builder.setOption(Options.INDEX_ON_THE_FLY);
			builder.setOption(Options.DO_NOT_WRITE_GENOTYPES);

			String refGenome = proj.REFERENCE_GENOME_FASTA_FILENAME.getValue();
			proj.getLog().reportTimeInfo("Using reference genome" + refGenome);
			SAMSequenceDictionary samSequenceDictionary = new ReferenceGenome(refGenome, proj.getLog()).getIndexedFastaSequenceFile().getSequenceDictionary();
			MarkerSet markerSet = proj.getMarkerSet();
			if (samSequenceDictionary.getSequenceIndex("chrXY") == -1 && markerSet.getIndicesByChr()[25].length > 0) {
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

	@Override
	public void close() {
		if (writer != null) {
			writer.close();
		}
	}
}
