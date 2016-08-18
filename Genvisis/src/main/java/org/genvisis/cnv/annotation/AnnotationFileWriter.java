package org.genvisis.cnv.annotation;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.WritingFilePrep;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.VCOps;

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

/**
 * @author lane0212 Basically we hijack a vcf file for storing our annotations...seems like most of
 *         the things we want to annotate will have a chr, start and stop
 */
public abstract class AnnotationFileWriter extends AnnotationFile implements WritingFilePrep {

  private VariantContextWriter writer;
  private final boolean overWriteExisting;
  private boolean additionMode;
  private String tmpFile;
  private CloseableIterator<VariantContext> additionReader;

  public AnnotationFileWriter(Project proj, AnalysisParams[] analysisParams,
                              Annotation[] annotations, String annotationFilename,
                              boolean overWriteExisting) {
    super(proj, annotationFilename);
    setAnnotations(annotations);
    setParams(analysisParams);

    this.overWriteExisting = overWriteExisting;
    additionMode = false;
    tmpFile = null;
    additionReader = null;
    init();
  }

  @Override
  public boolean overwriteExisting() {
    return overWriteExisting;
  }

  /**
   * @param locusAnnotation we convert this to a minimal {@link VariantContext} and write to the
   *        underlying {@link VariantContextWriter}
   * @param skipDefaultValue if true, annotations that are set to their default value will not be
   *        added to the {@link VariantContext}
   */
  public void write(LocusAnnotation locusAnnotation, boolean skipDefaultValue,
                    boolean skipAlleleCheck) {
    if (writer != null) {
      VariantContext vcAnno = locusAnnotation.annotateLocus(skipDefaultValue);
      if (vcAnno.getStart() <= 0 || vcAnno.getEnd() <= 0) {
        String error = "Entry " + vcAnno.toStringWithoutGenotypes()
                       + " had postion less than or equal to zero. Most readers will skip";
        proj.getLog().reportTimeError(error);
        throw new IllegalArgumentException(error);
      }
      if (additionReader != null) {
        if (!additionReader.hasNext()) {
          String error =
              "Mismatched number of entries in " + annotationFilename + " , cancelling addition...";
          proj.getLog().reportTimeError(error);
          throw new IllegalStateException(error);
        } else {
          VariantContext vcAdd = additionReader.next();
          // when we skipAlleleCheck, we cannot tell the stop position, based on alleles
          if ((skipAlleleCheck || vcAdd.hasSameAllelesAs(vcAnno))
              && ((skipAlleleCheck
                   && VCOps.getSegment(vcAdd).getStart() == VCOps.getSegment(vcAnno).getStart())
                  || VCOps.getSegment(vcAdd).equals(VCOps.getSegment(vcAnno)))
              && vcAdd.getID().equals(vcAnno.getID())) {
            VariantContextBuilder builder = new VariantContextBuilder(vcAdd);
            for (String att : vcAnno.getAttributes().keySet()) {
              builder.attribute(att, vcAnno.getAttributes().get(att));
            }
            vcAnno = builder.make();
          } else {
            String error = "Entries must be in the exact same order to be added together, ";
            error += vcAdd.getID() + "\t" + vcAdd.toStringWithoutGenotypes() + " from "
                     + annotationFilename;
            error += vcAnno.getID() + "\t" + vcAnno.toStringWithoutGenotypes() + " being added";
            error += "\nAllele Check = " + skipAlleleCheck;
            error += "\nSame alleles = " + vcAdd.hasSameAllelesAs(vcAnno);
            error += "\nSeg add = " + VCOps.getSegment(vcAdd).getUCSClocation();
            error += "\nSeg anno = " + VCOps.getSegment(vcAnno).getUCSClocation();
            error += "\nSeg equals = " + VCOps.getSegment(vcAdd).equals(VCOps.getSegment(vcAnno));
            error += "\nName  add = " + vcAdd.getID();
            error += "\nName  anno = " + vcAnno.getID();
            error += "\nName equals = " + vcAdd.getID().equals(vcAnno.getID());

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
      Files.copyFileUsingFileChannels(new File(annotationFilename),
                                      new File(annotationFilename + "." + ts + ".bak"),
                                      proj.getLog());
      Files.copyFileUsingFileChannels(new File(annotationFilename + ".tbi"),
                                      new File(annotationFilename + ".tbi." + ts + ".bak"),
                                      proj.getLog());
    }

    boolean valid = true;
    if (additionMode) {
      try {
        additionReader = new VCFFileReader(new File(annotationFilename), false).iterator();
        tmpFile = getTmpFile(annotationFilename);
      } catch (Exception e) {
        proj.getLog().reportTimeError("Trying to initialize addition mode, but "
                                      + annotationFilename + " did not pass vcf file checks");
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
        proj.getLog().reportTimeError("Could not find required file "
                                      + proj.REFERENCE_GENOME_FASTA_FILENAME.getValue());
      }
    }

    return valid;
  }

  @Override
  public void init() {
    if (validate()) {
      VCFHeader vcfHeader =
          additionMode ? new VCFFileReader(new File(annotationFilename), true).getFileHeader()
                       : new VCFHeader();// take previous header if needed
      for (int i = 0; i < annotations.length; i++) {
        if (!vcfHeader.hasInfoLine(annotations[i].getName())) {
          VCFInfoHeaderLine vHeaderLine = null;
          if (annotations[i].getCount() != null) {
            vHeaderLine =
                new VCFInfoHeaderLine(annotations[i].getName(), annotations[i].getCount(),
                                      annotations[i].getType(), annotations[i].getDescription());
          } else {
            vHeaderLine =
                new VCFInfoHeaderLine(annotations[i].getName(), annotations[i].getNumber(),
                                      annotations[i].getType(), annotations[i].getDescription());
          }

          vcfHeader.addMetaDataLine(vHeaderLine);
        } else {
          proj.getLog()
              .reportTimeWarning("Detected that info line " + annotations[i].getName()
                                 + " is already present, any new data added will overwrite previous");
        }
      }
      if (params != null) {
        for (AnalysisParams param : params) {
          if (vcfHeader.getOtherHeaderLine(param.getKey()) == null) {
            vcfHeader.addMetaDataLine(param.developHeaderLine());
          }
        }
      }
      VariantContextWriterBuilder builder =
          new VariantContextWriterBuilder().setOutputFile(additionMode ? tmpFile
                                                                       : annotationFilename);// tmp
                                                                                             // file
                                                                                             // if
                                                                                             // needed
      builder.clearOptions();
      builder.setOption(Options.INDEX_ON_THE_FLY);
      builder.setOption(Options.DO_NOT_WRITE_GENOTYPES);

      String refGenome = proj.REFERENCE_GENOME_FASTA_FILENAME.getValue();
      proj.getLog().reportTimeInfo("Using reference genome" + refGenome);

      SAMSequenceDictionary samSequenceDictionary =
          new ReferenceGenome(refGenome, proj.getLog()).getIndexedFastaSequenceFile()
                                                       .getSequenceDictionary();
      SAMSequenceDictionary upDatedSamSequenceDictionary =
          getUpdatedSamSequenceDictionary(proj, samSequenceDictionary);

      builder.setReferenceDictionary(upDatedSamSequenceDictionary);
      vcfHeader.setSequenceDictionary(upDatedSamSequenceDictionary);
      writer = builder.build();
      writer.writeHeader(vcfHeader);

    } else {
      proj.getLog().reportTimeError("Could not intialize annotation file " + annotationFilename);
    }
  }

  /**
   * If the {@link SAMSequenceDictionary} of the reference fasta does not contain a contig for a
   * project, it will cause an error that looks like:<br>
   * java.lang.NullPointerException at
   * htsjdk.tribble.index.tabix.TabixIndexCreator.advanceToReference(TabixIndexCreator.java:116) at
   * htsjdk.tribble.index.tabix.TabixIndexCreator.addFeature(TabixIndexCreator.java:96)
   */
  private static SAMSequenceDictionary getUpdatedSamSequenceDictionary(Project proj,
                                                                       SAMSequenceDictionary samSequenceDictionary) {
    MarkerSet markerSet = proj.getMarkerSet();
    List<SAMSequenceRecord> samSequenceRecords = samSequenceDictionary.getSequences();
    ArrayList<SAMSequenceRecord> updatedRecords = new ArrayList<SAMSequenceRecord>();
    SAMSequenceRecord mitoRecord = null;
    int currentIndex = 0;
    int[][] indicesByChr = markerSet.getIndicesByChr();
    if (samSequenceDictionary.getSequenceIndex("chr0") == -1 && indicesByChr[0].length > 0) {// not
                                                                                             // always
                                                                                             // present
                                                                                             // in
                                                                                             // ref
      int[] chr0Len = Array.subArray(markerSet.getPositions(), indicesByChr[0]);
      proj.getLog()
          .reportTimeInfo("Since the project contained markers designated as chr0, a chr0 contig is being added");
      if (Array.countIf(chr0Len, 0) > 0) {
        proj.getLog()
            .reportTimeWarning("VCF files cannot have positions of 0, positions will be updated to chr0:1");
      }
      SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord("chr0", Array.max(chr0Len) + 1);
      samSequenceRecord.setSequenceIndex(currentIndex);
      updatedRecords.add(samSequenceRecord);
      currentIndex++;
    }
    for (SAMSequenceRecord samSequenceRecord : samSequenceRecords) {
      if (!samSequenceRecord.getSequenceName().equals("chrM")) {

        int contig = Positions.chromosomeNumber(samSequenceRecord.getSequenceName());
        if (contig > 0) {
          samSequenceRecord.setSequenceIndex(currentIndex);

          int contigProjLength =
              Array.max(Array.subArray(markerSet.getPositions(), indicesByChr[contig]));
          if (samSequenceRecord.getSequenceLength() < contigProjLength) {
            proj.getLog()
                .reportTimeError(samSequenceRecord.getSequenceName() + " had length "
                                 + samSequenceRecord.getSequenceLength()
                                 + " but the project had a max length of " + contigProjLength
                                 + " ,please choose check your reference build, but will update for now");
            return null;
            // samSequenceRecord.setSequenceLength(contigProjLength);
          }
          currentIndex++;
          updatedRecords.add(samSequenceRecord);
        } else {
          proj.getLog().reportTimeWarning(samSequenceRecord.getSequenceName()
                                          + " is an invalid chromosome for Genvisis, skipping");
        }
      } else {
        mitoRecord = samSequenceRecord;
      }
    }

    if (samSequenceDictionary.getSequenceIndex("chrXY") == -1 && indicesByChr[25].length > 0) {// not
                                                                                               // always
                                                                                               // present
                                                                                               // in
                                                                                               // ref
      int[] chrxyLen = Array.subArray(markerSet.getPositions(), indicesByChr[25]);
      proj.getLog()
          .reportTimeInfo("Since the project contained markers designated as pseudo-autosomal, a chrXY contig is being added");
      SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord("chrXY", Array.max(chrxyLen) + 1);
      samSequenceRecord.setSequenceIndex(currentIndex);
      currentIndex++;
      updatedRecords.add(samSequenceRecord);
    }
    if (indicesByChr[26].length > 0) {// not always present in ref
      int[] chrMLen = Array.subArray(markerSet.getPositions(), indicesByChr[26]);
      proj.getLog()
          .reportTimeInfo("Since the project contained markers designated as mitochondrial, a chrM entry is being added");
      mitoRecord =
          mitoRecord != null ? mitoRecord : new SAMSequenceRecord("chrM", Array.max(chrMLen) + 1);
      mitoRecord.setSequenceIndex(currentIndex);
      updatedRecords.add(mitoRecord);
    }

    return new SAMSequenceDictionary(updatedRecords);
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
      proj.getLog().reportTimeInfo("Copying temporary file " + tmpFile + " to " + annotationFilename
                                   + ", " + tmpFile + " can be deleted on successful completion");
      Files.copyFileUsingFileChannels(new File(tmpFile), new File(annotationFilename),
                                      proj.getLog());
      Files.copyFileUsingFileChannels(new File(tmpFile + ".tbi"),
                                      new File(annotationFilename + ".tbi"), proj.getLog());
      new File(tmpFile).delete();
      new File(tmpFile + ".tbi").delete();
    }
  }
}
