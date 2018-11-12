package org.genvisis.seq.manage.mosdepth;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.AllelePair;
import org.genvisis.common.Files;
import org.genvisis.common.GenomicPosition;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.BEDFileReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class MosdepthPipeline {

  Logger log;

  String projDir;
  String propFileDir;
  String projName;
  String mosSrcDir;
  String mosSrcExt;
  double scaleFactor = 10;

  public void setProjectName(String projName2) {
    projName = projName2;
  }

  public void setProjectPropertiesDir(String propFileDir2) {
    propFileDir = propFileDir2;
  }

  public void setProjectDir(String projDir2) {
    projDir = projDir2;
  }

  Project proj;

  String useBed;
  String markerVCF;
  String genoVCF;

  public void setBinsToUseBED(String bedFile) {
    this.useBed = bedFile;
  }

  public void setSelectedMarkerVCF(String vcfFile) {
    this.markerVCF = vcfFile;
  }

  public void setGenotypeVCF(String vcfFile) {
    this.genoVCF = vcfFile;
  }

  public void setMosdepthDirectory(String dir, String ext) {
    this.mosSrcDir = org.genvisis.common.ext.verifyDirFormat(dir);
    // TODO if ext != ".norm" (or whatever the default is in MosdepthNormalizer), throw an error and return, with a suggestion to either set override flag (to allow raw) or to process with MosdepthNormalizer
    this.mosSrcExt = ext;
    this.mosdepthFiles = Files.list(dir, ext);
  }

  Segment[] useBins;
  VCFFileReader snpReader;
  VCFFileReader genoReader;

  String[] mosdepthFiles;

  MarkerDetailSet snpPosMDS;
  MarkerDetailSet binPosMDS;
  List<Marker> snpMarkers = new ArrayList<>();
  List<Marker> binMarkers = new ArrayList<>();

  void run() {
    if (log == null) {
      log = new Logger();
    }
    createProject();

    loadBins();

    loadSNPs();

    //    genoReader = new VCFFileReader(new File(genoVCF), true);

    // TODO setup mosdepth readers
    // TODO check bins validity across files

  }

  protected void createProject() {
    String propFile = ext.verifyDirFormat(propFileDir)
                      + ext.replaceWithLinuxSafeCharacters(projName) + ".properties";
    if (!Files.exists(propFile)) {
      Files.write((new Project()).PROJECT_NAME.getName() + "=" + projName, propFile);
      proj = new Project(propFile);
      proj.PROJECT_NAME.setValue(projName);
      proj.PROJECT_DIRECTORY.setValue(projDir);
      proj.SOURCE_DIRECTORY.setValue(mosSrcDir);
      proj.XY_SCALE_FACTOR.setValue(scaleFactor);
      proj.TARGET_MARKERS_FILENAMES.setValue(new String[] {});
      proj.SOURCE_FILENAME_EXTENSION.setValue(mosSrcExt);
      proj.ID_HEADER.setValue("NULL");
      proj.GENOME_BUILD_VERSION.setValue(GENOME_BUILD.HG19);

      proj.ARRAY_TYPE.setValue(ARRAY.NGS);

      proj.saveProperties();
      log.reportTime("Created project properties file: " + propFile);
    } else {
      log.reportTime("Project properties file already exists at " + propFile
                     + "; skipping creation.");
      proj = new Project(propFile);
    }
  }

  private void loadBins() {
    // 1-based indexing!
    BEDFileReader reader = new BEDFileReader(useBed, false);
    useBins = reader.loadAll(log).getStrictSegments();
    reader.close();
    log.reportTime("Loaded list of bins to be imported.");
  }

  private void loadSNPs() {
    snpReader = new VCFFileReader(new File(markerVCF), true);
    for (Segment seg : useBins) {
      List<VariantContext> markerVC = snpReader.query(seg.getChromosomeUCSC(), seg.getStart() - 1,
                                                      seg.getStop())
                                               .toList();
      if (markerVC.size() == 0) {
        log.reportTimeWarning("No snp found for bin " + seg.getUCSClocation()
                              + ". Bin will be skipped.");
        continue;
      }
      if (markerVC.size() > 1) {
        log.reportError("Multiple markers-to-use found for bin " + seg.getUCSClocation()
                        + ".  Bin will be skipped.");
        continue;
      }
      VariantContext chosenSnp = markerVC.get(0);

      Marker mSnpPos = new Marker(chosenSnp.getID(),
                                  new GenomicPosition(Positions.chromosomeNumber(chosenSnp.getContig()),
                                                      chosenSnp.getStart()),
                                  AllelePair.of(chosenSnp.getReference(),
                                                chosenSnp.getAlternateAlleles().get(0)));
      int stt = new Integer(chosenSnp.getAttribute("BINSTART").toString()).intValue();
      int stp = new Integer(chosenSnp.getAttribute("BINSTOP").toString()).intValue();
      int mid = stt + ((stp - stt) / 2) + 1;

      Marker mBinPos = new Marker(chosenSnp.getID(),
                                  new GenomicPosition(Positions.chromosomeNumber(chosenSnp.getContig()),
                                                      mid),
                                  AllelePair.of(chosenSnp.getReference(),
                                                chosenSnp.getAlternateAlleles().get(0)));

      snpMarkers.add(mSnpPos);
      binMarkers.add(mBinPos);
    }
    snpReader.close();

    MarkerDetailSet mdsSnp = new MarkerDetailSet(snpMarkers);
    mdsSnp.serialize(getMDSName(false));
    MarkerDetailSet mdsBin = new MarkerDetailSet(binMarkers);
    mdsBin.serialize(getMDSName(true));

    log.reportTime("Found " + snpMarkers.size() + " markers to be parsed.");
  }

  private String getMDSName(boolean binnedVersion) {
    return binnedVersion ? ext.rootOf(proj.MARKER_DETAILS_FILENAME.getValue(), false) + "_bins.ser"
                         : proj.MARKER_DETAILS_FILENAME.getValue();
  }

  public static void main(String[] args) {
    MosdepthPipeline mi = new MosdepthPipeline();
    mi.setProjectDir("G:\\bamTesting\\mosdepth\\project\\");
    mi.setProjectName("EwingsMosdepth");
    mi.setProjectPropertiesDir("D:\\projects\\");

    mi.setBinsToUseBED("G:\\bamTesting\\snpSelection\\ReferenceGenomeBins.bed");
    mi.setGenotypeVCF("G:\\bamTesting\\EwingWGS\\ES_recailbrated_snps_indels.vcf.gz");
    mi.setSelectedMarkerVCF("G:\\bamTesting\\snpSelection\\selected.vcf");
    mi.setMosdepthDirectory("G:\\bamTesting\\mosdepth\\00src\\", ".norm");
    mi.run();
  }

}
