package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.BoolRequirement;
import org.genvisis.cnv.workflow.Requirement.DirRequirement;
import org.genvisis.cnv.workflow.Requirement.FileRequirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.genvisis.seq.manage.mosdepth.CRAMSnpReader;
import org.genvisis.seq.manage.mosdepth.FASTAToBedConversion;
import org.genvisis.seq.manage.mosdepth.MosdepthImport;
import org.genvisis.seq.manage.mosdepth.MosdepthPipeline;
import org.genvisis.seq.manage.mosdepth.NGSBinSNPSelector;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;

public class MosdepthImportStep extends Step {

  public static final String NAME = "Mosdepth Import";
  public static final String DESC = "";

  public static MosdepthImportStep create(Project proj) {

    FileRequirement binBedReq = new Requirement.FileRequirement("binBed",
                                                                "BED File defining bins (used to select SNPs and run mosdepth)",
                                                                "Can be created with " + FASTAToBedConversion.class.getName(),
                                                                new File("bins.bed"));

    FileRequirement selectedSNPReq = new Requirement.FileRequirement("selectedSNP",
                                                                     "VCF file containing selected SNPs for each bin",
                                                                     "Can be created with " + NGSBinSNPSelector.class.getName(),
                                                                     new File("selected.vcf"));

    DirRequirement mosDirReq = new Requirement.DirRequirement("mosDir",
                                                              "Directory with mosdepth output (.bed files)",
                                                              new File("00mos"));

    DirRequirement cntsDirReq = new Requirement.DirRequirement("cntsDir",
                                                               "Directory with allele count output (.bed files)",
                                                               "Created with " + CRAMSnpReader.class.getName(),
                                                               new File("00cnt"));

    BoolRequirement dontImportGenosReq = new BoolRequirement("noGenos", "Do not import genotypes",
                                                             true);
    FileRequirement genotypeVCF = new FileRequirement("genoVCF",
                                                      "VCF file containing genotype values for selected SNPs",
                                                      new File(""));

    return new MosdepthImportStep(proj, mosDirReq, cntsDirReq, binBedReq, selectedSNPReq,
                                  dontImportGenosReq, genotypeVCF);
  }

  private final Project proj;
  private final DirRequirement mosDirReq;
  private final DirRequirement alleleCountReq;
  private final FileRequirement binBedReq;
  private final FileRequirement snpVCFReq;

  private MosdepthImportStep(Project proj, DirRequirement mosReq, DirRequirement alleleReq,
                             FileRequirement binBedReq, FileRequirement snpVCFReq,
                             BoolRequirement dontImportGenosReq, FileRequirement genotypeVCF) {
    super(NAME, DESC,
          RequirementSetBuilder.and().add(mosReq).add(alleleReq).add(binBedReq).add(snpVCFReq)
                               .add(RequirementSetBuilder.or().add(dontImportGenosReq)
                                                         .add(genotypeVCF)),
          EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
                     Requirement.Flag.MULTITHREADED));
    this.proj = proj;
    this.mosDirReq = mosReq;
    this.alleleCountReq = alleleReq;
    this.snpVCFReq = snpVCFReq;
    this.binBedReq = binBedReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // No-op
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    return Files.ensurePathExists(proj.MARKER_DATA_DIRECTORY.getValue())
           && Files.exists(proj.MARKER_DATA_DIRECTORY.getValue() + "outliers.ser");
  }

  @Override
  public String getCommandLine(Variables variables) {
    StringBuilder cmd = new StringBuilder(Files.getRunString());
    cmd.append(" ").append(MosdepthImport.class.getName());
    cmd.append(" propDir=").append(ext.parseDirectoryOfFile(proj.getPropertyFilename()));
    cmd.append(" projName=").append(proj.PROJECT_NAME.getValue());
    cmd.append(" projDir=").append(proj.PROJECT_DIRECTORY.getValue());
    cmd.append(" mosDir=").append(variables.get(mosDirReq).getPath());
    cmd.append(" cramReadsDir=").append(variables.get(alleleCountReq).getPath());
    cmd.append(" selectedVCF=").append(variables.get(snpVCFReq).getPath());
    cmd.append(" binBed=").append(variables.get(binBedReq).getPath());
    cmd.append(" build=").append(proj.GENOME_BUILD_VERSION.getValue());
    cmd.append(" mosExt=.bed");
    return cmd.toString();
  }

  @Override
  public void run(Variables variables) {
    MosdepthPipeline mi = new MosdepthPipeline();
    mi.setProjectDir(proj.PROJECT_DIRECTORY.getValue());
    mi.setProjectName(proj.PROJECT_NAME.getValue());
    mi.setProjectPropertiesDir(ext.parseDirectoryOfFile(proj.getPropertyFilename()));
    mi.setNumThreads(Runtime.getRuntime().availableProcessors());

    mi.setBinsToUseBED(variables.get(binBedReq).getPath());
    // mi.setGenotypeVCF("G:\\bamTesting\\EwingWGS\\ES_recalibrated_snps_indels.vcf.gz");
    mi.setSelectedMarkerVCF(variables.get(snpVCFReq).getPath());
    mi.setMosdepthDirectory(variables.get(mosDirReq).getPath(), ".bed");
    mi.setCRAMReadDirectory(variables.get(alleleCountReq).getPath());
    mi.setJobID(null);
    try {
      mi.run();
    } catch (IOException | Elision e) {
      e.printStackTrace();
    }
  }

}
