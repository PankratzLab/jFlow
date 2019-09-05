package org.genvisis.cnv.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Vector;
import java.util.stream.Collectors;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.genvisis.cnv.LaunchProperties;
import org.genvisis.cnv.Resources;
import org.genvisis.cnv.analysis.pca.PCImputeRace;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.gwas.pca.ancestry.AncestryPCA;
import org.genvisis.cnv.gwas.pca.ancestry.PlinkDataMatrixLoader;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.var.SampleData;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;
import org.pankratzlab.utils.filesys.SnpMarkerSet;
import org.pankratzlab.utils.gwas.RelationAncestryQc;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.ImmutableTable;
import com.google.common.collect.Maps;
import com.google.common.collect.MoreCollectors;
import com.google.common.collect.Sets;
import com.google.common.collect.Table;

public class Ancestry {

  private static final String PLINK_BIM_UNAMBIGUOUS_RENAMED_TXT = "plink.bim_unambiguous_renamed.txt";
  private static final String PLINK_BIM_UNAMBIGUOUS_TXT = "plink.bim_unambiguous.txt";
  public static final int DEFAULT_NUM_COMPONENTS_ANCESTRY = 10;
  public static final String DEFAULT_HAPMAP_PLINKROOT = "/home/pankrat2/shared/bin/HapMap/unambiguousHapMapFounders";
  public static final String RACE_IMPUTATIONS_FILENAME = "raceImputations.mds";
  public static final String RACE_FREQS_FILENAME = "freqsByRace.xln";

  private static final String HAPMAP_CLASS_NAME = "HapMap";
  private static final String PCA_OUTPUT_NAME = "combo.extrapolatedPCs.txt.gz";
  private static final String EIGENSTRAT_IID_LABEL = "IID";
  private static final String EIGENSTRAT_FID_LABEL = "FID";
  private static final String PCA_PC1_LABEL = "PC1";
  private static final String PCA_PC2_LABEL = "PC2";

  private static final Set<String> POSSIBLE_MISSNP_EXTS = ImmutableSet.of(".missnp",
                                                                          "-merge.missnp");

  private static enum HapMapPopulation {
    CEU, YRI, CHB, JPT;
  }

  private final String dir;
  private final Project proj;
  private final String plinkroot;
  private final String plinkExe;
  private final Logger log;

  public Ancestry(String dir, Project proj, String plinkExe) {
    this(dir, "plink", proj, plinkExe, proj.getLog());
  }

  public Ancestry(String dir, String plinkroot, Project proj, String plinkExe, Logger log) {
    super();
    this.dir = new File(ext.verifyDirFormat(dir)).getAbsolutePath() + File.separator;
    this.proj = proj;
    this.plinkroot = plinkroot;
    this.plinkExe = plinkExe;
    this.log = log;
  }

  public Ancestry(String dir, String dummyProjectPrefix, String plinkExe, Logger log) {
    this(dir, "plink", createDummyProject(dir, dummyProjectPrefix), plinkExe, log);
  }

  private static Project createDummyProject(String dir, String dummyProjectPrefix) {
    String projectName = dummyProjectPrefix + "_AncestryResults";

    final Project dummyProject;

    if (Files.exists(LaunchProperties.formProjectPropertiesFilename(projectName))) {
      dummyProject = new Project(LaunchProperties.formProjectPropertiesFilename(projectName));
      dummyProject.getLog()
                  .reportTimeWarning("Project with name " + projectName
                                     + " already exists, continuing with existing project");
    } else {
      dummyProject = Project.initializeProject(projectName, dir);
    }

    String plinkFamFile = dir + "plink.fam";
    String[][] plinkFam = HashVec.loadFileToStringMatrix(plinkFamFile, false, null);
    String[] dummySampleList = new String[plinkFam.length];
    String[][] dummyPedigree = new String[plinkFam.length][];
    for (int i = 0; i < plinkFam.length; i++) {
      String dummyDNA = plinkFam[i][0] + "_" + plinkFam[i][1];
      dummySampleList[i] = dummyDNA;
      dummyPedigree[i] = ArrayUtils.addStrToArray(dummyDNA, plinkFam[i]);
    }
    new SampleList(dummySampleList).serialize(dummyProject.SAMPLELIST_FILENAME.getValue());
    Files.writeMatrix(dummyPedigree, dummyProject.PEDIGREE_FILENAME.getValue(), "\t");
    try {
      SampleData.createSampleData(dummyProject.PEDIGREE_FILENAME.getValue(), null, dummyProject);
    } catch (Elision e) {
      dummyProject.getLog()
                  .reportError("Could not generate Sample Data for dummy project "
                               + dummyProject.getPropertyFilename()
                               + ", ancestry results cannot be visualized");
    }
    return dummyProject;
  }

  private void maybeAddHapMapToSampleData(SampleData.ClassHeader hapMapClassHeader,
                                          Table<String, String, HapMapPopulation> hapmaps) {
    if (!proj.getSampleData(false).hasClass(HAPMAP_CLASS_NAME)) {
      String[] hapMapColumnHeaders = new String[] {"FID", "IID", "DNA",
                                                   hapMapClassHeader.toString()};
      Map<String, Integer> hapMapCodes = hapMapClassHeader.getCodeOptions().inverse();
      String[][] hapMapColumnData = new String[hapmaps.size()][hapMapColumnHeaders.length];
      int i = 0;
      for (Table.Cell<String, String, HapMapPopulation> hapMapCell : hapmaps.cellSet()) {
        String fid = hapMapCell.getRowKey();
        String iid = hapMapCell.getColumnKey();
        HapMapPopulation pop = hapMapCell.getValue();
        hapMapColumnData[i++] = new String[] {fid, iid, iid,
                                              hapMapCodes.get(pop.toString()).toString()};
      }

      proj.getSampleData(false).addSamples(hapMapColumnData, hapMapColumnHeaders);
    }
  }

  public void runPipeline(String putativeWhitesFile, String hapMapPlinkRoot,
                          String snpRSIDLookupFile) {
    if (hapMapPlinkRoot == null) {
      hapMapPlinkRoot = DEFAULT_HAPMAP_PLINKROOT;
    }
    if (!Files.exists(dir + "homogeneity/" + MergeDatasets.CHI_SQUARE_DROPS_FILENAME)
        && Files.list(dir + "homogeneity/", ".Rout").length == 0) {
      log.report("Running homogeneity checks...");
      checkHomogeneity(putativeWhitesFile, hapMapPlinkRoot);
    }
    String homogeneityDrops = parseHomogeneity();
    mergeHapMap(hapMapPlinkRoot, homogeneityDrops, snpRSIDLookupFile);
    runPCA(DEFAULT_NUM_COMPONENTS_ANCESTRY);
    imputeRaceFromPCA();
  }

  private void runPCA(int numComps) {
    setupAncestry();
    AncestryPCA ancestryPCA = AncestryPCA.generatePCs(new PlinkDataMatrixLoader(dir,
                                                                                "unrelateds/plink",
                                                                                log),
                                                      numComps, log);
    ancestryPCA.getSvd().dumpLoadingsToText(dir + "combo", "MARKER", log);
    AncestryPCA.extrapolatePCs(ancestryPCA, new PlinkDataMatrixLoader(dir, "combo", log), log)
               .dumpToText(dir + PCA_OUTPUT_NAME, "FID\tIID", log);
  }

  private void checkHomogeneity(String putativeWhitesFile, String hapMapPlinkRoot) {
    String homoDir = dir + "homogeneity/";
    String homoProjDir = homoDir + ext.removeDirectoryInfo(plinkroot) + "/";
    String homoHapMapDir = homoDir + ext.removeDirectoryInfo(hapMapPlinkRoot) + "/";
    new File(homoProjDir).mkdirs();
    new File(homoHapMapDir).mkdirs();
    String cleanPutativeWhitesFile = validatePutativeWhites(dir + plinkroot + PSF.Plink.FAM,
                                                            putativeWhitesFile);
    CmdLine.runDefaults(plinkExe + " --bfile " + dir + plinkroot + " --keep "
                        + cleanPutativeWhitesFile + " --hardy --noweb", homoProjDir, log);
    CmdLine.runDefaults(plinkExe + " --bfile " + hapMapPlinkRoot + " --keep "
                        + ext.parseDirectoryOfFile(hapMapPlinkRoot)
                        + "CEUFounders.txt --hardy --noweb", homoHapMapDir, log);

    MergeDatasets.checkForHomogeneity(homoDir, null, null, "UNAFF", log);
  }

  private String validatePutativeWhites(String projFamFile, String putativeWhitesFile) {
    if (proj == null) return putativeWhitesFile;
    PlinkData.ExportIDScheme projIDScheme = PlinkData.detectExportIDScheme(proj, projFamFile);
    if (projIDScheme == null) {
      log.reportTimeWarning("Could not validate project fam file to project IDs, assuming FID/IID matches plink");
      return putativeWhitesFile;
    }
    String cleanPutativeWhitesFile = PlinkData.convertIDScheme(proj, putativeWhitesFile,
                                                               projIDScheme);
    if (cleanPutativeWhitesFile == null) {
      log.reportTimeWarning("Could not convert " + putativeWhitesFile
                            + " to match ID scheme, leaving as-is");
      return putativeWhitesFile;
    }
    return cleanPutativeWhitesFile;
  }

  private String parseHomogeneity() {
    int rOuts = Files.list(dir + "homogeneity/", ".Rout").length;
    if (rOuts == 0) {
      log.report("No Fisher's Exact results found, using Chi Square to choose homogeneous markers");
      return dir + "homogeneity/" + MergeDatasets.CHI_SQUARE_DROPS_FILENAME;
    }
    MergeDatasets.parseHomo(dir + "homogeneity/");
    return dir + "homogeneity/" + MergeDatasets.FISHER_OR_CHI_SQUARE_DROPS_FILENAME;
  }

  private void mergeHapMap(String hapMapPlinkRoot, String dropMarkersFile, String snpIDLookupFile) {
    if (!Files.exists(dir + RelationAncestryQc.UNRELATEDS_FILENAME)) {
      log.reportError("Error - need a file called " + RelationAncestryQc.UNRELATEDS_FILENAME
                      + " with FID and IID pairs before we can proceed");
      return;
    }

    BidiMap<String, String> lookup = null;
    if (snpIDLookupFile != null) {
      if (!Files.exists(snpIDLookupFile)) {
        log.reportError("SNP Name Lookup file was specified, but couldn't be found at: "
                        + snpIDLookupFile);
      } else {
        lookup = new DualHashBidiMap<String, String>(HashVec.loadFileColumnToMap(snpIDLookupFile, 0,
                                                                                 1, false, log));
      }
    }

    String srcData = plinkroot;
    if (lookup != null) {
      srcData = plinkroot + "_renamed";
      log.report(ext.getTime() + "]\tRenaming snps using lookup file: " + snpIDLookupFile);
      CmdLine.runDefaults(plinkExe + " --no-web --bfile " + plinkroot + " --update-name "
                          + snpIDLookupFile + " --allow-no-sex --make-bed --out " + srcData, dir,
                          log);
    }

    if (!Files.exists(dir + PLINK_BIM_UNAMBIGUOUS_TXT)) {
      log.report(ext.getTime() + "]\tGenerating list of unambiguous SNPs");
      new SnpMarkerSet(dir + srcData + ".bim", true,
                       log).listUnambiguousMarkers(dir + PLINK_BIM_UNAMBIGUOUS_TXT, dropMarkersFile,
                                                   true);
    }

    if (!Files.exists(dir + "unambiguousHapMap.bed")) {
      log.report(ext.getTime() + "]\tExtracting unambiguous SNPs for HapMap founders");
      CmdLine.runDefaults(plinkExe + " --bfile " + hapMapPlinkRoot + " --extract "
                          + PLINK_BIM_UNAMBIGUOUS_TXT
                          + " --make-bed --allow-no-sex --out unambiguousHapMap --noweb", dir, log);
    }

    if (!Files.exists(dir + "overlap.txt")) {
      log.report(ext.getTime() + "]\tGenerating list of overlapping SNPs");
      Files.writeIterable(HashVec.loadFileToVec(dir + "unambiguousHapMap.bim", false, new int[] {1},
                                                false),
                          dir + "overlap.txt");
    }
    POSSIBLE_MISSNP_EXTS.stream().map(ext -> dir + "combo" + ext).filter(Files::exists)
                        .map(File::new).forEach(File::delete);

    if (!Files.exists(dir + "unambiguous.bed")) {
      log.report(ext.getTime() + "]\tExtracting overlapping SNPs for study samples");
      CmdLine.runDefaults(plinkExe + " --bfile " + srcData
                          + " --extract overlap.txt --make-bed --allow-no-sex --out unambiguous --noweb",
                          dir, log);
    }
    log.report(ext.getTime() + "]\tMerging study data and HapMap data for overlapping SNPs");
    CmdLine.runDefaults(plinkExe
                        + " --bfile unambiguous --bmerge unambiguousHapMap.bed unambiguousHapMap.bim unambiguousHapMap.fam --make-bed --out combo --noweb",
                        dir, log);

    if (POSSIBLE_MISSNP_EXTS.stream().map(ext -> dir + "combo" + ext).anyMatch(Files::exists)) {
      POSSIBLE_MISSNP_EXTS.stream().map(ext -> dir + "combo.1" + ext).filter(Files::exists)
                          .map(File::new).forEach(File::delete);
      if (Files.exists(dir + "combo.1.missnp")) {
        new File(dir + "combo.1.missnp").delete();
      }
      log.report(ext.getTime() + "]\tChecking for flipped alleles");
      POSSIBLE_MISSNP_EXTS.stream().map(ext -> dir + "combo" + ext).filter(Files::exists)
                          .map(File::new).collect(MoreCollectors.onlyElement())
                          .renameTo(new File(dir + "combo.1.missnp"));
      CmdLine.runDefaults(plinkExe
                          + " --bfile unambiguous --flip combo.1.missnp --make-bed --allow-no-sex --out unambiguousFlipped --noweb",
                          dir, log);
      CmdLine.runDefaults(plinkExe
                          + " --bfile unambiguousFlipped --bmerge unambiguousHapMap.bed unambiguousHapMap.bim unambiguousHapMap.fam --make-bed --allow-no-sex --out combo --noweb",
                          dir, log);
      if (POSSIBLE_MISSNP_EXTS.stream().map(ext -> dir + "combo" + ext).anyMatch(Files::exists)) {
        if (Files.exists(dir + "combo.2.missnp")) {
          new File(dir + "combo.2.missnp").delete();
        }
        log.report(ext.getTime() + "]\tDropping SNPs that cannot be resolved by flipping alleles");
        POSSIBLE_MISSNP_EXTS.stream().map(ext -> dir + "combo" + ext).filter(Files::exists)
                            .map(File::new).collect(MoreCollectors.onlyElement())
                            .renameTo(new File(dir + "combo.2.missnp"));
        CmdLine.runDefaults(plinkExe
                            + " --bfile unambiguousFlipped --exclude combo.2.missnp --make-bed --allow-no-sex --out unambiguousFlippedDropped --noweb",
                            dir, log);
        CmdLine.runDefaults(plinkExe
                            + " --bfile unambiguousHapMap --exclude combo.2.missnp --make-bed --allow-no-sex --out unambiguousDroppedHapMap --noweb",
                            dir, log);
        CmdLine.runDefaults(plinkExe
                            + " --bfile unambiguousFlippedDropped --bmerge unambiguousDroppedHapMap.bed unambiguousDroppedHapMap.bim unambiguousDroppedHapMap.fam --make-bed --allow-no-sex --out combo --noweb",
                            dir, log);
      }
    }

    if (!Files.exists(dir + "finalSNPs.txt")) {
      log.report(ext.getTime() + "]\tWriting final list of SNPs to use");
      Files.writeIterable(HashVec.loadFileToVec(dir + "combo.bim", false, new int[] {1}, false),
                          dir + "finalSNPs.txt");
    }
  }

  private void setupAncestry() {
    Logger log;

    log = new Logger(dir + "ancestry.log");

    String unrelatedsDir = dir + "unrelateds/";

    if (!Files.exists(unrelatedsDir + RelationAncestryQc.UNRELATEDS_FILENAME)) {
      log.report(ext.getTime() + "]\tGenerating combined "
                 + RelationAncestryQc.UNRELATEDS_FILENAME);
      new File(unrelatedsDir).mkdir();
      Vector<String> unrelateds = HashVec.loadFileToVec(dir + "unambiguousHapMap.fam", false,
                                                        new int[] {0, 1}, false);
      unrelateds.addAll(HashVec.loadFileToVec(dir + RelationAncestryQc.UNRELATEDS_FILENAME, false,
                                              false, false));
      Files.writeIterable(unrelateds, unrelatedsDir + RelationAncestryQc.UNRELATEDS_FILENAME);
    }

    if (!Files.exists(unrelatedsDir + "plink.bed")) {
      log.report(ext.getTime() + "]\tGenerating PLINK files based on combined "
                 + RelationAncestryQc.UNRELATEDS_FILENAME);
      CmdLine.runDefaults(plinkExe + " --bfile ../combo --keep "
                          + RelationAncestryQc.UNRELATEDS_FILENAME
                          + " --make-bed --allow-no-sex --noweb", unrelatedsDir, log);
    }

  }

  private void imputeRaceFromPCA() {
    if (!Files.exists(dir + RACE_IMPUTATIONS_FILENAME)) {
      Table<String, String, HapMapPopulation> fidIidHapMapPopTable = parseHapMapAncestries();
      if (fidIidHapMapPopTable == null) return;

      Set<PCImputeRace.Sample> samples = Sets.newHashSet();
      Set<PCImputeRace.Sample> europeans = Sets.newHashSet();
      Set<PCImputeRace.Sample> africans = Sets.newHashSet();
      Set<PCImputeRace.Sample> asians = Sets.newHashSet();
      try (BufferedReader eigenReader = Files.getAppropriateReader(dir + PCA_OUTPUT_NAME)) {
        String headerLine = eigenReader.readLine();
        String delim = ext.determineDelimiter(headerLine);
        String[] header = headerLine.split("\t");
        Map<String, Integer> eigenstratIndices = ext.indexMap(new String[] {EIGENSTRAT_FID_LABEL,
                                                                            EIGENSTRAT_IID_LABEL,
                                                                            PCA_PC1_LABEL,
                                                                            PCA_PC2_LABEL},
                                                              header, true, false);

        while (eigenReader.ready()) {
          String[] line = eigenReader.readLine().split(delim);
          String fid = line[eigenstratIndices.get(EIGENSTRAT_FID_LABEL)];
          String iid = line[eigenstratIndices.get(EIGENSTRAT_IID_LABEL)];
          double pc1;
          double pc2;
          try {
            pc1 = Double.parseDouble(line[eigenstratIndices.get(PCA_PC1_LABEL)]);
          } catch (NumberFormatException nfe) {
            pc1 = Double.NaN;
          }
          try {
            pc2 = Double.parseDouble(line[eigenstratIndices.get(PCA_PC2_LABEL)]);
          } catch (NumberFormatException nfe) {
            pc2 = Double.NaN;
          }

          PCImputeRace.Sample sample = new PCImputeRace.Sample(fid, iid, pc1, pc2);
          samples.add(sample);
          HapMapPopulation hapMapPop = fidIidHapMapPopTable.get(fid, iid);
          if (hapMapPop != null) {
            switch (hapMapPop) {
              case CEU:
                europeans.add(sample);
                break;
              case YRI:
                africans.add(sample);
                break;
              case CHB:
              case JPT:
                asians.add(sample);
                break;
              default:
                log.reportError("Unexpected HapMap population ignored: " + hapMapPop);
                break;
            }
          }
        }
      } catch (FileNotFoundException e) {
        log.reportFileNotFound(dir + PCA_OUTPUT_NAME);
      } catch (IOException e) {
        log.reportIOException(dir + PCA_OUTPUT_NAME);
      }

      PCImputeRace pcir = new PCImputeRace(proj, samples, europeans, africans, asians, log);
      pcir.correctPCsToRace(dir + RACE_IMPUTATIONS_FILENAME);
    } else {
      log.reportTimeWarning("Skipping imputation - output already exists: "
                            + (dir + RACE_IMPUTATIONS_FILENAME));
    }

    if (!Files.exists(dir + RACE_FREQS_FILENAME)) {
      freqsByRace(dir + RACE_IMPUTATIONS_FILENAME, dir + RACE_FREQS_FILENAME);
    } else {
      log.reportTimeWarning("Skipping race freq calculation - output already exists: "
                            + (dir + RACE_FREQS_FILENAME));
    }

  }

  private Table<String, String, HapMapPopulation> parseHapMapAncestries() {
    String hapMapAncestries = Resources.hapMap(log).getHapMapAncestries().get();
    if (hapMapAncestries == null) {
      log.reportError("Cannot impute race without the HapMap ancestries resource");
      return null;
    }
    try (BufferedReader hapMapAncestryReader = Files.getAppropriateReader(hapMapAncestries)) {
      String headerLine = hapMapAncestryReader.readLine();
      String delimiter = ext.determineDelimiter(headerLine);
      String[] ancestriesHeader = headerLine.split(delimiter);
      int fidIndex = ext.indexOfStr("FID", ancestriesHeader, false, true);
      int iidIndex = ext.indexOfStr("IID", ancestriesHeader, false, true);
      int hapIndex = ext.indexOfStartsWith(SampleData.ClassHeader.generateClassHeaderLabel(HAPMAP_CLASS_NAME),
                                           ancestriesHeader, false);
      if (Collections.max(Arrays.asList(fidIndex, iidIndex, hapIndex)) < 0) {
        log.reportError("Cannot impute: malformed HapMap ancestries resource: " + hapMapAncestries);
        return null;
      }
      SampleData.ClassHeader hapMapClassHeader = new SampleData.ClassHeader(ancestriesHeader[hapIndex],
                                                                            log);
      Map<Integer, HapMapPopulation> hapMapCodeMap;
      try {
        hapMapCodeMap = hapMapClassHeader.getCodeOptions().entrySet().stream()
                                         .collect(Collectors.toMap(Map.Entry::getKey,
                                                                   e -> HapMapPopulation.valueOf(e.getValue())));
      } catch (IllegalArgumentException iae) {
        log.reportError("Unexpected HapMap population");
        log.reportException(iae);
        return null;
      }
      ImmutableTable.Builder<String, String, HapMapPopulation> fidIidHapMapPopTableBuilder = ImmutableTable.builder();
      while (hapMapAncestryReader.ready()) {
        String[] line = hapMapAncestryReader.readLine().split(delimiter);
        String fid = line[fidIndex];
        String iid = line[iidIndex];
        HapMapPopulation pop;
        try {
          pop = hapMapCodeMap.get(Integer.parseInt(line[hapIndex]));
        } catch (NumberFormatException nfe) {
          log.reportError("Non-numeric HapMap code: " + line[hapIndex]);
          return null;
        }
        fidIidHapMapPopTableBuilder.put(fid, iid, pop);
      }
      Table<String, String, HapMapPopulation> fidIidHapMapPopTable = fidIidHapMapPopTableBuilder.build();
      if (proj != null) {
        maybeAddHapMapToSampleData(hapMapClassHeader, fidIidHapMapPopTable);
      }
      return fidIidHapMapPopTable;
    } catch (FileNotFoundException e) {
      log.reportError("Cannot find " + hapMapAncestries + " to load HapMap ancestries");
      return null;
    } catch (IOException e) {
      log.reportError("Error reading " + hapMapAncestries + " to load HapMap ancestries");
      return null;
    }
  }

  private void freqsByRace(String resultFile, String outFile) {
    String overallFrqFile = ext.rootOf(resultFile, false) + "_all.frq";
    CmdLine.runDefaults(plinkExe + " --noweb --bfile " + plinkroot + " --freq" + " --out "
                        + ext.rootOf(overallFrqFile, false), dir);
    String[] header = Files.getHeaderOfFile(overallFrqFile, log);
    int key = ext.indexOfStr("SNP", header);
    int[] targets = new int[] {ext.indexOfStr("A1", header), ext.indexOfStr("A2", header),
                               ext.indexOfStr("MAF", header)};
    Hashtable<String, String> overallFreq = HashVec.loadFileToHashString(overallFrqFile, key,
                                                                         targets, "\t", true);

    Map<PCImputeRace.RACE, String> raceListFiles = PCImputeRace.raceListFilenames(resultFile);
    @SuppressWarnings("unchecked")
    Map<PCImputeRace.RACE, Map<String, String>> raceFreqs = Maps.newHashMap();

    PrintWriter writer = Files.getAppropriateWriter(outFile);
    writer.print("SNP\tA1\tA2\tOverall A1F (n=" + PCImputeRace.countFounders(dir + plinkroot, null)
                 + ")");

    for (PCImputeRace.RACE race : PCImputeRace.RACE.values()) {
      String raceListFile = raceListFiles.get(race);
      String raceFrqFile = ext.rootOf(raceListFile, false) + ".frq";
      CmdLine.runDefaults(plinkExe + " --noweb --bfile " + plinkroot + " --keep " + raceListFile
                          + " --freq" + " --out " + ext.rootOf(raceFrqFile, false), dir);
      if (Files.exists(raceFrqFile)) {
        header = Files.getHeaderOfFile(raceFrqFile, log);
        key = ext.indexOfStr("SNP", header);
        targets = new int[] {ext.indexOfStr("A1", header), ext.indexOfStr("A2", header),
                             ext.indexOfStr("MAF", header)};
        raceFreqs.put(race, HashVec.loadFileToHashString(raceFrqFile, key, targets, "\t", true));

        writer.print("\t" + race + " A1F (n="
                     + PCImputeRace.countFounders(dir + plinkroot, raceListFile) + ")");
      } else {
        log.reportTimeWarning("Could not calculate frequencies for " + race.getDescription()
                              + ", this occurs when no project samples were found for this race.");
        raceFreqs.put(race, new HashMap<String, String>());
      }
    }
    writer.println();

    for (Entry<String, String> overallEntry : overallFreq.entrySet()) {
      String marker = overallEntry.getKey();
      String[] data = overallEntry.getValue().split("\t");
      String A1 = data[0];
      String A2 = data[1];
      ;
      String overallMAF;
      try {
        overallMAF = Double.toString(ext.roundToSignificantFigures(Double.parseDouble(data[2]), 4));
      } catch (NumberFormatException nfe) {
        log.reportError("Invalid MAF (" + data[2] + ") for SNP '" + marker + "'");
        overallMAF = ".";
      }

      writer.print(marker + "\t" + A1 + "\t" + A2 + "\t" + overallMAF);

      for (PCImputeRace.RACE race : PCImputeRace.RACE.values()) {
        Map<String, String> raceFreq = raceFreqs.get(race);
        String a1f;
        if (marker == null || !raceFreq.containsKey(marker)) {
          a1f = ".";
        } else {
          String[] raceData = raceFreq.get(marker).split("\t");
          String raceA1 = raceData[0];
          String raceA2 = raceData[1];
          try {
            double raceMaf = Double.parseDouble(raceData[2]);
            if (A1.equals(raceA1) && A2.equals(raceA2)) {
              a1f = Double.toString(ext.roundToSignificantFigures(raceMaf, 4));
            } else if (A1.equals(raceA2) && A2.equals(raceA1)) {
              a1f = Double.toString(ext.roundToSignificantFigures(1.0 - raceMaf, 4));
            } else {
              log.reportError("Alleles for SNP '" + marker + "' and " + raceListFiles.get(race)
                              + " (" + raceA1 + ", " + raceA2 + ") do not match overall alleles ("
                              + A1 + ", " + A2 + " )");
              a1f = ".";
            }
          } catch (NumberFormatException nfe) {
            log.reportError("Invalid MAF (" + raceData[2] + ") for SNP '" + marker + "' and "
                            + raceListFiles.get(race));
            a1f = ".";
          }
        }
        writer.print("\t" + a1f);
      }
      writer.println();
    }
    writer.close();
  }

  public static void main(String[] args) {

    int numArgs = args.length;
    String dir = "./";
    String hapMapPlinkRoot = DEFAULT_HAPMAP_PLINKROOT;
    String putativeWhites = null;
    boolean checkHomo = false;
    boolean run = false;
    boolean runPipeline = false;
    boolean imputeRace = false;
    Project proj = null;
    String dummyProjectPrefix = null;
    String snpRSIDLookupFile = null;
    String logfile = null;
    String plinkExe = CLI.DEF_PLINK_EXE;
    Logger log;

    String usage = "\n" + "gwas.Ancestry requires 3+ arguments\n"
                   + "   (1) Run directory with plink.* files and "
                   + RelationAncestryQc.UNRELATEDS_FILENAME + " (i.e. dir=" + dir + " (default))\n"
                   + "   (2) PLINK root of Unambiguous HapMap Founders (i.e. hapMapPlinkRoot="
                   + hapMapPlinkRoot + " (default))\n" + "   (3) " + CLI.DESC_PLINK_EXE + " (i.e. "
                   + CLI.ARG_PLINK_EXE + "=" + plinkExe + " (default))"
                   + "   (4) Logfile (i.e. log=" + "ancestry.log" + " (default))\n" + "  AND\n"
                   + "   (4) Run full pipeline (i.e. -runPipeline (not the default, requires arguments for each step))\n"
                   + "  OR\n"
                   + "   (5) Check Homogeneity using Chi-Square (Generates PBS script to run Fisher's exact, if desired) (i.e. -checkHomo (not the default))\n"
                   + "   (6) File of FID/IID pairs of putative whites to use for finding homogenous markers by comparison to CEU (i.e. putativeWhites=whites.txt (required, no default))\n"
                   + "  OR\n"
                   + "   (7) Parse homogeneity checks and run Eigenstrat (i.e. -run (not the default))\n"
                   + "  OR\n" + "   (8) Impute race (i.e. -imputeRace (not the default))\n"
                   + "   (9) Project properties file (i.e. proj=example.properties (not the default))\n"
                   + "  OR\n"
                   + "   (10) Prefix for dummy project to allow viewing of results when run without a project (i.e. dummyProject=example (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = ext.parseStringArg(arg, "./");
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("hapMapPlinkRoot=")) {
        hapMapPlinkRoot = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("-runPipeline")) {
        runPipeline = true;
        numArgs--;
      } else if (arg.startsWith("-checkHomo")) {
        checkHomo = true;
        numArgs--;
      } else if (arg.startsWith("putativeWhites=")) {
        putativeWhites = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("-run")) {
        run = true;
        numArgs--;
      } else if (arg.startsWith("-imputeRace")) {
        imputeRace = true;
        numArgs--;
      } else if (arg.startsWith("proj")) {
        proj = new Project(ext.parseStringArg(arg, "./"));
        numArgs--;
      } else if (arg.startsWith("snpLookup")) {
        snpRSIDLookupFile = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith("dummyProject")) {
        dummyProjectPrefix = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith(CLI.ARG_PLINK_EXE)) {
        plinkExe = ext.parseStringArg(arg);
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    if (logfile == null) {
      if (proj == null) {
        log = new Logger(dir + "ancestry.log");
      } else {
        log = proj.getLog();
      }
    } else {
      log = new Logger(logfile);
    }
    final Ancestry ancestry;
    if (proj == null && dummyProjectPrefix != null) {
      ancestry = new Ancestry(dir, dummyProjectPrefix, plinkExe, log);
    } else {
      ancestry = new Ancestry(dir, "plink", proj, plinkExe, log);
    }
    try {
      if (runPipeline && putativeWhites != null) {
        ancestry.runPipeline(putativeWhites, hapMapPlinkRoot, snpRSIDLookupFile);
      } else if (checkHomo && putativeWhites != null) {
        ancestry.checkHomogeneity(putativeWhites, hapMapPlinkRoot);
      } else if (run) {
        String homogeneityDrops = ancestry.parseHomogeneity();
        ancestry.mergeHapMap(hapMapPlinkRoot, homogeneityDrops, snpRSIDLookupFile);
        ancestry.runPCA(DEFAULT_NUM_COMPONENTS_ANCESTRY);
      } else if (imputeRace) {
        ancestry.imputeRaceFromPCA();
      } else {
        System.err.println(usage);
        System.exit(1);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }

  }

}
