package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.stream.Collectors;
import javax.annotation.Nullable;
import org.genvisis.cnv.analysis.pca.PCImputeRace;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.filesys.SnpMarkerSet;
import com.google.common.collect.ImmutableTable;
import com.google.common.collect.Sets;
import com.google.common.collect.Table;

public class Ancestry {

  public static final String DEFAULT_HAPMAP_PLINKROOT = "/home/pankrat2/shared/bin/HapMap/unambiguousHapMapFounders";
  public static final String RACE_IMPUTATIONAS_FILENAME = "raceImputations.mds";
  public static final String RACE_FREQS_FILENAME = "freqsByRace.xln";

  private static final String HAPMAP_CLASS_NAME = "HapMap";
  private static final String EIGENSTRAT_OUTPUT_NAME = "combo_fancy_postnormed_eigens.xln";
  private static final String EIGENSTRAT_IID_LABEL = "IID";
  private static final String EIGENSTRAT_FID_LABEL = "FID";
  private static final String EIGENSTRAT_PC1_LABEL = "C1";
  private static final String EIGENSTRAT_PC2_LABEL = "C2";

  private static enum HapMapPopulation {
    CEU, YRI, CHB, JPT;
  }

  public static void runPipeline(String dir, String putativeWhitesFile, Project proj, Logger log) {
    runPipeline(dir, putativeWhitesFile, null, proj, log);
  }

  private static Project createDummyProject(String dir, String dummyProjectPrefix,
                                            String putativeWhitesFile, Logger log) {
    String projectName = dummyProjectPrefix + "_AncestryResults";
    String projFilename = MitoPipeline.initGenvisisProject() + projectName
                          + MitoPipeline.PROJECT_EXT;
    Files.write((new Project()).PROJECT_NAME.getName() + "=" + projectName, projFilename);
    Project dummyProject = new Project(projFilename, null, false);
    dummyProject.PROJECT_NAME.setValue(projectName);
    dummyProject.PROJECT_DIRECTORY.setValue(dir);
    dummyProject.saveProperties();

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
      log.reportError("Could not generate Sample Data for dummy project "
                      + dummyProject.getPropertyFilename()
                      + ", ancestry results cannot be visualized");
      return null;
    }
    return dummyProject;
  }

  private static void maybeAddHapMapToSampleData(Project proj,
                                                 SampleData.ClassHeader hapMapClassHeader,
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

  public static void runPipeline(String dir, String putativeWhitesFile, String hapMapPlinkRoot,
                                 Project proj, Logger log) {
    if (hapMapPlinkRoot == null) {
      hapMapPlinkRoot = DEFAULT_HAPMAP_PLINKROOT;
    }
    if (!Files.exists(dir + "homogeneity/" + MergeDatasets.CHI_SQUARE_DROPS_FILENAME)
        && Files.list(dir + "homogeneity/", ".Rout").length == 0) {
      log.report("Running homogeneity checks...");
      checkHomogeneity(dir, putativeWhitesFile, dir + "plink", hapMapPlinkRoot, proj, log);
    }
    String homogeneityDrops = parseHomogeneity(dir, log);
    mergeHapMap(dir, dir + "plink", hapMapPlinkRoot, homogeneityDrops, log);
    runEigenstrat(dir);
    imputeRace(dir, proj, log);
  }

  private static void checkHomogeneity(String dir, String putativeWhitesFile,
                                       String projectPlinkRoot, String hapMapPlinkRoot,
                                       @Nullable Project proj, Logger log) {
    String homoDir = dir + "homogeneity/";
    String homoProjDir = homoDir + ext.removeDirectoryInfo(projectPlinkRoot) + "/";
    String homoHapMapDir = homoDir + ext.removeDirectoryInfo(hapMapPlinkRoot) + "/";
    new File(homoProjDir).mkdirs();
    new File(homoHapMapDir).mkdirs();
    String cleanPutativeWhitesFile = validatePutativeWhites(proj, projectPlinkRoot + PSF.Plink.FAM,
                                                            putativeWhitesFile, log);
    CmdLine.runDefaults("plink2 --bfile " + projectPlinkRoot + " --keep " + cleanPutativeWhitesFile
                        + " --hardy", homoProjDir, log);
    CmdLine.runDefaults("plink2 --bfile " + hapMapPlinkRoot + " --keep "
                        + ext.parseDirectoryOfFile(hapMapPlinkRoot) + "CEUFounders.txt --hardy",
                        homoHapMapDir, log);

    MergeDatasets.checkForHomogeneity(homoDir, null, null, "UNAFF", log);
  }

  private static String validatePutativeWhites(@Nullable Project proj, String projFamFile,
                                               String putativeWhitesFile, Logger log) {
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

  private static String parseHomogeneity(String dir, Logger log) {
    int rOuts = Files.list(dir + "homogeneity/", ".Rout").length;
    if (rOuts == 0) {
      log.report("No Fisher's Exact results found, using Chi Square to choose homogeneous markers");
      return dir + "homogeneity/" + MergeDatasets.CHI_SQUARE_DROPS_FILENAME;
    }
    MergeDatasets.parseHomo(dir + "homogeneity/");
    return dir + "homogeneity/" + MergeDatasets.FISHER_OR_CHI_SQUARE_DROPS_FILENAME;
  }

  private static void mergeHapMap(String dir, String projectPlinkRoot, String hapMapPlinkRoot,
                                  String dropMarkersFile, Logger log) {
    if (!Files.exists(dir + RelationAncestryQc.UNRELATEDS_FILENAME)) {
      log.reportError("Error - need a file called " + RelationAncestryQc.UNRELATEDS_FILENAME
                      + " with FID and IID pairs before we can proceed");
      return;
    }
    if (!Files.exists(dir + "plink.bim_unambiguous.txt")) {
      log.report(ext.getTime() + "]\tGenerating list of unambiguous SNPs");
      new SnpMarkerSet(dir + "plink.bim", true,
                       log).listUnambiguousMarkers(dir + "plink.bim_unambiguous.txt",
                                                   dropMarkersFile, true);
    }

    if (!Files.exists(dir + "unambiguousHapMap.bed")) {
      log.report(ext.getTime() + "]\tExtracting unambiguous SNPs for HapMap founders");
      CmdLine.runDefaults("plink2 --bfile " + hapMapPlinkRoot
                          + " --extract plink.bim_unambiguous.txt --make-bed --out unambiguousHapMap --noweb",
                          dir, log);
    }

    if (!Files.exists(dir + "overlap.txt")) {
      log.report(ext.getTime() + "]\tGenerating list of overlapping SNPs");
      Files.writeIterable(HashVec.loadFileToVec(dir + "unambiguousHapMap.bim", false, new int[] {1},
                                                false),
                          dir + "overlap.txt");
    }

    if (Files.exists(dir + "combo.missnp")) {
      new File(dir + "combo.missnp").delete();
    }

    if (!Files.exists(dir + "unambiguous.bed")) {
      log.report(ext.getTime() + "]\tExtracting overlapping SNPs for study samples");
      CmdLine.runDefaults("plink2 --bfile plink --extract overlap.txt --make-bed --out unambiguous --noweb",
                          dir, log);
    }

    if (!Files.exists(dir + "combo.missnp")) {
      log.report(ext.getTime() + "]\tMerging study data and HapMap data for overlapping SNPs");
      CmdLine.runDefaults("plink --bfile unambiguous --bmerge unambiguousHapMap.bed unambiguousHapMap.bim unambiguousHapMap.fam --out combo --noweb",
                          dir, log);
    }

    if (Files.exists(dir + "combo.missnp")) {
      if (Files.exists(dir + "combo.1.missnp")) {
        new File(dir + "combo.1.missnp").delete();
      }
      log.report(ext.getTime() + "]\tChecking for flipped alleles");
      new File(dir + "combo.missnp").renameTo(new File(dir + "combo.1.missnp"));
      CmdLine.runDefaults("plink2 --bfile unambiguous --flip combo.1.missnp --make-bed --out unambiguousFlipped --noweb",
                          dir, log);
      CmdLine.runDefaults("plink --bfile unambiguousFlipped --bmerge unambiguousHapMap.bed unambiguousHapMap.bim unambiguousHapMap.fam --make-bed --out combo --noweb",
                          dir, log);
      if (Files.exists(dir + "combo.missnp")) {
        if (Files.exists(dir + "combo.2.missnp")) {
          new File(dir + "combo.2.missnp").delete();
        }
        log.report(ext.getTime() + "]\tDropping SNPs that cannot be resolved by flipping alleles");
        new File(dir + "combo.missnp").renameTo(new File(dir + "combo.2.missnp"));
        CmdLine.runDefaults("plink2 --bfile unambiguousFlipped --exclude combo.2.missnp --make-bed --out unambiguousFlippedDropped --noweb",
                            dir, log);
        CmdLine.runDefaults("plink2 --bfile unambiguousHapMap --exclude combo.2.missnp --make-bed --out unambiguousDroppedHapMap --noweb",
                            dir, log);
        CmdLine.runDefaults("plink --bfile unambiguousFlippedDropped --bmerge unambiguousDroppedHapMap.bed unambiguousDroppedHapMap.bim unambiguousDroppedHapMap.fam --make-bed --out combo --noweb",
                            dir, log);
      }
    }

    if (!Files.exists(dir + "finalSNPs.txt")) {
      log.report(ext.getTime() + "]\tWriting final list of SNPs to use");
      Files.writeIterable(HashVec.loadFileToVec(dir + "combo.bim", false, new int[] {1}, false),
                          dir + "finalSNPs.txt");
    }
  }

  private static void runEigenstrat(String dir) {
    Logger log;

    dir = ext.verifyDirFormat(dir);
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
      CmdLine.runDefaults("plink2 --bfile ../combo --keep " + RelationAncestryQc.UNRELATEDS_FILENAME
                          + " --make-bed --noweb", unrelatedsDir, log);
    }

    if (!Files.exists(unrelatedsDir + "master")) {
      log.report(ext.getTime() + "]\tCreating Eigenstrat");
      CmdLine.runDefaults(Files.getRunString() + " gwas.Eigenstrat source=plink -create",
                          unrelatedsDir, log);
    }

    if (!Files.exists(unrelatedsDir + "plink.pca.evec")) {
      log.report(ext.getTime() + "]\tRunning master");
      CmdLine.runDefaults("./master", unrelatedsDir, log);
      CmdLine.runDefaults(Files.getRunString() + " gwas.Eigenstrat convert=plink.pca.evec",
                          unrelatedsDir, log);
    }

    if (!Files.exists(dir + "convertf.par")) {
      log.report(ext.getTime() + "]\tRunning convertf");
      CmdLine.runDefaults(Files.getRunString() + " gwas.Eigenstrat source=combo -create", dir, log);
      CmdLine.runDefaults("convertf -p convertf.par", dir, log);
    }

    if (!Files.exists(dir + EIGENSTRAT_OUTPUT_NAME)) {
      CmdLine.runDefaults("plink2 --bfile unrelateds/plink --freq --out unrelateds/plink --noweb",
                          dir, log);
      CmdLine.runDefaults(Files.getRunString()
                          + " gwas.Eigenstrat source=unrelateds/plink target=combo -parse -eigenFormat",
                          dir, log);
    }
  }

  private static void imputeRace(String dir, @Nullable Project proj, Logger log) {
    if (!Files.exists(dir + RACE_IMPUTATIONAS_FILENAME)) {
      Table<String, String, HapMapPopulation> fidIidHapMapPopTable = parseHapMapAncestries(proj,
                                                                                           log);
      if (fidIidHapMapPopTable == null) return;

      Set<PCImputeRace.Sample> samples = Sets.newHashSet();
      Set<PCImputeRace.Sample> europeans = Sets.newHashSet();
      Set<PCImputeRace.Sample> africans = Sets.newHashSet();
      Set<PCImputeRace.Sample> asians = Sets.newHashSet();
      try (BufferedReader eigenReader = Files.getAppropriateReader(dir + EIGENSTRAT_OUTPUT_NAME)) {
        String headerLine = eigenReader.readLine();
        String delim = ext.determineDelimiter(headerLine);
        String[] header = headerLine.split("\t");
        Map<String, Integer> eigenstratIndices = ext.indexMap(new String[] {EIGENSTRAT_FID_LABEL,
                                                                            EIGENSTRAT_IID_LABEL,
                                                                            EIGENSTRAT_PC1_LABEL,
                                                                            EIGENSTRAT_PC2_LABEL},
                                                              header, true, false);

        while (eigenReader.ready()) {
          String[] line = eigenReader.readLine().split(delim);
          String fid = line[eigenstratIndices.get(EIGENSTRAT_FID_LABEL)];
          String iid = line[eigenstratIndices.get(EIGENSTRAT_IID_LABEL)];
          double pc1;
          double pc2;
          try {
            pc1 = Double.parseDouble(line[eigenstratIndices.get(EIGENSTRAT_PC1_LABEL)]);
          } catch (NumberFormatException nfe) {
            pc1 = Double.NaN;
          }
          try {
            pc2 = Double.parseDouble(line[eigenstratIndices.get(EIGENSTRAT_PC2_LABEL)]);
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
        log.reportFileNotFound(dir + EIGENSTRAT_OUTPUT_NAME);
      } catch (IOException e) {
        log.reportIOException(dir + EIGENSTRAT_OUTPUT_NAME);
      }

      PCImputeRace pcir = new PCImputeRace(proj, samples, europeans, africans, asians, log);
      pcir.correctPCsToRace(dir + RACE_IMPUTATIONAS_FILENAME);
    } else {
      log.reportTimeWarning("Skipping imputation - output already exists: "
                            + (dir + RACE_IMPUTATIONAS_FILENAME));
    }

    if (!Files.exists(dir + RACE_FREQS_FILENAME)) {
      PCImputeRace.freqsByRace(dir + RACE_IMPUTATIONAS_FILENAME, dir + "plink",
                               dir + RACE_FREQS_FILENAME, log);
    } else {
      log.reportTimeWarning("Skipping race freq calculation - output already exists: "
                            + (dir + RACE_FREQS_FILENAME));
    }

  }

  private static Table<String, String, HapMapPopulation> parseHapMapAncestries(@Nullable Project proj,
                                                                               Logger log) {
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
        maybeAddHapMapToSampleData(proj, hapMapClassHeader, fidIidHapMapPopTable);
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
    String logfile = null;
    Logger log;

    String usage = "\n" + "gwas.Ancestry requires 3+ arguments\n"
                   + "   (1) Run directory with plink.* files and "
                   + RelationAncestryQc.UNRELATEDS_FILENAME + " (i.e. dir=" + dir + " (default))\n"
                   + "   (2) PLINK root of Unambiguous HapMap Founders (i.e. hapMapPlinkRoot="
                   + hapMapPlinkRoot + " (default))\n" + "   (3) Logfile (i.e. log="
                   + "ancestry.log" + " (default))\n" + "  AND\n"
                   + "   (4) Run full pipeline (i.e. -runPipeline (not the default, requires arguments for each step))\n"
                   + "  OR\n"
                   + "   (5) Check Homogeneity using Chi-Square (Generates PBS script to run Fisher's exact, if desired) (i.e. -checkHomo (not the default))\n"
                   + "   (6) File of FID/IID pairs of putative whites to use for finding homogenous markers by comparison to CEU (i.e. putativeWhites=whites.txt (not the default))\n"
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
      } else if (arg.startsWith("dummyProject")) {
        dummyProjectPrefix = ext.parseStringArg(arg, null);
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
    if (proj == null && dummyProjectPrefix != null) {
      proj = createDummyProject(ext.verifyDirFormat(new File(dir).getAbsolutePath()),
                                dummyProjectPrefix, putativeWhites, log);
    }
    try {
      if (runPipeline) {
        runPipeline(dir, putativeWhites, hapMapPlinkRoot, proj, log);
      } else if (checkHomo) {
        checkHomogeneity(dir, putativeWhites, dir + "plink", hapMapPlinkRoot, proj, log);
      } else if (run) {
        String homogeneityDrops = parseHomogeneity(dir, log);
        mergeHapMap(dir, dir + "plink", hapMapPlinkRoot, homogeneityDrops, log);
        runEigenstrat(dir);
      } else if (imputeRace) {
        imputeRace(dir, proj, log);
      } else {
        System.err.println(usage);
        System.exit(1);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }

  }

}
