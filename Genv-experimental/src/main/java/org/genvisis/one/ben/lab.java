package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.CNVFilter;
import org.genvisis.common.CNVFilter.CNVFilterPass;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.DosageData;
import org.genvisis.filesys.Segment;
import org.genvisis.gwas.MergeExtractPipeline;

public class lab {


  public enum TEST {
                    TEST1, TEST2, TEST3;
  }

  public static void breakCentromeric() throws IOException {
    CNVFilter filter = new CNVFilter(null);
    filter.setBreakupCentromeres(true);
    filter.setCentromereBoundariesFromFile("D:/data/gedi_gwas/data/markers.bim");
    filter.computeCentromereMidPoints();

    CNVariant[] centromeric = CNVariant.loadPlinkFile("D:/SIDS and IQ/IQ/merged.cnv", false);

    PrintWriter writer = new PrintWriter(new FileWriter("D:/SIDS and IQ/IQ/merged_split.cnv"));
    writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));

    for (CNVariant cnv : centromeric) {
      CNVFilterPass fp = filter.getCNVFilterPass(cnv);
      if (fp.isCentromeric()) {
        CNVariant[] broken = filter.breakUpCentromere(fp, cnv);
        for (CNVariant newcnv : broken) {
          writer.println(newcnv.toPlinkFormat());
        }
      } else {
        writer.println(cnv.toPlinkFormat());
      }
    }
    writer.flush();
    writer.close();

  }

  private static void countCNVsForIndividuals(String indivFile, String cnvFile,
                                              String outFile) throws IOException {
    Hashtable<String, String> sampleKeyHash = new Hashtable<String, String>();
    BufferedReader reader = new BufferedReader(new FileReader(indivFile));
    String line = null;
    while ((line = reader.readLine()) != null) {
      String[] tmp = line.split("\t");
      sampleKeyHash.put(tmp[0] + "\t" + tmp[1], "");
    }
    reader.close();

    Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(cnvFile, sampleKeyHash, true, false);

    PrintWriter writer = Files.getAppropriateWriter(outFile);
    for (CNVariant cnv : cnvs) {
      writer.println(cnv.toPlinkFormat());
    }
    writer.flush();
    writer.close();
  }

  public static <T extends Enum<T>> void enumValues(Class<T> enumType) {
    for (T c : enumType.getEnumConstants()) {
      System.out.println(c.name());
    }
  }


  public static void filterCentromeric(String dir, String in, String out,
                                       String markerSetFilenameToBreakUpCentromeres, int build,
                                       Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    CNVariant cnv;
    Segment[] centromereMidpoints;
    int[][] centromereBoundaries;

    centromereBoundaries =
        Positions.determineCentromereBoundariesFromMarkerSet(markerSetFilenameToBreakUpCentromeres,
                                                             build, log);
    centromereMidpoints = Positions.computeCentromereMidpoints(centromereBoundaries);

    try {
      reader = new BufferedReader(new FileReader(dir + in));
      writer = new PrintWriter(new FileWriter(dir + out));
      writer.println(reader.readLine());
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        cnv = new CNVariant(line);
        if (cnv.overlaps(centromereMidpoints[cnv.getChr()])) {
          writer.println(Array.toStr(line));
        }
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + in + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + in + "\"");
      return;
    }
  }

  private static void idSwap(Project proj, String fileIn) throws IOException {
    BufferedReader reader = new BufferedReader(new FileReader(fileIn));
    String outFile = ext.rootOf(fileIn, false) + ".ids";
    PrintWriter writer = new PrintWriter(new FileWriter(outFile));

    SampleData sampleData = proj.getSampleData(0, false);

    while (reader.ready()) {
      String line = reader.readLine();
      String[] keyVal = sampleData.lookup(line);
      writer.println(Array.toStr(keyVal, "\t"));
    }
    writer.flush();
    writer.close();
    reader.close();
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    Project proj;
    String filename = "lab.dat";
    String logfile = null;
    String file = null;

    boolean test = true;
    if (test) {

      String outFile = "/scratch.global/cole0482/test.db.xln.gz";
      String mapOutFile = "/scratch.global/cole0482/mapOut.xln";


      MergeExtractPipeline pipeline = new MergeExtractPipeline();
      // pipeline.setMarkers(markersFile);
      pipeline.setRunDirectory("/scratch.global/cole0482/merge/", true);
      pipeline.setOutputFormat(DosageData.DATABASE_DOSE_FORMAT);
      pipeline.setOutputFiles(outFile, mapOutFile);
      pipeline.setRenameMarkers(true);
      // pipeline.addDataSource("/scratch.global/cole0482/merge/blacks/", "gwas.bed", "gwas.bim",
      // "gwas.fam");
      pipeline.addDataSource("exome", "/scratch.global/cole0482/merge/blacks/", "exome.bed",
                             "exome.bim", "exome.fam");
      pipeline.addDataSource("metab", "/scratch.global/cole0482/merge/blacks/", "metab.bed",
                             "metab.bim", "metab.fam");
      // add more;
      pipeline.run();


      // String doseFile1 =
      // "/home/pankarne/shared/ARIC_Genomics_Data/GWAS_Chip/1000G/ARIC.whites.impute2/chr3.90069244.95069244.impute2.gz";
      // String mapFile1 =
      // "/home/pankarne/shared/ARIC_Genomics_Data/GWAS_Chip/1000G/ARIC.whites.impute2/chr3.90069244.95069244.impute2_info";
      //
      // String doseFile2 =
      // "/home/pankarne/shared/ARIC_Genomics_Data/GWAS_Chip/1000G/ARIC.whites.impute2/chr3.95069244.100069244.impute2.gz";
      // String mapFile2 =
      // "/home/pankarne/shared/ARIC_Genomics_Data/GWAS_Chip/1000G/ARIC.whites.impute2/chr3.95069244.100069244.impute2_info";
      //
      // String idFile = "/home/pankarne/cole0482/EA.indiv.dup";
      // String outFile = "/scratch.global/cole0482/test.db.xln.gz";
      // String mapOutFile = "/scratch.global/cole0482/mapOut.xln";
      //
      // DosageData dd1 = new DosageData(doseFile1, idFile, mapFile1,
      // DosageData.IMPUTE2_DOSE_FORMAT, null, true, null);
      // DosageData dd2 = new DosageData(doseFile2, idFile, mapFile2,
      // DosageData.IMPUTE2_DOSE_FORMAT, null, true, null);
      // DosageData dd3 = DosageData.combine(dd1, dd2);
      // dd1 = null;
      // dd2 = null;
      // dd1 = DosageData.loadPlinkBinary(dir, plinkRoot);
      // dd2 = DosageData.combine(dd3, dd1);
      // dd1 = null;
      // dd3 = null;
      // dd1 = DosageData.loadPlinkBinary(dir2, plinkRoot2);
      // dd3 = DosageData.combine(dd2, dd1);
      // dd3.writeToFile(outFile, mapOutFile, null, DosageData.DATABASE_DOSE_FORMAT, null);
      // System.out.println("complete!");

      // createRandomSelectionFile();
      // try {
      // String checkfile = "D:/data/gedi_gwas/data/cluster.genome.gz";
      // BufferedReader reader = Files.getAppropriateReader(checkfile);
      // String line = reader.readLine();
      // int cnt = 1;
      // while ((line = reader.readLine()) != null) {
      // cnt++;
      // }
      // reader.close();
      // System.out.println("Read " + cnt + " lines");



      // String s = "10000000000";
      // long l = Long.parseLong(s);
      // System.out.println(l > Integer.MAX_VALUE);
      // int d = Integer.parseInt(s);

      //
      // Process p = Runtime.getRuntime().exec("where notepad.exe");
      // int waitCode = p.waitFor();
      // int outCode = p.exitValue();
      // System.out.println("wait: " + waitCode + "| out: " + outCode);
      // processData();
      // concatData();
      // ConditionalAnalysisPipeline.processOnly(args[0], args[1], args[2]);
      // } catch (IOException e) {
      // // TODO Auto-generated catch block
      // e.printStackTrace();
      // } /*catch (InterruptedException e) {
      // TODO Auto-generated catch block
      // e.printStackTrace();
      // }*/
      //
      // String dir = "F:/CARDIA processing/";
      //
      // String[] files = (new File(dir)).list(new FilenameFilter() {
      // @Override
      // public boolean accept(File dir, String name) {
      // return name.endsWith("positions.xln");
      // }
      // });
      //
      // for (String file2 : files) {
      // try {
      // processSNPPositions(dir + file2, dir + ext.rootOf(file2) + ".proc.xln");
      // } catch (IOException e) {
      // TODO Auto-generated catch block
      // e.printStackTrace();
      // }
      // }


      return;
    }


    String usage = "";

    if (numArgs == 0) {
      try {
        countCNVsForIndividuals("D:/data/ny_registry/new_york/stats/puv_ids.txt",
                                "D:/data/ny_registry/new_york/stats/recodedM.cnv",
                                "D:/data/ny_registry/new_york/stats/puv_cnvs.cnv");
        // testClipboard();
        // BufferedReader reader = new BufferedReader(new
        // FileReader("D:/ForestPlot/Hb_SingleSNP.csv"));
        // String line = reader.readLine();
        // do {
        // System.out.println(line);
        // } while((line = reader.readLine()) != null);

        // filterForMarkers("D:/height/scratch/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq_parsed.xln",
        // "D:/height/scratch/samples_logan.bim");
        // writeRegions();
        // stripMarkerNames();
        // writeExclusions();
        // concatFiles();
        // splitFile();
        // idSwap(new Project("D:/projects/gedi_gwas.properties", false),
        // "D:/data/gedi_gwas/overlap_ids.txt");
        // compareMarkers();
        // filterLDFiles(0.5);
        // formatLDResults();
        // filter();
        // breakCentromeric();
        // filterWrong();
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      // mockupGUI();
      return;
    }


    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("file=")) {
        file = arg.split("=")[1];
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
      proj = new Project(filename, logfile, false);
      if (file != null) {
        idSwap(proj, file);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}


