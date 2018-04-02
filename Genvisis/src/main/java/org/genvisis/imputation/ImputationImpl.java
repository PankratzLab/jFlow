package org.genvisis.imputation;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.CHROMOSOME;
import org.genvisis.cnv.manage.Resources.Resource;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.qsub.Qsub;

public interface ImputationImpl {

  class MiniMac {

    Logger log;
    String minimacPath;
    String outDir;
    int threads;

    HashMap<Integer, String> refMap = new HashMap<>();
    HashMap<Integer, String> hapsMap = new HashMap<>();

    public MiniMac(Project proj, int[] chrs, String hapsDir, String outDir) {
      this.log = proj.getLog();
      this.outDir = outDir;
      new File(outDir).mkdirs();
      this.threads = proj.NUM_THREADS.getValue();
      Resource rsc = Resources.miniMac(log).getMiniMac3();
      if (rsc == null) {
        log.reportError("Minimac not available!  Error fetching Minimac resource.");
        return;
      }
      if (!rsc.isAvailable(true)) {
        return;
      }
      minimacPath = rsc.getAbsolute();

      // set without knowing if files exist:
      for (int c : chrs) {
        hapsMap.put(c, hapsDir + ShapeIt.outBase + c + ".vcf");
      }
      // OR:
      // create based on existing output files:
      // String[] hapsFiles = (new File(hapsDir)).list(new FilenameFilter() {
      // @Override
      // public boolean accept(File dir, String name) {
      // return name.endsWith(".haps");
      // }
      // });
      // for (String h : hapsFiles) {
      // String chrStr = h.substring(ShapeIt.outBase.length(), h.lastIndexOf("."));
      // int chr = -1;
      // try {
      // chr = Integer.parseInt(chrStr);
      // hapsMap.put(chr, hapsDir + h);
      // } catch (NumberFormatException e) {
      // log.reportError("Couldn't parse chromosome number of haps file: " + h);
      // }
      // }

      for (Integer c : hapsMap.keySet()) {
        Resource refPanel = new Resources.Chr(proj.GENOME_BUILD_VERSION.getValue(),
                                              CHROMOSOME.valueOf("C" + c), log)
                                                                               .getG1Kphase3v5RefPanel();
        if (refPanel == null) {
          log.reportError("Reference panel not found for chr " + c + "!");
          continue;
        }
        refMap.put(c, refPanel.getAbsolute());
      }

    }

    public void createScripts() {
      String scriptFile = outDir + RUN;
      PrintWriter writer = Files.getAppropriateWriter(scriptFile);

      Set<Integer> chrs = new HashSet<>(hapsMap.keySet());

      for (Integer chr : chrs) {
        String shapeItCmd = TEMPLATE.replace(MIN, minimacPath).replace(REF, refMap.get(chr))
                                    .replace(HAP, hapsMap.get(chr)).replace(THD, threads + "")
                                    .replace(OUT, outBase + chr.intValue());
        writer.println(shapeItCmd);
      }

      writer.flush();
      writer.close();

      Files.chmod(scriptFile);
      String commands = new StringBuilder().append("cd ")
                                           .append(ext.parseDirectoryOfFile(scriptFile))
                                           .append("\n").append("./").append(RUN).toString();
      Qsub.qsub(ext.rootOf(new File(scriptFile).getAbsolutePath(), false) + ".qsub", commands,
                Files.PBS_MEM, Files.PBS_PROC, 1);
    }

    String outBase = "mini_chr";

    String RUN = "runMinimac.sh";
    String MIN = "$minimac";
    String REF = "$ref";
    String HAP = "$haps";
    String OUT = "$out";
    String THD = "$threads";
    String TEMPLATE = MIN + " --refHaps " + REF + " --haps " + HAP + " --prefix " + OUT + " --cpus "
                      + THD;

  }

  class ShapeIt {

    Logger log;
    String shapeItPath;
    int threads;
    String outDir;
    HashMap<Integer, String> mapMap = new HashMap<>();
    HashMap<Integer, String> plinkFileMap = new HashMap<>();
    HashMap<Integer, Boolean> plinkFileTypeMap = new HashMap<>();

    public ShapeIt(Project proj, int[] chrs, String plinkFileDir, String plinkChrFilePrefix,
                   String outDir) {
      this.log = proj.getLog();
      this.threads = proj.NUM_THREADS.getValue();
      this.outDir = outDir;
      new File(outDir).mkdirs();
      Resource rsc = Resources.shapeit(log).getShapeit();
      if (rsc == null) {
        log.reportError("ShapeIt not available!  Error fetching ShapeIt resource.");
        return;
      }
      if (!rsc.isAvailable(true)) {
        return;
      }
      shapeItPath = rsc.getAbsolute();

      for (int i = 0; i < chrs.length; i++) {
        Resource chrMap = new Resources.Chr(proj.GENOME_BUILD_VERSION.getValue(),
                                            CHROMOSOME.valueOf("C" + chrs[i]), log).getGeneticMap();
        if (chrMap == null) {
          log.reportError("Genome map not found for chr " + chrs[i] + "!");
          continue;
        }
        mapMap.put(chrs[i], chrMap.getAbsolute());
      }

      for (int i = 0; i < chrs.length; i++) {
        boolean bedSet = PSF.Plink.allFilesExist(plinkFileDir + plinkChrFilePrefix + chrs[i], true);
        boolean pedSet = PSF.Plink.allFilesExist(plinkFileDir + plinkChrFilePrefix + chrs[i], true);
        if (bedSet || pedSet) {
          plinkFileMap.put(chrs[i], plinkFileDir + plinkChrFilePrefix + chrs[i]);
          plinkFileTypeMap.put(chrs[i], bedSet);
        }
      }
    }

    public void createScripts() {
      String scriptFile = outDir + RUN;
      PrintWriter writer = Files.getAppropriateWriter(scriptFile);

      Set<Integer> chrs = new HashSet<>();
      chrs.addAll(mapMap.keySet());
      chrs.retainAll(plinkFileMap.keySet());

      for (Integer chr : chrs) {
        String shapeItCmd = TEMPLATE.replace(SHP, shapeItPath)
                                    .replace(TYP, plinkFileTypeMap.get(chr) ? T_B : T_P)
                                    .replace(PLK, plinkFileMap.get(chr))
                                    .replace(MAP, mapMap.get(chr)).replace(THD, threads + "")
                                    .replace(OUT, outBase + chr.intValue());
        writer.println(shapeItCmd);
      }

      writer.flush();
      writer.close();

      Files.chmod(scriptFile);
      String commands = new StringBuilder().append("cd ")
                                           .append(ext.parseDirectoryOfFile(scriptFile))
                                           .append("\n").append("./").append(RUN).toString();
      Qsub.qsub(ext.rootOf(new File(scriptFile).getAbsolutePath(), false) + ".qsub", commands,
                Files.PBS_MEM, Files.PBS_PROC, 1);

      createConvertScripts();
    }

    private void createConvertScripts() {
      String scriptFile = outDir + CON;
      PrintWriter writer = Files.getAppropriateWriter(scriptFile);

      Set<Integer> chrs = new HashSet<>();
      chrs.addAll(mapMap.keySet());
      chrs.retainAll(plinkFileMap.keySet());

      for (Integer chr : chrs) {
        String hapFile = outBase + chr;
        String outFile = outBase + chr + ".vcf";
        String convertCmd = CONVERT.replace(SHP, shapeItPath).replace(HAP, hapFile)
                                   .replace(OUT, outFile).replace(THD, threads + "");
        writer.println(convertCmd);
      }

      writer.flush();
      writer.close();

      Files.chmod(scriptFile);
      String commands = new StringBuilder().append("cd ")
                                           .append(ext.parseDirectoryOfFile(scriptFile))
                                           .append("\n").append("./").append(RUN).toString();
      Qsub.qsub(ext.rootOf(new File(scriptFile).getAbsolutePath(), false) + ".qsub", commands,
                Files.PBS_MEM, Files.PBS_PROC, 1);
    }

    static final String RUN = "runShapeIt.sh";
    static final String CON = "runShapeItConvert.sh";
    static final String SHP = "$shape";
    static final String PLK = "$plink";
    static final String T_B = "-B";
    static final String T_P = "-P";
    static final String TYP = "$type";
    static final String MAP = "$map";
    static final String THD = "$threads";
    static final String OUT = "$out";
    static final String HAP = "$haps";
    static final String outBase = "out_chr";

    static final String TEMPLATE = SHP + " " + TYP + " " + PLK + " -M " + MAP + " -T " + THD
                                   + " -O " + OUT + " --duohmm";
    static final String CONVERT = SHP + " -convert --input-haps " + "$haps" + " --output-vcf " + OUT
                                  + " -T " + THD;
  }

}
