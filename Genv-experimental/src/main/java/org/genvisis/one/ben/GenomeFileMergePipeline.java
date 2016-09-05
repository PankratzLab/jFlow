package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.io.UnsupportedEncodingException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.Map.Entry;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.gwas.Qc;

import scala.math.BigInt;

import com.google.common.collect.Lists;

public class GenomeFileMergePipeline {
  
  private static final String[] GENOME_COLUMNS = {"FID1", "IID1", "FID2", "IID2", "Z0", "Z1", "Z2", "PI_HAT"};
  
  private HashMap<String, String> idMap = new HashMap<String, String>();
  private ArrayList<Project> plinkProjects = new ArrayList<Project>();
  private ArrayList<Project> qcProjects = new ArrayList<Project>();
  private ArrayList<Project> projects = new ArrayList<Project>();
  private ArrayList<GenomeFile> files = new ArrayList<GenomeFileMergePipeline.GenomeFile>();

  private boolean runPlinkOrQCIfMissing = false;
  private String outputFile = "audit.xln";
  
  private Logger log = new Logger();
  
  private class CompactCharSequence implements CharSequence, Serializable {
    // from http://www.javamex.com/tutorials/memory/ascii_charsequence.shtml
    static final long serialVersionUID = 1L;

    private static final String ENCODING = "ISO-8859-1";
    private final int offset;
    private final int end;
    private final byte[] data;

    public CompactCharSequence(String str) {
      try {
        data = str.getBytes(ENCODING);
        offset = 0;
        end = data.length;
      } catch (UnsupportedEncodingException e) {
        throw new RuntimeException("Unexpected: " + ENCODING + " not supported!");
      }
    }

    private CompactCharSequence(byte[] data, int offset, int end) {
      this.data = data;
      this.offset = offset;
      this.end = end;
    }
    
    public char charAt(int index) {
      int ix = index+offset;
      if (ix >= end) {
        throw new StringIndexOutOfBoundsException("Invalid index " +
          index + " length " + length());
      }
      return (char) (data[ix] & 0xff);
    }

    public int length() {
      return end - offset;
    }
    
    public CharSequence subSequence(int start, int end) {
      if (start < 0 || end >= (this.end-offset)) {
        throw new IllegalArgumentException("Illegal range " +
          start + "-" + end + " for sequence of length " + length());
      }
      return new CompactCharSequence(data, start + offset, end + offset);
    }
    
    public String toString() {
      try {
        return new String(data, offset, end-offset, ENCODING);
      } catch (UnsupportedEncodingException e) {
        throw new RuntimeException("Unexpected: " + ENCODING + " not supported");
      }
    }
    
  }
  
  private static class GenomeFile {
    public GenomeFile(String name2, String file, Logger log) {
      this.name = name2;
      this.genomeFile = file;
      this.lineCount = Files.countLines(genomeFile, 1);
      log.report("Counted " + this.lineCount + " lines in " + genomeFile);
    }
    String name;
    String genomeFile;
    int lineCount;
  }

  public void loadIDLookupFile(String file, boolean ignoreFirstLine) {
    String[][] ids = HashVec.loadFileToStringMatrix(file, ignoreFirstLine, new int[] {0, 1}, false);
    for (int i = 0; i < ids.length; i++) {
      if (idMap.containsKey(ids[i][0])) {
        log.report("Warning - duplicate FID found on line " + i + ": " + ids[i][0]);
      }
      if (idMap.containsKey(ids[i][1])) {
        log.report("Warning - duplicate IID found on line " + i + ": " + ids[i][1]);
      }
      idMap.put(ids[i][0], ids[i][1]);
      idMap.put(ids[i][1], ids[i][0]);
    }
  }
  
  public void setRunPlinkOrQCIfMissing(boolean run) {
    this.runPlinkOrQCIfMissing = run;
  }
  
  public void addProject(String propertiesFile) {
    Project proj = new Project(propertiesFile, false);
    String plinkDir = proj.PROJECT_DIRECTORY.getValue() + "plink/";
    String genDir = plinkDir + Qc.GENOME_DIR;
    String genFile = genDir = "plink.genome";
    boolean plink = false;
    boolean qc = false;
    String msg;
    if (!Files.exists(plinkDir)) {
      msg = "no \"plink/\" directory in project directory " + proj.PROJECT_DIRECTORY.getValue();
      if (runPlinkOrQCIfMissing) {
        log.reportTime("Warning - " + msg);
        log.reportTime("PLINK files will be created and QC'd for project " + proj.PROJECT_NAME.getValue());
        plink = true;
      } else {
        log.reportError("Error - " + msg);
        return;
      }
    } else if (!Files.exists(genDir)) {
      msg = "no \"genome/\" directory in QC folders for PLINK data; looking for " + genDir;
      if (runPlinkOrQCIfMissing) {
        log.reportTime("Warning - " + msg);
        log.reportTime("PLINK files will be QC'd for project " + proj.PROJECT_NAME.getValue());
        qc = true;
      } else {
        log.reportError("Error - " + msg);
        return;
      }
    } else if (!Files.exists(genFile)) {
      msg = "no plink.genome file found for project " + proj.PROJECT_NAME.getValue();
      if (runPlinkOrQCIfMissing) {
        log.reportTime("Warning - " + msg);
        log.reportTime("PLINK files will be QC'd for project " + proj.PROJECT_NAME.getValue() + "; however, the output folders were found but the plink.genome file was missing, so other errors may be present.");
        qc = true;
      } else {
        log.reportError("Error - " + msg);
        return;
      }
    }
    this.projects.add(proj);
    if (plink) {
      plinkProjects.add(proj);
    }
    if (qc) {
      qcProjects.add(proj);
    }
  }
  
  public void addGenomeFile(String name, String genomeFile) {
    if (!Files.exists(genomeFile)) {
      log.reportTimeError("Error - genome file \"" + genomeFile + "\" not found!");
      return;
    }
    files.add(new GenomeFile(name, genomeFile, log));
  }
  
  public void setOutputFile(String outFile) {
    this.outputFile = outFile;
  }
  
  private static CharSequence[] charSequenceArray(int size, String initValue) {
    CharSequence[] array = new CharSequence[size];
    for (int i = 0; i < size; i++) {
      array[i] = initValue;
    }
    return array;
  }
  
  public void run() {
    BufferedReader reader;
    String line, outKey1, outKey2;
    String fid1, fid2, iid1, iid2;
    CompactCharSequence ibd0, ibd1, ibd2, piHt;
    String genFile;
    String[] parts;
    CharSequence[] outLine;
    int outLineCount, lineCount, count;
    int projInd0, projInd1, projInd2, projInd3;
    int[] factors;
    Project proj;
    HashMap<String, CharSequence[]> outLineMap;
    StringTokenizer st;
    ArrayList<String> pts = new ArrayList<String>(500);
    PrintWriter writer;
    
    // for all projects in plinkProjects, run PLINK export, add all to qcProjects
    
    // for all projects in qcProjects, run gwas.Qc
    
    projectsLoop : for (int p = projects.size() - 1; p >= 0; p++) {
      proj = projects.get(p);
      genFile = proj.PROJECT_DIRECTORY.getValue() + "plink/" + Qc.GENOME_DIR + "plink.genome";
      if (!Files.exists(genFile)) {
        log.reportTimeError("Error - plink.genome file missing for project " + proj.PROJECT_NAME.getValue() + "; data from project will not be included.");
        projects.remove(proj);
      } else {
        try {
          reader = Files.getAppropriateReader(genFile);
          line = reader.readLine();
          reader.close();
          factors = ext.indexFactors(GENOME_COLUMNS, line.trim().split("[\\s]+"), false, false);
          for (int i = 0; i < factors.length; i++) {
            if (factors[i] == -1) {
              log.reportTimeError("Error - column " + GENOME_COLUMNS[i] + " was missing from plink.genome file: " + genFile + " ; data from this file will not be included in final output.");
              projects.remove(proj);
              continue projectsLoop;
            }
          }
        } catch (IOException e) {
          log.reportException(e);
          projects.remove(proj);
        }
      }
    }
    
    for (int p = 0; p < projects.size(); p++) {
      String name = projects.get(p).PROJECT_NAME.getValue();
      String file = projects.get(p).PROJECT_DIRECTORY.getValue() + "plink/" + Qc.GENOME_DIR + "plink.genome";
      files.add(new GenomeFile(name, file, log));
    }
    
    System.gc();

    BigInteger max = new BigInteger("0");
    for (GenomeFile f : files) {
      max.add(new BigInteger("" + f.lineCount));
    }
    outLineCount = 4 + (4 * files.size());  // fid1 iid1 fid2 iid2 + 4 columns per projects (ibd0, ibd1, ibd2, pi_hat)
    outLineMap = new HashMap<String, CharSequence[]>(max.intValue(), 0.95f);
    
    try {
      for (int p = 0; p < files.size(); p++) {
        log.reportTime("Loading data from " + files.get(p).genomeFile);
        projInd0 = 4 + p * 4 + 0;
        projInd1 = projInd0 + 1;
        projInd2 = projInd1 + 1;
        projInd3 = projInd2 + 1;
        
        count = files.get(p).lineCount / 10;
        lineCount = 0;
        reader = Files.getAppropriateReader(files.get(p).genomeFile);
        line = reader.readLine();
        factors = ext.indexFactors(GENOME_COLUMNS, line.trim().split("[\\s]+"), false, false);
        // do id lookup, determine which column of 1/2 and 3/4 (if not both) are in lookup
        while((line = reader.readLine()) != null) {
          line = line.trim();
          if ("".equals(line)) continue;
          st = new StringTokenizer(line, " ", false);
          pts.clear();
          while(st.hasMoreTokens()) {
            pts.add(st.nextToken());
          }
          st = null;
          fid1 = pts.get(factors[0]).toString();    
          iid1 = pts.get(factors[1]).toString();  
          fid2 = pts.get(factors[2]).toString();  
          iid2 = pts.get(factors[3]).toString();  
          ibd0 = new CompactCharSequence(pts.get(factors[4]).toString().intern());  // intern because many of these are either "1" or "0"
          ibd1 = new CompactCharSequence(pts.get(factors[5]).toString().intern());  
          ibd2 = new CompactCharSequence(pts.get(factors[6]).toString().intern());  
          piHt = new CompactCharSequence(pts.get(factors[7]).toString().intern());

          if (idMap.containsKey(iid1) && !idMap.containsKey(fid1)) {
            fid1 = idMap.get(iid1);
          } else if (idMap.containsKey(fid1) && !idMap.containsKey(iid1)) {
            iid1 = idMap.get(fid1);
          }
          if (idMap.containsKey(iid2) && !idMap.containsKey(fid2)) {
            fid2 = idMap.get(iid2);
          } else if (idMap.containsKey(fid2) && !idMap.containsKey(iid2)) {
            iid2 = idMap.get(fid2);
          }
          
          outKey1 = fid1 + "|" + iid1 + "||" + fid2 + "|" + iid2;
          outKey2 = fid2 + "|" + iid2 + "||" + fid1 + "|" + iid1;
          
          if (outLineMap.containsKey(outKey1)) {
            outLine = outLineMap.get(outKey1);
          } else if (outLineMap.containsKey(outKey2)) {
            outLine = outLineMap.get(outKey2);
            // TODO will IBD/PIHAT need to be altered due to flipped ids?
          } else {
            outLine = charSequenceArray(outLineCount, ".");
            outLine[0] = fid1;
            outLine[1] = iid1;
            outLine[2] = fid2;
            outLine[3] = iid2;
            outLineMap.put(outKey1, outLine);
          }
          
          outLine[projInd0] = ibd0;
          outLine[projInd1] = ibd1;
          outLine[projInd2] = ibd2;
          outLine[projInd3] = piHt;
          
          lineCount++;
          if (lineCount % count == 0) {
            System.gc();
          }
        }
        System.gc();
        reader.close();
      }
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    
    // write outmap (options for exclusions?  split pops, ea/aa, etc?)
    log.reportTime("Writing output to " + outputFile);
    
    writer = Files.getAppropriateWriter(outputFile);
    
    outLine = charSequenceArray(outLineCount, "");
    for (int p = 0; p < files.size(); p++) {
      outLine[4 + (4 * p)] = files.get(p).name;
    }
    writer.println(Array.toStr(outLine, "\t"));
    outLine = charSequenceArray(outLineCount, "");
    outLine[0] = "FID1";
    outLine[1] = "IID1";
    outLine[2] = "FID2";
    outLine[3] = "IID2";
    for (int i = 4; i < outLine.length; i += 4) {
      outLine[i] = "P(IBD=0)";
      outLine[i+1] = "P(IBD=1)";
      outLine[i+2] = "P(IBD=2)";
      outLine[i+3] = "PI_HAT";
    }
    writer.println(Array.toStr(outLine, "\t"));
    
    for (Entry<String, CharSequence[]> lines : outLineMap.entrySet()) {
      writer.println(Array.toStr(lines.getValue(), "\t"));
    }
    
    writer.flush();
    writer.close();
    
  }
  
  public static void main(String[] args) {
    
    GenomeFileMergePipeline gfmp = new GenomeFileMergePipeline();
    gfmp.setRunPlinkOrQCIfMissing(false);
    gfmp.loadIDLookupFile("/scratch.global/cole0482/genomeFiles/ids.txt", false);
    gfmp.addGenomeFile("Exome_EA", "/scratch.global/cole0482/genomeFiles/exome_EA_plink.genome");
    gfmp.addGenomeFile("IBC_EA", "/scratch.global/cole0482/genomeFiles/IBC_EA_plink.genome");
    gfmp.addGenomeFile("WES", "/scratch.global/cole0482/genomeFiles/wes_plink.genome");
    gfmp.addGenomeFile("Exome_AA", "/scratch.global/cole0482/genomeFiles/exome_AA_plink.genome");
    gfmp.addGenomeFile("IBC_AA", "/scratch.global/cole0482/genomeFiles/IBC_AA_plink.genome");
    gfmp.setOutputFile("/scratch.global/cole0482/genomeFiles/tempAudit.xln");
    gfmp.run();
  
  }
  
  
}
