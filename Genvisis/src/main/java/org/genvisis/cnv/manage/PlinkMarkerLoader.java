package org.genvisis.cnv.manage;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Set;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.PlinkData.ExportIDScheme;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.GenomicPosition;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.utils.gwas.DosageData;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class PlinkMarkerLoader implements Runnable {

  public static void main(String[] args) {
    String[] markerNames = HashVec.loadFileToStringArray("D:/PlinkGeno/mkrs10000.txt", false, null,
                                                         false);
    String plinkFileRoot = "D:/PlinkGeno/plink";
    (new PlinkMarkerLoader(null, plinkFileRoot, markerNames)).run();
    System.out.println(ext.getTime() + "]\tFinished");
  }

  String fileRoot;
  String[] markerList;
  int[] markerPosInBim;
  volatile byte[][] genotypes;
  volatile boolean[] loaded;
  Logger log;
  HashMap<String, Integer> markerIndicesLookup;
  HashMap<String, Integer> dnaPlinkIndexLookup;

  boolean idListsDiffer;
  boolean initialized;
  boolean killed;
  boolean killComplete;
  Thread thread;

  public PlinkMarkerLoader(Project proj, String plinkFileRoot, String[] markers) {
    fileRoot = plinkFileRoot;
    markerList = markers;
    log = proj == null ? new Logger() : proj.getLog();

    if (markers == null) {
      log.reportError("The list of markers for MarkerDataLoader to load was null");
      killed = true;
      return;
    }

    if (markers.length == 0) {
      log.reportError("The list of markers for MarkerDataLoader to load was empty (n=0)");
      killed = true;
      return;
    }

    markerPosInBim = new int[markerList.length];
    genotypes = new byte[markerList.length][];
    loaded = ArrayUtils.booleanArray(markerList.length, false);

    lookupMarkerPositions(proj);
    lookupIDs(proj);
  }

  private void initiate() {
    initialized = true;
  }

  public boolean isKilled() {
    return killed;
  }

  private void registerThread(Thread thread) {
    this.thread = thread;
  }

  public void kill() {
    killed = true;
    if (!thread.isAlive()) {
      killComplete = true;
    }
  }

  public boolean killComplete() {
    return killComplete;
  }

  public Thread getThread() {
    return thread;
  }

  private void lookupMarkerPositions(Project proj) {
    Set<String> lookFor = Sets.newHashSet(markerList);
    markerIndicesLookup = new HashMap<>();
    HashMap<String, Integer> markerIndicesLookupTemp = new HashMap<>();
    String bim = fileRoot + ".bim";
    try (BufferedReader reader = Files.getAppropriateReader(bim)) {
      String temp;
      int cnt = 0;
      while ((temp = reader.readLine()) != null) {
        String[] line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
        String mkr = line[1];
        if (lookFor.remove(mkr)) {
          markerIndicesLookupTemp.put(mkr, cnt);
        }
        cnt++;
      }
    } catch (FileNotFoundException fnfe) {
      // TODO should KILL here
      log.reportException(fnfe);
    } catch (IOException ioe) {
      // TODO should KILL here
      log.reportException(ioe);
    }

    if (!lookFor.isEmpty()) {
      log.reportTimeWarning(lookFor.size() + " markers were not found in " + bim
                            + ", attempting to search by position");
      Map<String, GenomicPosition> markerPosMap = Maps.transformValues(proj.getMarkerSet()
                                                                           .getMarkerNameMap(),
                                                                       Marker::getGenomicPosition);
      Maps.asMap(lookFor, markerPosMap::get).values();
      Map<GenomicPosition, NavigableSet<Marker>> projMarkerPositions = proj.getMarkerSet()
                                                                           .getGenomicPositionMap();
      try (BufferedReader reader = Files.getAppropriateReader(bim)) {
        String temp;
        int cnt = 0;
        while ((temp = reader.readLine()) != null) {
          String[] line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
          String chr = line[PSF.Plink.BIM_CHR_INDEX];
          String pos = line[PSF.Plink.BIM_POS_INDEX];
          GenomicPosition markerPos = new GenomicPosition(Positions.chromosomeNumber(chr),
                                                          Integer.parseInt(pos));
          Set<Marker> matchedMarkers = projMarkerPositions.get(markerPos);
          if (matchedMarkers != null) {
            final int index = cnt;
            matchedMarkers.stream().map(Marker::getName).filter(lookFor::remove)
                          .forEach(mkr -> markerIndicesLookupTemp.put(mkr, index));
          }
          cnt++;
        }
        if (!lookFor.isEmpty()) log.reportTimeWarning(lookFor.size()
                                                      + " markers were still not found by position in "
                                                      + bim);
      } catch (FileNotFoundException fnfe) {
        // TODO should KILL here
        log.reportException(fnfe);
      } catch (IOException ioe) {
        // TODO should KILL here
        log.reportException(ioe);
      }
    }

    for (int i = 0; i < markerList.length; i++) {
      markerPosInBim[i] = markerIndicesLookupTemp.get(markerList[i]) == null ? -1
                                                                             : markerIndicesLookupTemp.get(markerList[i])
                                                                                                      .intValue();
      markerIndicesLookup.put(markerList[i], i);
    }

  }

  private void lookupIDs(Project proj) {
    BufferedReader reader;
    String[] line;

    String famFile = PSF.Plink.getFAM(fileRoot);

    ExportIDScheme plinkIDScheme = PlinkData.detectExportIDScheme(proj, famFile);
    dnaPlinkIndexLookup = new HashMap<>();

    try {
      reader = Files.getAppropriateReader(famFile);
      String temp = null;
      int cnt = 0;
      int unmatched = 0;
      while ((temp = reader.readLine()) != null) {
        line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
        String dna = plinkIDScheme.getProjDNA(proj, line[PSF.Plink.FAM_FID_INDEX],
                                              line[PSF.Plink.FAM_IID_INDEX]);
        if (dna == null)
          unmatched++;
        else if (dnaPlinkIndexLookup.putIfAbsent(dna, cnt) != null) {
          log.reportError("Duplicate sample in " + famFile + ": " + line[PSF.Plink.FAM_FID_INDEX]
                          + "\t" + line[PSF.Plink.FAM_IID_INDEX]);
        }
        cnt++;
      }
      if (unmatched > 0) log.reportTimeWarning(unmatched + " samples in " + famFile + " (of " + cnt
                                               + " total samples) did not match to a project sample");
      reader.close();
      reader = null;
    } catch (FileNotFoundException fnfe) {
      // TODO should KILL here
      log.reportException(fnfe);
    } catch (IOException ioe) {
      // TODO should KILL here
      log.reportException(ioe);
    }
  }

  @Override
  public void run() {
    RandomAccessFile in = null;

    HashMap<String, byte[]> mkrGenotypes = new HashMap<>();

    if (killed) {
      return;
    }

    initiate();
    int cnt = 0;
    long time = new Date().getTime();

    try {
      in = new RandomAccessFile(fileRoot + ".bed", "r");

      byte[] magicBytes = new byte[3];
      in.read(magicBytes);
      if (magicBytes[2] == 0) {
        log.reportError("Error - .bed file is sample-dominant.");
      } else {
        int famCnt = Files.countLines(fileRoot + ".fam", 0);
        int blockSize = (int) Math.ceil(famCnt / 4.0d);
        for (int i = 0; i < markerList.length; i++) {
          if (markerPosInBim[i] == -1) {
            // missing marker, not present in PLINK files
            mkrGenotypes.put(markerList[i], ArrayUtils.byteArray(famCnt, (byte) -1));
            cnt++;
            continue;
          }
          in.seek(3 + (long) markerPosInBim[i] * blockSize);

          byte[] markerBytes = new byte[blockSize];
          byte[] sampGeno = new byte[famCnt];
          in.read(markerBytes);
          for (int bitInd = 0; bitInd < markerBytes.length; bitInd++) {
            byte bedByte = markerBytes[bitInd];
            byte[] genos = DosageData.decodeBedByte(bedByte);

            for (int g = 0; g < genos.length; g++) {
              int idInd = bitInd * 4 + g;
              if (idInd >= famCnt) {
                break;
              }
              sampGeno[idInd] = genos[g];
            }
          }

          mkrGenotypes.put(markerList[i], sampGeno);
          cnt++;
        }
      }

      in.close();

      for (int i = 0; i < markerList.length; i++) {
        genotypes[i] = mkrGenotypes.get(markerList[i]);
        loaded[i] = true;
      }

    } catch (FileNotFoundException e) {
      log.reportException(e);
      // TODO should KILL here
    } catch (IOException e) {
      log.reportException(e);
      // TODO should KILL here
    } catch (Elision e) {
      log.reportException(e);
      // TODO should KILL here
    } finally {
      if (in != null) {
        try {
          in.close();
        } catch (Exception e) {
          e.printStackTrace();
        }
      }
    }

    int aa, ab, bb, m;
    aa = 0;
    ab = 0;
    bb = 0;
    m = 0;
    for (int i = 0; i < genotypes[0].length; i++) {
      byte g = genotypes[0][i];
      if (g == 0) {
        aa++;
      } else if (g == 1) {
        ab++;
      } else if (g == 2) {
        bb++;
      } else if (g == -1) {
        m++;
      }
    }
    System.out.println("PLINK GENOTYPES  --->  AA: " + aa + " AB: " + ab + " BB: " + bb + " Miss: "
                       + m);

    if (killed) {
      log.report("PlinkMarkerLoader killed");
      markerList = null;
      markerPosInBim = null;
      loaded = null;
      genotypes = null;
      System.gc();
      killComplete = true;
    } else {
      log.report("Independent thread has finished loading " + cnt + " markers in "
                 + ext.getTimeElapsed(time));
    }

  }

  /**
   * @param proj Project to match Sample and Marker IDs from
   * @param marker Marker to get genotype for
   * @param sample DNA of sample to get genotype for
   * @return genotype
   */
  public byte getGenotypeForIndi(Project proj, String marker, String sample) {
    Integer idIndex = dnaPlinkIndexLookup.get(sample);
    if (idIndex == null) {
      return (byte) -1;
    } else {
      int markerIndex = markerIndicesLookup.get(marker) == null ? -1
                                                                : markerIndicesLookup.get(marker);
      if (markerIndex == -1) {
        return (byte) -1;
      }
      while (!loaded[markerIndex]) {
        Thread.yield();
      }
      return genotypes[markerIndex][idIndex];
    }
  }

  public static PlinkMarkerLoader loadPlinkDataFromListInSeparateThread(Project proj,
                                                                        String plinkDirFileRoot,
                                                                        String[] markerList) {
    PlinkMarkerLoader plinkMarkerLoader;
    Thread thread;

    proj.getLog().report("PLINK marker data is loading in an independent thread.");
    plinkMarkerLoader = new PlinkMarkerLoader(proj, plinkDirFileRoot, markerList);
    if (plinkMarkerLoader.isKilled()) {
      return null;
    }
    plinkMarkerLoader.initiate();
    thread = new Thread(plinkMarkerLoader);
    thread.start();
    plinkMarkerLoader.registerThread(thread);

    return plinkMarkerLoader;
  }

}
