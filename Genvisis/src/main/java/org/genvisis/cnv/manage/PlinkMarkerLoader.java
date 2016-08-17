package org.genvisis.cnv.manage;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeMap;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class PlinkMarkerLoader implements Runnable {

  public static PlinkMarkerLoader loadPlinkDataFromListInSeparateThread(Project proj,
                                                                        String plinkDirFileRoot,
                                                                        String[] markerList,
                                                                        String[] sampleIDList) {
    PlinkMarkerLoader plinkMarkerLoader;
    Thread thread;
    // int amountToLoadAtOnceInMB;

    proj.getLog().report("PLINK marker data is loading in an independent thread.");
    // amountToLoadAtOnceInMB = proj.MAX_MEMORY_USED_TO_LOAD_MARKER_DATA.getValue();
    plinkMarkerLoader = new PlinkMarkerLoader(proj, plinkDirFileRoot, markerList, sampleIDList);// ,
                                                                                                // amountToLoadAtOnceInMB);
    if (plinkMarkerLoader.isKilled()) {
      return null;
    }
    plinkMarkerLoader.initiate();
    thread = new Thread(plinkMarkerLoader);
    thread.start();
    plinkMarkerLoader.registerThread(thread);

    return plinkMarkerLoader;
  }

  public static void main(String[] args) {
    String[] markerNames =
        HashVec.loadFileToStringArray("D:/PlinkGeno/mkrs10000.txt", false, null, false);
    String plinkFileRoot = "D:/PlinkGeno/plink";
    (new PlinkMarkerLoader(null, plinkFileRoot, markerNames, null)).run();
    System.out.println(ext.getTime() + "]\tFinished");
  }

  public static String[] sortMarkers(String[] markers, int[] positions) {
    final ArrayList<String> missing = new ArrayList<String>();
    final HashMap<String, Integer> unsortedMap = new HashMap<String, Integer>();
    for (int i = 0; i < markers.length; i++) {
      if (positions[i] == -1) {
        missing.add(markers[i]);
      } else {
        unsortedMap.put(markers[i], positions[i]);
      }
    }
    TreeMap<String, Integer> sortedMap = new TreeMap<String, Integer>(new Comparator<String>() {
      @Override
      public int compare(String o1, String o2) {
        return unsortedMap.get(o1).compareTo(unsortedMap.get(o2));
      }
    });
    sortedMap.putAll(unsortedMap);
    String[] sortedMarkers = new String[sortedMap.size()];
    Set<String> sortedKeys = sortedMap.keySet();
    sortedMarkers = sortedKeys.toArray(sortedMarkers);
    String[] allMarkers = new String[missing.size() + sortedMarkers.length];
    for (int i = 0; i < missing.size(); i++) {
      allMarkers[i] = missing.get(i);
    }
    for (int i = missing.size(); i < allMarkers.length; i++) {
      allMarkers[i] = sortedMarkers[i - missing.size()];
    }
    return allMarkers;
  }

  String fileRoot;
  String[] markerList;
  String[] famIDList;
  String[] idList;
  int[] markerPositions;
  volatile byte[][] genotypes;
  volatile boolean[] loaded;
  Logger log;

  HashMap<String, Integer> markerIndicesLookup;
  HashMap<String, Integer> idIndicesLookup;
  boolean idListsDiffer;
  boolean initialized;
  boolean killed;

  boolean killComplete;

  Thread thread;

  public PlinkMarkerLoader(Project proj, String plinkFileRoot, String[] markers,
                           String[] sampleIDsOrder) {
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

    markerPositions = new int[markerList.length];
    genotypes = new byte[markerList.length][];
    loaded = Array.booleanArray(markerList.length, false);

    lookupMarkerPositions();
    lookupIDs();
    idList = sampleIDsOrder == null ? famIDList : sampleIDsOrder;
    idIndicesLookup = new HashMap<String, Integer>();
    for (int i = 0; i < famIDList.length; i++) {
      idIndicesLookup.put(famIDList[i], i);
      idIndicesLookup.put(idList[i], i);
    }
  }

  public byte getGenotypeForIndi(String marker, String fidiid) {
    int idIndex = idIndicesLookup.get(fidiid) == null ? -1 : idIndicesLookup.get(fidiid);
    if (idIndex == -1) {
      return (byte) -1;
    } else {
      int markerIndex =
          markerIndicesLookup.get(marker) == null ? -1 : markerIndicesLookup.get(marker);
      if (markerIndex == -1) {
        return (byte) -1;
      }
      while (!loaded[markerIndex]) {
        Thread.yield();
      }
      // synchronized(genotypes) {
      return genotypes[markerIndex][idIndex];
      // }
    }
  }

  public byte[] getGenotypesForMarker(String marker) {
    int markerIndex =
        markerIndicesLookup.get(marker) == null ? -1 : markerIndicesLookup.get(marker);
    if (markerIndex == -1) {
      return Array.byteArray(idList.length, (byte) -1);
    }
    while (!loaded[markerIndex]) {
      // try {
      // Thread.sleep(200);
      // } catch (InterruptedException e) {
      // // TODO Auto-generated catch block
      // continue;
      // }
      Thread.yield();
    }
    // synchronized(genotypes) {
    return genotypes[markerIndex];
    // }
  }

  public Thread getThread() {
    return thread;
  }

  private void initiate() {
    initialized = true;
  }

  public boolean isKilled() {
    return killed;
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

  private void lookupIDs() {
    BufferedReader reader;
    String[] line;

    famIDList = new String[Files.countLines(fileRoot + ".fam", 0)];

    try {
      reader = new BufferedReader(new FileReader(fileRoot + ".fam"));
      String temp = null;
      int cnt = 0;
      while ((temp = reader.readLine()) != null) {
        line = temp.trim().split("[\\s]+");
        famIDList[cnt] = line[0] + "\t" + line[1];
        cnt++;
      }
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

  private void lookupMarkerPositions() {
    BufferedReader reader;
    String[] line;
    HashSet<String> lookFor = new HashSet<String>();
    for (String s : markerList) {
      lookFor.add(s);
    }
    markerIndicesLookup = new HashMap<String, Integer>();
    HashMap<String, Integer> markerIndicesLookupTemp = new HashMap<String, Integer>();
    int cnt = 0;
    try {
      reader = new BufferedReader(new FileReader(fileRoot + ".bim"));
      String temp = null;
      while ((temp = reader.readLine()) != null) {
        line = temp.trim().split("\\s+");
        String mkr = line[1];
        if (lookFor.contains(mkr)) {
          markerIndicesLookupTemp.put(mkr, cnt);
        }
        cnt++;
      }
      reader.close();

      for (int i = 0; i < markerList.length; i++) {
        markerPositions[i] =
            markerIndicesLookupTemp.get(markerList[i]) == null ? -1
                                                               : markerIndicesLookupTemp.get(markerList[i])
                                                                                        .intValue();
        markerIndicesLookup.put(markerList[i], i);

      }
    } catch (FileNotFoundException fnfe) {
      // TODO should KILL here
      log.reportException(fnfe);
    } catch (IOException ioe) {
      // TODO should KILL here
      log.reportException(ioe);
    }
  }

  private void registerThread(Thread thread) {
    this.thread = thread;
  }

  @Override
  public void run() {
    RandomAccessFile in = null;

    HashMap<String, byte[]> mkrGenotypes = new HashMap<String, byte[]>();

    if (killed) {
      return;
    }

    initiate();
    int cnt = 0;
    long time = new Date().getTime();

    try {
      // synchronized(genotypes) {
      in = new RandomAccessFile(fileRoot + ".bed", "r");

      byte[] magicBytes = new byte[3];
      in.read(magicBytes);
      if (magicBytes[2] == 0) {
        log.reportError("Error - .bed file is sample-dominant.");
      } else {
        int famCnt = Files.countLines(fileRoot + ".fam", 0);
        int blockSize = (int) (famCnt / 4.0d); // TODO check for non-completeness (i.e. N not evenly
                                               // divisible by 4)
        // String[] markersSorted = sortMarkers(markerList, markerPositions);
        // int[] markerIndices = sortMarkerIndices(markerPositions);
        // int[] markerIndices = Sort.putInOrder(markerPositions);

        for (int i = 0; i < markerPositions.length; i++) {
          if (markerPositions[i] == -1) {
            // missing marker, not present in PLINK files
            mkrGenotypes.put(markerList[i], Array.byteArray(idList.length, (byte) -1));
            cnt++;
            continue;
          }
          // in.seek(markerPositions[i] * blockSize);
          in.seek(3 + i * blockSize);

          byte[] markerBytes = new byte[blockSize];
          byte[] sampGeno = new byte[idList.length];
          in.read(markerBytes);
          for (int bitInd = 0; bitInd < markerBytes.length; bitInd++) {
            byte bedByte = markerBytes[bitInd];
            byte[] genotypes = PlinkData.decodeBedByte(bedByte);

            for (int g = 0; g < genotypes.length; g++) {
              int idInd =
                  idIndicesLookup.get(famIDList[bitInd * 4
                                                + g]) == null ? -1
                                                              : idIndicesLookup.get(famIDList[bitInd
                                                                                              * 4
                                                                                              + g]);
              if (idInd == -1 || idInd > sampGeno.length) {
                continue;
              }
              sampGeno[idInd] = genotypes[g];
            }
          }

          mkrGenotypes.put(markerList[i], sampGeno);
          cnt++;
        }
      }

      in.close();

      for (int i = 0; i < markerList.length; i++) {
        // synchronized(genotypes) {
        genotypes[i] = mkrGenotypes.get(markerList[i]);
        loaded[i] = true;
        // }
      }
      // }

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

    if (killed) {
      log.report("PlinkMarkerLoader killed");
      markerList = null;
      markerPositions = null;
      loaded = null;
      genotypes = null;
      System.gc();
      killComplete = true;
    } else {
      log.report("Independent thread has finished loading " + cnt + " markers in "
                 + ext.getTimeElapsed(time));
    }

  }



}
