// autogenerate the master list
package org.genvisis.park;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class PlateLists {
  public static boolean PROGENI_NOT_CARES = true;
  public static final String[] COLS = {"Plate", "Well", "DNA#", "FamInd"};
  public static final int NUM_ON_PLATE = 86;
  public static final int NUM_IN_ROW = 12;
  public static final boolean UPDATE = true;
  // public static final String NEW_FILENAME = "PROGENI_067b.txt";
  public static final String NEW_FILENAME = "PROGENI_068.txt";
  // public static final String NEW_FILENAME = "firstCARES.txt";
  // public static final String NEW_FILENAME = "C-008.txt";
  public static final String PLATE_LIST_DIRECTORY = tools.MASTER_PLATELIST_DIR;
  public static final String PROGENI = "PROGENI";
  public static final String CARES = "CARES";
  public static final String[][] LABELS = {{"probands.xln", "Probands", "Sporadics"},
                                           {"other_affecteds.xln", "OtherAffecteds",
                                            "ShouldNotBeOtherAffecteds"},
                                           {"unaffecteds.xln", "Unaffecteds", "Controls"},
                                           {"temps.xln", "TEMP", "CARES_TEMPS"}};
  public static final String[] IGNORE_CODES = new String[] {"CONF", "HS", "SWAPPED", "MZLINKAGE",
                                                            "DEL4LINKAGE", "POORYIELD", "MUT",
                                                            "PHENO"};
  public static final String[] CREATE_HEADER = {"DNA", "ID", "OldPlate", "OldWell", "NewPlate",
                                                "NewWell", "Batch", "Dx", "Note"};
  public static final String[] UPDATE_FROM_CREATE_HEADER = {"DNA", "ID", "NewPlate", "NewWell",
                                                            "Batch", "Dx", "Note"};
  public static final String[] UPDATE_HEADER =
                                             {"DNA", "ID", "Plate", "Well", "Batch", "Dx", "Note"};

  public static void create(String dir, String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, fams, inds, dnas;
    String trav, fam;
    Hashtable<String, Hashtable<String, Vector<DnaSample>>> hash =
                                                                 new Hashtable<String, Hashtable<String, Vector<DnaSample>>>();
    Hashtable<String, Vector<DnaSample>> family;
    Vector<DnaSample> v = new Vector<DnaSample>();
    int[] indices, keys;
    Vector<Vector<DnaSample>> categories;
    DnaSample dna;
    CheckIDsAgainstDNAs idCheck = new CheckIDsAgainstDNAs();
    Hashtable<String, String> diagnoses = tools.getExactPDdx(PROGENI_NOT_CARES ? 2 : 6);
    int proband, probandLevel;
    String dx;
    Vector<String> mzs;
    Hashtable<String, String> deldna, delind;

    mzs = new Vector<String>();
    deldna = new Hashtable<String, String>();
    delind = new Hashtable<String, String>();
    try {
      reader = tools.getNinfoReader(3);
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t");
        if (line[0].equals("MZ")) {
          mzs.add(tools.getUniqueID(line[3], line[4]));
        } else if (line[0].equals("DELIND")) {
          delind.put(tools.getUniqueID(line[1], line[2]), line[4]);
        } else if (line[0].equals("DELDNA")) {
          deldna.put(line[3], line[4]);
        } else if (ext.indexOfStr(line[0], IGNORE_CODES) >= 0) {

        } else {
          System.err.println("Error - unknown ninfo3 code: " + line[0]);
        }
      }
      reader.close();
    } catch (IOException e) {
      System.err.println("Error parsing ninfo3.dat");
      System.exit(2);
    }

    categories = new Vector<Vector<DnaSample>>();
    for (String[] element : LABELS) {
      categories.add(new Vector<DnaSample>());
    }
    try {
      reader = new BufferedReader(new FileReader(dir + filename));
      line = reader.readLine().trim().split("\t");
      indices = ext.indexFactors(COLS, line, false, true);
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t");
        dna =
            new DnaSample(line[indices[ext.indexOfStr("DNA#", COLS)]],
                          line[indices[ext.indexOfStr("FamInd", COLS)]],
                          line[indices[ext.indexOfStr("Plate", COLS)]],
                          line[indices[ext.indexOfStr("Well", COLS)]],
                          // ext.rootOf(filename),
                          ".", false);
        dna.setNote(".");

        trav = dna.getID();
        dna.setDx(diagnoses.get(trav));

        if (delind.containsKey(trav)) {
          dna.setNote(delind.get(trav));
          categories.elementAt(3).add(dna);
        } else if (deldna.contains(dna.getDNA())) {
          dna.setNote(deldna.get(trav));
          categories.elementAt(3).add(dna);
        } else if (mzs.contains(trav)) {
          dna.setNote("MZ twin pair");
          categories.elementAt(3).add(dna);
        } else {
          fam = trav.substring(0, 5);
          idCheck.checkUsingUniqueID(trav, dna.getDNA(), false);
          if (hash.containsKey(fam)) {
            family = hash.get(fam);
          } else {
            hash.put(fam, family = new Hashtable<String, Vector<DnaSample>>());
          }
          if (family.containsKey(trav)) {
            v = family.get(trav);
          } else {
            family.put(trav, v = new Vector<DnaSample>());
          }
          v.add(dna);
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + filename + "\"");
      System.exit(2);
    }

    try {
      fams = HashVec.getKeys(hash, true, true);
      for (String fam2 : fams) {
        family = hash.get(fam2);
        inds = HashVec.getKeys(family, true, true);
        proband = probandLevel = -1;
        for (int j = 0; j < inds.length; j++) {
          v = family.get(inds[j]);
          dnas = new String[v.size()];
          for (int k = 0; k < v.size(); k++) {
            dnas[k] = v.elementAt(k).getDNA();
          }
          keys = Sort.quicksort(dnas);
          for (int k = v.size() - 1; k >= 0; k--) {
            if (k != keys[0]) {
              if (v.elementAt(k).getDNA().equals(v.elementAt(keys[0]).getDNA())) {
                v.elementAt(k).setNote("Redraw");
              } else {
                v.elementAt(k).setNote("Duplicate of " + v.elementAt(keys[0]).getDNA());
              }
              categories.elementAt(3).add(v.remove(k));
            }
          }

          dx = diagnoses.get(inds[j]);
          if (dx == null) {
            System.err.println("Error - no diagnosis for " + inds[j]
                               + " (will be assigned to the unaffecteds plate)");
            dx = "Unknown";
          } else {
            if (dx.equals("CONF_PD")) {
              if (probandLevel < 3) {
                proband = j;
                probandLevel = 3;
              }
            } else if (dx.equals("VPD")) {
              if (probandLevel < 2) {
                proband = j;
                probandLevel = 2;
              }
            } else if (dx.equals("NVPD")) {
              if (probandLevel < 1) {
                proband = j;
                probandLevel = 1;
              }
            }
          }
          v.elementAt(0).setDx(dx);
        }

        if (proband >= 0) {
          categories.elementAt(0).add(family.remove(inds[proband]).elementAt(0));
        }
        inds = HashVec.getKeys(family, true, true);
        for (String ind : inds) {
          dx = diagnoses.get(ind);
          if (dx == null) {
            // System.err.println("Warning - assigning "+inds[j]+"
            // to unaffecteds plate");
            categories.elementAt(2).add(family.remove(ind).elementAt(0));
          } else if (ext.indexOfStr(dx, tools.PD_DX) >= 0) {
            categories.elementAt(1).add(family.remove(ind).elementAt(0));
          } else if (ext.indexOfStr(dx, tools.UNAFFECTED) >= 0) {
            categories.elementAt(2).add(family.remove(ind).elementAt(0));
          }
        }
      }
    } catch (Exception e) {
      System.err.println("Error writing to file");
    }

    try {
      writer = new PrintWriter(new FileWriter(dir + "MasterFile.xln"));
      writer.println(Array.toStr(CREATE_HEADER));
      for (int i = 0; i < LABELS.length; i++) {
        printDNAsToNewPlates(dir, LABELS[i][0], writer, categories.elementAt(i),
                             LABELS[i][PROGENI_NOT_CARES ? 1 : 2], 0);
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to master file");
      e.printStackTrace();
    }
  }

  public static void update(String dir, String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, fams, inds, dnas;
    String trav, fam;
    Hashtable<String, Hashtable<String, Vector<DnaSample>>> hash =
                                                                 new Hashtable<String, Hashtable<String, Vector<DnaSample>>>();
    Hashtable<String, Vector<DnaSample>> family;
    Vector<DnaSample> v = new Vector<DnaSample>();
    int[] indices, keys;
    Vector<Vector<DnaSample>> categories;
    DnaSample dna;
    CheckIDsAgainstDNAs idCheck = new CheckIDsAgainstDNAs();
    Hashtable<String, String> diagnoses = tools.getExactPDdx(PROGENI_NOT_CARES ? 2 : 6);
    int proband, probandLevel, prior, priorLevel;
    String dx;
    Vector<String> mzs;
    Hashtable<String, String> deldna, delind;
    int[] maxes;
    String plate;
    int plateNumber, seriesIndex;

    mzs = new Vector<String>();
    deldna = new Hashtable<String, String>();
    delind = new Hashtable<String, String>();
    try {
      reader = tools.getNinfoReader(3);
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t");
        if (line[0].equals("MZ")) {
          mzs.add(tools.getUniqueID(line[3], line[4]));
        } else if (line[0].equals("DELIND")) {
          delind.put(tools.getUniqueID(line[1], line[2]), line[4]);
        } else if (line[0].equals("DELDNA")) {
          deldna.put(line[3], line[4]);
        } else if (ext.indexOfStr(line[0], IGNORE_CODES) >= 0) {

        } else {
          System.err.println("Error - unknown ninfo3 code: " + line[0]);
        }
      }
      reader.close();
    } catch (IOException e) {
      System.err.println("Error parsing ninfo3.dat");
      System.exit(2);
    }

    categories = new Vector<Vector<DnaSample>>();
    for (String[] element : LABELS) {
      categories.add(new Vector<DnaSample>());
    }
    maxes = new int[LABELS.length];
    for (int i = 0; i < LABELS.length; i++) {
      try {
        reader = new BufferedReader(new FileReader(dir + (PROGENI_NOT_CARES ? PROGENI : CARES) + "/"
                                                   + LABELS[i][0]));
        line = reader.readLine().trim().split("\t");
        if (ext.eqArrays(line, CREATE_HEADER)) {
          indices = ext.indexFactors(UPDATE_FROM_CREATE_HEADER, line, false, true);
        } else {
          indices = ext.indexFactors(UPDATE_HEADER, line, false, true);
        }
        while (reader.ready()) {
          line = reader.readLine().trim().split("\t");
          dna = new DnaSample(line[indices[0]], line[indices[1]], line[indices[2]],
                              line[indices[3]], line[indices[4]], true);
          dna.setDx(line[indices[5]]);
          dna.setNote(line[indices[6]]);
          if (i == 0) {
            dna.setProband();
          }
          trav = dna.getID();

          if (delind.containsKey(trav)) {
            dna.setNote(delind.get(trav));
            categories.elementAt(3).add(dna);
          } else if (deldna.contains(dna.getDNA())) {
            dna.setNote(deldna.get(trav));
            categories.elementAt(3).add(dna);
          } else if (mzs.contains(trav)) {
            dna.setNote("MZ twin pair");
            categories.elementAt(3).add(dna);
          } else {
            fam = trav.substring(0, 5);
            idCheck.checkUsingUniqueID(trav, dna.getDNA(), false);
            if (hash.containsKey(fam)) {
              family = hash.get(fam);
            } else {
              hash.put(fam, family = new Hashtable<String, Vector<DnaSample>>());
            }
            if (family.containsKey(trav)) {
              v = family.get(trav);
            } else {
              family.put(trav, v = new Vector<DnaSample>());
            }
            v.add(dna);
            plate = dna.getPlate();
            try {
              plateNumber = Integer.parseInt(plate.substring(plate.length() - 2));
            } catch (Exception e) {
              // System.out.println("No plate number for:
              // "+plate);
              plateNumber = 1;
            }
            seriesIndex = translateFromPlateWell(plateNumber, dna.getWell());
            maxes[i] = Math.max(maxes[i], seriesIndex);
          }
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + LABELS[i][0] + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + LABELS[i][0] + "\"");
        System.exit(2);
      }
    }

    try {
      reader = new BufferedReader(new FileReader(dir + filename));
      line = reader.readLine().trim().split("\t");
      indices = ext.indexFactors(new String[] {"UniqueID", "DNA"}, line, false, true);
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t");
        dna = new DnaSample(line[indices[1]], line[indices[0]], ext.rootOf(filename));

        trav = dna.getID();

        if (delind.containsKey(trav)) {
          dna.setNote(delind.get(trav));
          categories.elementAt(3).add(dna);
        } else if (deldna.contains(dna.getDNA())) {
          dna.setNote(deldna.get(trav));
          categories.elementAt(3).add(dna);
        } else if (mzs.contains(trav)) {
          dna.setNote("MZ twin pair");
          categories.elementAt(3).add(dna);
        } else {
          fam = trav.substring(0, 5);
          idCheck.checkUsingUniqueID(trav, dna.getDNA(), false);
          if (hash.containsKey(fam)) {
            family = hash.get(fam);
          } else {
            hash.put(fam, family = new Hashtable<String, Vector<DnaSample>>());
          }
          if (family.containsKey(trav)) {
            v = family.get(trav);
          } else {
            family.put(trav, v = new Vector<DnaSample>());
          }
          v.add(dna);
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + filename + "\"");
      System.exit(2);
    }

    try {
      fams = HashVec.getKeys(hash, true, true);
      for (String fam2 : fams) {
        family = hash.get(fam2);
        inds = HashVec.getKeys(family, true, true);
        proband = probandLevel = prior = priorLevel = -1;
        for (int j = 0; j < inds.length; j++) {
          v = family.get(inds[j]);
          dnas = new String[v.size()];
          keys = null;
          for (int k = 0; k < v.size(); k++) {
            dnas[k] = v.elementAt(k).getDNA();
            if (v.elementAt(k).getProband()) {
              keys = new int[] {k};
            }
          }
          if (keys == null) {
            keys = Sort.quicksort(dnas);
          }
          for (int k = v.size() - 1; k >= 0; k--) {
            if (k != keys[0]) {
              if (v.elementAt(k).getDNA().equals(v.elementAt(keys[0]).getDNA())) {
                v.elementAt(k).setNote("Redraw");
              } else {
                v.elementAt(k).setNote("Duplicate of " + v.elementAt(keys[0]).getDNA());
              }
              categories.elementAt(3).add(v.remove(k));
            }
          }

          dx = diagnoses.get(inds[j]);
          if (dx == null) {
            System.err.println("Error - no diagnosis for " + inds[j]
                               + " (will be assigned to the unaffecteds plate)");
            dx = "Unknown";
          } else {
            if (v.elementAt(0).getProband()) {
              prior = j;
              priorLevel = dx.equals("CONF_PD") ? 3
                                                : (dx.equals("VPD") ? 2
                                                                    : (dx.equals("NVPD") ? 1 : 0));
            }
            if (dx.equals("CONF_PD")) {
              if (probandLevel < 3) {
                proband = j;
                probandLevel = 3;
              }
            } else if (dx.equals("VPD")) {
              if (probandLevel < 2) {
                proband = j;
                probandLevel = 2;
              }
            } else if (dx.equals("NVPD")) {
              if (probandLevel < 1) {
                proband = j;
                probandLevel = 1;
              }
            } else if (dx.equals("RPD")) {
              if (probandLevel < 0) {
                proband = j;
                probandLevel = 0;
              }
            }
          }
          // if (v.elementAt(0).getDx() == null) { // remove after
          // checking, so that this will update in the future
          v.elementAt(0).setDx(dx);
          // }
        }

        if (probandLevel > priorLevel && priorLevel != -1) {
          System.err.println("FYI: " + family.get(inds[prior]).elementAt(0).getID() + " ("
                             + diagnoses.get(family.get(inds[prior]).elementAt(0).getID())
                             + ") could be replaced by "
                             + family.get(inds[proband]).elementAt(0).getID() + " ("
                             + diagnoses.get(family.get(inds[proband]).elementAt(0).getID()) + ")");
        }
        if (prior >= 0) {
          categories.elementAt(0).add(family.remove(inds[prior]).elementAt(0));
        } else if (proband >= 0) {
          categories.elementAt(0).add(family.remove(inds[proband]).elementAt(0));
        }
        inds = HashVec.getKeys(family, true, true);
        for (String ind : inds) {
          dx = diagnoses.get(ind);
          if (dx == null) {
            // System.err.println("Warning - assigning "+inds[j]+"
            // to unaffecteds plate");
            categories.elementAt(2).add(family.remove(ind).elementAt(0));
          } else if (ext.indexOfStr(dx, tools.PD_DX) >= 0) {
            categories.elementAt(1).add(family.remove(ind).elementAt(0));
          } else if (ext.indexOfStr(dx, tools.UNAFFECTED) >= 0) {
            categories.elementAt(2).add(family.remove(ind).elementAt(0));
          } else {
            System.err.println("Error - failed to place " + ind + " ("
                               + family.get(ind).elementAt(0).getDNA() + ")");
            family.get(ind).elementAt(0).setNote("No valid phenotype for this indiviudal");
            // categories.elementAt(0).add(family.remove(inds[j]).elementAt(0));
          }
        }
      }
    } catch (Exception e) {
      System.err.println("Error writing to file");
      e.printStackTrace();
    }

    try {
      writer = new PrintWriter(new FileWriter(dir + "MasterFile.xln"));
      writer.println(Array.toStr(UPDATE_HEADER));
      for (int i = 0; i < LABELS.length; i++) {
        printDNAsToNewPlates(dir, LABELS[i][0], writer, categories.elementAt(i),
                             LABELS[i][PROGENI_NOT_CARES ? 1 : 2], maxes[i] + 1);
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to master file");
      e.printStackTrace();
    }
  }

  public static void printDNAsToNewPlates(String dir, String filename, PrintWriter master,
                                          Vector<DnaSample> v, String plateName, int startIndex) {
    PrintWriter writer;
    DnaSample dna;
    int count;

    try {
      writer = new PrintWriter(new FileWriter(dir + filename));
      writer.println(Array.toStr(startIndex == 0 ? CREATE_HEADER : UPDATE_HEADER));
      for (int i = 0; i < v.size(); i++) {
        dna = v.elementAt(i);
        if (dna.getPrior()) {
          writer.println(dna.getDNA() + "\t" + dna.getID()
                         + (startIndex == 0 ? "\tprior\tprior" : "") + "\t" + dna.getPlate() + "\t"
                         + dna.getWell() + "\t" + dna.getBatch() + "\t" + dna.getDx() + "\t"
                         + dna.getNote());
          master.println(dna.getDNA() + "\t" + dna.getID()
                         + (startIndex == 0 ? "\tprior\tprior" : "") + "\t" + dna.getPlate() + "\t"
                         + dna.getWell() + "\t" + dna.getBatch() + "\t" + dna.getDx() + "\t"
                         + dna.getNote());
        }
      }
      count = 0;
      for (int i = 0; i < v.size(); i++) {
        dna = v.elementAt(i);
        if (!dna.getPrior()) {
          writer.println(dna.getDNA() + "\t" + dna.getID()
                         + (startIndex == 0 ? "\t" + dna.getPlate() + "\t" + dna.getWell() : "")
                         + "\t"
                         + translateToPlateWell(plateName, startIndex + count,
                                                v.size() < NUM_ON_PLATE)
                         + "\t" + dna.getBatch() + "\t" + dna.getDx() + "\t" + dna.getNote());
          master.println(dna.getDNA() + "\t" + dna.getID()
                         + (startIndex == 0 ? "\t" + dna.getPlate() + "\t" + dna.getWell() : "")
                         + "\t"
                         + translateToPlateWell(plateName, startIndex + count,
                                                v.size() < NUM_ON_PLATE)
                         + "\t" + dna.getBatch() + "\t" + dna.getDx() + "\t" + dna.getNote());
          count++;
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + dir + filename);
    }
  }

  public static String translateToPlateWell(String plateName, int index, boolean singlePlate) {
    int plateNumber = index / NUM_ON_PLATE;
    int wellIndex = index - plateNumber * NUM_ON_PLATE;

    return plateName + (singlePlate ? "" : ext.formNum(plateNumber + 1, 2)) + "\t"
           + translateToWell(wellIndex);
  }

  public static int translateFromPlateWell(int plate, String well) {
    int index;

    index = (plate - 1) * NUM_ON_PLATE;
    index += ext.getIndexFromColumn(well.substring(0, 1)) * NUM_IN_ROW
             + Integer.parseInt(well.substring(1)) - 1;

    return index;
  }

  public static String translateToWell(int index) {
    int row = index / NUM_IN_ROW;
    int col = index - row * NUM_IN_ROW;

    return ext.getColumn(row) + (col + 1);
  }

  public static void batchCombine() {
    String str = "Send, {ENTER}" + "\n" + "Sleep, 3000" + "\n" + "" + "\n" + "loop, 4 {" + "\n"
                 + "Send, {SHIFTDOWN}{F11}{SHIFTUP}" + "\n" + "Sleep, 200" + "\n"
                 + "Send, {ALTDOWN}{TAB}{ALTUP}" + "\n" + "Sleep, 200" + "\n"
                 + "Send, {DOWN}{ENTER}" + "\n" + "Sleep, 1500" + "\n"
                 + "Send, {CTRLDOWN}a{CTRLUP}{CTRLDOWN}c{F4}{CTRLUP}" + "\n" + "Sleep, 300" + "\n"
                 + "Send, {ENTER}" + "\n" + "Sleep, 500" + "\n" + "Send, {CTRLDOWN}v{CTRLUP}" + "\n"
                 + "Sleep, 500" + "\n" + "Send, {ALTDOWN}{TAB}{ALTUP}" + "\n" + "Sleep, 500" + "\n"
                 + "Send, {F2}" + "\n" + "Sleep, 200" + "\n"
                 + "Send, {SHIFTDOWN}{LEFT}{LEFT}{LEFT}{LEFT}{SHIFTUP}" + "\n" + "Sleep, 100" + "\n"
                 + "Send, {CTRLDOWN}c{CTRLUP}" + "\n" + "Sleep, 100" + "\n" + "Send, {ESC}" + "\n"
                 + "Sleep, 200" + "\n" + "Send, {ALTDOWN}{TAB}{ALTUP}" + "\n" + "Sleep, 500" + "\n"
                 + "MouseClick, left,  102,  966" + "\n" + "Sleep, 100" + "\n"
                 + "MouseClick, left,  102,  966" + "\n" + "Sleep, 100" + "\n"
                 + "Send, {CTRLDOWN}v{CTRLUP}{ENTER}" + "\n" + "Sleep, 100" + "\n" + "}" + "\n";

    ext.setClipboard(str);
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = PLATE_LIST_DIRECTORY;
    String filename = NEW_FILENAME;
    boolean combine = true;

    String usage = "\\n" + "park.PlateLists requires 0-1 arguments\n"
                   + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    try {
      if (combine) {
        batchCombine();
      } else if (UPDATE) {
        update(dir, filename);
      } else {
        create(dir, filename);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}


class DnaSample {
  private final String dna;

  private final String id;

  private final String plate;

  private final String well;

  private final String batch;

  private String dx;

  private String note;

  private final boolean prior;

  private boolean proband;

  public DnaSample(String dna, String id, String batch) {
    this(dna, id, "new", "new", batch, false);
  }

  public DnaSample(String dna, String id, String plate, String well, String batch, boolean prior) {
    this.dna = dna;
    this.id = id;
    this.plate = plate;
    this.well = well;
    this.batch = batch;
    note = ".";
    this.prior = prior;
    proband = false;
  }

  public String getID() {
    return id;
  }

  public String getDNA() {
    return dna;
  }

  public String getPlate() {
    return plate;
  }

  public String getWell() {
    return well;
  }

  public String getBatch() {
    return batch;
  }

  public void setDx(String dx) {
    this.dx = dx;
  }

  public String getDx() {
    return dx;
  }

  public void setNote(String reason) {
    note = reason;
  }

  public String getNote() {
    return note;
  }

  public void setProband() {
    proband = true;
  }

  public boolean getProband() {
    return proband;
  }

  public boolean getPrior() {
    return prior;
  }
}
