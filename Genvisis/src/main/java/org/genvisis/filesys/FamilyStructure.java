package org.genvisis.filesys;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

public class FamilyStructure {

  public static final String[][] TYPICAL_HEADERS = {{"FID", "famid"}, {"IID", "id"}, {"fa"}, {"mo"},
                                                    {"sex"}};

  private static final byte MISSING_VALUE_BYTE = (byte) -9;

  public static final int FID_INDEX = 0;
  public static final int IID_INDEX = 1;
  public static final int FA_INDEX = 2;
  public static final int MO_INDEX = 3;
  public static final int SEX_INDEX = 4;
  public static final int AFF_INDEX = 5;
  public static final int DNA_INDEX = 6;
  public static final int MZ_TWIN_INDEX = 7;

  public static final String MISSING_ID_STR = "0";

  public ArrayList<String[]> cached_poPairsIDs = null;
  public ArrayList<int[]> cached_poPairsCompleteOnly = null;
  public HashMap<String, ArrayList<String>> cached_parentToChildrenMap = null;
  public boolean cached_parentMapIsCompleteOnly = false;
  public ArrayList<String[]> cached_all_trios = null;
  public ArrayList<int[]> cached_complete_trios = null;
  public ArrayList<String[]> cached_sib_pairs = null;
  public HashMap<String, Integer> cached_fidiidToIndexMap = null;

  public void clearCache() {
    cached_poPairsIDs = null;
    cached_poPairsCompleteOnly = null;
    cached_parentToChildrenMap = null;
    cached_parentMapIsCompleteOnly = false;
    cached_all_trios = null;
    cached_complete_trios = null;
    cached_sib_pairs = null;
    cached_fidiidToIndexMap = null;
  }

  protected String[][] ids;
  protected String[] iids;
  protected String[] fids;
  protected String[] fas;
  protected String[] mos;
  protected byte[] genders;
  protected byte[] affections;
  protected String[] dnas;
  protected String[] mzTwinIds;
  protected Logger log;

  public FamilyStructure(String[][] ids, byte[] genders, byte[] affections) {
    this(ids, genders, affections, null);
  }

  public FamilyStructure(String[][] ids, byte[] genders, byte[] affections, String[] dnas) {
    this(ids, genders, affections, dnas, new String[ids.length], new Logger());
  }

  public FamilyStructure(String[][] ids, byte[] genders, byte[] affections, String[] dnas,
                         String[] mzTwinIds, Logger log) {
    this.ids = ids;
    this.genders = genders;
    this.affections = affections;
    this.dnas = dnas;
    this.mzTwinIds = mzTwinIds;
    this.log = log;
  }

  public FamilyStructure(String filename) {
    this(filename, false);
  }

  public FamilyStructure(String filename, boolean loadDNAs) {
    this(filename, loadDNAs, new Logger());
  }

  public FamilyStructure(String filename, boolean loadDNAs, Logger log) {

    String[][] pedCols = Matrix.transpose(HashVec.loadFileToStringMatrix(filename, false, null));
    if (pedCols == null) {
      return;
    }
    int count = pedCols[0].length;

    ids = new String[count][];
    fids = pedCols[FID_INDEX];
    iids = pedCols[IID_INDEX];
    fas = pedCols[FA_INDEX];
    mos = pedCols[MO_INDEX];
    genders = new byte[count];
    affections = new byte[count];
    dnas = loadDNAs ? pedCols[DNA_INDEX] : null;
    mzTwinIds = new String[count];
    for (int i = 0; i < count; i++) {
      ids[i] = new String[] {fids[i], iids[i], fas[i], mos[i]};
      genders[i] = ext.isMissingValue(pedCols[SEX_INDEX][i]) ? FamilyStructure.MISSING_VALUE_BYTE
                                                             : Byte.parseByte(pedCols[SEX_INDEX][i]);
      affections[i] = ext.isMissingValue(pedCols[AFF_INDEX][i]) ? FamilyStructure.MISSING_VALUE_BYTE
                                                                : Byte.parseByte(pedCols[AFF_INDEX][i]);

      if (pedCols.length > MZ_TWIN_INDEX) {
        mzTwinIds[i] = ext.isMissingValue(pedCols[MZ_TWIN_INDEX][i]) ? null
                                                                     : pedCols[MZ_TWIN_INDEX][i];
      }
    }
  }

  public String[][] getIDs() {
    return ids;
  }

  /**
   * @param fidiid the FID and IID of the individual, separated by a tab
   * @return the index of the provided fidiid in the {@link FamilyStructure} or -1 for missing
   */
  public int getIndIndex(String fidiid) {
    Integer index = getfidiidToIndexMap().get(fidiid);
    return index == null ? -1 : index;
  }

  /**
   * @param fid
   * @param iid
   * @return the index of the provided fidiid in the {@link FamilyStructure} or -1 for missing
   */
  public int getIndIndex(String fid, String iid) {
    return getIndIndex(fid + "\t" + iid);
  }

  public HashMap<String, Integer> getfidiidToIndexMap() {
    if (cached_fidiidToIndexMap == null) {
      buildfidiidToIndexMap();
    }
    return cached_fidiidToIndexMap;
  }

  private synchronized void buildfidiidToIndexMap() {
    if (cached_fidiidToIndexMap == null) {
      HashMap<String, Integer> fidiidMap = new HashMap<String, Integer>();
      for (int i = 0; i < ids.length; i++) {
        if (fidiidMap.put(fids[i] + "\t" + iids[i], i) != null) {
          System.err.println("Warning - Pedigree contains non-unique FID/IID combinations!");
        }
      }
      cached_fidiidToIndexMap = fidiidMap;
    }
  }

  public String getiDNA(int i) {
    return dnas[i];
  }

  public String getFID(int index) {
    return ids[index][FID_INDEX];
  }

  public String getIID(int index) {
    return ids[index][IID_INDEX];
  }

  public String getFA(int index) {
    return ids[index][FA_INDEX];
  }

  public int getIndexOfFaInIDs(int indivIndex) {
    return MISSING_ID_STR.equals(fas[indivIndex]) ? -1
                                                  : getIndIndex(fids[indivIndex], fas[indivIndex]);
  }

  public int getIndexOfMoInIDs(int indivIndex) {
    return MISSING_ID_STR.equals(mos[indivIndex]) ? -1
                                                  : getIndIndex(fids[indivIndex], mos[indivIndex]);
  }

  public String getMO(int index) {
    return ids[index][MO_INDEX];
  }

  public byte[] getGenders() {
    return genders;
  }

  public byte getGender(int index) {
    return genders[index];
  }

  public byte[] getAffections() {
    return affections;
  }

  public void setAffections(byte[] affs) {
    affections = affs;
  }

  public String[] getDnas() {
    return dnas;
  }

  public String[] getMzTwinIds() {
    return mzTwinIds;
  }

  public String getMzTwinId(int index) {
    return mzTwinIds[index];
  }

  public String getIndividualHeader(int index, boolean displayDNA) {
    return ArrayUtils.toStr(ids[index]) + "\t" + genders[index] + "\t" + affections[index]
           + (displayDNA && dnas != null ? "\t" + dnas[index] : "");
  }

  public void writeToFile(String filename, boolean displayDNA) {
    PrintWriter writer;

    try {
      writer = Files.openAppropriateWriter(filename);
      for (int i = 0; i < ids.length; i++) {
        writer.println(getIndividualHeader(i, displayDNA));
      }
      writer.close();
    } catch (Exception e) {
      log.reportException(e);
    }
  }

  public static boolean likelyPedHeader(String[] line) {
    return ArrayUtils.countIf(ext.indexFactors(TYPICAL_HEADERS, line, false, true, false), -1) < 3;
  }
}
