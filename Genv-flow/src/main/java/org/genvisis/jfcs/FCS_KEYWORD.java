package org.genvisis.jfcs;

import java.util.HashSet;
import java.util.Set;
import java.util.regex.Pattern;

public enum FCS_KEYWORD {
  // required by spec:
  BEGINANALYSIS("$BEGINANALYSIS", Pattern.quote("$BEGINANALYSIS"), true),
  BEGINDATA("$BEGINDATA", Pattern.quote("$BEGINDATA"), true),
  BEGINSTEXT("$BEGINSTEXT", Pattern.quote("$BEGINSTEXT"), true),
  BYTEORD("$BYTEORD", Pattern.quote("$BYTEORD"), true),
  DATATYPE("$DATATYPE", Pattern.quote("$DATATYPE"), true),
  ENDANALYSIS("$ENDANALYSIS", Pattern.quote("$ENDANALYSIS"), true),
  ENDDATA("$ENDDATA", Pattern.quote("$ENDDATA"), true),
  ENDSTEXT("$ENDSTEXT", Pattern.quote("$ENDSTEXT"), true),
  MODE("$MODE", Pattern.quote("$MODE"), true),
  NEXTDATA("$NEXTDATA", Pattern.quote("$NEXTDATA"), true),
  PAR("$PAR", Pattern.quote("$PAR"), true),
  PnB("$PnB", "\\$P\\d+B", true),
  PnE("$PnE", "\\$P\\d+E", true),
  PnN("$PnN", "\\$P\\d+N", true),
  PnR("$PnR", "\\$P\\d+R", true),
  TOT("$TOT", "\\$TOT", true),

  ABRT("$ABRT", "\\$ABRT", false),
  BTIM("$BTIM", "\\$BTIM", false),
  CELLS("$CELLS", "\\$CELLS", false),
  COM("$COM", "\\$COM", false),
  CSMODE("$CSMODE", "\\$CSMODE", false),
  CSVBITS("$CSVBITS", "\\$CSVBITS", false),
  CSVnFLAG("$CSVnFLAG", "\\$CSV\\d+FLAG", false),
  CYT("$CYT", "\\$CYT", false),
  CYTSN("$CYTSN", "\\$CYTSN", false),
  DATE("$DATE", "\\$DATE", false),
  ETIM("$ETIM", "\\$ETIM", false),
  EXP("$EXP", "\\$EXP", false),
  FIL("$FIL", "\\$FIL", false),
  GATE("$GATE", "\\$GATE", false),
  GATING("$GATING", "\\$GATING", false),
  GnE("$GnE", "\\$G\\d+E", false),
  GnF("$GnF", "\\$G\\d+F", false),
  GnN("$GnN", "\\$G\\d+N", false),
  GnP("$GnP", "\\$G\\d+P", false),
  GnR("$GnR", "\\$G\\d+R", false),
  GnS("$GnS", "\\$G\\d+S", false),
  GnT("$GnT", "\\$G\\d+T", false),
  GnV("$GnV", "\\$G\\d+V", false),
  INST("$INST", "\\$INST", false),
  LAST_MODIFIED("$LAST_MODIFIED", "\\$LAST_MODIFIED", false),
  LAST_MODIFIER("$LAST_MODIFIER", "\\$LAST_MODIFIER", false),
  LOST("$LOST", "\\$LOST", false),
  OP("$OP", "\\$OP", false),
  ORIGINALITY("$ORIGINALITY", "\\$ORIGINALITY", false),
  PKn("$PKn", "\\$PK\\d+", false),
  PKNn("$PKNn", "\\$PKN\\d+", false),
  PLATEID("$PLATEID", "\\$PLATEID", false),
  PLATENAME("$PLATENAME", "\\$PLATENAME", false),
  PnCALIBRATION("$PnCALIBRATION", "\\$P\\d+CALIBRATION", false),
  PnD("$PnD", "\\$P\\d+D", false),
  PnF("$PnF", "\\$P\\d+F", false),
  PnG("$PnG", "\\$P\\d+G", false),
  PnL("$PnL", "\\$P\\d+L", false),
  PnO("$PnO", "\\$P\\d+O", false),
  PnP("$PnP", "\\$P\\d+P", false),
  PnS("$PnS", "\\$P\\d+S", false),
  PnT("$PnT", "\\$P\\d+T", false),
  PnV("$PnV", "\\$P\\d+V", false),
  PROJ("$PROJ", "\\$PROJ", false),
  RnI("$RnI", "\\$R\\d+I", false),
  RnW("$RnW", "\\$R\\d+W", false),
  SMNO("$SMNO", "\\$SMNO", false),
  SPILLOVER("$SPILLOVER", "\\$SPILLOVER", false),
  SRC("$SRC", "\\$SRC", false),
  SYS("$SYS", "\\$SYS", false),
  TIMESTEP("$TIMESTEP", "\\$TIMESTEP", false),
  TR("$TR", "\\$TR", false),
  VOL("$VOL", "\\$VOL", false),
  WELLID("$WELLID", "\\$WELLID", false),
  // custom
  SPILL("SPILL", "\\SPILL", false);

  FCS_KEYWORD(String key, String pattern, boolean required) {
    this.required = required;
    this.keyword = key;
    this.patternStr = pattern;
    this.pattern = Pattern.compile(patternStr);
    this.isInd = key.contains("n");
  }

  static Set<FCS_KEYWORD> getRequiredKeywords() {
    HashSet<FCS_KEYWORD> keys = new HashSet<>();
    for (FCS_KEYWORD k : values()) {
      if (k.required) {
        keys.add(k);
      }
    }
    return keys;
  }

  boolean required;
  String keyword;
  String patternStr;
  boolean isInd;
  Pattern pattern;

  String getKey() {
    return keyword;
  }

  String getKey(int index) {
    return isInd ? keyword.replace("n", Integer.toString(index)) : keyword;
  }

  static FCS_KEYWORD findKey(String poss) {
    try {
      // check if valid
      FCS_KEYWORD word = FCS_KEYWORD.valueOf(poss.startsWith("$") ? poss.substring(1) : poss);
      return word;
    } catch (IllegalArgumentException e) {
      for (FCS_KEYWORD word : values()) {
        if (!word.isInd) {
          continue;
        } else if (word.pattern.matcher(poss).matches()) {
          return word;
        }
      }
      return null;
    }
  }

  public static boolean isKey(String poss) {
    return findKey(poss) != null;
  }
}
