package org.genvisis.seq.manage;

import static org.junit.Assert.assertEquals;
import java.util.Map.Entry;
import org.genvisis.seq.manage.StrandOps.AlleleMatch;
import org.genvisis.seq.manage.StrandOps.AlleleStatus;
import org.junit.Test;
import org.pankratzlab.common.ArrayUtils;
import com.google.common.collect.ImmutableMap;

public class TestStrandOps {

  @Test
  public void testDetermineStrandConfig() {
    ImmutableMap.Builder<String[][], StrandOps.AlleleConfig> testMapBld = ImmutableMap.builder();
    /*
     * @format:off
     */
    // missing:
    testMapBld.put(new String[][]{{"NA", "NA"}, {"A", "T"}}, new StrandOps.AlleleConfig(AlleleStatus.NULL, AlleleStatus.NULL, AlleleMatch.UNKNOWN, AlleleMatch.UNKNOWN));
    testMapBld.put(new String[][]{{"N", "N"}, {"A", "T"}}, new StrandOps.AlleleConfig(AlleleStatus.NULL, AlleleStatus.NULL, AlleleMatch.UNKNOWN, AlleleMatch.UNKNOWN));
    testMapBld.put(new String[][]{{"-", "-"}, {"A", "T"}}, new StrandOps.AlleleConfig(AlleleStatus.NULL, AlleleStatus.NULL, AlleleMatch.UNKNOWN, AlleleMatch.UNKNOWN));
    // valid:
    testMapBld.put(new String[][]{{"A", "G"}, {"A", "G"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF1, AlleleMatch.REF2));
    testMapBld.put(new String[][]{{"A", "G"}, {"G", "A"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF2, AlleleMatch.REF1));
    testMapBld.put(new String[][]{{"A", "G"}, {"T", "C"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF1_FLIP, AlleleMatch.REF2_FLIP));
    testMapBld.put(new String[][]{{"A", "G"}, {"C", "T"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF2_FLIP, AlleleMatch.REF1_FLIP));
    // ambiguous:
    testMapBld.put(new String[][]{{"A", "T"}, {"A", "T"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.AMBIGUOUS, AlleleMatch.AMBIGUOUS));
    testMapBld.put(new String[][]{{"T", "A"}, {"A", "T"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.AMBIGUOUS, AlleleMatch.AMBIGUOUS));
    // technically valid
    testMapBld.put(new String[][]{{"A", "G"}, {"A", "N"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF1, AlleleMatch.UNKNOWN_REF2));
    testMapBld.put(new String[][]{{"T", "G"}, {"A", "N"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF1_FLIP, AlleleMatch.UNKNOWN_REF2));
    testMapBld.put(new String[][]{{"T", "N"}, {"A", "N"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.NULL, AlleleMatch.REF1_FLIP, AlleleMatch.UNKNOWN));
    testMapBld.put(new String[][]{{"N", "T"}, {"A", "N"}}, new StrandOps.AlleleConfig(AlleleStatus.NULL, AlleleStatus.VALID, AlleleMatch.UNKNOWN, AlleleMatch.REF1_FLIP));
    testMapBld.put(new String[][]{{"N", "A"}, {"A", "N"}}, new StrandOps.AlleleConfig(AlleleStatus.NULL, AlleleStatus.VALID, AlleleMatch.UNKNOWN, AlleleMatch.REF1));
    testMapBld.put(new String[][]{{"N", "T"}, {"A", "N"}}, new StrandOps.AlleleConfig(AlleleStatus.NULL, AlleleStatus.VALID, AlleleMatch.UNKNOWN, AlleleMatch.REF1_FLIP));
    // invalid
    testMapBld.put(new String[][]{{"A", "G"}, {"T", "G"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF1_FLIP, AlleleMatch.REF2));
    testMapBld.put(new String[][]{{"A", "G"}, {"C", "G"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.NO_MATCH, AlleleMatch.AMBIGUOUS));
    testMapBld.put(new String[][]{{"A", "C"}, {"C", "G"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.NO_MATCH, AlleleMatch.AMBIGUOUS));
    testMapBld.put(new String[][]{{"A", "C"}, {"T", "C"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF1_FLIP, AlleleMatch.REF2));
    // invalid ref
    testMapBld.put(new String[][]{{"A", "T"}, {"T", "T"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.INVALID_REF, AlleleMatch.INVALID_REF));
    testMapBld.put(new String[][]{{"A", "T"}, {"A", "A"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.INVALID_REF, AlleleMatch.INVALID_REF));
    testMapBld.put(new String[][]{{"A", "T"}, {"T", "Q"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.INVALID_REF, AlleleMatch.INVALID_REF));
    testMapBld.put(new String[][]{{"A", "C"}, {"T", "Q"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.INVALID_REF, AlleleMatch.INVALID_REF));
    testMapBld.put(new String[][]{{"A", "Q"}, {"T", "Q"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.INVALID, AlleleMatch.INVALID_REF, AlleleMatch.INVALID_REF));
    testMapBld.put(new String[][]{{"A", "V"}, {"T", "Q"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.INVALID, AlleleMatch.INVALID_REF, AlleleMatch.INVALID_REF));
    testMapBld.put(new String[][]{{"J", "A"}, {"T", "T"}}, new StrandOps.AlleleConfig(AlleleStatus.INVALID, AlleleStatus.VALID, AlleleMatch.INVALID_REF, AlleleMatch.INVALID_REF));
    testMapBld.put(new String[][]{{"J", "A"}, {"J", "A"}}, new StrandOps.AlleleConfig(AlleleStatus.INVALID, AlleleStatus.VALID, AlleleMatch.INVALID_REF, AlleleMatch.INVALID_REF));
    // invalid input
    testMapBld.put(new String[][]{{"J", "A"}, {"N", "A"}}, new StrandOps.AlleleConfig(AlleleStatus.INVALID, AlleleStatus.VALID, AlleleMatch.UNKNOWN, AlleleMatch.REF2));
    testMapBld.put(new String[][]{{"F", "R"}, {"T", "G"}}, new StrandOps.AlleleConfig(AlleleStatus.INVALID, AlleleStatus.INVALID, AlleleMatch.UNKNOWN, AlleleMatch.UNKNOWN));
    testMapBld.put(new String[][]{{"H", "M"}, {"T", "G"}}, new StrandOps.AlleleConfig(AlleleStatus.INVALID, AlleleStatus.INVALID, AlleleMatch.UNKNOWN, AlleleMatch.UNKNOWN));
    // I/D
    testMapBld.put(new String[][]{{"I", "D"}, {"I", "D"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF1, AlleleMatch.REF2));
    testMapBld.put(new String[][]{{"I", "D"}, {"D", "I"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF2, AlleleMatch.REF1));
    testMapBld.put(new String[][]{{"I", "D"}, {"I", "I"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.INVALID_REF, AlleleMatch.INVALID_REF));
    testMapBld.put(new String[][]{{"I", "D"}, {"D", "D"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.INVALID_REF, AlleleMatch.INVALID_REF));
    testMapBld.put(new String[][]{{"A", "AGG"}, {"A", "AGG"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF1, AlleleMatch.REF2));
    testMapBld.put(new String[][]{{"A", "AGG"}, {"AGG", "A"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF2, AlleleMatch.REF1));
    testMapBld.put(new String[][]{{"A", "AGG"}, {"T", "TCC"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF1_FLIP, AlleleMatch.REF2_FLIP));
    testMapBld.put(new String[][]{{"A", "AGG"}, {"TCC", "T"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF2_FLIP, AlleleMatch.REF1_FLIP));
    // multi-allele N's
    testMapBld.put(new String[][]{{"A", "AGN"}, {"A", "AGN"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF1, AlleleMatch.UNKNOWN_REF2));
    testMapBld.put(new String[][]{{"A", "AGT"}, {"A", "AGN"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF1, AlleleMatch.UNKNOWN_REF2));
    testMapBld.put(new String[][]{{"A", "AGN"}, {"A", "AGT"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF1, AlleleMatch.UNKNOWN_REF2));
    testMapBld.put(new String[][]{{"A", "ANN"}, {"A", "AGT"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF1, AlleleMatch.UNKNOWN_REF2));
    testMapBld.put(new String[][]{{"A", "AGN"}, {"TCC", "T"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.VALID, AlleleMatch.REF2_FLIP, AlleleMatch.UNKNOWN_REF1));
    // multi-allele invalid
    testMapBld.put(new String[][]{{"A", "N"}, {"TCC", "T"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.NULL, AlleleMatch.REF2_FLIP, AlleleMatch.UNKNOWN_REF2));
    testMapBld.put(new String[][]{{"N", "A"}, {"TCC", "T"}}, new StrandOps.AlleleConfig(AlleleStatus.NULL, AlleleStatus.VALID, AlleleMatch.UNKNOWN_REF2, AlleleMatch.REF2_FLIP));
    testMapBld.put(new String[][]{{"A", "N"}, {"G", "T"}}, new StrandOps.AlleleConfig(AlleleStatus.VALID, AlleleStatus.NULL, AlleleMatch.REF2_FLIP, AlleleMatch.UNKNOWN));    
    
    ImmutableMap<String[][], StrandOps.AlleleConfig> testMap = testMapBld.build();
    for (Entry<String[][],  StrandOps.AlleleConfig> ent : testMap.entrySet()) {
      System.out.println(ArrayUtils.toStr(ent.getKey()[0], ",") + " = " + ArrayUtils.toStr(ent.getKey()[1], ","));
      StrandOps.AlleleConfig conf = StrandOps.strandConfig(ent.getKey()[0], ent.getKey()[1]);
      assertEquals(ent.getValue().allele1, conf.allele1);
      assertEquals(ent.getValue().allele2, conf.allele2);
      assertEquals(ent.getValue().match1, conf.match1);
      assertEquals(ent.getValue().match2, conf.match2);
    }
  }

}
