package org.genvisis.one.JL.topMed;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;

import com.google.common.collect.ImmutableSet;

public class ParseDupGroups {

  public static void main(String[] args) {

    String dupFile = "/Volumes/Beta2/NGS/topmed/Duplicate_list_PEIZO_genotypes.dups.txt";

    String[][] dups = HashVec.loadFileToStringMatrix(dupFile, true, new int[] {0, 1}, "\t", 100,
                                                     false);

    Map<String, Set<String>> duplicates = new HashMap<>();
    for (String[] dup : dups) {
      Set<String> dupeSet = new HashSet<>();
      dupeSet.add(dup[0]);
      dupeSet.add(dup[1]);
      // Merge with any pre-existing sets
      if (duplicates.containsKey(dup[0])) {
        dupeSet.addAll(duplicates.remove(dup[0]));
      }
      if (duplicates.containsKey(dup[1])) {
        dupeSet.addAll(duplicates.remove(dup[1]));
      }
      // Update the mapping for all duplicates in this set
      for (String duplicate : dupeSet) {
        duplicates.put(duplicate, dupeSet);
      }
    }
    for (String key : duplicates.keySet()) {
      System.out.println(duplicates.get(key).size());
    }

    Collection<Set<String>> uniqueDuplicateSets = ImmutableSet.copyOf(duplicates.values());
    System.out.println(uniqueDuplicateSets.size());
    try {
      PrintWriter writerSet = Files.openAppropriateWriter(dupFile + "_duplicatesSet.dat");
      writerSet.println("SAMPLE\tSET\tSET_SIZE");

      PrintWriter writer = Files.openAppropriateWriter(dupFile + "_duplicates.dat");
      int set = 1;
      for (Set<String> dupeSet : uniqueDuplicateSets) {
        for (String fidIid : dupeSet) {
          writerSet.println(fidIid + "\t" + set + "\t" + dupeSet.size());
          writer.print(fidIid);
          writer.print("\t");

        }
        writer.print("\t");
        writer.print(set);

        writer.println();
        set++;
      }
      writerSet.close();
      writer.close();
    } catch (IOException e) {
      System.err.println("Error writing to duplicateSets.txt file");
      e.printStackTrace();
    }
  }

}
