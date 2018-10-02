package org.genvisis.one.JL;

import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

/**
 * @author Kitty One off for processing .metal files to .beta files
 */
public class processMetal {

  public static void main(String[] args) {
    String[] conv = new String[] {"rsID", "alt", "ref", "beta", "p"};
    String[] required = new String[] {"MarkerName", "Allele1", "Allele2", "Effect", "P.value"};
    String dir = "/Volumes/Beta/data/wbcGwasMetal/oldWbcMeta/";
    String[] metals = Files.listFullPaths(dir, ".metal");
    Logger log = new Logger();
    for (String metal : metals) {
      String[] header = Files.getHeaderOfFile(metal, log);
      int[] indices = ext.indexFactors(required, header, true);
      String[][] file = HashVec.loadFileToStringMatrix(metal, false, indices);
      file[0] = conv;
      Files.writeMatrix(file, ext.rootOf(metal, false) + ".beta", "\t");
    }

  }

}
