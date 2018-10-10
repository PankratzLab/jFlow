package org.pankratzlab.common;


public class VCFUtils {

    public static final String GZ = ".gz";
    public static final String GZ_INDEX = ".tbi";

    public static final String VCF = ".vcf";
    public static final String VCF_INDEX = ".idx";

    /**
     * @param filename
     * @return
     */
    public static final String getVcfIndex(String filename) {
      if (filename.endsWith(GZ)) {
        return filename + GZ_INDEX;
      } else {
        return filename + VCF_INDEX;
      }
    }
}
