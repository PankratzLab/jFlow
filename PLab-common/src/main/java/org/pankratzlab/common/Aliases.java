package org.pankratzlab.common;

public class Aliases {

  public static final String[] MARKER_NAMES = {"MarkerName", "Marker", "Name", "SNP", "SNP Name",
                                               "SNP_Name", "SNPID", "SNP.id", "rs_id", "rsID",
                                               "Probe Set ID", "ProbeSetName", "ProbeSet",
                                               "Variant", "VariantName", "AnalysisUnit",
                                               "Variant_ID", "SingleVariant", "BurdenTests", "RSID",
                                               "BinName", "dbSNP ID"};
  public static final String[] GENE_UNITS = {"Gene", "SKATgene"};

  public static final String[] REGION = {"Region", "UCSC", "Band", "Arm"};
  public static final String[] CHRS = {"Chr", "Chromosome", "CHROM", "#CHROM"};
  public static final String[] POSITIONS = {"Position", "position", "pos", "Pos", "POS", "BP",
                                            "MapInfo", "PositionOfFirstMarkerInGene",
                                            "Physical Position"};
  public static final String[] POSITIONS_START = {"Start", "Begin"};
  public static final String[] POSITIONS_STOP = {"Stop", "End", "Stop Position"};
  public static final String[] CENTIMORGANS = {"centiMorgans", "cM"};

  public static final String[][] ALLELES = {{"coded_all", "A1", "Al1", "Allele1", "ALT",
                                             "Effect_allele", "EA", "Coded.Allele"},
                                            {"noncoded_all", "A2", "Al2", "Allele2", "REF", "OTHER",
                                             "Reference_allele", "NEA", "Other_allele",
                                             "NonCoded.Allele", "NON_EFFECT_ALLELE"},};
  public static final String[] REF_ALLELES = {"REF", "Ref", "ref", "RefAllele", "refAllele", "a0",
                                              "Allele1"};
  public static final String[] ALT_ALLELES = {"ALT", "Alt", "alt", "AltAllele", "altAllele", "a1",
                                              "Allele2"};
  public static final String[] EFFECTS = {"beta", "beta_SNP_add", "Effect", "BETA_SNP"};
  public static final String[] STD_ERRS = {"se", "StdErr", "sebeta_SNP_add", "ROBUST_SE_SNP"};

  public static final String[] PVALUES = {"pval", "P", "p-val", "p-value", "Pvalue", "mbpval",
                                          "minPval", "ROBUST_PVAL_SNP", "chi_P_2df", "EMP1"};
  public static final String[] NS = {"N", "NMISS", "sampleN", "n_total", "ESampleSize"};

  public static final String[] ALLELE_FREQS = {"freq", "AlleleFreq", "A1Freq", "AF", "AAF", "MAF",
                                               "sampleMAF", "Effect_allele_frequency",
                                               "effect_allele_freq", "EAF", "FRQ", "Coded.Af",
                                               "AF_coded_all"};

  public static final String[] IMPUTATION_EFFICIENCY = {"imp_info", "rsq", "Qual"};

  public static final String[] REFERENCE_FOLDERS = {"C:/bin/NCBI/", "G:/My Drive/statgen/NCBI/",
                                                    "N:/statgen/NCBI/",
                                                    "/panfs/roc/groups/5/pankrat2/public/bin/NCBI/",
                                                    "/home/npankrat/NCBI/",
                                                    "/panfs/roc/groups/5/pankrat2/public/bin/"};

  public static final String[] INDIVIDUAL_ID = {"ID", "IID", "I_ID", "IndID", "Ind_ID"};
  public static final String[] FAMILY_ID = {"Family ID", "FamID", "FID", "F_ID"};
  public static final String[] DNA = {"DNA/Sample", "DNA", "DNA#", "Sample", "LabID", "SampleID",
                                      "Sample ID"};

  /**
   * Searches all of the reference directories to see if it contains the specified file
   *
   * @param filename the filename to search for
   * @param verbose whether to report an error if none of the locations exists or if none of the
   *          locations contains the specified file
   * @param log a Logger to report any errors if verbose
   * @return the full path to the file of interest if it exists in one of the directories, otherwise
   *         null
   */
  public static String getPathToFileInReferenceDirectory(String filename, boolean verbose,
                                                         Logger log) {
    // TODO allow customizable paths from some future Genvisis global properties file
    return Files.firstPathToFileThatExists(REFERENCE_FOLDERS, filename, verbose, false, log);
  }

}
