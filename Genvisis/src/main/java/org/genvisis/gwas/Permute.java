package org.genvisis.gwas;

public class Permute {


  public static void permute(String phenoFile) {

  }



  public static void main(String[] args) {
    int numArgs = args.length;
    String pheno = "pheno.dat";

    String usage = "\n" + "gwas.Permute requires 0-1 arguments\n" + "   (1) filename (i.e. pheno="
                   + pheno + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        pheno = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      permute(pheno);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
