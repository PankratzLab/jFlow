package org.genvisis.park;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.Files;

public class batch {
  public static int[] maxes = {281, 253, 216, 202, 193, 180, 174, 153, 151, 168, 144, 164, 104, 127,
      106, 115, 125, 126, 89, 95, 36, 47, 176};

  public static void main(String[] args) throws IOException {
    if (args.length == 0) {
      System.out.println("USAGE: java park.batch [case] [run] [bs]");
      System.out
          .println("(optional) to run the batch file add \"run\" to the end of the command line");
      System.out.println("(optional) number of batch files (i.e. nb=3 to create 3 files)");
      System.out.println("(optional) trait to run (i.e. trait=AOO (default))");
      System.out.println("case  0 = create Filesystem");
      System.out.println("case  1 = create .pre files");
      System.out.println("case  2 = setup for Relpair");
      System.out.println("case  3 = create loc files and run MAPMAKER/SIBS");
      System.out.println("case  4 = create loc files and run MAPMAKER/SIBS using 2pt");
      System.out.println("case  5 = create files and run ASPEX");
      System.out.println("case  6 = create files and run Genehunter 2");
      System.out.println("case  7 = edit map files and run MLINK");
      System.out.println("case  8 = create files and run Allegro");
      System.out
          .println("case  9 = create files and run Autosomal Dominant model [w/ and w/o het]");
      System.out
          .println("case 10 = create files and run Autosomal Recessive model  [w/ and w/o het]");
      System.out.println("case 11 = create loc files and run MAPMAKER/SIBS with aoo");
      System.out.println("case 13 = create files and run Solar with aoo");
      System.out.println("case 14 = run Solar with a different phenotype (file=pheno.dat)");
      System.out.println("case 16 = create files and run all Genehunter 2.1 analyses with aoo");
      System.out.println("case 17 = create files and run Merlin");
      System.out.println("case 18 = simulate null hypthesis for MAPMAKER/SIBS");
      System.out.println("case 19 = batch for Ordered Subset Analysis using MAPMAKER/SIBS");
      System.out.println("case 20 = batch for Ordered Subset Analysis using GenehunterPlus");
      System.out.println();
    } else {
      try {
        new batch(args);
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
  }

  public int[] chrs =
      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};

  public batch(String[] arguments) throws IOException {
    PrintWriter writer;
    PrintWriter[] writers = null;
    String[] line;
    String chrome;
    int step = Integer.valueOf(arguments[0]).intValue();
    boolean run = false;
    int start = 1, stop = 23;
    BufferedReader reader;
    String filename = "solar.ptypes", trait = "AOO";
    String classpath = "";

    int numBatches = step == 13 || step == 18 ? 4 : 1;

    for (int i = 1; i < arguments.length; i++) {
      if (arguments[i].equals("run")) {
        run = true;
      } else if (arguments[i].startsWith("trait=")) {
        trait = arguments[i].split("=")[1];
      } else if (arguments[i].startsWith("nb=")) {
        numBatches = Integer.parseInt(arguments[i].split("=")[1]);
      } else {
        try {
          if ((step == 13 || step == 14) && arguments[i].startsWith("file=")) {
            filename = arguments[i].split("=")[1];
            if (!new File(filename).exists()) {
              System.err.println("Error - could not find " + filename + " in current directory");
              System.exit(2);
            }
            try {
              reader = new BufferedReader(new FileReader(filename));
              line = reader.readLine().split(",");
              if (line.length < 3) {
                System.err.println("Error - " + filename + " is not comma delimited");
                System.exit(11);
              } else {
                trait = line[2];
              }
              reader.close();
            } catch (IOException ioe) {
              System.err.println("Error parsing " + filename);
            }
          } else if (start == 1 && stop == 23) {
            start = stop = Integer.valueOf(arguments[i]).intValue();
          } else if (start == stop) {
            stop = Integer.valueOf(arguments[i]).intValue();
          } else {
            System.err.println("Error: way too many arguments");
          }
        } catch (Exception e) {
          System.err.println("Error: invalid arguments (valid = command [chromosome #] [run]");
        }
      }
    }

    writers = new PrintWriter[numBatches];
    System.out.println(new File(".").getAbsolutePath());
    classpath = "-classpath /home/npankrat/" + org.genvisis.common.PSF.Java.GENVISIS + "";
    if (new File(".").getAbsolutePath().contains("bc2/pankratz")) {
      classpath = "-classpath /home/bc2/pankratz/" + org.genvisis.common.PSF.Java.GENVISIS + "";
    }

    for (int i = 0; i < numBatches; i++) {
      writers[i] = new PrintWriter(new FileWriter(i == 0 && numBatches == 1
          ? (Files.isWindows() ? "batch.bat" : "batch") : "batch." + (i + 1)));
      writers[i].println("#/bin/sh\n");
      if (new File(".").getAbsolutePath().contains("bc2/pankratz")) {
        writers[i].println("module load java\n");
      }
      if (step == 13 || step == 14 || step == 18) {
        writers[i].println("sleep 20\n");
      }
    }

    for (int i = start; i <= stop; i++) {
      chrome = (i < 10) ? "0" + i : "" + i;

      writer = writers[i % numBatches];

      switch (step) {
        case 0:
          writer.println("java " + classpath + " link.bat.Filesystem " + i);
          writer.println("echo completed filesystem for chromosome " + i);
          writer.println();
          break;
        case 1:
          writer.println("java " + classpath + " link.LinkageFormat chr=" + i);
          writer.println("echo made pre file for chromosome " + i);
          writer.println();
          break;
        case 2:
          // if (i != 23) {
          // writer.println("echo -e
          // \"re_chrom"+chrome+".pre\\nchrom"+chrome+".ped\\nn\\nn\\nn\\nn\\n0\\nn\\n\"
          // | /usr/local/bin/makeped");
          // writer.println("java "+classpath+" park.bat.borel "+i);
          // writer.println();
          // }
          if (i == start) {
            writer.println("java " + classpath + " park.bat.relpair");
            writer.println();
          }
          break;
        case 3:
          writer.println("java " + classpath + " park.bat.dat2loc map" + chrome + ".dat");
          if (i != 23) {
            writer.println("echo -e \"pairs\\n3\\nload map" + chrome + ".loc\\nprep re_chrom"
                + chrome + ".pre\\nn\\nscan\\ninfo\\ninfo" + chrome + ".out\\ninfo" + chrome
                + ".ps\\nestimate\\ny\\nchrom" + chrome + "-mls.out\\nchrom" + chrome
                + "-share.ps\\nchrom" + chrome + "-mls.ps\\nestimate\\nn\\nchrom" + chrome
                + "-Dv-mls.out\\nchrom" + chrome + "-Dv-share.ps\\nchrom" + chrome
                + "-Dv-mls.ps\\nquit\\n\" | /software/bin/sibs > chrom" + chrome + ".log");
          } else {
            writer.println("echo -e \"sex on\\npairs\\n3\\nload map" + chrome
                + ".loc\\nprep re_chrom" + chrome + ".pre\\nn\\nscan\\ninfo\\ninfo" + chrome
                + ".out\\ninfo" + chrome + ".ps\\nestimate\\nchrom" + chrome + "-mls.out\\nchrom"
                + chrome + "-multi-mls.ps\\nchrom" + chrome
                + "-mls.ps\\nquit\\n\" | /software/bin/sibs > chrom" + chrome + ".log");
          }
          writer.println();
          break;
        case 4:
          writer.println("java " + classpath + " park.bat.dat2loc map" + chrome + ".dat");
          if (i != 23) {
            writer.println("echo -e \"pairs\\n3\\nload map" + chrome + ".loc\\nprep re_chrom"
                + chrome + ".pre\\nn\\nsingle point on\\nscan\\nestimate\\ny\\nchrom" + chrome
                + "-2pt.out\\nestimate\\nn\\nchrom" + chrome
                + "-Dv-2pt.out\\nquit\\n\" | /software/bin/sibs > chrom" + chrome + ".log");
          } else {
            writer
                .println("echo -e \"sex on\\npairs\\n3\\nload map" + chrome + ".loc\\nprep re_chrom"
                    + chrome + ".pre\\nn\\nsingle point on\\nscan\\nestimate\\nchrom" + chrome
                    + "-2pt.out\\nquit\\n\" | /software/bin/sibs > chrom" + chrome + ".log");
          }
          writer.println();
          break;
        case 5:
          writer.println("java " + classpath + " park.bat.createAspex " + i);
          // if (i == 23) {
          // writer.println("cat chrom23h.param sex_linked > temp");
          // writer.println("mv temp chrom23h.param");
          // writer.println("cat chrom23hd.param sex_linked > temp");
          // writer.println("mv temp chrom23hd.param");
          // writer.println("cat chrom23k.param sex_linked > temp");
          // writer.println("mv temp chrom23k.param");
          // writer.println("cat chrom23kd.param sex_linked > temp");
          // writer.println("mv temp chrom23kd.param");
          // writer.println("rm sex_linked");
          // }
          writer.println("/opt/local/software/aspex-2.3-solaris/sib_phase -v -f chrom" + chrome
              + "k-mpt.param chrom" + chrome + ".pre > chrom" + chrome + "-K-mpt.out");
          writer.println("/opt/local/software/aspex-2.3-solaris/sib_phase -f chrom" + chrome
              + "k-2pt.param chrom" + chrome + ".pre > chrom" + chrome + "-K-2pt.out");
          writer.println();
          break;
        case 6:
          writer.println("echo -e \"analysis npl\\npost on\\nload map" + chrome
              + ".dat\\nscan re_chrom" + chrome + ".pre\\ntotal\\nchrom" + chrome
              + "-npl.ps\\nchrom" + chrome + "-lod.ps\\nchrom" + chrome
              + "-info.ps\\nquit\\n\" | /opt/local/software/gh2/gh.sol > chrom" + chrome + ".out");
          writer.println();
          break;
        case 7:
          writer.println("makeped re_chrom" + chrome + ".pre re_chrom" + chrome + ".ped n");
          writer.println("cp re_chrom" + chrome + ".ped pedin.dat");
          writer.println("cp map" + chrome + ".dat datain.dat");
          writer.println("java " + classpath + " park.bat.makeMap4MLINK 0.005 0.80 0.80 0.03");

          writer.println();
          break;
        case 8:
          writer.println("java " + classpath + " link.Allegro chr=" + i + " model=NPL");
          writer.println("allegro useful" + chrome + ".opt");
          // writer.println("find . -name \"lqhat*\" -exec rm {} \\;");
          writer.println();
          break;
        case 9:
          writer.println("java " + classpath + " link.Allegro chr=" + i + " model=DOM");
          writer.println("allegro useful" + chrome + ".opt");
          // writer.println("find . -name \"lqhat*\" -exec rm {} \\;");
          writer.println();
          break;
        case 10:
          writer.println("java " + classpath + " link.Allegro chr=" + i + " model=REC");
          writer.println("allegro useful" + chrome + ".opt");
          // writer.println("find . -name \"lqhat*\" -exec rm {} \\;");
          writer.println();
          break;
        case 11:
          writer.println("java " + classpath + " park.bat.createLinkage " + i + " -pheno");
          writer.println("java " + classpath + " park.bat.dat2loc map" + chrome + ".dat");
          if (i != 23) {
            writer.println("echo -e \"pairs\\n3\\nload map" + chrome + ".loc\\nprep re_chrom"
                + chrome + ".pre\\ny\\npheno.sibs\\nscan\\np\\nml\\nchrom" + chrome
                + "-ml.out\\nchrom" + chrome + "-ml.ps\\nhase\\nchrom" + chrome
                + "-trad-he.out\\nchrom" + chrome + "-em-he.out\\nchrom" + chrome
                + "-he.ps\\nquit\\n\" | /software/bin/sibs > chrom" + chrome + ".log");
          } else {
            writer.println("echo -e \"sex on\\npairs\\n3\\nload map" + chrome
                + ".loc\\nprep re_chrom" + chrome + ".pre\\ny\\npheno.sibs\\nscan\\np\\ninfo\\ninfo"
                + chrome + ".out\\ninfo" + chrome + ".ps\\nml\\nchrom" + chrome + "-ml.out\\nchrom"
                + chrome + "-ml.ps\\nhase\\nchrom" + chrome + "-trad-he.out\\nchrom" + chrome
                + "-em-he.out\\nchrom" + chrome + "-he.ps\\nquit\\n\" | /software/bin/sibs > chrom"
                + chrome + ".log");
          }
          writer.println();
          break;
        case 12:

        case 13:
          if (!chrome.equals("23")) {
            writer.println("mkdir chrom" + chrome);
            writer.println(
                "java " + classpath + " park.bat.createSolar chr=" + i + " trait=" + trait);
            // writer.println("java "+classpath+" park.zeroByInd
            // filter.pre "+i);
            writer.println("cp solar.fam chrom" + chrome);
            writer.println("mv solar.map." + i + " chrom" + chrome);
            writer.println("mv solar.freqs." + i + " chrom" + chrome);
            writer.println("mv solar.gtypes." + i + " chrom" + chrome);
            writer.println("cp " + filename + " chrom" + chrome);
            writer.println("cd chrom" + chrome);
            writer.println("echo -e \"load pedigree solar.fam\\nload freq solar.freqs." + i
                + "\\nload marker solar.gtypes." + i
                + "\\nibddir .\\nverbosity min\\nibd\\nload map solar.map." + i
                + "\\nibddir .\\nmibddir .\\nmibd 0 " + maxes[i - 1] + " 1\\nmibddir .\\nautomodel "
                + filename + " " + trait + "\\npolygenic -screen\\nmibddir .\\nchromosome " + i
                + "\\ninterval 1\\nmultipoint -overwrite\\nquit\\n\" | solar > solar.log");
            writer.println("cd ..");
            writer.println();
          }
          break;
        case 14:
          if (!chrome.equals("23")) {
            writer.println("cp " + filename + " chrom" + chrome);
            writer.println("cd chrom" + chrome);
            writer.println("echo -e \"ibddir .\\nmibddir .\\nautomodel " + filename + " " + trait
                + "\\npolygenic -screen\\nmibddir .\\nchromosome " + i
                + "\\ninterval 1\\nmultipoint -overwrite\\nquit\\n\" | solar > solar.log");
            writer.println("cd ..");
            writer.println();
          }
          break;
        case 15:
          writer.println("java " + classpath + " park.bat.dat2loc map" + chrome + ".dat");
          if (i != 23) {
            writer.println("echo -e \"pairs\\n2\\nload map" + chrome + ".loc\\nprep re_chrom"
                + chrome + ".pre\\nn\\nscan\\ninfo\\ninfo" + chrome + ".out\\ninfo" + chrome
                + ".ps\\nestimate\\ny\\nchrom" + chrome + "-mls.out\\nchrom" + chrome
                + "-share.ps\\nchrom" + chrome + "-mls.ps\\nestimate\\nn\\nchrom" + chrome
                + "-Dv-mls.out\\nchrom" + chrome + "-Dv-share.ps\\nchrom" + chrome
                + "-Dv-mls.ps\\nquit\\n\" | /software/bin/sibs > chrom" + chrome + ".log");
          } else {
            writer.println("echo -e \"sex on\\npairs\\n3\\nload map" + chrome
                + ".loc\\nprep re_chrom" + chrome + ".pre\\nn\\nscan\\ninfo\\ninfo" + chrome
                + ".out\\ninfo" + chrome + ".ps\\nestimate\\nchrom" + chrome + "-mls.out\\nchrom"
                + chrome + "-multi-mls.ps\\nchrom" + chrome
                + "-mls.ps\\nquit\\n\" | /software/bin/sibs > chrom" + chrome + ".log");
          }
          writer.println();
          break;
        case 16:
          writer.println("java " + classpath + " park.bat.addPheno2gh " + i + " 8");
          writer.println(
              "echo -e \"load gh_map" + chrome + ".dat\\nincrement distance 1\\nscan gh_chrom"
                  + chrome + ".pre\\nvariance\\nn\\nn\\nvc" + chrome + ".out\\ncorr" + chrome
                  + ".out\\nn\\nml\\nml" + chrome + ".out\\nhase\\ntrad-he" + chrome
                  + ".out\\nem-he" + chrome + ".out\\nnon\\nnp.out\\nno dom\\nml_ndv" + chrome
                  + ".out\\nquit\\n\" | /home/npankrat/bin/gh2 > /dev/null");
          // writer.println("echo -e \"load
          // gh_map"+chrome+".dat\\nincrement distance 1\\nscan
          // gh_chrom"+chrome+".pre\\nvariance\\ny\\nn\\nvc_y_n"+chrome+".out\\ncorr_y_n"+chrome+".out\\nn\\nvariance\\nn\\ny\\nvc_n_y"+chrome+".out\\ncorr_n_y"+chrome+".out\\nn\\nquit\\n\"
          // | /home/npankrat/bin/gh2 > /dev/null");

          writer.println();
          break;
        case 17:
          writer.println("java " + classpath + " link.bat.createMerlin " + i);
          if (i != 23) {
            writer.println("merlin -d merlin_data." + chrome + " -p re_chrom" + chrome
                + ".pre -m merlin_map." + chrome + " --npl --grid 1 > chrom" + chrome + ".out");
          } else {
            writer.println("minx -d merlin_data." + chrome + " -p re_chrom" + chrome
                + ".pre -m merlin_map." + chrome + " --npl --grid 1 > chrom" + chrome + ".out");
          }
          writer.println();
          break;
        case 18:
          writer.println("mkdir sims" + chrome);
          writer.println("cp re_chrom" + chrome + ".pre sims" + chrome);
          writer.println("mv map" + chrome + ".dat sims" + chrome);
          writer.println("cd sims" + chrome);
          writer.println("jcp empiricalPvalues " + i + " 2000");
          writer.println("cd ..");
          writer.println();
          break;
        case 19:
          writer
              .println("java " + classpath + " park.bat.osa file=trait.dat chr=" + i + " dv=true");
          writer.println("cd chrom" + chrome);
          writer.println("sleep 5");
          writer.println("./batch");
          writer.println("cd ..");
          writer.println("java " + classpath + " park.bat.gatherOSA " + i);
          writer.println();
          break;
        case 20:
          writer.println("java " + classpath + " park.bat.osaHauser trait.dat " + i);
          writer.println("cp trait.dat chrom" + chrome + "/");
          writer.println("cp map" + chrome + ".dat chrom" + chrome + "/");
          writer.println("cp re_chrom" + chrome + ".pre chrom" + chrome + "/");
          writer.println("cd chrom" + chrome);
          writer.println("echo -e \"map kos\\nload map" + chrome + ".dat\\nscan re_chrom" + chrome
              + ".pre\\nquit\\n\" | /software/bin/ghp > /dev/null");
          writer.println("kaclod 200");
          writer.println("osa2.1 osa" + chrome + ".dat > osa" + chrome + ".log");
          writer.println("cd ..");
          writer.println();
          break;
        default:
          writer.println();
          break;
      }
    }

    if (step == 20) {
      writers[0].println("java " + classpath + " park.bat.gatherOSAhauser");
    }

    for (int i = 0; i < numBatches; i++) {
      writers[i].close();
      Files.chmod((i == 0 && numBatches == 1 ? "batch" : "batch." + (i + 1)));
    }
    if (numBatches == 1 && run) {
      Runtime.getRuntime().exec("./batch > log.log");
    }

  }
}
