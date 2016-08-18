// -Xms1024M -Xmx1024M
package org.genvisis.park;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Hashtable;

import org.genvisis.common.Internat;
import org.genvisis.common.ext;

public class temp {
  public temp(String[] args) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line;
    String temp, trav, filename, dir;
    Hashtable<String, String> hash;
    // Vector<String> v = new Vector<String>();
    int count;
    long time;


    for (int i = 219; i <= 219; i++) {
      Internat.downloadFile("http://www.eperc.mcw.edu/FileLibrary/User/jrehm/fastfactpdfs/Concept"
                            + ext.formNum(i, 3) + "2.pdf", "Concept" + ext.formNum(i, 3) + ".pdf");
    }

    System.exit(1);

    dir = filename = trav = "";
    time = new Date().getTime();
    hash = new Hashtable<String, String>();
    count = hash.size();
    System.out.println(count);
    try {
      reader = new BufferedReader(new FileReader(dir + filename));
      writer = new PrintWriter(new FileWriter(trav));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        temp = line[0];
        count = temp.length();
      }
      writer.close();
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + filename + "\"");
      System.exit(2);
    }
    System.out.println(" in " + ext.getTimeElapsed(time));
  }

  public static void main(String[] args) throws IOException {
    try {
      new temp(args);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}

// for moving MACH files between systems...
// for (int chr = 1; chr<=22; chr++) {
// System.out.println("gunzip -c chr"+chr+".tar.gz | tar xvf -");
// }
//
// System.exit(1);

// example of each different Sort order...
// int[] positions = new int[] {5, 3, 4, 1, 6};
// int[] keys = Sort.quicksort(positions);
// int[] traditional = new int[positions.length];
// int[] reverse = new int[positions.length];
//
// for (int i = 0; i<positions.length; i++) {
// traditional[i] = positions[keys[i]];
// reverse[keys[i]] = positions[i];
// }
//
// System.out.println(Array.toStr(traditional));
// System.out.println(Array.toStr(reverse));
//
// System.exit(1);

// might yet be useful with PD GWAS consortium
// line = new String[] {"CIDR", "Fung", "Sing550", "LEAPS", "Miami", "NGRC"};
// String[] suffi = new String[] {"final", "pheno", "pheno", "pheno", "final", "final"};
//
// for (int i = 0; i<line.length; i++) {
// System.out.println("cd imputation"+line[i]+"/");
// System.out.println("jcp gwas.Probabel pheno="+line[i]+"_Aff_"+suffi[i]+".dat -parse");
// System.out.println("gzip "+line[i]+"_Aff_"+suffi[i]+"_add.xln");
// System.out.println("mv "+line[i]+"_Aff_"+suffi[i]+"_add.xln.gz ..");
// System.out.println("jcp gwas.Probabel pheno="+line[i]+"_AAO_"+suffi[i]+".dat -parse");
// System.out.println("gzip "+line[i]+"_AAO_"+suffi[i]+"_add.xln");
// System.out.println("mv "+line[i]+"_AAO_"+suffi[i]+"_add.xln.gz ..");
// System.out.println("cd ..");
//
// System.out.println();
// }
//
// System.exit(1);



// prove the same birthday in a room thing
// int num = 23, reps = 10000000, count = 0;
// int[] done = new int[num];
// boolean match;
//
// for (int i = 0; i<reps; i++) {
// match = false;
// for (int j = 0; j<num; j++) {
// done[j] = (int)(Math.random()*365);
// for (int k = 0; k<j; k++) {
// if (done[j] == done[k]) {
// match = true;
// }
// }
// }
// if (match) {
// count++;
// }
// }
// System.out.println("Same birthday among "+num+" people, "+reps+" reps, has a
// probability of "+ext.formDeci((double)count/(double)reps*100, 1, true)+"%");
