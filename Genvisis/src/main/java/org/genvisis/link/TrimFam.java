package org.genvisis.link;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class TrimFam {
  public static final int SCORE_99_NAMING_SCHEME = 0;
  public static final int ALPHA_NAMING_SCHEME = 1;
  public static final int NUMBER_NAMING_SCHEME = 2;
  public static final int DEFAULT_NAMING_SCHEME = NUMBER_NAMING_SCHEME;
  public static final String[] NAMING_SCHEME_KEY = {"99", "alpha", "num"};

  private Hashtable<String, Person> hash;
  private Vector<String> extended;
  private Vector<String> nukular;
  private Vector<String> unused;
  private boolean rename;
  private int namingScheme;
  private boolean hasSubfamilies;
  private int countGenoAndPheno;
  private int countAffectedPOpairs;
  private int favGen;
  private Logger log;
  private boolean error;

  public class Person extends Object {
    public String famid;
    public String id;
    public String recoded_id;
    public String dad;
    public String mom;
    public String sex;
    public boolean genotyped;
    public boolean phenotyped;
    public Person mother;
    public Person father;
    public Vector<Person> children;
    public boolean written;
    public boolean necessary;
    public int subFamily;
    public int generation;
    public double priorityScore;

    public Person(String newfamid, String newid, String dadsid, String momsid) {
      famid = newfamid;
      id = newid;
      dad = dadsid.equals("0") ? "root" : dadsid;
      mom = momsid.equals("0") ? "root" : momsid;
      recoded_id = isRoot() ? "0" : "-1";

      sex = "0";
      genotyped = false;
      phenotyped = false;
      generation = -99;

      mother = null;
      father = null;
      children = new Vector<Person>();

      written = false;
      necessary = false;

      subFamily = 0;
    }

    public boolean isRoot() {
      return id.equals("root");
    }

    public boolean isMale() {
      return (sex.equals("M") || sex.equals("1"));
    }

    public boolean isFemale() {
      return (sex.equals("F") || sex.equals("2"));
    }

    public boolean defineNecessaryPath() {
      if (necessary) {
        return true;
      }

      for (int i = 0; i < children.size(); i++) {
        if (children.elementAt(i).defineNecessaryPath()) {
          necessary = true;
        }
      }

      if (genotyped && phenotyped) {
        necessary = true;
        countGenoAndPheno++;
        if ((father.genotyped && father.phenotyped) || (mother.genotyped && mother.phenotyped)) {
          countAffectedPOpairs++;
        }
      }

      return necessary;
    }

    public boolean trimUselessDescendents() {
      for (int i = 0; i < children.size(); i++) {
        if (children.elementAt(i).trimUselessDescendents()) {
          necessary = true;
        }
      }

      if (!isRoot()) {
        if (genotyped && father.necessary && !father.genotyped) {
          necessary = true;
        }

        if (genotyped && mother.necessary && !mother.genotyped) {
          necessary = true;
        }

        if (!necessary) {
          hash.remove(id);
        }
      }

      return necessary;
    }

    public void spreadSubFamID(int subFam) {
      if (subFamily > 0) {
        return;
      } else {
        subFamily = subFam;

        for (int i = 0; i < children.size(); i++) {
          children.elementAt(i).spreadSubFamID(subFam);
        }

        if (!dad.equals("root")) {
          father.spreadSubFamID(subFam);
          mother.spreadSubFamID(subFam);
        }
        return;
      }
    }

    public int spreadGeneration(int gen, boolean commit) {
      int min;

      min = gen;
      if (generation == -99 || (commit && generation <= -98)) {
        if (commit) {
          generation = gen;
        } else {
          generation = -98;
        }

        for (int i = 0; i < children.size(); i++) {
          min = children.elementAt(i).spreadGeneration(gen - 1, commit);
        }

        if (!dad.equals("root")) {
          father.spreadGeneration(gen + 1, commit);
        }
        if (!mom.equals("root")) {
          mother.spreadGeneration(gen + 1, commit);
        }
      }

      return Math.min(min, gen);
    }

    // would checking for presence speed this up?
    public void removeParentalLines() {
      if (isRoot() || !hash.containsKey(id)) {
        return;
      }

      father.removeParentalLines();
      mother.removeParentalLines();

      for (int i = 0; i < children.size(); i++) {
        children.elementAt(i).removeAllDescendents();
      }

      if (necessary) {
        dad = "root";
        mom = "root";
      } else {
        hash.remove(id);
      }
    }

    public void removeAllDescendents() {
      for (int i = 0; i < children.size(); i++) {
        children.elementAt(i).removeAllDescendents();
      }
      if (!necessary) {
        hash.remove(id);
      }
    }

    public void createExtendedOutput(boolean listGeneration) {
      Person kid;
      String temp = "";

      while (!children.isEmpty()) {
        kid = children.lastElement();
        kid.createExtendedOutput(listGeneration);
        children.removeElementAt(children.size() - 1);
      }
      if (necessary && !isRoot() && !written) {
        if (hasSubfamilies) {
          switch (namingScheme) {
            case SCORE_99_NAMING_SCHEME:
              for (int i = 1; i < subFamily; i++) {
                temp += "99";
              }
              temp += famid;
              break;
            case ALPHA_NAMING_SCHEME:
              temp += famid + "_" + ext.getExcelColumn(subFamily - 1);
              break;
            case NUMBER_NAMING_SCHEME:
              temp += famid + "_" + subFamily;
              break;
            default:
              return;
          }
        } else {
          temp += famid;
        }
        temp += "\t" + (genotyped || !rename ? id : recoded_id);
        temp += "\t"
                + (father.isRoot() ? "0"
                                   : (father.genotyped || !rename ? father.id : father.recoded_id))
                + "\t"
                + (mother.isRoot() ? "0"
                                   : (mother.genotyped || !rename ? mother.id : mother.recoded_id))
                + "\t" + sex + "\t" + (phenotyped ? "2" : "0") + "\t"
                + (genotyped ? "1\t1" : "0\t0");

        if (listGeneration) {
          temp += "\tgen" + generation;
        }

        written = true;
        extended.add(temp);
      }
    }

    public void createNuclearOutput() {
      Person purse, kid;
      String temp = "", maxtMom = "-1";
      Vector<String> v;
      Hashtable<String, Vector<String>> poss = new Hashtable<String, Vector<String>>();
      int maxt = 0;
      String[] keys;

      for (int i = 0; i < children.size(); i++) {
        kid = children.elementAt(i);

        if (kid.genotyped) {
          temp = kid.mom;

          if (poss.containsKey(temp)) {
            v = poss.get(temp);
          } else {
            v = new Vector<String>();
            poss.put(temp, v);
          }

          v.add(famid + "\t" + kid.id + "\t"
                + (kid.father.genotyped || !rename ? kid.father.id : kid.father.recoded_id) + "\t"
                + (kid.mother.genotyped || !rename ? kid.mother.id : kid.mother.recoded_id) + "\t"
                + kid.sex + "\t" + (kid.phenotyped ? "2" : "0") + "\t"
                + (kid.genotyped ? "1\t1" : "0\t0"));
        }
      }

      keys = HashVec.getKeys(poss, false, false);
      for (String key : keys) {
        v = poss.get(key);
        if (v.size() > maxt) {
          maxt = v.size();
          maxtMom = key;
        }
      }
      if (maxt == 1) {
        nukular.removeAllElements();
        return;
      }

      purse = hash.get(maxtMom);
      v = poss.get(maxtMom);

      nukular.add(famid + "\t" + (genotyped || !rename ? id : recoded_id) + "\t0\t0\t" + sex);
      nukular.add(famid + "\t" + (purse.genotyped || !rename ? purse.id : purse.recoded_id)
                  + "\t0\t0\t" + purse.sex);

      for (int i = 0; i < v.size(); i++) {
        nukular.add(v.elementAt(i));
      }
    }

  }

  public TrimFam(Vector<String> preinfo) {
    this(preinfo, true, false, true, DEFAULT_NAMING_SCHEME, 0, false, false, new Logger());
  }

  public TrimFam(Vector<String> preinfo, boolean deleteSinglets, boolean unrelatedsOnly,
                 boolean canRename, int scheme, int favorGeneration, boolean listGeneration,
                 boolean allowMissingIndividuals, Logger logger) {
    Hashtable<String, String> genoed;
    String[] line, keys;
    String famid;
    Person purse;
    int numSubFams;
    String mostKids = null;

    hash = new Hashtable<String, Person>();
    extended = new Vector<String>();
    nukular = new Vector<String>();
    unused = new Vector<String>();
    rename = canRename;
    namingScheme = scheme;
    hasSubfamilies = false;
    countGenoAndPheno = 0;
    favGen = favorGeneration;
    log = logger;
    error = false;

    if (preinfo.size() == 0) {
      return;
    }

    famid = null;

    genoed = new Hashtable<String, String>();
    for (int i = 0; i < preinfo.size(); i++) {
      line = preinfo.elementAt(i).trim().split("[\\s]+");
      if (line.length < 7) {
        log.reportError("Error - there must be at least seven columns in order to run TrimFam (6 standards plus the VIP column)");
        error = true;
        return;
      }
      if (famid == null) {
        famid = line[0];
        hash.put("root", new Person(famid, "root", "-1", "-2"));
      } else if (!line[0].equals(famid)) {
        log.reportError("Error - TrimFam is now only set up to trim one family at a time. Use/write a wrapper if need be.");
        error = true;
        return;
      }
      purse = new Person(line[0], line[1], line[2], line[3]);
      purse.sex = line[4];
      purse.phenotyped = Integer.parseInt(line[5]) > 0;
      if (line[6].equals("1")) {
        purse.genotyped = true;
        genoed.put(line[1], "unused");
      } else if (!line[6].equals("0")) {
        log.reportError("Error - individual " + line[0] + "-" + line[1]
                        + " was not properly declared as genotyped (1) or not genotyped (0): '"
                        + line[6] + "'");
      }
      if (line.length > 7) {
        purse.priorityScore = Double.parseDouble(line[7]);
      } else {
        purse.priorityScore = purse.genotyped ? 1 : -1;
      }
      hash.put(line[1], purse);
    }

    if (hash.size() > 1) {
      if (!link(allowMissingIndividuals)) {
        return;
      }
      if (!link(allowMissingIndividuals)) {
        return;
      }
      hash.get("root").defineNecessaryPath();
      hash.get("root").trimUselessDescendents();
    }

    if (deleteSinglets) {
      if (countGenoAndPheno < 2) {
        hash.clear();
      }

      if (countGenoAndPheno == 2 && countAffectedPOpairs == 1) {
        hash.clear();
      }
    }

    if (hash.size() > 1) {
      link(allowMissingIndividuals);
      trimUselessFounders();
    }

    if (hash.size() > 1) {
      link(allowMissingIndividuals);
      trimSinglets(deleteSinglets);
    }

    if (deleteSinglets && hash.size() == 3 + 1) {
      hash.clear();
    }

    if (unrelatedsOnly) {
      codeGeneration();
      link(allowMissingIndividuals);
      delineateUniqueUnrelateds(allowMissingIndividuals);
      // I thought this was necessary to eliminate ungenotyped selections, but it appears to be
      // redundant
      // link();
      // hash.get("root").trimUselessDescendents();
      // link();
      // trimUselessFounders();
    } else if (hash.size() > 1) {
      mostKids = delineateFatherOfLargestSibship();
      if (rename) {
        recodeFromLargestSibship(mostKids);
      }
    }


    if (hash.size() > 1) {
      link(allowMissingIndividuals);
      if (!unrelatedsOnly) {
        numSubFams = detectSplitFams();
        if (numSubFams > 1) {
          log.report("Family " + famid + " needs to be split into " + numSubFams
                     + " sub families.");
        }
        hash.get(mostKids).createNuclearOutput();
      }
      hash.get("root").createExtendedOutput(listGeneration);

      // mark those that are genotyped and made it through to the last cut
      keys = HashVec.getKeys(hash, false, false);
      for (String key : keys) {
        if (genoed.containsKey(key)) {
          genoed.put(key, "used");
        }

      }
    }

    keys = HashVec.getKeys(genoed, false, false);
    for (String key : keys) {
      if (genoed.get(key).equals("unused")) {
        unused.add(key);
      }

    }
  }

  public boolean hadError() {
    return error;
  }

  public Vector<String> getExtendedFamilyInformation() {
    Vector<String> output = new Vector<String>();
    int[] keys = reorder(extended);

    for (int i = 0; i < extended.size(); i++) {
      output.add(extended.elementAt(keys[i]));
    }

    return output;
  }

  public Vector<String> getNuclearFamilyInformation() {
    Vector<String> output = new Vector<String>();
    int[] keys = reorder(nukular);

    for (int i = 0; i < nukular.size(); i++) {
      output.add(nukular.elementAt(keys[i]));
    }

    return output;
  }

  public boolean hasUnused() {
    return (!unused.isEmpty());
  }

  public Vector<String> getUnused() {
    return unused;
  }

  public boolean link(boolean allowMissingIndividuals) {
    Person purse;
    String[] keys;

    keys = HashVec.getKeys(hash);
    for (String key : keys) {
      hash.get(key).children.removeAllElements();
    }

    for (int i = 0; i < keys.length; i++) {
      purse = hash.get(keys[i]);
      if (!keys[i].equals("root")) {
        purse.father = hash.get(purse.dad);
        purse.mother = hash.get(purse.mom);
        if (purse.father == null) {
          if (allowMissingIndividuals) {
            purse.father = new Person(purse.famid, purse.dad, "root", "root");
            purse.father.sex = "1";
            purse.father.phenotyped = false;
            purse.father.genotyped = false;
            purse.father.necessary = true;
            purse.father.priorityScore = -1;
            hash.put(purse.dad, purse.father);
          } else {
            log.reportError("Link error - Pedigree '" + purse.famid
                            + "' does not contain all the indivudals referenced as parents (father: "
                            + purse.dad + ")");
            error = true;
          }
        }
        if (purse.mother == null) {
          if (allowMissingIndividuals) {
            purse.mother = new Person(purse.famid, purse.mom, "root", "root");
            purse.mother.sex = "2";
            purse.mother.phenotyped = false;
            purse.mother.genotyped = false;
            purse.mother.necessary = true;
            purse.mother.priorityScore = -1;
            hash.put(purse.mom, purse.mother);
          } else {
            log.reportError("Link error - Pedigree '" + purse.famid
                            + "' does not contain all the indivudals referenced as parents (mother: "
                            + purse.mom + ")");
            error = true;
          }
        }
        if (!error) {
          if (purse.father.isFemale()) {
            log.reportError("Error - '" + purse.id + "' has a female father");
            error = true;
          }
          if (purse.mother.isMale()) {
            log.reportError("Error - '" + purse.id + "' has a male mother");
            error = true;
          }

          purse.father.children.add(purse);
          if (!purse.dad.equals("root") || !purse.mom.equals("root")) {
            purse.mother.children.add(purse);
          }
        }
      }
    }

    return !error;
  }

  public void trimUselessFounders() {
    Person source, purse;

    source = hash.get("root");
    for (int i = 0; i < source.children.size(); i++) {
      purse = source.children.elementAt(i);
      if (purse.children.size() == 1 && purse.children.firstElement().father.children.size() == 1
          && purse.children.firstElement().mother.children.size() == 1
          && !purse.children.firstElement().father.genotyped
          && !purse.children.firstElement().mother.genotyped
          && purse.children.firstElement().mother.dad.equals("root")
          && purse.children.firstElement().father.dad.equals("root")) {

        purse.children.firstElement().dad = "root";
        purse.children.firstElement().mom = "root";
        hash.remove(purse.id);
        source.children.add(purse.children.firstElement());
      }
    }
  }

  public void trimSinglets(boolean deleteGenotypedSinglets) {
    Person source, purse;

    source = hash.get("root");
    for (int i = 0; i < source.children.size(); i++) {
      purse = source.children.elementAt(i);
      if (purse.children.size() == 0 && (deleteGenotypedSinglets || !purse.genotyped)) {
        hash.remove(purse.id);
      }
    }
  }

  // keeps the first highest ranked valid individual in each founder branch, deletes all individuals
  // related to that individual, and repeats until there is only one
  public void delineateUniqueUnrelateds(boolean allowMissingIndividuals) {
    Person source, purse;
    Hashtable<String, Vector<String>> famIndHash;
    String[] inds, subfams;
    double[] scores;
    boolean done;

    inds = HashVec.getKeys(hash);
    for (int i = 0; i < inds.length; i++) {
      hash.get(inds[i]).necessary = false;
    }

    for (String ind : inds) {
      purse = hash.get(ind);
      purse.priorityScore = purse.priorityScore * 1000 + favGen * purse.generation + 100;
    }

    done = false;
    source = hash.get("root");
    while (!done) {
      done = true;
      famIndHash = new Hashtable<String, Vector<String>>();
      inds = HashVec.getKeys(hash);
      for (String ind : inds) {
        purse = hash.get(ind);
        if (!purse.necessary && purse.genotyped) {
          HashVec.addToHashVec(famIndHash, purse.subFamily + "", ind, false);
        }
      }
      subfams = HashVec.getKeys(famIndHash);
      for (String subfam : subfams) {
        inds = Array.toStringArray(famIndHash.get(subfam));
        scores = new double[inds.length];
        for (int j = 0; j < scores.length; j++) {
          scores[j] = hash.get(inds[j]).priorityScore;
          // if (inds[j].equals("12748")) {
          // System.out.println("oi!");
          // }
        }
        purse = hash.get(inds[Sort.quicksort(scores)[scores.length - 1]]);
        purse.necessary = true;
        purse.removeParentalLines();
        link(allowMissingIndividuals);
        done = false;
      }
    }
    for (int i = 0; i < source.children.size(); i++) {
      purse = source.children.elementAt(i);
      if (purse.children.size() == 1 && purse.children.firstElement().father.children.size() == 1
          && purse.children.firstElement().mother.children.size() == 1
          && !purse.children.firstElement().father.genotyped
          && !purse.children.firstElement().mother.genotyped
          && purse.children.firstElement().mother.dad.equals("root")
          && purse.children.firstElement().father.dad.equals("root")) {

        purse.children.firstElement().dad = "root";
        purse.children.firstElement().mom = "root";
        hash.remove(purse.id);
        source.children.add(purse.children.firstElement());
      }
    }
    hash.get("root").createExtendedOutput(false);
  }

  public String delineateFatherOfLargestSibship() {
    Person purse;
    int mostKids, count;
    String withMostKids = "root";
    String[] keys;

    mostKids = 0;
    keys = HashVec.getKeys(hash, false, false);
    for (String key : keys) {
      count = 0;
      purse = hash.get(key);
      if (purse.isMale() && purse.children.size() >= mostKids) {
        for (int j = 0; j < purse.children.size(); j++) {
          if (purse.children.elementAt(j).genotyped) {
            count++;
          }
          if (purse.children.elementAt(j).id.equals("1")) {
            count++;
          }
        }
        if (count > mostKids) {
          mostKids = count;
          withMostKids = key;
        }
      }
    }

    return withMostKids;
  }

  public void recodeFromLargestSibship(String fatherOfLargestSibship) {
    Person purse, daddy;
    int count;
    Vector<String> usedRecodedIDs;
    String[] keys;

    // make parents of largest sibship 101 and 102
    daddy = hash.get(fatherOfLargestSibship);
    usedRecodedIDs = new Vector<String>();
    if (!daddy.genotyped && !daddy.isRoot()) {
      daddy.recoded_id = "101";
      usedRecodedIDs.add("101");
    }
    for (int i = 0; i < daddy.children.size(); i++) {
      purse = daddy.children.elementAt(i);
      if (purse.mother.recoded_id.equals("-1")) {
        count = 102;
        while (usedRecodedIDs.contains(count + "")) {
          count += 2;
        }
        purse.mother.recoded_id = count + "";
        usedRecodedIDs.add(count + "");
      }
    }

    // go up and down tree, assigning generation numbers
    if (!daddy.isRoot()) {
      daddy.spreadGeneration(1, true);
    }

    keys = HashVec.getKeys(hash);
    for (int i = 0; i < keys.length; i++) {
      purse = hash.get(keys[i]);
      if (!keys[i].equals("root") && purse.recoded_id.equals("-1")) {
        count = (purse.genotyped ? 900 : (purse.generation > 1 ? purse.generation * 100 : 100))
                + (purse.isMale() ? 1 : 2);
        while (usedRecodedIDs.contains(count + "")) {
          count += 2;
        }
        purse.recoded_id = count + "";
        usedRecodedIDs.add(count + "");
      }
    }
  }

  public void codeGeneration() {
    Person purse, root;
    int min;

    // go up and down tree, assigning generation numbers
    root = hash.get("root");
    for (int i = 0; i < root.children.size(); i++) {
      purse = root.children.elementAt(i);
      if (purse.generation == -99) {
        min = purse.spreadGeneration(1, false);
        min = purse.spreadGeneration(1 - min, true);
      }

    }
  }

  public static int[] reorder(Vector<String> original) {
    Vector<String> unordered = new Vector<String>();
    Vector<String> pool;
    Vector<String> subset = new Vector<String>();
    int[] order = new int[original.size()], keys, numArray;
    int count = 0, maxID = 0, id;
    String[] line;

    for (int i = 0; i < original.size(); i++) {
      line = original.elementAt(i).trim().split("[\\s]+");
      try {
        id = Integer.parseInt(line[1]);
      } catch (NumberFormatException nfe) {
        return Array.arrayOfIndices(original.size());
      }
      if (id > maxID) {
        maxID = id;
      }
      unordered.add(id + "");
    }
    pool = HashVec.cloneVectorString(unordered);

    for (int threshold = ((maxID / 100) * 100); threshold >= 0; threshold = threshold - 100) {
      for (int i = 0; i < pool.size(); i++) {
        if (Integer.parseInt(pool.elementAt(i)) > threshold) {
          subset.add(pool.elementAt(i));
        }
      }
      if (subset.size() > 0) {
        numArray = new int[subset.size()];
        for (int i = 0; i < subset.size(); i++) {
          numArray[i] = Integer.parseInt(subset.elementAt(i));
        }
        keys = Sort.quicksort(numArray);
        for (int i = 0; i < subset.size(); i++) {
          order[count++] = unordered.indexOf(subset.elementAt(keys[i]));
          pool.remove(subset.elementAt(i));
        }
        subset.removeAllElements();
      }
    }

    return order;
  }

  public int detectSplitFams() {
    Person zeroman = hash.get("root");
    Person founder = zeroman.children.elementAt(0);
    int subFam = 1;

    founder.spreadSubFamID(subFam);

    for (int i = 0; i < zeroman.children.size(); i++) {
      if (zeroman.children.elementAt(i).subFamily == 0) {
        subFam++;
        hasSubfamilies = true;
        zeroman.children.elementAt(i).spreadSubFamID(subFam);
      }
    }
    return subFam;
  }

  public static void demo() throws IOException {
    Vector<String> pre = new Vector<String>();
    Vector<String> v;

    pre.add("70504	101	201	202	1	0	0");
    pre.add("70504	102	0	0	2	0	0");
    pre.add("70504	103	0	0	1	0	0");
    pre.add("70504	104	201	202	2	0	0");
    pre.add("70504	201	0	0	1	0	0");
    pre.add("70504	202	0	0	2	0	0");
    pre.add("70504	301	0	0	1	0	0");
    pre.add("70504	304	0	0	2	0	0");
    pre.add("70504	1	101	102	1	2	1");
    pre.add("70504	6	1	304	1	2	1");
    pre.add("70504	18	1	304	1	2	1");
    pre.add("70504	24	101	102	2	2	1");
    pre.add("70504	26	301	24	1	2	1");
    pre.add("70504	27	301	24	2	2	1");
    pre.add("70504	46	301	24	1	2	1");
    pre.add("70504	71	103	104	2	2	1");

    try {
      v = new TrimFam(pre, false, false, false, DEFAULT_NAMING_SCHEME, 0, false, false,
                      new Logger()).getExtendedFamilyInformation();
      for (int i = 0; i < v.size(); i++) {
        System.out.println(v.elementAt(i));
      }

    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  // priority scores are not required, default scores will favor those genotyped over those
  // ungenotyped, need a separate file to favor affected over unaffected
  public static void trimAllFams(String pedigreeFile, String outfile,
                                 boolean renameUngenotypedFounders, boolean deleteSinglets,
                                 boolean unrelatedsOnly, String priorityScoresFile,
                                 int namingScheme, int favorGeneration, boolean listGeneration,
                                 boolean allowMissingIndividuals, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Hashtable<String, String> priorityScores;
    Hashtable<String, Vector<String>> preHash;
    Vector<String> fams, v;
    int count;
    long time;
    TrimFam tfam;
    String score;
    int countMissing;

    time = new Date().getTime();

    if (priorityScoresFile == null) {
      priorityScores = null;
    } else {
      priorityScores =
                     HashVec.loadFileToHashString(priorityScoresFile, new int[] {0, 1},
                                                  new int[] {2}, false, "\t", false, false, false);
    }
    fams = new Vector<String>();
    preHash = new Hashtable<String, Vector<String>>();
    countMissing = 0;
    try {
      reader = new BufferedReader(new FileReader(pedigreeFile));

      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        count = 0;
        for (int i = 6; i < line.length; i++) {
          if (line[i].equals("0")) {
            count++;
          }
        }
        if (priorityScores == null) {
          score = "";
        } else if (priorityScores.containsKey(line[0] + "\t" + line[1])) {
          score = "\t" + priorityScores.get(line[0] + "\t" + line[1]);
        } else if (count < line.length - 6) {
          if (countMissing <= 20) {
            log.reportError("Error - no valid value in hash for individual " + line[0] + "-"
                            + line[1]);
          } else {
            log.reportError("Error - no valid value in hash for individual " + line[0] + "-"
                            + line[1], true, false);
          }
          score = "\t-1";
          countMissing++;
        } else {
          score = "\t-1";
        }
        HashVec.addToHashVec(preHash,
                             line[0],
                             line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t"
                                      + line[4] + "\t" + line[5] + "\t"
                                      + (count == line.length - 6 ? "0" : "1") + score,
                             false);
        HashVec.addIfAbsent(line[0], fams);
      }
      if (countMissing > 0) {
        log.reportError("Just so you know... " + countMissing
                        + " indiviudal(s) are missing priority scores");
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + pedigreeFile + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + pedigreeFile + "\"");
      return;
    }

    try {
      writer = new PrintWriter(new FileWriter(outfile));
      count = 0;
      for (int i = 0; i < fams.size(); i++) {
        tfam = new TrimFam(preHash.get(fams.elementAt(i)), deleteSinglets, unrelatedsOnly,
                           renameUngenotypedFounders, namingScheme, favorGeneration, listGeneration,
                           allowMissingIndividuals, log);
        if (tfam.hadError()) {
          log.reportError("Aborting trimFam due to error");
          writer.close();
          return;
        }
        v = tfam.getExtendedFamilyInformation();
        for (int j = 0; j < v.size(); j++) {
          writer.println(v.elementAt(j));
        }
        writer.flush();
        if (v.size() > 0) {
          count++;
        }
      }
      writer.close();
      log.report(count + " of the original " + fams.size() + " families remain");
    } catch (Exception e) {
      log.reportException(e);
    }

    log.report("Finished in " + ext.getTimeElapsed(time));
  }

  public static int determineScheme(String scheme, Logger log) {
    if (scheme.startsWith("99")) {
      return SCORE_99_NAMING_SCHEME;
    } else if (scheme.toLowerCase().startsWith("alpha")) {
      return ALPHA_NAMING_SCHEME;
    } else if (scheme.toLowerCase().startsWith("num")) {
      return NUMBER_NAMING_SCHEME;
    } else {
      log.reportError("Invalid scheme type; choose from: 99s, alpha, or numeric");
      return -1;
    }
  }

  public static void fromParameters(String filename, Logger log) {
    Vector<String> paramV;

    paramV = Files.parseControlFile(filename, "trimFam",
                                    new String[] {"fams.pre", "out=trimmed.pre",
                                                  "renameUngenotypedFounders=false",
                                                  "deleteSinglets=false", "unrelatedsOnly=false",
                                                  "scores=null", "scheme=num", "favorGeneration=0",
                                                  "allowMissing=false",
                                                  "# add -h as an argument if you want a full description"},
                                    log);
    if (paramV != null) {
      paramV.addElement("logfile=" + log.getFilename());
      main(Array.toStringArray(paramV));
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String filename = "pedigrees.pre";
    // String filename = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\LOAD\\finalPedigree_geno.pre";
    // String filename = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\LOAD\\finalPedigree_testSplit.pre";
    // String filename = "C:\\Documents and Settings\\npankrat\\My Documents\\hearing\\struct.dat";
    // String filename = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\hearing\\00src\\PCA_data\\struct_PCA1_pruned.pre";
    // String filename = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\BOSS\\master_pedigree.txt";
    // String filename = "D:\\\\LOAD\\finalPedigree_geno.pre";
    // String filename = "D:\\BOSS\\Height\\plinkFilesCaucasianUnrelated\\master_pedigree.txt";
    String filename = "D:\\BOSS\\ancestry\\corrected_master_pedigree.txt";

    String outfile = null;

    // boolean deleteSinglets = true;
    boolean deleteSinglets = false;
    boolean unrelatedsOnly = false;
    String priorityScoresFile = null;
    String logfile = null;
    String scheme = null;
    int namingScheme = DEFAULT_NAMING_SCHEME;
    Logger log;
    int favorGeneration = 0;
    boolean renameUngenotypedFounders = false;
    boolean listGeneration = false;
    boolean allowMissingIndividuals = false;

    String usage = "\n" + "link.TrimFam requires 0-1 arguments\n" + "   (1) filename (i.e. file="
                   + filename + " (default))\n"
                   + "   (2) output filename (i.e. out=[input file]-trimmed.pre (default))\n"
                   + "   (3) delete singlets (i.e. deleteSinglets=" + deleteSinglets
                   + " (default))\n"
                   + "   (4) only keep unrelated indiviudals (i.e. unrelatedsOnly=" + unrelatedsOnly
                   + " (default))\n" + "   (5) priority scores filename (i.e. score="
                   + priorityScoresFile + " (default))\n"
                   + "   (6) naming scheme for split fams (i.e. scheme="
                   + NAMING_SCHEME_KEY[namingScheme] + " (default))\n"
                   + "   (7) favor older or earlier generations if unrelatedsOnly (i.e. favorGeneration="
                   + favorGeneration + " (default; +1 for older, -1 for younger, 0 for neither))\n"
                   + "   (8) rename ungenotyped founders 101, 102, etc. (i.e. renameUngenotypedFounders="
                   + renameUngenotypedFounders + " (default))\n"
                   + "   (9) list generation in final file (i.e. listGeneration=" + listGeneration
                   + " (default))\n"
                   + "   (10) if an individual is missing, create them and code as a founder of appropriateSex (i.e. allowMissing="
                   + allowMissingIndividuals + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out=")) {
        outfile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("renameUngenotypedFounders=")) {
        renameUngenotypedFounders = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("deleteSinglets=")) {
        deleteSinglets = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("unrelatedsOnly=")) {
        unrelatedsOnly = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("scores=")) {
        priorityScoresFile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("scheme=")) {
        scheme = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("logfile=")) {
        logfile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("favorGeneration=")) {
        favorGeneration = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("listGeneration=")) {
        listGeneration = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("allowMissing=")) {
        allowMissingIndividuals = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.contains("=")) {
        System.err.println("Error - unknown paramter: " + arg);
      } else {
        filename = arg;
        numArgs--;
      }
    }

    log = new Logger(logfile);
    if (logfile == null) {
      log.report("No filename for a logfile was delineated, reporting only to stdout");
    }
    if (scheme != null) {
      namingScheme = determineScheme(scheme, log);
      log.report("Changing scheme to " + NAMING_SCHEME_KEY[namingScheme]);
    } else {
      log.report("Using default scheme (" + NAMING_SCHEME_KEY[namingScheme] + ")");
    }
    if (outfile == null) {
      outfile = filename + "-trimmed.pre";
    }
    if (unrelatedsOnly && deleteSinglets) {
      log.reportError("Error - delineating unrelated indiviudals and then deleting singlets as well will result in a null set");
      return;
    }
    if (numArgs != 0) {
      log.reportError(usage);
      return;
    }
    try {
      // demo();
      trimAllFams(filename, outfile, renameUngenotypedFounders, deleteSinglets, unrelatedsOnly,
                  priorityScoresFile, namingScheme, favorGeneration, listGeneration,
                  allowMissingIndividuals, log);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
