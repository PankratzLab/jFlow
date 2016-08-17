package org.genvisis.seq.pathway;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Serializable;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.GeneTrack;

public class Pathways implements Serializable {
  private static class KeggPathwayProducer extends AbstractProducer<Pathway> {
    private final GeneTrack geneTrack;
    private final BufferedReader pathwayReader;
    private final Logger log;

    public KeggPathwayProducer(GeneTrack geneTrack, BufferedReader pathwayReader, Logger log) {
      super();
      this.geneTrack = geneTrack;
      this.pathwayReader = pathwayReader;
      this.log = log;
    }

    @Override
    public boolean hasNext() {
      boolean ready = false;

      try {
        ready = pathwayReader.ready();
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }

      return ready;

    }

    @Override
    public Callable<Pathway> next() {
      String[] line;
      try {
        line = pathwayReader.readLine().trim().split("\t");
        return new KeggPathwayWorker(line, geneTrack, log);

      } catch (IOException e) {
        log.reportException(e);
        e.printStackTrace();
      }
      // TODO Auto-generated method stub
      return null;
    }

    @Override
    public void shutdown() {
      try {
        pathwayReader.close();
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      // TODO Auto-generated method stub

    }

  }
  private static class KeggPathwayWorker implements Callable<Pathway> {

    private final String[] pathwayLine;
    private final GeneTrack geneTrack;
    private final Logger log;

    public KeggPathwayWorker(String[] pathwayLine, GeneTrack geneTrack, Logger log) {
      super();
      this.pathwayLine = pathwayLine;
      this.geneTrack = geneTrack;
      this.log = log;
    }

    @Override
    public Pathway call() throws Exception {
      return matchPathwayGenes(geneTrack, pathwayLine, log);
    }

  }

  /**
   * 
   */
  private static final long serialVersionUID = 1L;

  private static final String HUMAN = "hsa";
  private static final String RS = "rs:";
  private static final String KEGG_HUMAN_PATHWAY_LIST = "http://rest.kegg.jp/list/pathway/" + HUMAN;
  private static final String KEGG_HUMAN_PATHWAY_GENE_LINK =
      "http://rest.kegg.jp/link/" + HUMAN + "/";
  private static final String REF_SEQ_GENE_LINK = "http://rest.genome.jp/link/refnuc/";

  private static final String PATH = "path";

  private static Hashtable<String, ArrayList<String>> getNCBILookupForPathway(String pathway,
                                                                              Logger log) {
    ArrayList<String> keggGenes = new ArrayList<String>();

    try {
      BufferedReader in =
          new BufferedReader(new InputStreamReader(new URL(KEGG_HUMAN_PATHWAY_GENE_LINK + "/"
                                                           + pathway).openStream()));
      try {
        Thread.sleep(100);
      } catch (InterruptedException ie) {
      }
      while (in.ready()) {
        String[] line = in.readLine().trim().split("\t");
        String gene = line[1];
        keggGenes.add(gene);

      }
      in.close();

    } catch (MalformedURLException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    Hashtable<String, ArrayList<String>> lookup = new Hashtable<String, ArrayList<String>>();
    for (int i = 0; i < keggGenes.size(); i++) {
      lookup.put(keggGenes.get(i), new ArrayList<String>());

    }

    ArrayList<String[]> splits =
        Array.splitUpArray(keggGenes.toArray(new String[keggGenes.size()]),
                           Math.round((float) keggGenes.size() / 50) + 1, log);
    for (int i = 0; i < splits.size(); i++) {

      String q = Array.toStr(splits.get(i), "+");
      String url = REF_SEQ_GENE_LINK + q;
      try {
        URLConnection uc = new URL(url).openConnection();
        uc.connect();
        BufferedReader in = new BufferedReader(new InputStreamReader(uc.getInputStream()));
        try {
          Thread.sleep(100);
        } catch (InterruptedException ie) {
        }
        while (in.ready()) {
          String[] line = in.readLine().trim().split("\t");
          String gene = line[0];
          lookup.get(gene).add(line[1].replaceAll(RS, ""));
        }
        in.close();
      } catch (MalformedURLException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      // log.reportTimeInfo("Using url " + url);

    }
    // if (keggGenes.size() != lookup.size()) {
    // log.reportTimeWarning("Could not look up all refseqIds for pathway " + pathway + ", " +
    // keggGenes.size() + " genes and found " + lookup.size());
    // for (int j = 0; j < keggGenes.size(); j++) {
    // if (!lookup.containsKey(keggGenes.get(j).replaceAll(HUMAN + ":", ""))) {
    // log.reportTimeError("Missing " + keggGenes.get(j));
    // }
    // }
    // }

    return lookup;
  }

  public static Pathways load(String filename) {
    return (Pathways) SerializedFiles.readSerial(filename, false, false);
  }

  public static void main(String[] args) {
    test();
  }

  private static Pathway matchPathwayGenes(GeneTrack geneTrack, String[] line, Logger log) {

    HashSet<String> pathGenes = new HashSet<String>();
    String path = line[0].replaceAll(PATH + ":", "");
    String pathName = ext.replaceWithLinuxSafeCharacters(line[1], true);
    Hashtable<String, ArrayList<String>> lookup = getNCBILookupForPathway(path, log);
    GeneData[][] geneDatas = geneTrack.getGenes();
    for (String keggGene : lookup.keySet()) {
      ArrayList<String> tmp = lookup.get(keggGene);
      boolean found = false;
      for (GeneData[] geneData : geneDatas) {
        if (!found) {
          for (int j = 0; j < geneData.length; j++) {
            String[] ncbis = geneData[j].getNcbiAssessionNumbers();
            for (int k = 0; k < tmp.size(); k++) {
              if (ext.indexOfStr(tmp.get(k), ncbis) >= 0) {
                found = true;
                pathGenes.add(geneData[j].getGeneName());
              }
            }
          }
        }
      }
    }
    if (pathGenes.size() != lookup.size()) {
      log.reportTimeWarning("Could not match all entries for pathway " + path
                            + " , setting incomplete flag");
      // return new Pathway(path + ":Invalid", new GeneData[] {}, false, log);
    }
    // else {
    ArrayList<GeneData> genes = new ArrayList<GeneData>();
    for (String gene : pathGenes) {
      GeneData[] tmp = geneTrack.lookupAllGeneData(gene);
      for (GeneData element : tmp) {
        genes.add(element);
      }
    }
    return new Pathway(pathName, genes.toArray(new GeneData[genes.size()]), true,
                       pathGenes.size() == lookup.size(), log);
  }

  public static void test() {
    String geneTrack = "N:/statgen/NCBI/RefSeq_hg19.gtrack";
    GeneTrack geneTrack2 = GeneTrack.load(geneTrack, false);
    Logger log = new Logger();

    String keggSer = "D:/data/Project_Tsai_Project_021/KeggTest/kegg.ser";

    Pathways way = new Pathways(log);
    // if (!Files.exists(keggSer)) {
    way.downloadKeggData(geneTrack2, 3, log);
    way.serialize(keggSer);
    new GenomeRegions(geneTrack2, way, log).serialize(ext.addToRoot(keggSer, ".genomeRegions"));
    // } else {
    // log.reportTimeInfo("Loading pathways from " + keggSer);
    // way = load(keggSer);
    // }
  }

  private Pathway[] pathways;
  // private Logger log;

  public Pathways(Logger log) {
    super();
    // this.log = log;
  }

  public Pathways(Pathway[] pathways, Logger log) {
    super();
    this.pathways = pathways;
    // this.log = log;
  }

  private void downloadKeggData(GeneTrack geneTrack, int numthreads, Logger log) {
    ArrayList<Pathway> paths = new ArrayList<Pathway>();
    BufferedReader in;
    try {
      int index = 0;
      in = new BufferedReader(new InputStreamReader(new URL(KEGG_HUMAN_PATHWAY_LIST).openStream()));
      KeggPathwayProducer producer = new KeggPathwayProducer(geneTrack, in, log);
      WorkerTrain<Pathway> train = new WorkerTrain<Pathway>(producer, numthreads, numthreads, log);
      while (train.hasNext()) {
        index++;
        if (index % 20 == 0) {
          log.reportTimeInfo("Parsed " + index + " pathways");
          // this.pathways = paths.toArray(new Pathway[paths.size()]);
          // train.shutdown();
          // return;
        }
        paths.add(train.next());
      }
    } catch (MalformedURLException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    pathways = paths.toArray(new Pathway[paths.size()]);
  }

  public Pathway getPathway(String name) {
    for (Pathway pathway : pathways) {
      if (pathway.getPathwayName().equals(name)) {
        return pathway;
      }
    }
    return null;
  }

  public Pathway[] getPathways() {
    return pathways;
  }

  public Pathway[] getPathwaysFor(GeneData geneData) {
    ArrayList<Pathway> tmp = new ArrayList<Pathway>();
    for (Pathway pathway : pathways) {
      if (pathway.containsGene(geneData.getGeneName())) {
        tmp.add(pathway);
      }
    }
    return tmp.toArray(new Pathway[tmp.size()]);
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

}
