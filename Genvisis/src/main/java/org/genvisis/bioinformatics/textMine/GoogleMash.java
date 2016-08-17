package org.genvisis.bioinformatics.textMine;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Serializable;
import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.concurrent.Callable;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;

/**
 * @author lane0212
 *
 *         Mash some google searches together and see how many hits come out...
 */
public class GoogleMash implements Serializable {

  private static class GQuery implements Serializable {
    /**
     * 
     */
    private static final long serialVersionUID = 1L;
    private final QUERY_TYPE type;
    private final String base;
    private final String extra;
    private int hits;
    private final Logger log;

    public GQuery(QUERY_TYPE type, String base, String extra, Logger log) {
      super();
      this.type = type;
      this.base = base;
      this.extra = extra;
      this.log = log;
    }

    private String finalize(String result) {
      result = result.replaceAll(" results.*", "");
      result = result.replaceAll(",", "");
      return result;
    }

    public int getHits() {
      return hits;
    }

    private boolean hasCount(String spec, String line, String patternCount, String patternRegex) {
      boolean hasCount = false;
      if (line.contains(patternCount)) {
        String result = line.replaceAll(patternRegex, "");
        result = finalize(result);
        log.reportTimeInfo(spec + "\n" + result);
        hasCount = true;
        try {
          int tmpHits = Integer.parseInt(result);
          setHits(tmpHits);
        } catch (NumberFormatException nfe) {
          log.reportTimeError("invalid number " + result);
        }
      }
      return hasCount;
    }

    public void query() {
      try {
        String spec =
            type.getAddress() + URLEncoder.encode(base + " + \"" + extra + "\"", CHAR_SET);
        URL url = new URL(spec);
        log.reportTimeInfo(spec);
        try {
          double rand = Math.random();
          int randomSleep = (int) (rand * 10000);
          Thread.sleep(randomSleep + 1000);
        } catch (InterruptedException ie) {
        }
        final URLConnection connection = url.openConnection();
        connection.setConnectTimeout(60000);
        connection.setReadTimeout(60000);
        connection.addRequestProperty("User-Agent", "Chrome/12.0.742.5");
        BufferedReader in = new BufferedReader(new InputStreamReader(connection.getInputStream()));
        try {
          double rand = Math.random();

          int randomSleep = (int) (rand + 1 * 5000);
          Thread.sleep(1000 + randomSleep);
        } catch (InterruptedException ie) {
        }
        readme: while (in.ready()) {
          String line = in.readLine().trim();
          switch (type) {
            case GOOGLE_REGULAR:
              if (hasCount(spec, line, "id=\"resultStats\">", ".*id=\"resultStats\">About ")) {
                break readme;
              }
              break;
            case GOOGLE_SCHOLAR:
              if (hasCount(spec, line, "id=\"gs_ab_md\">About ", ".*id=\"gs_ab_md\">About ")) {
                break readme;
              }
              break;
            default:
              break;

          }
        }
        in.close();
      } catch (MalformedURLException e) {
        log.reportException(e);
        e.printStackTrace();
      } catch (UnsupportedEncodingException e) {
        log.reportException(e);
        e.printStackTrace();
      } catch (IOException e) {
        log.reportException(e);
        e.printStackTrace();
      }
    }

    public void setHits(int hits) {
      this.hits = hits;
    }
  }
  public enum QUERY_TYPE {
                          GOOGLE_SCHOLAR(DEFAULT_GOOGLE_SCHOLAR), GOOGLE_REGULAR(DEFAULT_GOOGLE);

    private String address;

    private QUERY_TYPE(String address) {
      this.address = address;
    }

    public String getAddress() {
      return address;
    }

    public void setAddress(String address) {
      this.address = address;
    }

  }
  private static class QueryProducer extends AbstractProducer<GQuery> {
    private final GQuery[] queries;
    private int index;

    public QueryProducer(GQuery[] queries) {
      super();
      this.queries = queries;
      index = 0;
    }

    @Override
    public boolean hasNext() {
      return index < queries.length;
    }

    @Override
    public Callable<GQuery> next() {
      final GQuery cQuery = queries[index];
      Callable<GQuery> worker = new Callable<GoogleMash.GQuery>() {

        @Override
        public GQuery call() throws Exception {
          cQuery.query();
          return cQuery;
        }
      };
      index++;
      return worker;
    }
  }

  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  // public static String GOOGLE_API =
  // "http://ajax.googleapis.com/ajax/services/search/web?v=1.0&q=";
  public static String DEFAULT_GOOGLE = "https://www.google.com/search?q=";
  public static String DEFAULT_GOOGLE_SCHOLAR = "https://scholar.google.com/scholar?hl=en&q=";

  public static String CHAR_SET = "UTF-8";

  public static void count(String baseFile, String queryFile, int numThreads, Logger log) {
    String[] bases = HashVec.loadFileToStringArray(baseFile, false, new int[] {0}, true);
    String[] extras = HashVec.loadFileToStringArray(queryFile, false, new int[] {0}, true);
    log.reportTimeInfo(bases.length + " bases with " + extras.length + " extras");
    GoogleMash[] googleMashs = new GoogleMash[bases.length];

    for (int i = 0; i < QUERY_TYPE.values().length; i++) {
      String output = ext.addToRoot(baseFile, "." + QUERY_TYPE.values()[i] + ".gQuery");

      try {
        PrintWriter writer = new PrintWriter(new FileWriter(output));
        writer.print("BASE_QUERY");
        for (String extra : extras) {
          writer.print("\t" + extra);
        }
        writer.println();
        for (int k = 0; k < bases.length; k++) {
          String serFile = ext.addToRoot(baseFile, "." + bases[k] + ".ser");
          if (!Files.exists(serFile)) {
            GQuery[] gQueries = new GQuery[extras.length];
            for (int k2 = 0; k2 < extras.length; k2++) {
              gQueries[k2] = new GQuery(QUERY_TYPE.values()[i], bases[k], extras[k2], log);
            }
            googleMashs[k] = new GoogleMash(gQueries, log);
            googleMashs[k].queryAll(numThreads);
          } else {
            log.reportFileExists(serFile);
            googleMashs[k] = load(serFile);
          }
          writer.print(bases[k]);
          for (int j = 0; j < googleMashs[k].getQueries().length; j++) {
            writer.print("\t" + googleMashs[k].getQueries()[j].getHits());
          }
          writer.println();
          writer.flush();
        }
        writer.close();
      } catch (Exception e) {
        log.reportError("Error writing to " + output);
        log.reportException(e);
      }
    }
  }

  public static GoogleMash load(String filename) {
    return (GoogleMash) SerializedFiles.readSerial(filename, false, false);
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String baseFile = "GoogleMash.genes";
    String queryFile = null;
    String logfile = null;
    int numThreads = 1;
    Logger log;

    String usage = "\n" + "bioinformatics.textMine.GoogleMash requires 0-1 arguments\n";
    usage += "   (1) full path to base query file, like a file of genes (i.e. baseFile=" + baseFile
             + " (default))\n" + "";
    usage += "   (2) full path to query extras file, like a file of phenotypes  (i.e. queryFile="
             + queryFile + " (default))\n" + "";
    usage += PSF.Ext.getNumThreadsCommand(3, numThreads);
    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("baseFile=")) {
        baseFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("queryFile=")) {
        queryFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
        numThreads = ext.parseIntArg(arg);
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
      log = new Logger(logfile);
      count(baseFile, queryFile, numThreads, log);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private GQuery[] queries;

  private final Logger log;

  public GoogleMash(GQuery[] queries, Logger log) {
    super();
    this.queries = queries;
    this.log = log;

  }

  public GQuery[] getQueries() {
    return queries;
  }

  public void queryAll(int numThreads) {
    GQuery[] tmp = new GQuery[queries.length];
    QueryProducer producer = new QueryProducer(queries);
    WorkerTrain<GQuery> train =
        new WorkerTrain<GoogleMash.GQuery>(producer, numThreads, numThreads, log);
    int index = 0;
    while (train.hasNext()) {
      tmp[index] = train.next();
      index++;
    }
    queries = tmp;
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

}
