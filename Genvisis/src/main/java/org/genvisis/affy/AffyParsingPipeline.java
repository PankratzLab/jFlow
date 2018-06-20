package org.genvisis.affy;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Set;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import com.google.common.collect.Sets;

public class AffyParsingPipeline {

  Project proj;
  private String callFile;
  private String confFile;
  private String sigFile;

  private String delim = "\t";

  private BufferedReader confReader;
  private BufferedReader callReader;
  private BufferedReader sigReader;

  private String[] confHeader;
  private String[] callHeader;
  private String[] sigHeader;

  int numSamples;
  long fingerprint;

  public void setProject(Project proj) {
    this.proj = proj;
  }

  public void setCallFile(String callFile) {
    this.callFile = callFile;
  }

  public void setConfFile(String confFile) {
    this.confFile = confFile;
  }

  public void setSigFile(String quantNormFile) {
    this.sigFile = quantNormFile;
  }

  public void run() {
    System.out.println("Conf: " + confFile);
    System.out.println("Call: " + callFile);
    System.out.println("Sig: " + sigFile);

    try {
      initReaders();
    } catch (IOException e) {
      // TODO handle gracefully
    }

    // TODO tweak this number:
    int numMarkersPerFile = 1000;
    // setup marker files (determine memory buffer for MarkerData)
    MarkerData md = null;

    // setup MarkerLookup
    // load and write buffer of MarkerData
  }

  @SuppressWarnings("deprecation")
  void initReaders() throws IOException {
    confReader = Files.getAppropriateReader(confFile);
    callReader = Files.getAppropriateReader(confFile);
    sigReader = Files.getAppropriateReader(sigFile);

    String line;

    line = null;
    while ((line = confReader.readLine()) != null && line.charAt(0) == '#') {
      // seek to header
    }
    confHeader = line.trim().split(delim, -1);

    line = null;
    while ((line = callReader.readLine()) != null && line.charAt(0) == '#') {
      // seek to header
    }
    callHeader = line.trim().split(delim, -1);

    line = null;
    while ((line = sigReader.readLine()) != null && line.charAt(0) == '#') {
      // seek to header
    }
    sigHeader = line.trim().split(delim, -1);

    ensureHeaderContentsSame();
    numSamples = confHeader.length - 1;
    fingerprint = org.genvisis.cnv.filesys.MarkerSet.fingerprint(ArrayUtils.subArray(confHeader,
                                                                                     1));
  }

  private void ensureHeaderContentsSame() {
    if (confHeader.length != callHeader.length || confHeader.length != sigHeader.length
        || callHeader.length != sigHeader.length) {
      throw new IllegalStateException("File headers have different lengths: [Conf: "
                                      + confHeader.length + "], [Call: " + callHeader.length
                                      + "], [Sig: " + sigHeader.length + "].");
    }

    Set<String> confs = Sets.newHashSet(confHeader);
    Set<String> calls = Sets.newHashSet(callHeader);
    Set<String> sigs = Sets.newHashSet(sigHeader);

    boolean diff = Sets.symmetricDifference(confs, calls).size() > 0
                   || Sets.symmetricDifference(confs, sigs).size() > 0
                   || Sets.symmetricDifference(calls, sigs).size() > 0;
    if (diff) {
      throw new IllegalStateException("File headers have different values.");
    }
  }

  MarkerData parseLine() throws IOException {
    String confLine = confReader.readLine();
    String callLine = callReader.readLine();
    String sigLine = sigReader.readLine();
    String[] confs = confLine.trim().split(delim, -1);
    String[] calls = callLine.trim().split(delim, -1);
    String[] sigsA = sigLine.trim().split(delim, -1);
    String[] sigsB = sigLine.trim().split(delim, -1);
    ensureSame(confs[0], calls[0], sigsA[0], sigsB[0]);

    String mkr = calls[0];

    float[] gcs = new float[numSamples];
    float[] xRaws = null;
    float[] yRaws = null;
    float[] xs = new float[numSamples];
    float[] ys = new float[numSamples];
    float[] thetas = null;
    float[] rs = null;
    float[] bafs = null;
    float[] lrrs = null;
    byte[] abGenos = new byte[numSamples];
    byte[] forwardGenos = null;

    for (int i = 0; i < numSamples; i++) {
      abGenos[i] = (byte) Integer.parseInt(calls[i]);
      gcs[i] = Float.parseFloat(confs[i]);
      xs[i] = (float) AffySNP6Tables.power2(sigsA[i]);
      ys[i] = (float) AffySNP6Tables.power2(sigsB[i]);
    }
    MarkerData md = new MarkerData(mkr, (byte) 0, 0, fingerprint, gcs, xRaws, yRaws, xs, ys, thetas,
                                   rs, bafs, lrrs, abGenos, forwardGenos);
    // TODO compute centroids
    // TODO recomputeLRR_BAFs
    return md;
  }

  private void ensureSame(String conf, String call, String sigA, String sigB) {
    if (conf.equals(call) && sigA.startsWith(conf) && sigB.startsWith(conf)) {
      return;
    }
    throw new IllegalStateException("Affymetrix parsing became desynchronised - expected these to all be the same: [Call: "
                                    + call + "], [Conf: " + conf + "], [SigA: " + sigA + ", SigB: "
                                    + sigB + "].");
  }

}
