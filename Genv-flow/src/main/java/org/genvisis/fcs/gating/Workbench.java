package org.genvisis.fcs.gating;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import org.genvisis.fcs.AbstractPanel2.AxisTransform;

public class Workbench {

  public static class SampleNode {

    public String id;
    public String fcsFile;
    public Gating gating;
    public HashMap<String, AxisTransform> savedTransforms;
    public org.w3c.dom.Element sampleNode;
    public org.w3c.dom.Document doc;
    public String wspFile;

    public Gating getGating() {
      return gating;
    }

    @Override
    public int hashCode() {
      final int prime = 31;
      int result = 1;
      result = prime * result + ((fcsFile == null) ? 0 : fcsFile.hashCode());
      result = prime * result + ((id == null) ? 0 : id.hashCode());
      result = prime * result + ((wspFile == null) ? 0 : wspFile.hashCode());
      return result;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (getClass() != obj.getClass()) return false;
      SampleNode other = (SampleNode) obj;
      if (fcsFile == null) {
        if (other.fcsFile != null) return false;
      } else if (!fcsFile.equals(other.fcsFile)) return false;
      if (id == null) {
        if (other.id != null) return false;
      } else if (!id.equals(other.id)) return false;
      if (wspFile == null) {
        if (other.wspFile != null) return false;
      } else if (!wspFile.equals(other.wspFile)) return false;
      return true;
    }
  }

  Gating templateGating;
  Set<SampleNode> samples;

  public Workbench() {
    templateGating = new Gating();
    samples = new HashSet<>();
  }

  public Set<String> getAllSamples() {
    return this.samples.stream().map(sn -> sn.id).collect(Collectors.toSet());
  }

  public String addNewSample(String fcsFile, boolean applyTemplate) {
    SampleNode sn = new SampleNode();
    try {
      sn.fcsFile = URLDecoder.decode(fcsFile, "utf-8");
    } catch (UnsupportedEncodingException e) {
      System.err.println("Error - " + e.getMessage());
      sn.fcsFile = fcsFile;
    }
    // sn.id = ext.removeDirectoryInfo(sn.fcsFile);// getNewSampleID();
    sn.id = getNewSampleID();
    if (applyTemplate) {
      sn.gating = templateGating.copy(fcsFile);
    } else {
      sn.gating = new Gating();
    }
    samples.add(sn);
    return sn.id;
  }

  private String getNewSampleID() {
    int id = samples.size() + 1;
    boolean done = false;
    while (!done) {
      if (getSample(Integer.toString(id)) != null) {
        id++;
        continue;
      }
      done = true;
    }
    return Integer.toString(id);
  }

  public SampleNode getSample(String currentSampleID) {
    Optional<SampleNode> samp =
        samples
            .stream()
            .filter(sn -> sn.id.equals(currentSampleID) || sn.fcsFile.equals(currentSampleID))
            .findFirst();
    return samp.orElse(null);
  }

  public void clearGating(String currentSampleID) {
    // TODO should gate IDs be deleted here? also, an UNDO might be nice...
    setGatingForSample(currentSampleID, new Gating());
  }

  public void setGatingForSample(String currentSampleID, Gating gateStrat) {
    getSample(currentSampleID).gating = gateStrat;
  }

  public boolean containsSampleFile(String filename) {
    return getSampleID(filename) != null;
  }

  public String getSampleID(String filename) {
    for (SampleNode sn : samples) {
      if (sn.fcsFile.equals(filename)) {
        return sn.id;
      }
      try {
        String f1;
        String f2;
        f1 = URLDecoder.decode(new File(sn.fcsFile).getCanonicalPath(), "UTF-8");
        f2 = URLDecoder.decode(new File(filename).getCanonicalPath(), "UTF-8");
        if (f1.equals(f2)) {
          return sn.id;
        }
        f1 = URLDecoder.decode(new File(sn.fcsFile).getName(), "UTF-8");
        f2 = URLDecoder.decode(new File(filename).getName(), "UTF-8");
        if (f1.equals(f2)) {
          return sn.id;
        }
      } catch (IOException e) {
      }
      if ((new File(sn.fcsFile).getName()).equals(new File(filename).getName())) {
        return sn.id;
      }
    }
    return null;
  }
}
