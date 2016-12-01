package org.genvisis.one.ben.fcs.gating;

import java.io.File;
import java.io.IOError;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.net.URLEncoder;
import java.util.HashMap;

import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.AbstractPanel2.AxisTransform;

public class Workbench {
  
  public static class SampleNode {
    public String id;
    public String fcsFile;
    public Gating gating;
    public HashMap<String, AxisTransform> savedTransforms;
    public org.w3c.dom.Element sampleNode;
    public org.w3c.dom.Document doc;
    
    public Gating getGating() {
      return gating;
    }
  }

  Gating templateGating;
  HashMap<String, SampleNode> samples;
  
  public Workbench() {
    templateGating = new Gating();
    samples = new HashMap<String, Workbench.SampleNode>();
  }
  
  public String addNewSample(String fcsFile, boolean applyTemplate) {
    SampleNode sn = new SampleNode();
    try {
      sn.fcsFile = URLDecoder.decode(fcsFile, "utf-8");
    } catch (UnsupportedEncodingException e) {
      System.err.println("Error - " + e.getMessage());
      sn.fcsFile = fcsFile;
    }
    sn.id = ext.removeDirectoryInfo(sn.fcsFile);//getNewSampleID();
    if (applyTemplate) {
      sn.gating = templateGating.copy(fcsFile);
    } else {
      sn.gating = new Gating();
    }
    samples.put(sn.id, sn);
    return sn.id;
  }
  
  private String getNewSampleID() {
    int id = samples.size() + 1;
    boolean done = false;
    while (!done) {
      if (samples.containsKey(Integer.toString(id))) {
        id++;
        continue;
      }
      done = true;
    }
    return Integer.toString(id);
  }
  
  public SampleNode getSample(String currentSampleID) {
    return samples.get(currentSampleID);
  }

  public void clearGating(String currentSampleID) {
    // TODO should gate IDs be deleted here?  also, an UNDO might be nice...
    setGatingForSample(currentSampleID, new Gating());
  }

  public void setGatingForSample(String currentSampleID, Gating gateStrat) {
    samples.get(currentSampleID).gating = gateStrat;
  }

  public boolean containsSampleFile(String filename) {
    for (SampleNode sn : samples.values()) {
      if (sn.fcsFile.equals(filename)) {
        return true;
      }
      try {
        String f1 = URLDecoder.decode(new File(sn.fcsFile).getCanonicalPath(), "UTF-8");
        String f2 = URLDecoder.decode(new File(filename).getCanonicalPath(), "UTF-8");
        if (f1.equals(f2)) {
          return true;
        }
      } catch (IOException e) {}
      if ((new File(sn.fcsFile)).equals(new File(filename))) {
        return true;
      }
    }
    return false;
  }

  public String getSampleID(String filename) {
    for (SampleNode sn : samples.values()) {
      if (sn.fcsFile.equals(filename)) {
        return sn.id;
      }
      try {
        String f1 = URLDecoder.decode(new File(sn.fcsFile).getCanonicalPath(), "UTF-8");
        String f2 = URLDecoder.decode(new File(filename).getCanonicalPath(), "UTF-8");
        if (f1.equals(f2)) {
          return sn.id;
        }
      } catch (IOException e) {}
      if ((new File(sn.fcsFile)).equals(new File(filename))) {
        return sn.id;
      }
    }
    return null;
  }

}
