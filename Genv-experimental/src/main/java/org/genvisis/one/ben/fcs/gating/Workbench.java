package org.genvisis.one.ben.fcs.gating;

import java.util.HashMap;

public class Workbench {
  
  public static class SampleNode {
    String id;
    String fcsFile;
    Gating gating;
    org.w3c.dom.Element sampleNode;
    org.w3c.dom.Document doc;
    
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
    sn.id = getNewSampleID();
    sn.fcsFile = fcsFile;
    if (applyTemplate) {
      sn.gating = templateGating.copy(fcsFile);
    } else {
      sn.gating = new Gating();
    }
    samples.put(sn.id, sn);
    return sn.id;
  }
  
  private String getNewSampleID() {
    int id = samples.size();
    boolean done = false;
    notDone : while (!done) {
      if (samples.containsKey("" + id)) {
        id++;
        continue notDone;
      }
      done = true;
    }
    return id + "";
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
    }
    return false;
  }

  public String getSampleID(String filename) {
    for (SampleNode sn : samples.values()) {
      if (sn.fcsFile.equals(filename)) {
        return sn.id;
      }
    }
    return null;
  }

}
