package org.genvisis.cnv.plots.data;

public class DataPipe {

  Pipe start;

  public String pipe(String value) {
    Pipe p = start;
    String v = value;
    boolean done = false;
    while (!done) {
      try {
        v = p.pipeValue(v);
        if (p.hasNextPipe()) {
          p = p.getNextPipe();
        } else {
          done = true;
        }
      } catch (RejectedValueException e) {
        done = true;
        v = null;
      }
    }
    return v;
  }

  public void addPipe(Pipe pipe) {
    if (this.start == null) {
      this.start = pipe;
      return;
    }
    Pipe p = this.start;
    while (p.hasNextPipe()) {
      p = p.getNextPipe();
    }
    if (p != null) {
      p.setNextPipe(pipe);
    }
  }

}
