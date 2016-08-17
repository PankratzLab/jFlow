package org.genvisis.cnv.filesys;

public interface ReadingFilePrep {
  public void close();

  public void init();

  public boolean validate();
}
