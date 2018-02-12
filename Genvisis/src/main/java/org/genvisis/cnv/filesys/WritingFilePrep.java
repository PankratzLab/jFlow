package org.genvisis.cnv.filesys;

public interface WritingFilePrep {

  public boolean overwriteExisting();

  public void init();

  public boolean validate();

  public void close();

}
