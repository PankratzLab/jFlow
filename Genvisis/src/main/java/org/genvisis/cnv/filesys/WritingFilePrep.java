package org.genvisis.cnv.filesys;

public interface WritingFilePrep {

  public void close();

  public void init();

  public boolean overwriteExisting();

  public boolean validate();

}
