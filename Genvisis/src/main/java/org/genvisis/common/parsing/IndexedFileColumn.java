package org.genvisis.common.parsing;

/**
 * An interface to define a common API for {@link FileColumn}s that point to a single index
 *
 * @param <T>
 */
public interface IndexedFileColumn<T> extends FileColumn<T> {

  /**
   * @return the index this {@link FileColumn} is pointing to
   */
  public int getIndex();
}
