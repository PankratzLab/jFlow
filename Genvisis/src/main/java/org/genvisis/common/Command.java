package org.genvisis.common;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import com.google.common.collect.ImmutableList;

public class Command {

  public static class Builder {

    private final String[] elements;
    private Collection<String> necessaryInputFiles = ImmutableList.of();
    private Collection<String> expectedOutputFiles = ImmutableList.of();
    private String dir = "";

    /**
     * @param elements the elements, spaces will be inserted between each element
     */
    public Builder(List<String> elements) {
      this.elements = elements.toArray(new String[elements.size()]);
    }

    /**
     * @param elements the elements, spaces will be inserted between each element
     */
    public Builder(String... elements) {
      this.elements = elements;
    }

    public Command build() {
      return new Command(this);
    }

    /**
     * @param necessaryInputFiles check these files for existence, will fail if they do not all
     *          exist prior to running
     * @return this {@link Builder}
     */
    public Builder necessaryInputFiles(Collection<String> necessaryInputFiles) {
      this.necessaryInputFiles = Collections.unmodifiableCollection(necessaryInputFiles);
      return this;
    }

    /**
     * @param expectedOutputFiles check these files for existence, will fail if they do not all
     *          exist after running
     * @return this {@link Builder}
     */
    public Builder expectedOutputFiles(Collection<String> expectedOutputFiles) {
      this.expectedOutputFiles = Collections.unmodifiableCollection(expectedOutputFiles);
      return this;
    }

    /**
     * @param necessaryInputFiles check these files for existence, will fail if they do not all
     *          exist prior to running
     * @return this {@link Builder}
     */
    public Builder necessaryInputFiles(String... necessaryInputFiles) {
      return necessaryInputFiles(Arrays.asList(necessaryInputFiles));
    }

    /**
     * @param expectedOutputFiles check these files for existence, will fail if they do not all
     *          exist after running
     * @return this {@link Builder}
     */
    public Builder expectedOutputFiles(String... expectedOutputFiles) {
      return expectedOutputFiles(Arrays.asList(expectedOutputFiles));
    }

    /**
     * @param dir directory to run {@code Command} in
     * @return this {@link Builder}
     */
    public Builder dir(String dir) {
      this.dir = dir;
      return this;
    }

  }

  private final String[] elements;
  private final Collection<String> necessaryInputFiles;
  private final Collection<String> expectedOutputFiles;
  private final String dir;

  private Command(Builder builder) {
    this.elements = builder.elements;
    this.necessaryInputFiles = builder.necessaryInputFiles;
    this.expectedOutputFiles = builder.expectedOutputFiles;
    this.dir = builder.dir;
  }

  /**
   * @return the elements of the command
   */
  public String[] getElements() {
    return elements;
  }

  /**
   * @return the necessary input files required to run the {@link Command}
   */
  public Collection<String> getNecessaryInputFiles() {
    return necessaryInputFiles;
  }

  /**
   * @return the expected output files required to indicate success in running the {@link Command}
   */
  public Collection<String> getExpectedOutputFiles() {
    return expectedOutputFiles;
  }

  /**
   * @return the directory to run the {@link Command} in
   */
  public String getDir() {
    return dir;
  }

  /**
   * A simple {@link Command} that does not include any settings other than the elements of the
   * command to run
   * 
   * @param elements
   * @return
   */
  public static Command basic(List<String> elements) {
    return new Builder(elements).build();
  }

  /**
   * A simple {@link Command} that does not include any settings other than the elements of the
   * command to run
   * 
   * @param elements
   * @return
   */
  public static Command basic(String... elements) {
    return new Builder(elements).build();
  }

  /**
   * A convenience for {@link Builder#Builder(List)}
   * 
   * @param elements the elements, spaces will be inserted between each element
   */
  public static Builder builder(List<String> elements) {
    return new Builder(elements);
  }

  /**
   * A convenience for {@link Builder#Builder(String...)}
   * 
   * @param elements the elements, spaces will be inserted between each element
   */
  public static Builder builder(String... elements) {
    return new Builder(elements);
  }

}
