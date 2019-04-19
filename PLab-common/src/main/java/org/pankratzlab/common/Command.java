package org.pankratzlab.common;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import javax.annotation.Nonnull;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

public class Command {

  public static class Builder {

    private String[] elements;
    @Nonnull
    private Collection<String> necessaryInputFiles = ImmutableList.of();
    @Nonnull
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
    public Builder necessaryInputFiles(@Nonnull Collection<String> necessaryInputFiles) {
      this.necessaryInputFiles = Collections.unmodifiableCollection(necessaryInputFiles);
      return this;
    }

    /**
     * @param expectedOutputFiles check these files for existence, will fail if they do not all
     *          exist after running
     * @return this {@link Builder}
     */
    public Builder expectedOutputFiles(@Nonnull Collection<String> expectedOutputFiles) {
      this.expectedOutputFiles = Collections.unmodifiableCollection(expectedOutputFiles);
      return this;
    }

    /**
     * @param necessaryInputFiles check these files for existence, will fail if they do not all
     *          exist prior to running
     * @return this {@link Builder}
     */
    public Builder necessaryInputFiles(@Nonnull String... necessaryInputFiles) {
      return necessaryInputFiles(Arrays.asList(necessaryInputFiles));
    }

    /**
     * @param expectedOutputFiles check these files for existence, will fail if they do not all
     *          exist after running
     * @return this {@link Builder}
     */
    public Builder expectedOutputFiles(@Nonnull String... expectedOutputFiles) {
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

    /**
     * Runs {@link #batch(String, Logger)} without logging
     */
    public Builder batch(String batFile) {
      return batch(batFile, null);
    }

    /**
     * Setup a batch file to run instead of running the command directly
     * 
     * @param batFile where the the command will be written
     * @param log {@link Logger} to log actual command run, null to not log
     * @return this {@link Builder}
     */
    public Builder batch(String batFile, Logger log) {
      String commandText = Joiner.on(' ').join(elements);
      if (log != null) {
        log.report(ext.getTime() + " Info - running command " + commandText + "\nUsing file "
                   + batFile);
      }
      Files.write(commandText, batFile);
      Files.chmod(batFile);
      elements = new String[] {batFile};
      return this;
    }

  }

  private final String[] elements;
  @Nonnull
  private final Collection<String> necessaryInputFiles;
  @Nonnull
  private final Collection<String> expectedOutputFiles;
  private final String dir;

  private Command(Builder builder) {
    this.elements = builder.elements;
    this.necessaryInputFiles = buildFiles(builder.dir, builder.necessaryInputFiles);
    this.expectedOutputFiles = buildFiles(builder.dir, builder.expectedOutputFiles);
    this.dir = builder.dir;
  }

  private Collection<String> buildFiles(String dir, Collection<String> files) {
    List<String> full = Lists.newArrayList();
    for (String f : files) {
      full.add(Files.isRelativePath(f) ? dir + f : f);
    }
    return full;
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
  public @Nonnull Collection<String> getNecessaryInputFiles() {
    return necessaryInputFiles;
  }

  /**
   * @return the expected output files required to indicate success in running the {@link Command}
   */
  public @Nonnull Collection<String> getExpectedOutputFiles() {
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
