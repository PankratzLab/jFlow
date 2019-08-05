package org.genvisis.cnv.workflow;

import java.io.File;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.cnv.Resources.Resource;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.stats.Maths;
import org.pankratzlab.utils.gwas.MarkerQC;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;

public abstract class Requirement<T> {

  public static class StepRequirement extends Requirement<Step> {

    private final Step requiredStep;

    public StepRequirement(Step requiredStep) {
      super(stepReqMessage(requiredStep));
      this.requiredStep = requiredStep;
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      return (stepSelections.contains(requiredStep)
              && requiredStep.hasRequirements(stepSelections, variables))
             || requiredStep.checkIfOutputExists(variables.get(requiredStep));
    }

    public Step getRequiredStep() {
      return requiredStep;
    }

    private static String stepReqMessage(Step requiredStep) {
      String msg = "[" + (requiredStep == null ? "" : requiredStep.getName())
                   + "] step must have been run already or must be selected";
      return msg;
    }

    @Override
    public Step parseValue(String raw) {
      throw new RuntimeException("Error - step parsing isn't implemented.");
    }

  }

  public abstract static class ParsedRequirement<T> extends Requirement<T> {

    public ParsedRequirement(String key, String description, RequirementInputType type) {
      super(description);
      this.key = key;
      this.type = type;
    }

    public ParsedRequirement(String key, String description, Requirement.RequirementInputType type,
                             T defaultValue) {
      super(description, defaultValue);
      this.key = key;
      this.type = type;
    }

    private final String key;

    private final Requirement.RequirementInputType type;

    public Requirement.RequirementInputType getType() {
      return type;
    }

    public String getKey() {
      return key;
    }

    @Override
    public int hashCode() {
      final int prime = 31;
      int result = super.hashCode();
      result = prime * result + ((key == null) ? 0 : key.hashCode());
      result = prime * result + ((type == null) ? 0 : type.hashCode());
      return result;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (!super.equals(obj)) return false;
      if (!(obj instanceof ParsedRequirement)) return false;
      ParsedRequirement other = (ParsedRequirement) obj;
      if (key == null) {
        if (other.key != null) return false;
      } else if (!key.equals(other.key)) return false;
      if (type != other.type) return false;
      return true;
    }

  }

  public static class FileRequirement extends ParsedRequirement<File> {

    public FileRequirement(String key, String description, File defaultValue) {
      super(key, description, Requirement.RequirementInputType.FILE, defaultValue);
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      return valid(arg);
    }

    protected boolean valid(String arg) {
      return Files.exists(arg) && !Files.isDirectory(arg);
    }

    @Override
    public File parseValue(String raw) {
      return new File(raw);
    }

  }

  public static class DirRequirement extends ParsedRequirement<File> {

    public DirRequirement(String key, String description, File defaultValue) {
      super(key, description, Requirement.RequirementInputType.DIR, defaultValue);
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      return valid(arg);
    }

    protected boolean valid(String arg) {
      return Files.exists(arg) && Files.isDirectory(arg);
    }

    @Override
    public File parseValue(String raw) {
      return new File(raw);
    }

  }

  public static class OptionalDirRequirement extends DirRequirement {

    /**
     * @param description
     * @param defaultValue
     */
    public OptionalDirRequirement(String key, String description, File defaultValue) {
      super(key, description, defaultValue);
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      return "".equals(arg) || super.checkRequirement(arg, stepSelections, variables);
    }

  }

  public static class OutputFileRequirement extends FileRequirement {

    public OutputFileRequirement(String key, String description, File defaultValue) {
      super(key, description, defaultValue);
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      return !Files.exists(arg);
    }
  }

  public static class OptionalFileRequirement extends FileRequirement {

    public OptionalFileRequirement(String key, String description, File defaultValue) {
      super(key, description, defaultValue);
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      return "".equals(arg) || super.checkRequirement(arg, stepSelections, variables);
    }
  }

  public static class ResourceRequirement extends Requirement<Resource> {

    private final Resource resource;

    /**
     * @param resourceDescription a description of the {@link Resource} required. This will be used
     *          as part of the full description describing the download/use of the resource
     * @param resource the {@link Resource} required
     */
    public ResourceRequirement(String resourceDescription, Resource resource) {
      super(generateRequirementDescription(resourceDescription, resource));
      this.resource = resource;
    }

    private static String generateRequirementDescription(String resourceDescription,
                                                         Resource resource) {
      String requirementDescription;
      if (resource.isLocallyAvailable()) {
        requirementDescription = "Use locally available " + resourceDescription;
      } else {
        requirementDescription = "Download remotely available " + resourceDescription;
      }
      return requirementDescription;
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      return true;
    }

    public Resource getResource() {
      return resource;
    }

    @Override
    public Resource parseValue(String raw) {
      if (!resource.getName().equals(raw)) {
        System.err.println("Requested Resource " + raw + "; got Resource " + resource.getName());
      }
      return resource;
    }

  }

  public static class BoolRequirement extends ParsedRequirement<Boolean> {

    public BoolRequirement(String key, String description, boolean defaultValue) {
      super(key, description, Requirement.RequirementInputType.BOOL, defaultValue);
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      return Boolean.parseBoolean(arg);
    }

    @Override
    public Boolean parseValue(String raw) {
      return Boolean.parseBoolean(raw);
    }

  }

  public static class OptionalBoolRequirement extends BoolRequirement {

    public OptionalBoolRequirement(String key, String description, boolean defaultValue) {
      super(key, description, defaultValue);
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      return true;
    }
  }

  public static class DoubleRequirement extends ParsedRequirement<Double> {

    private final double min;
    private final double max;

    public DoubleRequirement(String key, String description, double defaultValue, double min,
                             double max) {
      super(key, description, Requirement.RequirementInputType.NUMBER, defaultValue);
      this.min = min;
      this.max = max;
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      double value;
      try {
        value = Double.parseDouble(arg);
      } catch (NumberFormatException e) {
        return false;
      }
      return value >= min && value <= max;
    }

    @Override
    public Double parseValue(String raw) {
      return Double.parseDouble(raw);
    }

  }

  public static class IntRequirement extends ParsedRequirement<Integer> {

    private final int min;
    private final int max;

    public IntRequirement(String key, String description, int defaultValue, int min, int max) {
      super(key, description, Requirement.RequirementInputType.NUMBER, defaultValue);
      this.min = min;
      this.max = max;
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      int value;
      try {
        value = Integer.parseInt(arg);
      } catch (NumberFormatException e) {
        return false;
      }
      return value >= min && value <= max;
    }

    @Override
    public Integer parseValue(String raw) {
      return Integer.parseInt(raw);
    }

  }

  public static class PosIntRequirement extends IntRequirement {

    public PosIntRequirement(String key, String description, int defaultValue) {
      super(key, description, defaultValue, 1, Integer.MAX_VALUE);
    }

  }

  public static class ListSelectionRequirement extends ParsedRequirement<Collection<String>> {

    private static final char SELECTION_LIST_DELIM = ',';
    private static final Joiner SELECTION_LIST_JOINER = Joiner.on(SELECTION_LIST_DELIM);
    private static final Splitter SELECTION_LIST_SPLITTER = Splitter.on(SELECTION_LIST_DELIM);

    private final Collection<String> options;
    private final boolean allowNone;

    public ListSelectionRequirement(String key, String description, Collection<String> options,
                                    Collection<String> defaultOptions, boolean allowNone) {
      super(key, description, Requirement.RequirementInputType.LISTSELECTION, defaultOptions);
      if (!options.containsAll(defaultOptions)) throw new IllegalArgumentException("All defaultOptions are not in options");
      this.options = options;
      this.allowNone = allowNone;
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      return allowNone || !arg.isEmpty();
    }

    public Collection<String> getOptions() {
      return options;
    }

    public static String createArgValString(Iterable<?> selections) {
      return SELECTION_LIST_JOINER.join(selections);
    }

    public static List<String> parseArgValString(String arg) {
      return SELECTION_LIST_SPLITTER.splitToList(arg);
    }

    @Override
    public Collection<String> parseValue(String raw) {
      return parseArgValString(raw);
    }

  }

  public static class EnumRequirement<Y extends Enum<Y>> extends ParsedRequirement<Y> {

    public EnumRequirement(String key, String description, Y defaultValue) {
      super(key, description, RequirementInputType.ENUM, defaultValue);
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      return true;
    }

    @Override
    public Y parseValue(String raw) {
      Y v = this.getDefaultValue();
      return Enum.valueOf(v.getDeclaringClass(), raw);
    }

  }

  public static class ThresholdRequirement extends ParsedRequirement<String> {

    public ThresholdRequirement(String key, String description, String defaultValue) {
      super(key, description, RequirementInputType.STRING, defaultValue);
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      Maths.COMPARISON op = MarkerQC.findOperator(arg);
      if (op == null) return false;
      try {
        Double.parseDouble(arg.substring(op.getSymbol().length()));
      } catch (NumberFormatException nfe) {
        return false;
      }
      return true;
    }

    @Override
    public String parseValue(String raw) {
      return raw;
    }

  }

  public static class ProgramRequirement extends ParsedRequirement<String> {

    public ProgramRequirement(String key, String description, String defaultValue) {
      super(key, description, Requirement.RequirementInputType.FILE, defaultValue);
    }

    @Override
    public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                    Map<Step, Variables> variables) {
      return Files.programExists(arg);
    }

    @Override
    public String parseValue(String raw) {
      return raw;
    }

  }

  public enum RequirementInputType {
    FILE, DIR, STRING, NUMBER, BOOL, ENUM, LISTSELECTION
  }

  public enum Flag {
    MEMORY, RUNTIME, MULTITHREADED
  }

  private final String description;
  private final T defaultValue;

  /**
   * @param description
   */
  public Requirement(String description) {
    this(description, null);
  }

  /**
   * @param description
   * @param defaultValue
   */
  public Requirement(String description, T defaultValue) {
    super();
    this.description = description;
    this.defaultValue = defaultValue != null ? defaultValue : null;
  }

  public static int checkIntArgOrNeg1(String val) {
    int valInt = -1;
    try {
      valInt = Integer.parseInt(val);
    } catch (NumberFormatException e) {
      // leave as -1
    }
    return valInt;
  }

  public static double checkDoubleArgOrNeg1(String val) {
    double valDou = -1.0;
    try {
      valDou = Double.parseDouble(val);
    } catch (NumberFormatException e) {
      // leave as -1.0
    }
    return valDou;
  }

  public String getDescription() {
    return description;
  }

  public T getDefaultValue() {
    return defaultValue;
  }

  public abstract boolean checkRequirement(String arg, Set<Step> stepSelections,
                                           Map<Step, Variables> variables);

  public abstract T parseValue(String raw);

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((defaultValue == null) ? 0 : defaultValue.hashCode());
    result = prime * result + ((description == null) ? 0 : description.hashCode());
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (getClass() != obj.getClass()) return false;
    Requirement<?> other = (Requirement<?>) obj;
    if (defaultValue == null) {
      if (other.defaultValue != null) return false;
    } else if (!defaultValue.equals(other.defaultValue)) return false;
    if (description == null) {
      if (other.description != null) return false;
    } else if (!description.equals(other.description)) return false;
    return true;
  }

}
