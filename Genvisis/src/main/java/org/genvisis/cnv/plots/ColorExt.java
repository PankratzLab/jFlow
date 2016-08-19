package org.genvisis.cnv.plots;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Serializable;
import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Aliases;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Grafik;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;


public class ColorExt {
  public static Color generateRandomColor(Color mix) {
    Random random = new Random();
    int red = random.nextInt(256);
    int green = random.nextInt(256);
    int blue = random.nextInt(256);

    // mix the color
    if (mix != null) {
      red = (red + mix.getRed()) / 2;
      green = (green + mix.getGreen()) / 2;
      blue = (blue + mix.getBlue()) / 2;
    }

    Color color = new Color(red, green, blue);
    return color;
  }

  public static Color[] generatePastelShades() {
    Color one = generateRandomColor(new Color(255, 255, 255));
    Color two = generateRandomColor(one);
    return new Color[] {one, two};
  }

  public static Color[] generatRGBScale(int numSteps) {
    Color[] colors = new Color[numSteps];
    for (int i = 0; i < colors.length; i++) {
      int[] rgb = Grafik.getHeatmapColor((double) i / numSteps);
      colors[i] = new Color(rgb[0], rgb[1], rgb[2]);
    }
    return colors;
  }

  /**
   * This method tries to assign a color scale to a set of data, using the data's empirical
   * cumulative distribution function (CDF) to bin the data into discrete color bins. The goal with
   * using the CDF was to try to allow the color scale to highlight differences between the majority
   * of data points as opposed to simply scaling between the max and min values of the data range.
   * Maybe it accomplishes that goal, maybe it doesn't
   *
   * @param numColors corresponds to the number of different colors that will be generated and
   *        assigned to the data
   * @param data the data to assign colors to
   * @return an array of {@link Color}, one for each data point.
   */
  public static Color[] assignColorsForData(int numColors, double[] data) {
    NormalDistribution nd = new NormalDistribution(Array.mean(data, true),
                                                   Math.pow(Array.stdev(data, true), 2));
    Color[] colorsAvailable = ColorExt.generatRGBScale(numColors);
    Color[] colorsforData = new Color[data.length];
    for (int i = 0; i < data.length; i++) {
      colorsforData[i] = assignColor(data[i], nd, colorsAvailable);
    }
    return colorsforData;
  }

  /**
   * Based on a given {@link NormalDistribution}, assign a color to a given data point, selecting
   * from available colors
   */
  private static Color assignColor(double data, NormalDistribution nd, Color[] colorsAvailable) {
    int colorIndex = (int) Math.round(nd.cumulativeProbability(data) * colorsAvailable.length - 1);
    colorIndex = Math.max(0, colorIndex);
    colorIndex = Math.min(colorsAvailable.length, colorIndex);
    Color assigned = colorsAvailable[colorIndex];
    return assigned;
  }

  public static class ColorItem<E> implements Serializable {
    /**
     *
     */
    private static final long serialVersionUID = 1L;
    private final E key;
    private final Color color;

    public ColorItem(E key, Color color) {
      super();
      this.key = key;
      this.color = color;
    }

    public E getKey() {
      return key;
    }

    public Color getColor() {
      return color;
    }

  }

  private static HashMap<String, Color> parseColors(String colorHeader) {
    HashMap<String, Color> colorMap = new HashMap<String, Color>();
    String[] pts = colorHeader.split(";");
    for (int i = 1; i < pts.length; i++) {
      String code = pts[i].split("=")[0];
      String col = pts[i].split("=")[1].toLowerCase();
      Color color;
      try {
        Field field = Class.forName("java.awt.Color").getField(col);
        color = (Color) field.get(null);
      } catch (Exception e) {
        color = null; // Not defined
      }
      colorMap.put(code, color);
    }
    return colorMap;
  }

  public static HashMap<String, Color> getCols(Project proj, String file) {
    String[] header = Files.getHeaderOfFile(file, proj.getLog());
    int classIndex = ext.indexOfStartsWith("CLASS=MARKER_COLOR", header, false);
    HashMap<String, Color> cols = new HashMap<String, Color>();
    if (classIndex >= 0) {
      cols = parseColors(header[classIndex]);
    }
    return cols;

  }

  public static ColorManager<String> getColorManager(Project proj, String file) {
    String ser = ext.rootOf(file, false) + ".ser";
    if (Files.exists(ser)) {
      return ColorManager.loadStringSerial(ser, proj.getLog());
    }
    if (Files.exists(file)) {
      String[] header = Files.getHeaderOfFile(file, proj.getLog());
      int markerIndex = ext.indexFactors(new String[][] {Aliases.MARKER_NAMES}, header, true, true,
                                         false, proj.getLog(), false)[0];
      int classIndex = ext.indexOfStartsWith("CLASS=MARKER_COLOR", header, false);
      if (markerIndex < 0 || classIndex < 0) {
        if (markerIndex < 0) {
          proj.getLog().reportTimeError("Could not find any of the the following in the header of "
                                        + file + "\n" + Array.toStr(Aliases.MARKER_NAMES, "\n"));
        } else {
          proj.getLog()
              .reportTimeError("Could not find CLASS=MARKER_COLOR  in the header of " + file);

        }
        return null;

      } else {
        HashMap<String, Color> cols = parseColors(header[classIndex]);
        Hashtable<String, ColorItem<String>> manager =
                                                     new Hashtable<String, ColorExt.ColorItem<String>>();
        Hashtable<String, String> lookup = new Hashtable<String, String>();
        Hashtable<String, Integer> indices = proj.getMarkerIndices();
        try {

          for (String key : cols.keySet()) {
            manager.put(key, new ColorItem<String>(key, cols.get(key)));
          }

          BufferedReader reader = Files.getAppropriateReader(file);
          reader.readLine();
          while (reader.ready()) {

            String[] line = reader.readLine().trim().split("\t");
            if (!indices.containsKey(line[markerIndex])) {
              proj.getLog().reportTimeWarning("Did not detect marker " + line[markerIndex]
                                              + " in the project");
            }
            lookup.put(line[markerIndex], line[classIndex]);
          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          proj.getLog().reportError("Error: file \"" + file + "\" not found in current directory");
          return null;
        } catch (IOException ioe) {
          proj.getLog().reportError("Error reading file \"" + file + "\"");
          return null;
        }
        ColorManager<String> markerColorManager = new ColorManager<String>(lookup, manager) {

          /**
           *
           */
          private static final long serialVersionUID = 1L;
        };

        markerColorManager.writeSerial(ser);
        return markerColorManager;
      }
    } else {
      proj.getLog().reportFileNotFound(file);
    }

    return null;
  }

  public static abstract class ColorManager<E> implements Serializable {

    /**
     *
     */
    private static final long serialVersionUID = 1L;
    private final HashSet<E> toUse;// color categories in play
    private final Hashtable<E, E> lookup;// items associated with category
    // (marker->PoorQualityCategory)
    private final Hashtable<E, ColorItem<E>> manager; // categories associated
    // with color item
    // (PoorQualityCategory->blue)

    public ColorManager(Hashtable<E, E> lookup, Hashtable<E, ColorItem<E>> manager) {
      this.lookup = lookup;
      this.manager = manager;
      this.toUse = new HashSet<E>();
      for (E e : manager.keySet()) {
        toUse.add(e);
      }
    }

    public Set<E> getToUse() {
      return toUse;
    }

    public Hashtable<E, E> getLookup() {
      return lookup;
    }

    public Hashtable<E, ColorItem<E>> getManager() {
      return manager;
    }

    public boolean hasColorFor(E var) {
      return lookup.containsKey(var);
    }

    public ColorItem<E> getColorItemForVar(E var) {

      if (lookup.containsKey(var)) {
        E looked = lookup.get(var);
        if (toUse.contains(looked)) {
          return getColorFor(looked);
        } else {
          return null;
        }
      } else {
        return null;
      }
    }

    private ColorItem<E> getColorFor(E e) {
      if (manager.containsKey(e)) {
        return manager.get(e);
      } else {
        return null;
      }
    }

    public void writeSerial(String fileName) {
      SerializedFiles.writeSerial(this, fileName, true);
    }

    @SuppressWarnings("unchecked")
    private static ColorManager<String> loadStringSerial(String fileName, Logger log) {
      return (ColorManager<String>) SerializedFiles.readSerial(fileName, false, log, false, true);
    }

  }

}
