package cnv.plots;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Random;
import java.util.Set;

import common.Aliases;
import common.Array;
import common.Files;
import common.Grafik;
import common.ext;
import cnv.filesys.Project;

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
		return new Color[] { one, two };
	}

	public static Color[] generatRGBScale(int numSteps) {
		Color[] colors = new Color[numSteps];
		for (int i = 0; i < colors.length; i++) {
			int[] rgb = Grafik.getHeatmapColor((double) i / numSteps);
			colors[i] = new Color(rgb[0], rgb[1], rgb[2]);
		}
		return colors;
	}

	public static class ColorItem<E> {
		private E key;
		private Color color;

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
		if (Files.exists(file)) {
			String[] header = Files.getHeaderOfFile(file, proj.getLog());
			int markerIndex = ext.indexFactors(new String[][] { Aliases.MARKER_NAMES }, header, true, true, false, proj.getLog(), false)[0];
			int classIndex = ext.indexOfStartsWith("CLASS=MARKER_COLOR", header, false);
			if (markerIndex < 0 || classIndex < 0) {
				if (markerIndex < 0) {
					proj.getLog().reportTimeError("Could not find any of the the following in the header of " + file + "\n" + Array.toStr(Aliases.MARKER_NAMES, "\n"));
				} else {
					proj.getLog().reportTimeError("Could not find CLASS=MARKER_COLOR  in the header of " + file);

				}
				return null;

			} else {
				HashMap<String, Color> cols = parseColors(header[classIndex]);
				Hashtable<String, ColorItem<String>> manager = new Hashtable<String, ColorExt.ColorItem<String>>();
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
							proj.getLog().reportTimeWarning("Did not detect marker " + line[markerIndex] + " in the project");
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
				};
				return markerColorManager;
			}
		} else {
			proj.getLog().reportFileNotFound(file);
		}

		return null;
	}

	public static abstract class ColorManager<E> {

		private Set<E> toUse;// color categories in play
		private Hashtable<E, E> lookup;// items associated with category (marker->PoorQualityCategory)
		private Hashtable<E, ColorItem<E>> manager; // categories associated with color item (PoorQualityCategory->blue)

		public ColorManager(Hashtable<E, E> lookup, Hashtable<E, ColorItem<E>> manager) {
			this.lookup = lookup;
			this.manager = manager;
			this.toUse = manager.keySet();
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
	}

}
