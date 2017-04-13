package org.genvisis.cnv;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;
import java.util.Scanner;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.SwingConstants;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class LaunchProperties {
	public static final long serialVersionUID = 1L;

	private static String defaultDir = System.getProperty("user.home")
																		 + System.getProperty("file.separator") + ".genvisis"
																		 + System.getProperty("file.separator");

	private static String propertiesFile = defaultDir + "launch.properties";

	private static String customPropertiesDir = defaultDir;

	public static class LaunchKey {
		private final String def;
		private final boolean isDir;
		private final String keyName;

		public LaunchKey(String defaultValue, boolean isDir, String keyName) {
			def = defaultValue;
			this.isDir = isDir;
			this.keyName = keyName;
		}

		@Override
		public String toString() {
			return keyName;
		}

		/**
		 * @return Default value to use for this property if it is not present in the properties file
		 */
		public String defaultValue() {
			if (isDir && Files.isRelativePath(def))
				return directoryOfLaunchProperties() + def;
			return def;
		}

		/**
		 * @return true iff this property is a directory path that should be created on disk if it
		 *         doesn't already exist
		 */
		public boolean isDir() {
			return isDir;
		}
	}

	public static class DefaultLaunchKeys {
		public static LaunchKey PROJECTS_DIR = new LaunchKey("projects/", true,
																												 "PROJECTS_DIR");
		public static LaunchKey DEBUG_PROJECT_FILENAME = new LaunchKey("DEBUG_PROJECT", false,
																																	 "DEBUG_PROJECT_FILENAME");
		public static LaunchKey LAST_PROJECT_OPENED = new LaunchKey(Project.EXAMPLE_PROJ
																																+ ".properties", false,
																																"LAST_PROJECT_OPENED");
		public static LaunchKey RESOURCES_DIR = new LaunchKey("resources/", true,
																													"RESOURCES_DIR");

		public static LaunchKey[] values() {
			LaunchKey[] v = {PROJECTS_DIR, DEBUG_PROJECT_FILENAME, LAST_PROJECT_OPENED, RESOURCES_DIR};
			return v;
		}
	}


	private static Properties props = null;

	/**
	 * Ensure all properties are loaded, filling in default values if needed.
	 */
	private static synchronized void init() {
		if (props == null) {
			props = new Properties();
			try {
				File propFile = new File(propertiesFile);

				// if the file does not exist, this may be a first start up, so prompt the user for a path.
				if (!propFile.exists()) {
					// check for a custom path
					if (!hasCustomPath()) {
						// check for a file local to the .jar
						File localProps = new File("launch.properties");
						if (localProps.exists()) {
							propFile = localProps;
							customPropertiesDir = ".";
							updatePropertiesFile("launch.properties");
						} else {
							Files.ensurePathExists(defaultDir);
							JFileChooser jfc = new JFileChooser(defaultDir);
							jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
							jfc.setDialogTitle("Choose launch properties directory:");
							jfc.setMultiSelectionEnabled(false);
							int resp = jfc.showDialog(null, "Select");
							if (resp == JFileChooser.APPROVE_OPTION) {
								String newPath = jfc.getSelectedFile().getAbsolutePath();

								customPropertiesDir = newPath + System.getProperty("file.separator");
								updatePropertiesFile(customPropertiesDir + "launch.properties");
								propFile = new File(propertiesFile);
							}
						}
					} else {
						propertiesFile = customPropertiesDir + "launch.properties";
						propFile = new File(propertiesFile);
					}
					saveCustomPath();
				}

				if (!propFile.exists() && !propFile.createNewFile()) {
					System.err.println("Failed to create launch properties: " + propertiesFile
														 + ". Genvisis can not continue, please check error logs.");
					System.exit(-1);
				}

				InputStream is = new FileInputStream(propFile);
				props.load(is);
				is.close();

				for (LaunchKey k : DefaultLaunchKeys.values()) {
					if (!props.containsKey(k.toString())) {
						props.put(k.toString(), k.defaultValue());
					}
					String val = props.getProperty(k.toString());
					if (k.isDir() && !Files.exists(val)) {
						new Logger().reportTimeWarning("Did not detect launch property: " + k.toString()
																					 + ". Creating default: " + val);
						new File(val).mkdirs();
					}
				}
				save();
			} catch (Exception e) {
				System.err.println("Failed to load launch properties: " + propertiesFile
													 + ". Genvisis can not continue, please check error logs.");
				System.exit(-1);
			}
		}
	}

	/**
	 * Write the current properties to disk
	 */
	private static synchronized void save() {
		FileOutputStream out;

		try {
			out = new FileOutputStream(propertiesFile);
			props.store(out, null);
			out.close();
		} catch (Exception e) {
			System.err.println("Failed to save launch properties: " + propertiesFile);
		}
	}

	private static void saveCustomPath() {
		BufferedWriter out;

		try {
			out = new BufferedWriter(new FileWriter(defaultDir + "custom_path"));
			out.write(customPropertiesDir);
			out.close();
		} catch (Exception e) {
			System.err.println("Failed to save custom launch properties: " + customPropertiesDir);
		}
	}

	/**
	 * @return The value for this key in the launch properties file
	 */
	public static String get(LaunchKey key) {
		if (props == null) {
			init();
		}
		String val = props.getProperty(key.toString());
		if (key.isDir()) {
			val = ext.verifyDirFormat(val);
		}
		return val;
	}

	/**
	 * Set the value of the given key in the launch properties file
	 */
	public static void put(LaunchKey key, String value) {
		if (props == null) {
			init();
		}
		props.put(key.toString(), value);
		save();
	}

	public static synchronized String propertiesFile() {
		return propertiesFile;
	}

	public static synchronized void updatePropertiesFile(String newLoc) {
		System.err.println("Warning: changing genvisis properties file from: " + propertiesFile
											 + " to: " + newLoc + ". This could fundamentally alter behavior.");
		propertiesFile = newLoc;
	}

	/**
	 * @return A platform-independent fully-qualified path to the launch properties file
	 */
	public static String directoryOfLaunchProperties() {
		String path = null;
		try {
			path = ext.parseDirectoryOfFile(new File(propertiesFile).getCanonicalPath());
		} catch (IOException ioe) {
			path = "";
		}
		return path;
	}

	/**
	 * @return A list of the short name (no directory or {@code .properties}) for each known project.
	 */
	public static String[] getListOfProjectNames() {
		String[] projects = getListOfProjectProperties();
		String[] projectNames = new String[projects == null ? 0 : projects.length];
		for (int i = 0; i < projectNames.length; i++) {
			projectNames[i] = ext.rootOf(projects[i], true);
		}
		return projectNames;
	}

	/**
	 * @return An array of all known {@code *.properties} files.
	 */
	public static String[] getListOfProjectProperties() {
		return Files.list(get(DefaultLaunchKeys.PROJECTS_DIR), ".properties", false);
	}

	/**
	 * Open dialog to edit launch properties
	 */
	public static void openEditor() {
		// TODO
		// handy quick implementation from http://stackoverflow.com/a/790224/1027800
		// I would prefer to use ProjectPropertiesEditor styling but it would take some effort to
		// decouple
		JButton projects = new EditFileButton(DefaultLaunchKeys.PROJECTS_DIR);
		JButton resources = new EditFileButton(DefaultLaunchKeys.RESOURCES_DIR);
		final JComponent[] inputs = new JComponent[] {new JLabel(DefaultLaunchKeys.PROJECTS_DIR.toString()),
																									projects,
																									new JLabel(DefaultLaunchKeys.RESOURCES_DIR.toString()),
																									resources};
		int res = JOptionPane.showConfirmDialog(null, inputs, "Genvisis preferences",
																						JOptionPane.OK_CANCEL_OPTION,
																						JOptionPane.PLAIN_MESSAGE);
		if (JOptionPane.OK_OPTION == res) {
			put(DefaultLaunchKeys.PROJECTS_DIR, projects.getText());
			put(DefaultLaunchKeys.RESOURCES_DIR, resources.getText());
		}
	}

	/**
	 * Convenience {@link JButton} for editing specific {@link LaunchKey}s.
	 */
	private static class EditFileButton extends JButton {
		private static final long serialVersionUID = 1L;

		public EditFileButton(LaunchKey k) {
			super(get(k));
			setPreferredSize(new Dimension(350, 20));
			addActionListener(new EditFileProperty());
			setFocusPainted(false);
			setHorizontalAlignment(SwingConstants.LEFT);
			setFont(new Font("Arial", 0, 12));
		}
	}

	/**
	 * Convenience {@link ActionListener} for updating file chooser buttons
	 */
	private static class EditFileProperty implements ActionListener {

		@Override
		public void actionPerformed(ActionEvent e) {
			JButton button = (JButton) e.getSource();
			JFileChooser fc = new JFileChooser(button.getText());
			fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			if (JFileChooser.APPROVE_OPTION == fc.showOpenDialog(null)) {
				button.setText(fc.getSelectedFile().getAbsolutePath());
			}
		}

	}

	private static boolean hasCustomPath() {
		String path = defaultDir + "custom_path";

		try {
			File customPath = new File(path);
			if (customPath.exists()) {
				Scanner in = new Scanner(customPath);
				String p = in.next();
				in.close();
				customPropertiesDir = p;
				p += "launch.properties";
				File c = new File(p);
				return c.exists();
			}
		} catch (IOException e) {
			System.err.println("Failed to load custom launch properties: " + propertiesFile
												 + ". Genvisis can not continue, please check error logs.");
			System.exit(-1);
		}
		return false;
	}
}
