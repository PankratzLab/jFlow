package org.genvisis.cnv;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class LaunchProperties {
	public static final long serialVersionUID = 1L;

	private static String propertiesFile = "launch.properties";

	/**
	 * enum of all available launch properties.
	 */
	public enum LaunchKey {
		PROJECTS_DIR("projects/", true), 
		DEBUG_PROJECT_FILENAME("DEBUG_PROJECT", false), 
		LAST_PROJECT_OPENED(Project.EXAMPLE_PROJ + ".properties", false), 
		RESOURCES_DIR("resources/", true);

		private final String def;
		private final boolean isDir;

		private LaunchKey(String defaultValue, boolean isDir) {
			def = defaultValue;
			this.isDir = isDir;
		}

		/**
		 * @return Default value to use for this property if it is not present in the properties file
		 */
		public String defaultValue() {
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

	private static Properties props = null;

	/**
	 * Ensure all properties are loaded, filling in default values if needed.
	 */
	private static synchronized void init() {
		if (props == null) {
			props = new Properties();
			try {
				File propFile = new File(propertiesFile);
				if (!propFile.exists() && !propFile.createNewFile()) {
					System.err.println("Failed to create launch properties: " + propertiesFile
					                   + ". Genvisis can not continue, please check error logs.");
					System.exit(-1);
				}
				InputStream is = new FileInputStream(propFile);
				props.load(is);
				is.close();

				for (LaunchKey k : LaunchKey.values()) {
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
		                   + " to: " + newLoc + ". This could will fundamentally alter behavior.");
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
		String[] projectNames = new String[projects.length];
		for (int i = 0; i < projectNames.length; i++) {
			projectNames[i] = ext.rootOf(projects[i], true);
		}
		return projectNames;
	}

	/**
	 * @return An array of all known {@code *.properties} files.
	 */
	public static String[] getListOfProjectProperties() {
		return Files.list(get(LaunchKey.PROJECTS_DIR), ".properties", false);
	}

	/**
	 * Open dialog to edit launch properties
	 */
	public static void openEditor() {
		// TODO
		// handy quick implementation from http://stackoverflow.com/a/790224/1027800
		// I would prefer to use ProjectPropertiesEditor styling but it would take some effort to
		// decouple
		JButton projects = new EditFileButton(LaunchKey.PROJECTS_DIR);
		JButton resources = new EditFileButton(LaunchKey.RESOURCES_DIR);
		final JComponent[] inputs = new JComponent[] {new JLabel(LaunchKey.PROJECTS_DIR.toString()),
		                                              projects,
		                                              new JLabel(LaunchKey.RESOURCES_DIR.toString()),
		                                              resources};
		int res =
		        JOptionPane.showConfirmDialog(null, inputs, "Genvisis preferences",
		                                      JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE);
		if (JOptionPane.OK_OPTION == res) {
			put(LaunchKey.PROJECTS_DIR, projects.getText());
			put(LaunchKey.RESOURCES_DIR, resources.getText());
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
}
