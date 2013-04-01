package cnv.gui;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import cnv.filesys.Project;
import cnv.var.Region;

import common.Grafik;

public class RegionNavigator extends JPanel implements ActionListener {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private JTextField textField;
	JButton firstButton, leftButton, rightButton, lastButton;
	JLabel location;
	String[] regionsList; // List of the region files
	Project proj;
	Vector<Region> regions = new Vector<Region>();
	int regionIndex = 0;
	int lastRegionIndex = 0;

	public static final String DEFAULT_LOCATION = "chr6:161,624,000-163,776,000"; // PARK2 region

	/**
	 * Create the panel.
	 */
	public RegionNavigator(Project proj) {
		this.proj = proj;
		initButtons();

		// Parse the files and set up the regions Vector
		// This will be full of regions from the file, or the DEFAULT_LOCATION
		loadRegions();

		// for (int i = 0; i < regions.size(); i++) {
		// System.out.println("=(" + i + ")= " + regions.get(i).getRegion() + "\t" + regions.get(i).getLabel());
		// }

		// Set region to the first region in the list
		setRegion(regionIndex);
	}

	/**
	 * Initialize the graphics and add listeners
	 */
	public void initButtons() {
		firstButton = new JButton(Grafik.getImageIcon("images/firstLast/dFirst.gif", true));
		firstButton.addActionListener(this);
		add(firstButton);

		leftButton = new JButton(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
		leftButton.addActionListener(this);
		add(leftButton);

		textField = new JTextField();
		textField.setFont(new Font("Tahoma", Font.PLAIN, 14));
		textField.addActionListener(this);
		add(textField);
		textField.setColumns(20);

		rightButton = new JButton(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
		rightButton.addActionListener(this);
		add(rightButton);

		lastButton = new JButton(Grafik.getImageIcon("images/firstLast/dLast.gif", true));
		lastButton.addActionListener(this);
		add(lastButton);

		location = new JLabel();
		add(location);
	}

	/**
	 * 
	 * @return A Region object representing the current region
	 */
	public Region getRegion() {
		return regions.get(regionIndex);
	}

	/**
	 * Updates the region (text field, button tooltips, etc.) to reflect the new current position
	 * 
	 * @param index
	 *            The index in the regions Vector to which this should be set
	 */
	public void setRegion(int index) {
		regionIndex = index;

		// Update the text field
		textField.setText(regions.get(regionIndex).getRegion());

		// Update the region indicator
		location.setText("Region " + (regionIndex + 1) + " of " + regions.size());

		// Pass along the property change
		firePropertyChange("location", regions.get(lastRegionIndex), regions.get(regionIndex));

		// Set the tooltip text on the buttons to match the region label if any
		firstButton.setToolTipText(regions.get(0).getLabel());
		if (regionIndex >= 1) {
			leftButton.setToolTipText(regions.get(regionIndex - 1).getLabel());
		} else {
			leftButton.setToolTipText(regions.get(0).getLabel());
		}
		if (index < (regions.size() - 1)) {
			rightButton.setToolTipText(regions.get(regionIndex + 1).getLabel());
		} else {
			rightButton.setToolTipText(regions.get(regions.size() - 1).getLabel());
		}
		lastButton.setToolTipText(regions.get(regions.size() - 1).getLabel());
	}

	@Override
	/**
	 * Set the region based on which buttons are pressed
	 */
	public void actionPerformed(ActionEvent arg0) {
		JComponent source = (JComponent) arg0.getSource();
		lastRegionIndex = regionIndex;
		if (source.equals(firstButton)) {
			regionIndex = 0;
		} else if (source.equals(leftButton)) {
			if (regionIndex > 0) {
				regionIndex--;
			}
		} else if (source.equals(rightButton)) {
			if (regionIndex < (regions.size() - 1)) {
				regionIndex++;
			}
		} else if (source.equals(lastButton)) {
			regionIndex = regions.size() - 1;
		}
		setRegion(regionIndex);
		// else if (source.equals(textField)) {
		// // TODO: Provide facility to export the list of regions, and expand the list of regions any time someone enters one manually
		// System.out.println("Changed region to " + textField.getText());
		// setRegion(new Region(textField.getText()));
		// // regions.add(currentRegion);
		// }
	}

	/**
	 * Load the regions from the specified files in the project
	 */
	public void loadRegions() {
		BufferedReader reader;
		// Get a list of the regions
		regionsList = proj.getFilenames(Project.REGION_LIST_FILENAMES);
		regions = new Vector<Region>();

		try {
			if (regionsList.length == 0) {
				System.out.println("No regions file defined, using default of " + DEFAULT_LOCATION);
			} else {
				// Read each file line by line, format is:
				// chromosome:start-end/tlabel
				for (int i = 0; i < regionsList.length; i++) {
					System.out.println("Parsing file " + regionsList[i]);
					// TODO Regex to ensure the line is formatted correctly
					reader = new BufferedReader(new FileReader(regionsList[i]));
					String line = null;
					while ((line = reader.readLine()) != null) {
						if (line != "") {
							Region myRegion = new Region(line);
							regions.add(myRegion);
						}
					}
				}
			}
		} catch (Exception ex) {
			System.out.println(ex.getStackTrace());
		}

		if (regions.size() == 0) {
			// The file was invalid or didn't contain regions, create a default region
			System.out.println("Setting default location");
			regions.add(new Region(DEFAULT_LOCATION));
		}
	}
}
