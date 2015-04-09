package cnv.gui;

import java.awt.Desktop;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import cnv.filesys.Project;
import cnv.manage.UCSCtrack;
import cnv.plots.CompPlot;
import cnv.var.Region;

import common.Grafik;
import common.Positions;

public class RegionNavigator extends JPanel implements ActionListener {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private JTextField textField;
	JButton firstButton, leftButton, rightButton, lastButton; // Navigation buttons
	JButton UCSCButton; // Launches a browser instance
	JButton BEDButton; // Generates, compresses, and uploads a BED file base on the selected CNV file
	JButton LRRButton; // Launches a widget to compute median LRR values across regions of interest
	JLabel location;
	String[] regionsList; // List of the region files
	Project proj;
	CompPlot plot;
	Vector<Region> regions = new Vector<Region>();
	int regionIndex = 0;
	int lastRegionIndex = 0;

	public static final String DEFAULT_LOCATION = "chr6:161,624,000-163,776,000"; // PARK2 region

	/**
	 * Create the panel.
	 */
	public RegionNavigator(CompPlot cp) {
		this.proj = cp.getProject();
		plot = cp;
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

		UCSCButton = new JButton("To UCSC");
		if (Desktop.isDesktopSupported()) {
			UCSCButton.setToolTipText("View this location on UCSC in a browser");
			UCSCButton.addActionListener(this);
			UCSCButton.setEnabled(true);
		} else {
			UCSCButton.setToolTipText("Browser operations are not supported");
			UCSCButton.setEnabled(false);
		}
		add(UCSCButton);

		BEDButton = new JButton("Upload to UCSC");
		if (Desktop.isDesktopSupported()) {
			BEDButton.setToolTipText("Generate and upload a .BED file to UCSC");
			BEDButton.addActionListener(this);
			BEDButton.setEnabled(true);
		} else {
			BEDButton.setToolTipText("Browser operations are not supported");
			BEDButton.setEnabled(false);
		}
		add(BEDButton);
		
		LRRButton = new JButton("Median LRR");
		LRRButton.setToolTipText("Compute median Log R Ratios for a region");
		LRRButton.addActionListener(this);
		add(LRRButton);
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

	public void setLocation(int[] location) {
		textField.setText(Positions.getUCSCformat(location));
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
		} else if (source.equals(UCSCButton)) {
			Desktop desktop = Desktop.getDesktop();
			String URL = Positions.getUCSClink(Positions.parseUCSClocation(textField.getText()));

			// UCSC uses chrX and chrY instead of 23 and 24
			URL = URL.replaceAll("chr23", "chrX");
			URL = URL.replaceAll("chr24", "chrY");
			try {
				URI uri = new URI(URL);
				System.out.println("Browsing to " + URL);
				desktop.browse(uri);
			} catch (URISyntaxException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		} else if (source.equals(BEDButton)) {
			// Figure out which files are selected
			// Only allow upload if one file is selected (JDialog warning if multiples)
			ArrayList<String> files = plot.getFilterFiles();
			if (files.size() != 1) {
				JOptionPane.showMessageDialog(null, "One and only one file must be selected before a .BED File can be generated", "Error", JOptionPane.ERROR_MESSAGE);
			} else {
				// Find the full path to the selected file
				String[] filePaths = proj.getFilenames(proj.CNV_FILENAMES);
				String compressedFile = "";
				for (String file : filePaths) {
					if (file.endsWith(files.get(0))) {
						System.out.println("File path is " + file);
						compressedFile = file + ".bed..gz";
						// Generate BED file with:
						UCSCtrack.makeTrack(file, file, proj.getLog());
						break;
					}
				}

				// Direct the user to the BED upload page at UCSC Genome Browser
				Desktop desktop = Desktop.getDesktop();
				String URL = Positions.getUCSCUploadLink(Positions.parseUCSClocation(textField.getText()), compressedFile);

				// UCSC uses chrX and chrY instead of 23 and 24
				URL = URL.replaceAll("chr23", "chrX");
				URL = URL.replaceAll("chr24", "chrY");
				try {
					URI uri = new URI(URL);
					System.out.println("Browsing to " + URL);
					desktop.browse(uri);
				} catch (URISyntaxException e) {
					e.printStackTrace();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		if (lastRegionIndex != regionIndex) {
			setRegion(regionIndex);
		} else if (source.equals(textField)) {
			// TODO: Provide facility to export the list of regions, and expand the list of regions any time someone enters one manually
			int[] newLocation = Positions.parseUCSClocation(textField.getText());
			if ((newLocation[0] < 0) || (newLocation[1] < 0) || (newLocation[2] < 0)) {
				JOptionPane.showMessageDialog(this.getParent(), "Invalid UCSC location - " + textField.getText());
			} else {
				setLocation(newLocation);

				// Pass along the property change
				firePropertyChange("location", regions.get(lastRegionIndex), new Region(newLocation));
			}
		} else if (source.equals(LRRButton)) {
			new Thread(new LRRComp(proj, textField.getText())).start();
		}
	}

	/**
	 * Load the regions from the specified files in the project
	 */
	public void loadRegions() {
		BufferedReader reader;
		// Get a list of the regions
		regionsList = proj.getFilenames(proj.REGION_LIST_FILENAMES);
		regions = new Vector<Region>();

		try {
			if (regionsList.length == 0) {
				System.out.println("No regions file defined, using default of " + DEFAULT_LOCATION);
			} else {
				// Read each file line by line, format is:
				// chromosome:start-end/tlabel
				for (int i = 0; i < regionsList.length; i++) {
					// System.out.println("Parsing file " + regionsList[i]);
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
			// System.out.println("Setting default location");
			regions.add(new Region(DEFAULT_LOCATION));
		}
	}
}
