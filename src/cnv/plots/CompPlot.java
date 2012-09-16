/**
 * 
 */
package cnv.plots;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JFrame;

import cnv.filesys.Project;
import cnv.var.CNVariant;
import cnv.var.CNVariantHash;

import common.Positions;

/**
 * @author Michael Vieths
 * 
 */
public class CompPlot extends JFrame implements ActionListener {
	public static final long serialVersionUID = 1L;

	// public static final String DEFAULT_LOCATION = "chr17:55,609,472-55,824,368"; // USP32
	public static final String DEFAULT_LOCATION = "chr17:15,609,472-40,824,368"; // USP32

	Project proj;
	String[] files;

	public CompPlot(Project proj) {
		this.proj = proj;

		// Position
		// Starting position is a defined constant

		// Get a list of the .cnv files
		files = proj.getFilenames(Project.CNV_FILENAMES);
		System.out.println("There are " + files.length + " files");
		for (int i = 0; i < files.length; i++) {
			System.out.println("  " + files[i]);
		}

		// Read the data from the CNV files
		CNVariantHash[] cnvhs = new CNVariantHash[files.length];
		CNVariant[][] variants = new CNVariant[files.length][];

		// Parse out the location chromosome/start base/end base
		int[] location = Positions.parseUCSClocation(DEFAULT_LOCATION);
		for (int i = 0; i < cnvhs.length; i++) {
			// Load the CNVs out of the files
			cnvhs[i] = CNVariantHash.load(files[i], CNVariantHash.CONSTRUCT_ALL, false);
			// All CNVs are loaded into an array, which is then stored in an array by file
			// Format is:
			// file1: CNV1..CNVN
			// fileN: CNV1..CNVN
			// All CNVs in each file, when they are rendered, will be a single color
			// Any CNVs where CNVariant.getCN() > 2 will be a different shade of that color
			variants[i] = cnvhs[i].getAllInRegion((byte) location[0], location[1], location[2]);
		}

		for (int i = 0; i < cnvhs.length; i++) {
			System.out.println("=== variant[" + i + "]");
			for (int j = 0; j < variants[i].length; j++) {
				CNVariant x = variants[i][j];
				System.out.println("  Variant" + j + " at " + x.getUCSClocation() + " has " + x.getCN() + " copies");
			}
		}

		// setupGUI();
	}

	private void setupGUI() {
		// Set the default window size
		setSize(1000, 720);

		// Close this window but not the entire application on close
		// setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		// Close the whole thing for debugging purposes
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		// Create the panel
		CompPanel compPanel = new CompPanel(this);
		add(compPanel);
		// Set the panel visible
		setVisible(true);
	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		// TODO Auto-generated method stub

	}

}
