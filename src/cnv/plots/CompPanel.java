/**
 * 
 */
package cnv.plots;

import java.awt.BorderLayout;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JPanel;

import cnv.gui.ChromosomeViewer;
import cnv.gui.CompConfig;
import cnv.gui.FileNavigator;
import cnv.gui.RegionNavigator;

import common.Positions;

/**
 * @author Michael Vieths
 * 
 */
public class CompPanel extends AbstractPanel implements MouseListener, MouseMotionListener, KeyListener, ActionListener {
	public static final long serialVersionUID = 1L;
	private int[] location = new int[3];

	RegionNavigator regionNavigator;
	FileNavigator fileNavigator;
	CompConfig compConfig;
	ChromosomeViewer chromosomeViewer;

	public CompPanel(CompPlot sp) {
		location = Positions.parseUCSClocation(sp.getGeneLocation());

		setLayout(new BorderLayout());

		JPanel topPanel = new JPanel();
		topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.Y_AXIS));

		regionNavigator = new RegionNavigator(sp.getGeneLocation());
		topPanel.add(regionNavigator);

		fileNavigator = new FileNavigator(sp.getFiles());
		topPanel.add(fileNavigator);

		add(topPanel, BorderLayout.PAGE_START);

		chromosomeViewer = new ChromosomeViewer(location[0], location[1], location[2], sp.track);
		add(chromosomeViewer, BorderLayout.LINE_START);

		compConfig = new CompConfig();
		add(compConfig, BorderLayout.LINE_END);

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see cnv.plots.AbstractPanel#generatePoints()
	 */
	@Override
	void generatePoints() {
		// TODO Auto-generated method stub

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see cnv.plots.AbstractPanel#assignAxisLabels()
	 */
	@Override
	void assignAxisLabels() {
		// TODO Auto-generated method stub

	}

	@Override
	void highlightPoints() {
		// TODO Auto-generated method stub

	}

	@Override
	public void keyPressed(KeyEvent arg0) {
		// TODO Auto-generated method stub

	}

	@Override
	public void keyReleased(KeyEvent arg0) {
		if (arg0.getSource().equals(regionNavigator)) {
			// Update the location
			if (arg0.getID() == KeyEvent.VK_ENTER) {
				location = Positions.parseUCSClocation(regionNavigator.getRegion());
				System.out.println("Location is now chr = " + location[0] + " start = " + location[1] + " stop = " + location[2]);
				chromosomeViewer.updateView(location[0], location[1], location[2]);
			}
		}
	}

	@Override
	public void keyTyped(KeyEvent arg0) {
		// TODO Auto-generated method stub

	}

	public void actionPerformed(ActionEvent arg0) {
		JComponent source = (JComponent) arg0.getSource();
		System.out.println("CompPanel received an " + arg0.getActionCommand() + " event from " + source);
	}

	public void paintComponent(Graphics g) {
		location = Positions.parseUCSClocation(regionNavigator.getRegion());
		System.out.println("Location is now chr = " + location[0] + " start = " + location[1] + " stop = " + location[2]);
		chromosomeViewer.updateView(location[0], location[1], location[2]);
	}
}
