/**
 * 
 */
package cnv.plots;

import java.awt.GridLayout;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

import cnv.var.CNVariant;

import common.Positions;

import filesys.GeneTrack;

/**
 * @author Michael Vieths
 * 
 */
public class CompPanel extends AbstractPanel implements MouseListener, MouseMotionListener {
	public static final long serialVersionUID = 1L;
	private int stop, start;
	private byte chr;
	private CNVariant[][] cnvs;
	private String[] cnvLabels;
	public static final String DEFAULT_LOCATION = "chr17:55,609,472-55,824,368"; // USP32

	private GeneTrack track;

	public CompPanel(CompPlot sp) {
		setLayout(new GridLayout(3, 1));

		chr = (byte) Positions.parseUCSClocation(DEFAULT_LOCATION)[0];

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
}
