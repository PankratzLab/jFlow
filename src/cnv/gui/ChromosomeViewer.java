package cnv.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;

import javax.swing.JPanel;

import filesys.GeneData;
import filesys.GeneTrack;
import filesys.Segment;

public class ChromosomeViewer extends JPanel implements MouseMotionListener {
	public static final long serialVersionUID = 8L;

	public static final int WIDTH_BUFFER = 25;

	private int chromosome = 1;
	private int startPosition = 0;
	private int stopPosition = 0;
	private GeneTrack track = null; // Represents where the genes are
	float scalingFactor;
	private int height = 10;
	private ArrayList<GeneRect> geneRects;
	private GeneData[] oldGenes;

	public ChromosomeViewer(int chr, int start, int stop, GeneTrack gt) {
		chromosome = chr;
		startPosition = start;
		stopPosition = stop;
		track = gt;
		geneRects = new ArrayList<GeneRect>();
		oldGenes = new GeneData[0];
		addMouseMotionListener(this);

		setMaximumSize(new Dimension(getWidth(), 20));
		setMinimumSize(new Dimension(getWidth(), 20));
	}

	/**
	 * Moves the chromosome viewer to the specified position and updates the view
	 * 
	 * @param chr
	 * @param start
	 * @param stop
	 */
	public void updateView(int chr, int start, int stop) {
		chromosome = chr;
		startPosition = start;
		stopPosition = stop;

		int window = stopPosition - startPosition;
		scalingFactor = (float) getWidth() / window;

		repaint();
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		GeneData[] genes;
		int[][] exons;
		Vector<Segment> v = new Vector<Segment>();
		int width, begin, end;
		int lastNameEnd = 0;

		if (track == null) {
			System.out.println("Track is null");
		}

		try {
			if (track != null && stopPosition - startPosition < 10000000) {
				genes = track.getBetween(chromosome, startPosition, stopPosition, 30);
				if (!Arrays.equals(genes, oldGenes)) {
					// We've changed the list of genes, so update the rectangles list as we create them
					geneRects.clear();
					oldGenes = genes;
				}
				g.setColor(Color.BLACK);
				for (int i = 0; i < genes.length; i++) {
					begin = getX(genes[i].getStart());
					end = getX(genes[i].getStop());
					width = end - begin;
					if (width < 1) {
						width = 1;
					}

					// Maintain a list of the rectangles so we can check them on mouseover
					GeneRect rectangle = new GeneRect(genes[i].getGeneName(), new Rectangle(begin, 0, width, height));
					geneRects.add(rectangle);

					g.drawRoundRect(begin, 0, width, height, 2, 2);
					v.add(new Segment(begin, end));
					exons = genes[i].getExonBoundaries();
					for (int j = 0; j < exons.length; j++) {

						int exBegin = getX(exons[j][0]);
						int exEnd = getX(exons[j][1]);
						int exWidth = exEnd - exBegin;
						if (j == 0 || j == exons.length - 1) {
							g.fillRoundRect(exBegin, 0, exWidth, height, 2, 2);
						} else {
							g.fillRect(exBegin, 0, exWidth, height);
						}

					}
				}
				Segment.mergeOverlapsAndSort(v);
				g.setFont(new Font("Arial", 0, 14));

				for (int i = 0; i < genes.length; i++) {
					begin = getX(genes[i].getStart());
					width = g.getFontMetrics(g.getFont()).stringWidth(genes[i].getGeneName());

					// Keep track of the end of the name we rendered so we don't overwrite them
					if (begin > lastNameEnd) {
						g.drawString(genes[i].getGeneName(), begin, (height * 2) + 2);
						lastNameEnd = begin + width;
					}
				}

			}
		} catch (Exception ex) {
			System.out.println("Missing data in gene track");
		}
	}

	/**
	 * Calculate the X coordinate of the given base position
	 * 
	 * @param pos
	 * @return
	 */
	public int getX(int pos) {
		return (int) ((pos - startPosition) * scalingFactor);
	}

	@Override
	public void mouseDragged(MouseEvent e) {
		// Nothing to do on mouse drag
	}

	/**
	 * Creates a tooltip indicating which genes are included in the position being moused over
	 */
	public void mouseMoved(MouseEvent e) {
		// Create a mouseover tooltip with the name of any genes within this range
		setToolTipText(null);
		ArrayList<String> genes = new ArrayList<String>();

		for (GeneRect geneRect : geneRects) {
			if (geneRect.contains(e.getPoint())) {
				String name = geneRect.getName();
				if (!genes.contains(name))
					genes.add(name);
			}
		}

		String toolTipText = "";
		for (int i = 0; i < genes.size(); i++) {
			if (i == genes.size() - 1) {
				toolTipText += genes.get(i);
			} else {
				toolTipText += genes.get(i) + ", ";
			}
		}
		setToolTipText(toolTipText);

	}
}

/**
 * A class to keep track of bounding boxes for genes
 * 
 * @author Michael Vieths
 * 
 */
class GeneRect {
	private String geneName;
	private Rectangle geneRect;

	public GeneRect(String gene, Rectangle rect) {
		geneName = gene;
		geneRect = rect;
	}

	/**
	 * Returns whether the given point falls in this gene
	 * 
	 * @param point
	 * @return
	 */
	public boolean contains(Point point) {
		return geneRect.contains(point);
	}

	/**
	 * Returns the name of this gene
	 * 
	 * @return
	 */
	public String getName() {
		return geneName;
	}
}
