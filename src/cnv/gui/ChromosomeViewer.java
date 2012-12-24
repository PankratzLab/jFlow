package cnv.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.util.Vector;

import javax.swing.JPanel;

import filesys.GeneData;
import filesys.GeneTrack;
import filesys.Segment;

public class ChromosomeViewer extends JPanel {
	public static final long serialVersionUID = 8L;

	public static final int WIDTH_BUFFER = 25;

	private int chromosome = 1;
	private int startPosition = 0;
	private int stopPosition = 0;
	private GeneTrack track = null; // Represents where the genes are

	public ChromosomeViewer(int chr, int start, int stop, GeneTrack gt) {
		chromosome = chr;
		startPosition = start;
		stopPosition = stop;
		track = gt;

		// setMaximumSize(new Dimension(getWidth(), 20));
	}

	public void setChr(int chr) {
		chromosome = chr;
	}

	public void setStart(int start) {
		startPosition = start;
	}

	public void setStop(int stop) {
		stopPosition = stop;
	}

	public void setTrack(GeneTrack track) {
		this.track = track;
	}

	public void updateView(int chr, int start, int stop) {
		chromosome = chr;
		startPosition = start;
		stopPosition = stop;
		repaint();
	}

	public void paintComponent(Graphics g) {
		GeneData[] genes;
		int[][] exons;
		Vector<Segment> v = new Vector<Segment>();
		Vector<Segment> segVector = new Vector<Segment>();
		Segment[] segs, segArray;
		int width, begin, end;

		System.out.println("Start = " + startPosition + "  Stop = " + stopPosition + " Difference is " + (stopPosition - startPosition));
		if (track == null) {
			System.out.println("Track is null");
		}

		// // if (track != null && stopPosition - startPosition < 10 000 000) {
		// if (track != null) {
		// genes = track.getBetween(chromosome, startPosition, stopPosition, 30);
		// System.out.println("There are " + genes.length + " genes in this segment");
		//
		// g.setColor(Color.BLACK);
		// // Draw the chromosomes
		// for (int i = 0; i < genes.length; i++) {
		// int geneLength = genes[i].getStop() - genes[i].getStart();
		// System.out.println("Gene length is " + geneLength + " before getXY()");
		// begin = getX(genes[i].getStart());
		// end = getX(genes[i].getStop());
		//
		// geneLength = end - begin;
		// System.out.println("Gene length is " + geneLength + " after getXY()");
		//
		// System.out.println("Gene " + i + " starts at " + begin + " and ends at " + end);
		//
		// g.drawRoundRect(begin, 100, geneLength, 10, 2, 2);
		//
		// segVector.add(new Segment(begin, end));
		//
		// exons = genes[i].getExonBoundaries();
		// for (int j = 0; j < exons.length; j++) {
		// begin = getX(exons[j][0]);
		// end = getX(exons[j][1]);
		//
		// if (j == 0 || j == exons.length - 1) {
		// g.fillRoundRect(begin, 100, geneLength + 1, 10, 2, 2);
		// } else {
		// g.fillRect(begin, 100, geneLength + 1, 10);
		// }
		//
		// }
		// }
		//
		// Segment.mergeOverlapsAndSort(segVector);
		// segArray = Segment.toArray(segVector);
		// g.setFont(new Font("Arial", 0, 14));
		//
		// for (int i = 0; i < genes.length; i++) {
		// begin = getX(genes[i].getStart());
		// width = g.getFontMetrics(g.getFont()).stringWidth(genes[i].getGeneName());
		// if (!Segment.overlapsAny(new Segment(begin - width - 5, begin - 1), segArray)) {
		// g.drawString(genes[i].getGeneName(), begin - width - 3, 0 * 15 + 10);
		// }
		// }
		//
		// }

		if (track != null && stopPosition - startPosition < 10000000) {
			genes = track.getBetween(chromosome, startPosition, stopPosition, 30);
			// System.out.println(ext.getUCSCformat(new int[] {chr, start, stop}));
			g.setColor(Color.BLACK);
			for (int i = 0; i < genes.length; i++) {
				begin = getX(genes[i].getStart());
				end = getX(genes[i].getStop());
				g.drawRoundRect(begin, 0 * 15, end - begin, 10, 2, 2);
				v.add(new Segment(begin, end));
				exons = genes[i].getExonBoundaries();
				for (int j = 0; j < exons.length; j++) {
					begin = getX(exons[j][0]);
					end = getX(exons[j][1]);
					if (j == 0 || j == exons.length - 1) {
						g.fillRoundRect(begin, 0 * 15, end - begin + 1, 10, 2, 2);
					} else {
						g.fillRect(begin, 0 * 15, end - begin + 1, 10);
					}

				}
				// System.out.println(genes[i].getGeneName()+"\t"+genes[i].getStart()+"\t"+genes[i].getStop());
			}
			// System.out.println();
			Segment.mergeOverlapsAndSort(v);
			segs = Segment.toArray(v);
			g.setFont(new Font("Arial", 0, 14));

			for (int i = 0; i < genes.length; i++) {
				begin = getX(genes[i].getStart());
				width = g.getFontMetrics(g.getFont()).stringWidth(genes[i].getGeneName());
				if (!Segment.overlapsAny(new Segment(begin - width - 5, begin - 1), segs)) {
					g.drawString(genes[i].getGeneName(), begin - width - 3, 0 * 15 + 10);
				}
			}

		}

		// currentView = new Segment(chromosome, startPosition, stopPosition);
		// int firstBegin;
		// for (int i = 0; i<cnvs.length; i++) {
		// source = i;
		// firstBegin = Integer.MAX_VALUE;
		// for (int j = 0; j<cnvs[i].length; j++) {
		// if (cnvs[i][j].overlaps(currentView)) {
		// begin = getX(cnvs[i][j].getStart());
		// if (begin < firstBegin) {
		// firstBegin = begin;
		// }
		// end = getX(cnvs[i][j].getStop());
		// g.setColor(colorScheme[source][cnvs[i][j].getCN()<2?0:1]);
		// g.fillRoundRect(begin, (source+2)*15, end-begin+1, 10, 2, 2);
		// // g.drawString(ext.rootOf(cnvLabels[source]), begin+2, (source+1)*15+10);
		// }
		// }
		// if (firstBegin != Integer.MAX_VALUE) {
		// g.drawString(ext.rootOf(cnvLabels[source]), firstBegin-Grafik.getTextWidth(cnvLabels[source], g)-3, (source+2)*15+10);
		// }
		// }
	}

	public int getX(int pos) {
		return (int) ((double) (pos - startPosition) / (double) (stopPosition - startPosition) * (double) (getWidth() - 2 * WIDTH_BUFFER)) + WIDTH_BUFFER;
	}
}
