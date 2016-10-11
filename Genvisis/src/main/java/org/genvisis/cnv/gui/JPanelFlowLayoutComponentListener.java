package org.genvisis.cnv.gui;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;

import javax.swing.JPanel;

public class JPanelFlowLayoutComponentListener implements ComponentListener {
	@Override
	public void componentHidden(ComponentEvent e) {}

	@Override
	public void componentMoved(ComponentEvent e) {}

	@Override
	public void componentResized(ComponentEvent e) {
		// System.out.println("We've been resized!");
		figureOutOptimalPreferredSize(e);
	}

	@Override
	public void componentShown(ComponentEvent e) {
		// System.out.println("We're being shown!");
		figureOutOptimalPreferredSize(e);
	}

	public void figureOutOptimalPreferredSize(ComponentEvent e) {
		Component[] comps = ((JPanel) e.getComponent()).getComponents();
		int actualWidth, width, baseHeight, rows, height;


		baseHeight = 0;
		for (Component comp : comps) {
			height = 5 + comps[comps.length - 1].getHeight() + 5;
			if (height > baseHeight) {
				baseHeight = height;
			}
		}

		actualWidth = e.getComponent().getWidth();
		width = 5;
		rows = 1;
		for (Component comp : comps) {
			width += comp.getWidth() + 5;
			if (width > actualWidth) {
				rows++;
				width = 5 + comp.getWidth() + 5;
			}
		}

		e.getComponent().setPreferredSize(new Dimension(0, rows * baseHeight));
	}

}

//
// public void componentResized(ComponentEvent e) {
// Component[] comps, levelComps;
// int actualWidth, width, baseHeight, rows, height, totalHeight;
//
// levelComps = ((JPanel)e.getComponent()).getComponents();
// if (levelComps.length != levels) {
// System.err.println("Error - problem with the levels");
// }
// totalHeight = 0;
// for (int i = 0; i<levelComps.length; i++) {
// comps = ((JPanel)levelComps[i]).getComponents();
//
// baseHeight = 0;
// for (int j = 0; j<comps.length; j++) {
// height = 5 + comps[comps.length-1].getHeight() + 5;
// if (height > baseHeight) {
// baseHeight = height;
// }
// }
//
// actualWidth = e.getComponent().getWidth();
// width = 5;
// rows = 1;
// for (int j = 0; j<comps.length; j++) {
// width += comps[j].getWidth()+5;
// if (width > actualWidth) {
// rows++;
// width = 5+comps[j].getWidth()+5;
// }
// }
// totalHeight += rows*baseHeight;
// levelComps[i].setPreferredSize(new Dimension(0, rows*baseHeight));
// }
//
// e.getComponent().setPreferredSize(new Dimension(0, totalHeight));
// }
