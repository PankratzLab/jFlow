package org.genvisis.cnv.gui;

import java.awt.event.ActionEvent;

import javax.swing.AbstractAction;

import org.genvisis.cnv.plots.ScatterPlot;

public class AnnotationAction extends AbstractAction {
	public static final long serialVersionUID = 1L;
	
	private ScatterPlot sp;
	private char c;
	boolean add;
	
	public AnnotationAction(ScatterPlot sp, char c, boolean addNotRemove) {
		this.sp = sp;
		this.c = c;
		this.add = addNotRemove;
	}

	public void actionPerformed(ActionEvent e) {
		if (add) {
			sp.checkAnnotationBox(c);
		} else {
			sp.uncheckAnnotationBox(c);
		}
	}
}
