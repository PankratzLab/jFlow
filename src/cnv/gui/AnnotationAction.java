package cnv.gui;

import java.awt.event.ActionEvent;

import javax.swing.AbstractAction;

import cnv.plots.ScatterPlot;

public class AnnotationAction extends AbstractAction {
	public static final long serialVersionUID = 1L;
	
	private ScatterPlot sp;
	private char c;
	
	public AnnotationAction(ScatterPlot sp, char c) {
		this.sp = sp;
		this.c = c;
	}

	public void actionPerformed(ActionEvent e) {
		sp.toggleAnnotationBox(c);
	}
}
