package cnv.plots;

import java.awt.BorderLayout;
import java.util.ArrayList;

import javax.swing.JFrame;

import common.Logger;

public class HistogramPlot extends JFrame {
	
	private static final boolean TESTING = true;
	
	private Logger log;
	private HistogramPanel histPanel;
	private ArrayList<HistData> data = new ArrayList<HistData>();
	
	
	
	public HistogramPlot(Logger log) {
		super("Genvisis - Histogram Plot");
//		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		this.log = log;
		setup();
	}

	private void setup() {
		setLayout(new BorderLayout());
		histPanel = new HistogramPanel(this, log);
		add(histPanel, BorderLayout.CENTER);
		setBounds(20, 20, 1000, 600);
//		pack();
		setVisible(true);
		if (TESTING) {
			setTestData();
		}
	}
	
	protected ArrayList<HistData> getData() {
		return this.data;
	}
	
	private void setTestData() {
		this.data.add(new HistData(-0.003, 2));
		this.data.add(new HistData(-0.0028, 3));
		this.data.add(new HistData(-0.0026, 6));
		this.data.add(new HistData(-0.0024, 14));
		this.data.add(new HistData(-0.0022, 20));
		this.data.add(new HistData(-0.002, 43));
		this.data.add(new HistData(-0.0018, 60));
		this.data.add(new HistData(-0.0016, 81));
		this.data.add(new HistData(-0.0014, 95));
		this.data.add(new HistData(-0.0012, 105));
		this.data.add(new HistData(-0.001, 62));
		this.data.add(new HistData(-0.0008, 58));
		this.data.add(new HistData(-0.0006, 16));
		this.data.add(new HistData(-0.0004, 11));
		this.data.add(new HistData(-0.0002, 7));
		this.data.add(new HistData(0, 3));
		this.data.add(new HistData(0.0002, 2));
		this.data.add(new HistData(0, 0));
	}
	
	public static void main(String[] args) {
		new HistogramPlot(null);
	}
	
}
