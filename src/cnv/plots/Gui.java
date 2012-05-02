package cnv.plots;

import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;

public class Gui extends JPanel {

	JPanel gui;
	
	public Gui() {
//		menuBar = new JMenuBar();
//		this.setJMenuBar(menuBar);
//		menu = new JMenu("File");
//		menuBar.add(menu);
//		menu.add(new JMenuItem("Open"));
//		menu.add(new JMenuItem("Save"));
//		menu.add(new JMenuItem("Close"));
//		menu.add(new JMenuItem("Exit"));
//		menu = new JMenu("Edit");
//		menu.add(new JMenuItem("Cut"));
//		menu.add(new JMenuItem("Copy"));
//		menu.add(new JMenuItem("Paste"));
//		menu.add(new JMenuItem("Paste Image"));
//		menu.add(new JMenuItem("Exit"));
		
		JPanel gui = new JPanel();
		gui.add(new JLabel("File"));
		gui.add(new JLabel("Edit"));
		gui.add(new JLabel("Help"));
	}
	
	public static void main(String[] args) {
		
	}

}
