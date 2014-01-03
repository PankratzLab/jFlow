package cnv.gui;

import java.awt.Container;
import javax.swing.JFrame;

public class GuiManager {

	public static void disposeOfParentFrame(Container obj) {
		Container parent;
		boolean done;
		
		done = false;
		parent = obj;
		while (!done) {
//			System.out.println("parent of "+parent+" is "+parent.getParent());
			if (parent instanceof JFrame) {
				((JFrame)parent).dispose();
				done = true;
			}
			parent = parent.getParent();
		}
	}
}
