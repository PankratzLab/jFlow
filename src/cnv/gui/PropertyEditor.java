package cnv.gui;

import java.awt.BorderLayout;
//import java.awt.Scrollbar;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

//import javax.swing.AbstractButton;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
//import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.SpringLayout;

import cnv.filesys.Project;

import common.HashVec;

public class PropertyEditor implements ActionListener {

	public PropertyEditor(Project proj) {
    	JFrame editer;
//    	JPanel content;
    	BorderLayout mainFrameLayout;
    	JButton saveButton;
    	
    	mainFrameLayout = new BorderLayout();
    	editer = new JFrame("Properties of project '"+proj.getProjectDir()+"'");
    	editer.setBounds(250,50,850,600);
    	editer.setLayout(mainFrameLayout);
//    	editer.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		editer.pack();
		editer.setBounds(20, 20, 1000, 600);
//    	editer.addWindowListener(frame);
//		editer.setExtendedState(frame.getExtendedState()|JFrame.MAXIMIZED_BOTH);
//    	editer.add(new JScrollPane(editer.getContentPane(), JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS));
    	editer.setVisible(true);
//    	content = new JPanel();
//    	content.setLayout(contentLayout);
//    	content.setSize(600, 550);
//    	editer.setContentPane(content);
    	
    	editer.add(leftPanel(), BorderLayout.WEST);
    	editer.add(rightPanel(proj), BorderLayout.EAST);
    	saveButton = new JButton("Save");
    	saveButton.addActionListener(this);
    	editer.add(saveButton, BorderLayout.SOUTH);
    }

	private JPanel leftPanel() {
		JPanel result;
		result = new JPanel();
		result.add(new JTree());
		return result;
	}

	private JPanel rightPanel(Project proj) {
//	private JScrollPane rightPanel(Project proj) {
		JPanel result;
//		JScrollPane scroll;
    	JLabel propertyName;
    	SpringLayout contentLayout;
    	JComponent propertyContent;

    	result = new JPanel();
//    	scroll = new JScrollPane(result);
    	contentLayout = new SpringLayout();
    	result.setLayout(contentLayout);
    	String[] propertyKeys = HashVec.getKeys(proj);
    	for (int i=0; i<propertyKeys.length; i++) {
    		propertyName = new JLabel(propertyKeys[i], JLabel.TRAILING);
    		result.add(propertyName);
    		if (((String) proj.get(propertyKeys[i])).contains("TRUE") || ((String) proj.get(propertyKeys[i])).contains("FALSE")) {
    			propertyContent = new JCheckBox();
        		((JCheckBox) propertyContent).setSize(4,4);
        		((JCheckBox) propertyContent).setSelected(((String) proj.get(propertyKeys[i])).contains("TRUE"));
    		} else if (propertyKeys[i].contains("FILENAME") || propertyKeys[i].contains("DIRECTORY")) {
//    			fileOpenButton = new JButton("Open");
    			propertyContent = new JTextField(50);
        		((JTextField) propertyContent).setText((String) proj.get(propertyKeys[i]));
        		new JButton();
    		} else {
    			propertyContent = new JTextField(50);
        		((JTextField) propertyContent).setText((String) proj.get(propertyKeys[i]));
    		}
    		propertyName.setLabelFor(propertyContent);
    		result.add(propertyContent);
    		contentLayout.putConstraint(SpringLayout.EAST, propertyName, 250, SpringLayout.WEST, result.getRootPane());
    		contentLayout.putConstraint(SpringLayout.NORTH, propertyName, i*20, SpringLayout.NORTH, result.getRootPane());
    		contentLayout.putConstraint(SpringLayout.WEST, propertyContent, 260, SpringLayout.WEST, result.getRootPane());
    		contentLayout.putConstraint(SpringLayout.NORTH, propertyContent, i*20, SpringLayout.NORTH, result.getRootPane());
    	}
		return result;
//    	return scroll;
	}





//    public propertyEditor(Project proj) {
//    	JFrame editer;
////    	JPanel content;
//    	SpringLayout contentLayout;
//    	JLabel propertyName;
////    	JTextField propertyContent;
//    	JComponent propertyContent;
////    	JButton fileOpenButton;
//    	JButton saveButton;
//    	
//    	contentLayout = new SpringLayout();
//    	editer = new JFrame("Properties of project '"+proj.getProjectDir()+"'");
//    	editer.setBounds(250,50,850,600);
//    	editer.setLayout(contentLayout);
////    	editer.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//		editer.pack();
//		editer.setBounds(20, 20, 1000, 600);
////    	editer.addWindowListener(frame);
////		editer.setExtendedState(frame.getExtendedState()|JFrame.MAXIMIZED_BOTH);
//    	editer.setVisible(true);
////    	content = new JPanel();
////    	content.setLayout(contentLayout);
////    	content.setSize(600, 550);
////    	editer.setContentPane(content);
//    	
//    	String[] propertyKeys = HashVec.getKeys(proj);
//    	for (int i=0; i<propertyKeys.length; i++) {
//    		propertyName = new JLabel(propertyKeys[i], JLabel.TRAILING);
//    		editer.add(propertyName);
//    		if (((String) proj.get(propertyKeys[i])).contains("TRUE") || ((String) proj.get(propertyKeys[i])).contains("FALSE")) {
//    			propertyContent = new JCheckBox();
//        		((JCheckBox) propertyContent).setSize(4,4);
//        		((JCheckBox) propertyContent).setSelected(((String) proj.get(propertyKeys[i])).contains("TRUE"));
//    		} else if (propertyKeys[i].contains("FILENAME") || propertyKeys[i].contains("DIRECTORY")) {
////    			fileOpenButton = new JButton("Open");
//    			propertyContent = new JTextField(50);
//        		((JTextField) propertyContent).setText((String) proj.get(propertyKeys[i]));
//        		new JButton();
//    		} else {
//    			propertyContent = new JTextField(50);
//        		((JTextField) propertyContent).setText((String) proj.get(propertyKeys[i]));
//    		}
//    		propertyName.setLabelFor(propertyContent);
//    		editer.add(propertyContent);
//    		contentLayout.putConstraint(SpringLayout.EAST, propertyName, 250, SpringLayout.WEST, editer.getContentPane());
//    		contentLayout.putConstraint(SpringLayout.NORTH, propertyName, i*20, SpringLayout.NORTH, editer.getContentPane());
//    		contentLayout.putConstraint(SpringLayout.WEST, propertyContent, 260, SpringLayout.WEST, editer.getContentPane());
//    		contentLayout.putConstraint(SpringLayout.NORTH, propertyContent, i*20, SpringLayout.NORTH, editer.getContentPane());
//    	}
//    	saveButton = new JButton("Save");
//    	saveButton.addActionListener(this);
//    	editer.add(saveButton);
////		contentLayout.putConstraint(SpringLayout.NORTH, propertyContent, 20, SpringLayout.NORTH, editer.getContentPane());
//    }

	@Override
	public void actionPerformed(ActionEvent e) {
		// TODO Auto-generated method stub
		String command = e.getActionCommand();

		if (command.equals("Save")) {
			System.out.println("You just clicked Save");
		}
	}

}
