package cnv.gui;

import java.awt.BorderLayout;
//import java.awt.Scrollbar;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;

//import javax.swing.AbstractButton;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
//import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.SpringLayout;

import cnv.LaunchProperties;
import cnv.filesys.Project;

import common.Array;
import common.HashVec;

public class PropertyEditor implements ActionListener, WindowListener {
	private static final String SAVE = "Save";
	private static final String ADD = "Add New Property";

	Project proj;
//	private JLabel[] propertyNames;
	String[] propertyKeys;
    private JComponent[] propertyContents;
    private JFrame editor;
    JButton saveButton;

	public PropertyEditor(Project proj) {
		this.proj = proj;
//    	JPanel content;
    	BorderLayout mainFrameLayout;
    	mainFrameLayout = new BorderLayout();
    	this.editor = new JFrame("Properties of project '"+proj.getProjectDir()+"'");
    	this.editor.setBounds(250,50,850,600);
    	this.editor.setLayout(mainFrameLayout);
//    	this.editor.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	this.editor.pack();
    	this.editor.setBounds(20, 20, 1000, 600);
//    	this.editor.addWindowListener(frame);
//		this.editor.setExtendedState(frame.getExtendedState()|JFrame.MAXIMIZED_BOTH);
//    	this.editor.add(new JScrollPane(editer.getContentPane(), JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS));
    	this.editor.setVisible(true);
//    	content = new JPanel();
//    	content.setLayout(contentLayout);
//    	content.setSize(600, 550);
//    	editer.setContentPane(content);

    	initializePropertyList(proj);

//    	this.editor.add(leftPanel(), BorderLayout.WEST);
    	this.editor.add(rightPanel(), BorderLayout.EAST);
    	saveButton = new JButton(SAVE);
    	saveButton.addActionListener(this);
    	this.editor.add(saveButton, BorderLayout.SOUTH);
    }

	private JPanel leftPanel() {
		JPanel result;
		result = new JPanel();
		result.add(new JTree());
		return result;
	}

	private void initializePropertyList (Project proj) {
    	propertyKeys = HashVec.getKeys(proj);
//		propertyNames = new JLabel[propertyKeys.length];
		propertyContents = new JComponent[propertyKeys.length];

		for (int i=0; i<propertyKeys.length; i++) {
//    		propertyNames[i] = new JLabel(propertyKeys[i], JLabel.TRAILING);

    		if (((String) proj.get(propertyKeys[i])).contains("TRUE") || ((String) proj.get(propertyKeys[i])).contains("FALSE")) {
    			propertyContents[i] = new JCheckBox();
        		((JCheckBox) propertyContents[i]).setSize(4,4);
        		((JCheckBox) propertyContents[i]).setSelected(((String) proj.get(propertyKeys[i])).contains("TRUE"));

    		} else if (propertyKeys[i].contains("FILENAME") || propertyKeys[i].contains("DIRECTORY")) {
//    			fileOpenButton = new JButton("Open");
    			propertyContents[i] = new JTextField(50);
        		((JTextField) propertyContents[i]).setText((String) proj.get(propertyKeys[i]));

    		} else {
    			propertyContents[i] = new JTextField(50);
        		((JTextField) propertyContents[i]).setText((String) proj.get(propertyKeys[i]));
    		}
    	}
	}

//	private JPanel rightPanel(Project proj) {
	private JScrollPane rightPanel() {
		JPanel result;
		JScrollPane scrollPane;
    	SpringLayout contentLayout;
    	JLabel propertyName;

    	result = new JPanel();
    	scrollPane = new JScrollPane(result);
    	contentLayout = new SpringLayout();
    	result.setLayout(contentLayout);

    	for (int i=0; i<propertyKeys.length; i++) {
//    		result.add(propertyNames[i]);
//    		propertyNames[i].setLabelFor(propertyContents[i]);
    		propertyName = new JLabel(propertyKeys[i], JLabel.TRAILING);
    		result.add(propertyName);
    		propertyName.setLabelFor(propertyContents[i]);
    		result.add(propertyContents[i]);
//    		contentLayout.putConstraint(SpringLayout.EAST, propertyNames[i], 250, SpringLayout.WEST, result.getRootPane());
//    		contentLayout.putConstraint(SpringLayout.NORTH, propertyNames[i], i*20, SpringLayout.NORTH, result.getRootPane());
    		contentLayout.putConstraint(SpringLayout.EAST, propertyName, 250, SpringLayout.WEST, result.getRootPane());
    		contentLayout.putConstraint(SpringLayout.NORTH, propertyName, i*20, SpringLayout.NORTH, result.getRootPane());
    		contentLayout.putConstraint(SpringLayout.WEST, propertyContents[i], 260, SpringLayout.WEST, result.getRootPane());
    		contentLayout.putConstraint(SpringLayout.NORTH, propertyContents[i], i*20, SpringLayout.NORTH, result.getRootPane());
    	}

    	return scrollPane;
	}


//	private JPanel rightPanel(Project proj) {
////	private JScrollPane rightPanel(Project proj) {
//		JPanel result;
////		JScrollPane scroll;
//    	JLabel propertyName;
//    	SpringLayout contentLayout;
//    	JComponent propertyContent;
//
//    	result = new JPanel();
////    	scroll = new JScrollPane(result);
//    	contentLayout = new SpringLayout();
//    	result.setLayout(contentLayout);
//    	String[] propertyKeys = HashVec.getKeys(proj);
//    	for (int i=0; i<propertyKeys.length; i++) {
//    		propertyName = new JLabel(propertyKeys[i], JLabel.TRAILING);
//    		result.add(propertyName);
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
//    		result.add(propertyContent);
//    		contentLayout.putConstraint(SpringLayout.EAST, propertyName, 250, SpringLayout.WEST, result.getRootPane());
//    		contentLayout.putConstraint(SpringLayout.NORTH, propertyName, i*20, SpringLayout.NORTH, result.getRootPane());
//    		contentLayout.putConstraint(SpringLayout.WEST, propertyContent, 260, SpringLayout.WEST, result.getRootPane());
//    		contentLayout.putConstraint(SpringLayout.NORTH, propertyContent, i*20, SpringLayout.NORTH, result.getRootPane());
//    	}
//		return result;
////    	return scroll;
//	}





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
		String command = e.getActionCommand();

		if (command.equals(SAVE)) {
			PrintWriter out;
			String propertiesFilePath = proj.getPropertyFilename();
			
			new File(propertiesFilePath).renameTo(new File(propertiesFilePath + ".bak"));
			try {
				out = new PrintWriter(new FileOutputStream(propertiesFilePath));
				for (int i = 0; i < propertyKeys.length; i++) {
//					if (proj.get(propertyContents[i]).equals(propertyContents[i])) {
//						;
//					}
					out.println(propertyKeys[i] + "\t" + propertyContents[i]);
				}
			} catch (FileNotFoundException e1) {
				e1.printStackTrace();
			}
		} else if (command.equals(ADD)) {
			JTextField propertyName;
			
			propertyName = new JTextField();
			propertyName.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					JTextField propertyValue;

					propertyValue = new JTextField();
					editor.add(propertyValue);

					//TODO
					propertyKeys = Array.merge(propertyKeys, new String[] {e.getActionCommand()}, 0);
				}
				
			});
			editor.add(propertyName);
		}
	}

	@Override
	public void windowActivated(WindowEvent arg0) {
	}
	
	@Override
	public void windowClosed(WindowEvent arg0) {
	}
	
	@Override
	public void windowClosing(WindowEvent arg0) {
//		saveButton.firePropertyChange(propertyName, oldValue, newValue);
	}
	
	@Override
	public void windowDeactivated(WindowEvent arg0) {
	}
	
	@Override
	public void windowDeiconified(WindowEvent arg0) {
	}
	
	@Override
	public void windowIconified(WindowEvent arg0) {
	}
	
	@Override
	public void windowOpened(WindowEvent arg0) {
	}

}
