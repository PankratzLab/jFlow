package org.genvisis.bot;

import javax.swing.*;

import org.genvisis.common.*;

import java.awt.*;
import java.awt.event.*;

public class HelperCurrent extends JFrame {
	private static final long serialVersionUID = 1L;
	
	public static final String NATURE_GENETICS = "Author filler";
	public static final String IGV_TRIOS = "IGV Trios";

	public static final String JFRAME_TITLE = "Helper";
	public static final String[] DIRS = {"C:/Users/"};
	public static final String[] OPTIONS = {IGV_TRIOS, NATURE_GENETICS};
	public static final String LAUNCH_OPTION = "Launch option";
	public static final String BUTTON1 = "Next Author";
	public static final String BUTTON2 = "Next IGV";
	public static final String BUTTON3 = "Button 3";
	
	private JComboBox<String> optionsBox;
	private JTextArea text;
	
	private AhkBot bot;
	private Logger log;
	private boolean done;
	
	private NatureGenetics ng;
	private IGVTrios igvTrios;
	
	public HelperCurrent() {
		log = new Logger();

		try {
			bot = new AhkBot();
		} catch (AWTException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
//		ng = new NatureGenetics(bot, log);
		igvTrios = new IGVTrios(bot, log);
	}

	public boolean done() {
		return done;
	}
	
	public static String getTime() {
		return ext.replaceAllWith(ext.getDate()+"_"+ext.getTime(), ":", ".");
	}
	
    public void addComponentsToPane(final Container pane) {
        JPanel panel;
		GridLayout layout;
		JButton button;
		JScrollPane scrollPane;

        panel = new JPanel();
        panel.setBackground(Color.WHITE);
        text = new JTextArea();
        text.setFont(new Font("Arial", 0, 12));
        text.setLineWrap(true);
        log.linkTextArea(text);

        scrollPane = new JScrollPane(text); 
        scrollPane.setPreferredSize(new Dimension(350, 160));
        panel.add(scrollPane);
        text.setText(ext.getTime()+"\tOi!\r\n");

        pane.add(panel);

		optionsBox = new JComboBox<String>();
		optionsBox.setModel(new DefaultComboBoxModel<String>(OPTIONS));
		optionsBox.setSelectedIndex(0);
				
        panel = new JPanel();
        layout = new GridLayout(3, 2);
        panel.setLayout(layout);
        layout.setHgap(10);
        layout.setVgap(10);
        panel.setBackground(Color.WHITE);
        
        panel.add(button = new JButton(BUTTON1));
        button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (ng != null) {
					ng.next();
				}
			}
		});
        
        panel.add(button = new JButton(BUTTON2));
        button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				igvTrios.next();
			}
		});

        
        panel.add(optionsBox);
        
        panel.add(button = new JButton(LAUNCH_OPTION));
        button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (OPTIONS[optionsBox.getSelectedIndex()].equals(NATURE_GENETICS)) {
					ng = new NatureGenetics(bot, log);
				} else if (OPTIONS[optionsBox.getSelectedIndex()].equals(IGV_TRIOS)) {
					igvTrios = new IGVTrios(bot, log);
				}

			}
		});

        panel.add(button = new JButton(BUTTON3));
        button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

			}
		});
        
        panel.add(button = new JButton("TEST"));
        button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

			}
		});
        
        pane.add(panel, BorderLayout.SOUTH);
    }        

    private static void createAndShowGUI() {
    	JFrame frame;
    	HelperCurrent nature;
    	
    	nature = new HelperCurrent();

    	frame = new JFrame(JFRAME_TITLE);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        nature.addComponentsToPane(frame.getContentPane());
        frame.setAlwaysOnTop(true);
        frame.setLocationByPlatform(true);
        frame.pack();
//		frame.setLocation(800,640);
		frame.setLocation((int)Toolkit.getDefaultToolkit().getScreenSize().getWidth()-380,640);
        frame.setVisible(true);
    }

	public static void main(String[] args) {
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
//        		new Thread(new bot.AppKiller(determineDir(DIRS)+"plug")).start();
            	createAndShowGUI();
            }
        });
	}
}