package one.ben.fcs;

import javax.swing.JPanel;
import net.miginfocom.swing.MigLayout;
import javax.swing.JTextField;
import javax.swing.JLabel;
import javax.swing.JButton;
import javax.swing.JTable;
import java.awt.GridLayout;

public class DataControlPanel extends JPanel {
    private JTextField textField;
    private JTable table;
    private JPanel panel;
    private JButton btnLoadFile;
    private JButton btnMoveUp;
    private JButton btnMoveDown;
    private JButton btnRemove;

    /**
     * Create the panel.
     */
    public DataControlPanel() {
        setLayout(new MigLayout("", "[][grow][]", "[][grow][][]"));
        
        JLabel lblFileDirectory = new JLabel("File Directory:");
        add(lblFileDirectory, "cell 0 0,alignx trailing");
        
        textField = new JTextField();
        add(textField, "cell 1 0,growx");
        textField.setColumns(10);
        
        JButton button = new JButton(">");
        add(button, "cell 2 0");
        
        table = new JTable();
        add(table, "cell 0 1 3 1,grow");
        
        panel = new JPanel();
        add(panel, "cell 0 2 3 1,grow");
        panel.setLayout(new MigLayout("insets 0", "[80px:80px,grow][80px:80px,grow]", "[][][]"));
        
        btnLoadFile = new JButton("Load");
        panel.add(btnLoadFile, "cell 0 0 2 1,alignx center");
        
        btnMoveUp = new JButton("Move Up");
        panel.add(btnMoveUp, "cell 0 1,growx");
        
        btnMoveDown = new JButton("Move Down");
        panel.add(btnMoveDown, "cell 1 1,growx");
        
        btnRemove = new JButton("Remove");
        panel.add(btnRemove, "cell 0 2 2 1,alignx center");

    }

}
