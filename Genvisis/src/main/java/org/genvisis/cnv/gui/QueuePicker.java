package org.genvisis.cnv.gui;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.HashMap;
import java.util.List;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EmptyBorder;

import org.genvisis.common.Grafik;
import org.genvisis.qsub.JobQueue;
import org.genvisis.qsub.QueueProperties;

import net.miginfocom.swing.MigLayout;

public class QueuePicker extends JDialog {

	private final JPanel contentPanel = new JPanel();

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		try {
			QueuePicker dialog = new QueuePicker("");
			dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
			dialog.setVisible(true);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Create the dialog.
	 */
	public QueuePicker(String qsubFilenameSuggestion) {
		setTitle("QSUB Details");
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(new MigLayout("", "[][grow]", "[][][][][][]"));
		{
			JLabel lblQsubFilename = new JLabel("Qsub Filename:");
			contentPanel.add(lblQsubFilename, "cell 0 0,alignx trailing");
		}
		{
			txtFldFilename = new JTextField(qsubFilenameSuggestion);
			contentPanel.add(txtFldFilename, "cell 1 0,growx");
			txtFldFilename.setColumns(10);
		}
		{
			JSeparator separator = new JSeparator();
			contentPanel.add(separator, "cell 0 1 2 1,growx");
		}
		{
			JLabel lblPbsQueue = new JLabel("PBS Queue:");
			contentPanel.add(lblPbsQueue, "cell 0 2,alignx left");
		}
		{
			comboQueue = new JComboBox();
			comboQueue.addItemListener(new ItemListener() {
				public void itemStateChanged(ItemEvent ie) {
					if (ie.getStateChange() == ItemEvent.SELECTED) {
						queueSelected((String) ie.getItem());
					}
				}
			});
			comboQueue.setEditable(true);
			contentPanel.add(comboQueue, "cell 1 2,growx");
		}
		{
			JLabel lblProcessors = new JLabel("Processors:");
			contentPanel.add(lblProcessors, "cell 0 3");
		}
		{
			spinProc = new JSpinner();
			spinProc.setModel(new SpinnerNumberModel(1, 1, 999, 1));
			spinProc.setEditor(new JSpinner.NumberEditor(spinProc, "00"));
			contentPanel.add(spinProc, "flowx,cell 1 3,alignx left");
		}
		{
			JLabel lblMemory = new JLabel("Memory:");
			contentPanel.add(lblMemory, "cell 0 4");
		}
		{
			spinMem = new JSpinner();
			spinMem.setModel(new SpinnerNumberModel(1, 1, 999, 1));
			spinMem.setEditor(new JSpinner.NumberEditor(spinMem, "00"));
			contentPanel.add(spinMem, "flowx,cell 1 4");
		}
		{
			JLabel lblWallTime = new JLabel("Wall time:");
			contentPanel.add(lblWallTime, "cell 0 5");
		}
		{
			JSeparator separator = new JSeparator();
			contentPanel.add(separator, "cell 1 3,growx");
		}
		{
			JLabel lblQueueMinmax = new JLabel("Queue min/max:");
			contentPanel.add(lblQueueMinmax, "cell 1 3");
		}
		{
			comboMemUnit = new JComboBox();
			comboMemUnit.addItemListener(new ItemListener() {
				public void itemStateChanged(ItemEvent e) {
					updateMemDefault();
				}
			});
			comboMemUnit.setModel(new DefaultComboBoxModel(MEM_UNITS));
			comboMemUnit.setSelectedIndex(1);
			contentPanel.add(comboMemUnit, "cell 1 4");
		}
		{
			JSeparator separator = new JSeparator();
			contentPanel.add(separator, "cell 1 4,growx");
		}
		{
			JLabel label = new JLabel("Queue min/max:");
			contentPanel.add(label, "cell 1 4");
		}
		{
			spinWallDay = new JSpinner();
			spinWallDay.setModel(new SpinnerNumberModel(0, 0, 999, 1));
			spinWallDay.setEditor(new JSpinner.NumberEditor(spinWallDay, "00"));
			contentPanel.add(spinWallDay, "flowx,cell 1 5");
		}
		{
			JLabel label = new JLabel(":");
			contentPanel.add(label, "cell 1 5");
		}
		{
			spinWallHours = new JSpinner();
			spinWallHours.setModel(new SpinnerNumberModel(0, 0, 24, 1));
			spinWallHours.setEditor(new JSpinner.NumberEditor(spinWallHours, "00"));
			contentPanel.add(spinWallHours, "cell 1 5");
		}
		{
			JLabel label = new JLabel(":");
			contentPanel.add(label, "cell 1 5");
		}
		{
			JLabel lblWallSec = new JLabel("00");
			contentPanel.add(lblWallSec, "cell 1 5");
		}
		{
			JSeparator separator = new JSeparator();
			contentPanel.add(separator, "cell 1 5,growx");
		}
		{
			JLabel label = new JLabel("Queue min/max:");
			contentPanel.add(label, "cell 1 5");
		}
		{
			lblProcMinMax = new JLabel("?? / ??");
			contentPanel.add(lblProcMinMax, "cell 1 3");
		}
		{
			lblMemMinMax = new JLabel("?? / ??");
			contentPanel.add(lblMemMinMax, "cell 1 4");
		}
		{
			lblWalltimeMinMax = new JLabel("?? / ??");
			contentPanel.add(lblWalltimeMinMax, "cell 1 5");
		}
		{
			JPanel buttonPane = new JPanel();
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			buttonPane.setLayout(new MigLayout("ins 5", "[57px][][grow][47px][65px]", "[23px]"));
			{
				JButton okButton = new JButton("OK");
				okButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent arg0) {
						close(true);
					}
				});
				{
					JLabel label = Grafik.getToolTipIconLabel(TOOLTIP);
					buttonPane.add(label, "cell 0 0");
				}
				{
					chckbxSetAsDefaults = new JCheckBox("Set as queue defaults");
					buttonPane.add(chckbxSetAsDefaults, "cell 2 0,alignx right");
				}
				buttonPane.add(okButton, "cell 3 0,alignx left,aligny top");
				getRootPane().setDefaultButton(okButton);
			}
			{
				JButton cancelButton = new JButton("Cancel");
				cancelButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent arg0) {
						close(false);
					}
				});
				buttonPane.add(cancelButton, "cell 4 0,alignx left,aligny top");
			}
		}
		pack();
	}
	
	private void close(boolean ok) {
		cancelled = !ok;
		if (ok) {
			checkSettings();
		}
		setVisible(false);
	}
	
	private void checkSettings() {
		if (chckbxSetAsDefaults.isSelected()) {
			String qName = (String) comboQueue.getSelectedItem();
			JobQueue jq;
			if (premadeQueues.containsKey(qName)) {
  			jq = premadeQueues.get(qName);
  			if (jq == null) {
  				return;
  			}
			} else {
				jq = QueueProperties.createNewQueue(qName, true); 
			}
			jq.setDefaultMem(getMemoryInB());
			jq.setDefaultWalltime(getWalltimeHours());
			jq.setDefaultProcCnt(getProcessors());
			
			QueueProperties.setDefaultQueueName(jq.getName());
			QueueProperties.save(QueueProperties.PROPERTIES_FILE);
		}
	}
	
	protected void customQueue() {
		//
	}

	private void updateMemDefault() {
		if (settingLimits) return;
		JobQueue jq = premadeQueues.get(comboQueue.getSelectedItem());
		if (jq == null) return;
		
		long m;
		m = jq.getDefaultMem();
		
		long mMin = jq.getMinMem();
		if (mMin < 0) {
			mMin = 0;
		}
		long mMax = jq.getMaxMem();
		if (mMin < 0) {
			mMax = 99;
		}
		
		int trans = 1;
		String memUnit = (String) comboMemUnit.getSelectedItem();
		if (memUnit.equals(TB)) {
			trans = 1024 * 1024 * 1024;
		} else if (memUnit.equals(GB)) {
			trans = 1024 * 1024;
		} else if (memUnit.equals(MB)) {
			trans = 1024;
		} else if (memUnit.equals(B)) {
			// do nothing
		}
		spinMem.setModel(new SpinnerNumberModel((int) m / trans, (int) mMin / trans, (int) mMax / trans, 1));
		
		StringBuilder sb = new StringBuilder();
		sb.append(jq.getMinMem() >= 0 ? memTransform(jq.getMinMem()) : "??")
			.append(" / ")
			.append(jq.getMaxMem() >= 0 ? memTransform(jq.getMaxMem()) : "??");
		lblMemMinMax.setText(sb.toString());
	}
	
	protected void queueSelected(String queueName) {
		JobQueue jq = premadeQueues.get(queueName);
		if (jq == null) {
			customQueue();
			return;
		}
		setLimits(jq);
	}
	
	private long memTransform(long m) {
		String memUnit = (String) comboMemUnit.getSelectedItem();
		long mem = m;
		if (memUnit.equals(TB)) {
			mem = mem / 1024 / 1024 / 1024;
		} else if (memUnit.equals(GB)) {
			mem = mem / 1024 / 1024;
		} else if (memUnit.equals(MB)) {
			mem = mem / 1024;
		} else if (memUnit.equals(B)) {
			// do nothing
		}
		return mem;
	}
	
	volatile boolean settingLimits = false;
	private void setLimits(JobQueue jq) {
		settingLimits = true;
		StringBuilder sb = new StringBuilder();
		
		sb.append(jq.getMinProc() >= 0 ? jq.getMinProc() : "??")
			.append(" / ")
			.append(jq.getMaxProc() >= 0 ? jq.getMaxProc() : "??");
		lblProcMinMax.setText(sb.toString());
		sb = new StringBuilder();
		
		sb.append(jq.getMinMem() >= 0 ? memTransform(jq.getMinMem()) : "??")
			.append(" / ")
			.append(jq.getMaxMem() >= 0 ? memTransform(jq.getMaxMem()) : "??");
		lblMemMinMax.setText(sb.toString());
		sb = new StringBuilder();
		
		sb.append(jq.getMinWalltime() >= 0 ? jq.getMinWalltime() + "hr(s)" : "??")
			.append(" / ")
			.append(jq.getMaxWalltime() >= 0 ? jq.getMaxWalltime() + "hrs(s)" : "??");
		lblWalltimeMinMax.setText(sb.toString());
		
		long m;
		int p, wH = 0, wD = 0;
		m = jq.getDefaultMem();
		p = jq.getDefaultProc();
		wH = jq.getDefaultWalltime();
		if (wH > 24) {
			wD = wH / 24;
			wH = wH % 24;
		}
		
		int pMin = jq.getMinProc();
		if (pMin < 0) pMin = 0;
		int pMax = jq.getMaxProc();
		if (pMax < 0) {
			pMax = 99;
		}
		
		spinProc.setModel(new SpinnerNumberModel(p, pMin, pMax, 1));
		
		if (m < 0) {
			m = 0;
		} else {
			m = m / 1024;
		}
		long mMin = jq.getMinMem();
		if (mMin < 0) {
			mMin = 0;
		} else {
			mMin = mMin / 1024;
		}
		long mMax = jq.getMaxMem();
		if (mMin < 0) {
			mMax = 99;
		} else {
			mMax = mMax / 1024;
		}
		
		spinMem.setModel(new SpinnerNumberModel((int) m, (int) mMin, (int) mMax, 1));
		
		comboMemUnit.setSelectedItem(MB);
		
		spinProc.setValue(p);
		
		spinWallDay.setValue(wD);
		spinWallHours.setValue(wH);
		settingLimits = false;
	}
	
	private static final String TOOLTIP = "<html>"
																					+"<p>Enter desired walltime, memory, and processors.</p>"
																					+"<br/>"
																					+"<p>Queue name is not required unless you wish to save these values as defaults.</p>"
																					+"<br/>"
																					+"<p>Walltime is <code>days::hours::minutes</code></p>"
																					+"<br/>"
																					+"<p>\"Set as queue defaults\" will not affect qsub filename.</p>"
																					+"</html>";
	
	HashMap<String, JobQueue> premadeQueues = new HashMap<String, JobQueue>();
	private JLabel lblProcMinMax;
	private JLabel lblMemMinMax;
	private JLabel lblWalltimeMinMax;
	
	private final String TB = "tb";
	private final String GB = "gb";
	private final String MB = "mb";
	private final String B = "b";
	
	private final String[] MEM_UNITS = new String[] {TB, GB, MB, B};
	
	private boolean cancelled = false;
	private JSpinner spinProc;
	private JSpinner spinMem;
	private JComboBox comboMemUnit;
	private JSpinner spinWallDay;
	private JSpinner spinWallHours;
	private JComboBox comboQueue;

	private JCheckBox chckbxSetAsDefaults;
	private JTextField txtFldFilename;
	
	public boolean wasCancelled() {
		return cancelled;
	}
	
	public int getProcessors() {
		int proc = ((Number) spinProc.getValue()).intValue();
		return proc;
	}
	
	public int getWalltimeHours() {
		int days = ((Number) spinWallDay.getValue()).intValue();
		int hrs = ((Number) spinWallHours.getValue()).intValue();
		hrs += (days * 24);
		return hrs;
	}
	
	public long getMemoryInB() {
		String memUnit = (String) comboMemUnit.getSelectedItem();
		long mem = ((Number) spinMem.getValue()).longValue();
		if (memUnit.equals(TB)) {
			mem = mem * 1024 * 1024 * 1024;
		} else if (memUnit.equals(GB)) {
			mem = mem * 1024 * 1024;
		} else if (memUnit.equals(MB)) {
			mem = mem * 1024;
		} else if (memUnit.equals(B)) {
			// do nothing
		}
		return mem;
	}
	
	public int getMemoryInMb() {
		return (int) (getMemoryInB() / 1024);
	}

	public String getFilename() {
		return txtFldFilename.getText();
	}
	
	public void populate(List<JobQueue> queues) {
		String[] names = new String[queues.size() + 1];
		for (int i = 0; i < queues.size(); i++) {
			JobQueue jq = queues.get(i);
			names[i] = jq.getName();
			premadeQueues.put(jq.getName(), jq);
		}
		names[queues.size()] = "";
		comboQueue.setModel(new DefaultComboBoxModel(names));
		comboQueue.setSelectedItem(QueueProperties.getDefaultQueueName());
	}
	
	
}
