package org.genvisis.cnv.gui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.util.HashMap;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import net.miginfocom.swing.MigLayout;

import javax.swing.JLabel;
import javax.swing.JComboBox;
import javax.swing.JSpinner;
import javax.swing.JSeparator;
import javax.swing.SwingConstants;
import javax.swing.DefaultComboBoxModel;
import javax.swing.SpinnerNumberModel;

import org.genvisis.cnv.manage.GenvisisWorkflow.FLAG;
import org.genvisis.common.QueueControl;
import org.genvisis.common.QueueControl.JobQueue;

import java.awt.event.ItemListener;
import java.awt.event.ItemEvent;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class QueuePicker extends JDialog {

	private final JPanel contentPanel = new JPanel();

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		try {
			QueuePicker dialog = new QueuePicker();
			dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
			dialog.setVisible(true);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Create the dialog.
	 */
	public QueuePicker() {
		setTitle("Queue Details");
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(new MigLayout("", "[][grow]", "[][][][]"));
		{
			JLabel lblPbsQueue = new JLabel("PBS Queue:");
			contentPanel.add(lblPbsQueue, "cell 0 0,alignx trailing");
		}
		{
			JComboBox comboQueue = new JComboBox();
			comboQueue.addItemListener(new ItemListener() {
				public void itemStateChanged(ItemEvent ie) {
					if (ie.getStateChange() == ItemEvent.SELECTED) {
						queueSelected((String) ie.getItem());
					}
				}
			});
			comboQueue.setEditable(true);
			contentPanel.add(comboQueue, "cell 1 0,growx");
		}
		{
			JLabel lblProcessors = new JLabel("Processors:");
			contentPanel.add(lblProcessors, "cell 0 1");
		}
		{
			JSpinner spinProc = new JSpinner();
			spinProc.setModel(new SpinnerNumberModel(1, 1, 999, 1));
			spinProc.setEditor(new JSpinner.NumberEditor(spinProc, "00"));
			contentPanel.add(spinProc, "flowx,cell 1 1,alignx left");
		}
		{
			JLabel lblMemory = new JLabel("Memory:");
			contentPanel.add(lblMemory, "cell 0 2");
		}
		{
			JSpinner spinMem = new JSpinner();
			spinMem.setModel(new SpinnerNumberModel(1, 1, 999, 1));
			spinMem.setEditor(new JSpinner.NumberEditor(spinMem, "00"));
			contentPanel.add(spinMem, "flowx,cell 1 2");
		}
		{
			JLabel lblWallTime = new JLabel("Wall time:");
			contentPanel.add(lblWallTime, "cell 0 3");
		}
		{
			JSeparator separator = new JSeparator();
			contentPanel.add(separator, "cell 1 1,growx");
		}
		{
			JLabel lblQueueMinmax = new JLabel("Queue min/max:");
			contentPanel.add(lblQueueMinmax, "cell 1 1");
		}
		{
			JComboBox comboMemUnit = new JComboBox();
			comboMemUnit.setModel(new DefaultComboBoxModel(new String[] {"tb", "gb", "mb", "b"}));
			comboMemUnit.setSelectedIndex(1);
			contentPanel.add(comboMemUnit, "cell 1 2");
		}
		{
			JSeparator separator = new JSeparator();
			contentPanel.add(separator, "cell 1 2,growx");
		}
		{
			JLabel label = new JLabel("Queue min/max:");
			contentPanel.add(label, "cell 1 2");
		}
		{
			JSpinner spinWallDay = new JSpinner();
			spinWallDay.setModel(new SpinnerNumberModel(0, 0, 999, 1));
			spinWallDay.setEditor(new JSpinner.NumberEditor(spinWallDay, "00"));
			contentPanel.add(spinWallDay, "flowx,cell 1 3");
		}
		{
			JLabel label = new JLabel(":");
			contentPanel.add(label, "cell 1 3");
		}
		{
			JSpinner spinWallHours = new JSpinner();
			spinWallHours.setModel(new SpinnerNumberModel(0, 0, 99, 1));
			spinWallHours.setEditor(new JSpinner.NumberEditor(spinWallHours, "00"));
			contentPanel.add(spinWallHours, "cell 1 3");
		}
		{
			JLabel label = new JLabel(":");
			contentPanel.add(label, "cell 1 3");
		}
		{
			JLabel lblWallSec = new JLabel("00");
			contentPanel.add(lblWallSec, "cell 1 3");
		}
		{
			JSeparator separator = new JSeparator();
			contentPanel.add(separator, "cell 1 3,growx");
		}
		{
			JLabel label = new JLabel("Queue min/max:");
			contentPanel.add(label, "cell 1 3");
		}
		{
			lblProcMinMax = new JLabel("?? / ??");
			contentPanel.add(lblProcMinMax, "cell 1 1");
		}
		{
			lblMemMinMax = new JLabel("?? / ??");
			contentPanel.add(lblMemMinMax, "cell 1 2");
		}
		{
			lblWalltimeMinMax = new JLabel("?? / ??");
			contentPanel.add(lblWalltimeMinMax, "cell 1 3");
		}
		{
			JPanel buttonPane = new JPanel();
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			buttonPane.setLayout(new MigLayout("ins 5", "[57px][][grow][47px][65px]", "[23px]"));
			{
				JButton btnSave = new JButton("Save");
				btnSave.setEnabled(false);
				buttonPane.add(btnSave, "cell 0 0,alignx left,aligny top");
			}
			{
				JButton btnSetDefault = new JButton("Set Default");
				buttonPane.add(btnSetDefault, "cell 1 0");
			}
			{
				JButton okButton = new JButton("OK");
				okButton.setActionCommand("OK");
				buttonPane.add(okButton, "cell 3 0,alignx left,aligny top");
				getRootPane().setDefaultButton(okButton);
			}
			{
				JButton cancelButton = new JButton("Cancel");
				cancelButton.setActionCommand("Cancel");
				buttonPane.add(cancelButton, "cell 4 0,alignx left,aligny top");
			}
		}
		pack();
	}
	
	protected void customQueue() {
		System.out.println("Creating custom queue");
	}

	protected void queueSelected(String queueName) {
		JobQueue jq = premadeQueues.get(queueName);
		if (jq == null) {
			customQueue();
			return;
		}
		setLimits(jq);
	}
	
	private void setLimits(JobQueue jq) {
		StringBuilder sb = new StringBuilder();
		
		sb.append(jq.getMinProc() > 0 ? jq.getMinProc() : "??")
			.append(" / ")
			.append(jq.getMaxProc() > 0 ? jq.getMaxProc() : "??");
		lblProcMinMax.setText(sb.toString());
		sb = new StringBuilder();
		
		sb.append(jq.getMinMem() > 0 ? jq.getMinProc() : "??")
			.append(" / ")
			.append(jq.getMaxProc() > 0 ? jq.getMaxProc() : "??");
		lblMemMinMax.setText(sb.toString());
		sb = new StringBuilder();
		
		sb.append(jq.getMinWalltime() > 0 ? jq.getMinWalltime() + "hr(s)" : "??")
			.append(" / ")
			.append(jq.getMinWalltime() > 0 ? jq.getMinWalltime() + "hrs(s)" : "??");
		lblWalltimeMinMax.setText(sb.toString());
	}
	
	HashMap<String, JobQueue> premadeQueues = new HashMap<String, QueueControl.JobQueue>();
	private JLabel lblProcMinMax;
	private JLabel lblMemMinMax;
	private JLabel lblWalltimeMinMax;

	protected void populate(List<JobQueue> queues) {
		for (JobQueue jq : queues) {
			premadeQueues.put(jq.getName(), jq);
		}
	}
	
	
	

}
