package org.genvisis.cnv.workflow;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.GenvisisWorkflowGUI;
import org.genvisis.common.gui.Task;

public class StepTask extends Task<Void, Void> {
	GenvisisWorkflowGUI gui;
	protected Throwable failureException;
	private ArrayList<String> failReasons = new ArrayList<String>();
	protected Step.FINAL_CODE returnCode = Step.FINAL_CODE.CANCELLED;
	private boolean failed = false;
	private Project proj;
	private Step step;
	private Set<Step> selectedSteps;
	private Map<Step, Map<Requirement, String>> variables;
	private Thread bgThread;

	public StepTask(GenvisisWorkflowGUI gui, Step step, Project proj, Set<Step> selectedSteps,
									Map<Step, Map<Requirement, String>> variables) {
		this(gui, step, proj, selectedSteps, variables, 0);
	}

	public StepTask(GenvisisWorkflowGUI gui, Step step, Project proj,
									Set<Step> selectedSteps, Map<Step, Map<Requirement, String>> variables,
									int numUpdates) {
		super(step.getName(), numUpdates);
		this.proj = proj;
		this.gui = gui;
		this.step = step;
		this.selectedSteps = selectedSteps;
		this.variables = variables;
	}

	@Override
	protected Void doInBackground() throws Exception {
		this.bgThread = Thread.currentThread();
		Exception e = null;
		Step.FINAL_CODE code = Step.FINAL_CODE.COMPLETE;
		try {
			this.step.setNecessaryPreRunProperties(proj, variables);
			this.step.run(proj, variables);
		} catch (RuntimeException e1) {
			if (e1.getCause() instanceof InterruptedException) {
				code = Step.FINAL_CODE.CANCELLED;
			}
			e = e1;
		} catch (Exception e1) {
			e = e1;
		}
		if (code != Step.FINAL_CODE.CANCELLED
				&& (e != null || getFailed() || !this.step.checkIfOutputExists(variables))) {
			code = Step.FINAL_CODE.FAILED;
		}
		failureException = e;
		returnCode = code;
		return null;
	}

	@SuppressWarnings("deprecation")
	// uses Thread.stop();
	@Override
	protected void done() {
		if (this.isCancelled()) {
			bgThread.stop();
			returnCode = Step.FINAL_CODE.CANCELLED;
		}
		gui.nextStep(this, returnCode, selectedSteps, variables);
	}

	public void setFailed(String reason) {
		failed = true;
		failReasons.add(reason);
	}

	public List<String> getFailureMessages() {
		return failReasons;
	}

	public Throwable getFailureException() {
		return failureException;
	}

	public Step getStep() {
		return this.step;
	}

	public boolean getFailed() {
		return failed;
	}

}
