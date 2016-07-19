package org.genvisis.dead;
/**
to do list:

change file format to record when the last time a word was seen (i.e. 12415:XOOOO-2 finalized 2 days after date || alternatively, change date to today when finalized)
customize how far back words will come (i.e. don't show words I've seen today, yesterday, t - 2, etc.)

refresh one after changing anything

get rid of about panel?

come up with a master way of dealing with wrongHistory

make a global enablePane, disablePane to be used with word adder and stats

resizable window
ran - change innards to match
 - automatically keep lower bound when let go

add words internally
 - automatically fill in last known group heading

pick which group headings to include in the hash
 - might want to temporarily disable update when doing click all, click none
 - next to each heading might want to put error rate or % until "mastery"

 type it baby!
 - cut and paste

language
 - English/Turkish/Mix of the 2

automatic scrolling
 - with and without answers
 - sliding bar for speed
 - space for stop

for multiple definitions, pick one, delimit with &&

reaction time - speed baby! // this idea is so gay, I won't even try

*/

import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

//import borland.jbcl.layout.*;
import javax.swing.*;
import javax.swing.event.*;
//import javax.swing.table.*;

import org.genvisis.dead.Incipient.Word.TimePoint;


public class Incipient implements ActionListener, KeyListener, MenuListener {
    public class Word {
		public class TimePoint {
			private int when;
			private String performance;

			public TimePoint(int date, String anxiety) {
				when = date;
				performance = anxiety;
			}
			public boolean sameDay(int asToday) {
				return (when == asToday);
			}
			public boolean caughtUp() {
				return caughtUp(0);
			}
			public boolean caughtUp(int beginAt) {
				int reset = (currentQ2A?RESET_Q_SCORE:RESET_A_SCORE);
				int score = beginAt;

				for (int i=0; i<performance.length(); i++) {
					if (performance.charAt(i) == 'O') {
						score -= CORCT_INC;
					} else {
						if (score < reset-WRONG_INC) {
							score = reset;
						} else {
							score += WRONG_INC;
						}
						if (score < 1) {
							score = 1;
						}
						if (score > MAXSCORE) {
							score = MAXSCORE;
						}
					}
				}

				return (score <= 1);
			}
			
			public int countWrong() {
				int count;
				
				count = 0;
				for (int i = 0; i<performance.length(); i++) {
					if (performance.charAt(i) == 'X') {
						count++;
					}
                }
				
				return count;
			}
			
			public String toString() {
				return (when +":"+ performance);
			}
		}

		private String AnswerString;
		private String QuestionString;
		private Vector<TimePoint> qHistory;
		private Vector<TimePoint> aHistory;
		private String category;
		private int q_score;
		private int a_score;
		private int hashNum;

		public Word(String answer, String question, String whichCategory, int incdHashNum) {
			QuestionString = question;
			AnswerString = answer;
			category = whichCategory;
			hashNum = incdHashNum;
			qHistory = new Vector<TimePoint>();
			aHistory = new Vector<TimePoint>();
			q_score = RESET_Q_SCORE;
			a_score = RESET_A_SCORE;
		}

		public void addTimePoint(int when, boolean wasCorrect, boolean question2answer) {
			TimePoint t = null;
			Vector<TimePoint> history;

			if (question2answer) {
				history = qHistory;
			} else {
				history = aHistory;
			}

			if (history.size() > 0) {
				t = (TimePoint)history.elementAt(history.size()-1);
			}
			if (t != null && (!t.caughtUp() || t.sameDay(today))) {
				t.performance = t.performance + (wasCorrect?"O":"X");
			} else {
				history.addElement(new TimePoint(today, wasCorrect?"O":"X"));
			}
		}

		public void addTimeSeries(int when, String evals, boolean question2answer) {
			if (question2answer) {
				qHistory.addElement(new TimePoint(when, evals));
			} else {
				aHistory.addElement(new TimePoint(when, evals));
			}
		}

		public void calcScores() {
			int score, reset;
			TimePoint t;
			int numWrongRecently;
			Vector<TimePoint> history;

			for (int qORa=0; qORa<2; qORa++) {
				if (qORa==0) {
					history = qHistory;
					reset = RESET_Q_SCORE;
				} else {
					history = aHistory;
					reset = RESET_A_SCORE;
				}

				score = reset;
				numWrongRecently = 0;

				while (history.size() > MAXHIST) {
					history.removeElementAt(0);
				}
				for (int i=0; i<history.size(); i++) {
					t = (TimePoint)history.elementAt(i);
					for (int j=0; j<t.performance.length(); j++) {
                        if (t.performance.charAt(j) == 'O') {
							score -= CORCT_INC;
						} else {
							numWrongRecently++;
							if (score < reset-WRONG_INC) {
								score = reset;
							} else {
								score += WRONG_INC;
                            }
							if (score > MAXSCORE) {
								score = MAXSCORE;
							}
						}
						if (score < 1) {
							score = 1;
						}
					}
				}

				if (menuTaggedOption.getState() && score != 1) {
					score = 1;
				} else {
					// Time Modifier for Score (1 point for every RFRSH_INC days since it has been seen)
					if (menuRefreshOption.getState() && history.size() > 0) {
						score += (today - ((TimePoint)history.elementAt(history.size()-1)).when)/RFRSH_INC;
					}

					// Error Modifier for Score (number of times wrong during last MAXHIST time points added to score)
					if (menuHighErrOption.getState() && history.size() > 0 && !(((TimePoint)history.elementAt(history.size()-1)).sameDay(today) || !((TimePoint)history.elementAt(history.size()-1)).caughtUp())) {
						score += numWrongRecently;
					}

					// Latest Additions Modifier for Score (MAXHIST minus the number of times seen)
					if (menuLatestOption.getState() && (history.size() == 0 || !(((TimePoint)history.elementAt(history.size()-1)).sameDay(today) || !((TimePoint)history.elementAt(history.size()-1)).caughtUp()))) {
						score += MAXHIST - history.size();
					}

					// Perfection Modifier for Score (subtracts 2 points if it has a perfect 10 history)
					if (history.size() == 10 && numWrongRecently == 0) {
						score -= 2;
					}
				}

				if (score < 1) {
					score = 1;
				}
				if (score > MAXSCORE) {
					score = MAXSCORE;
				}

				if (qORa==0) {
					q_score = score;
				} else {
					a_score = score;
				}
			}
		}


		private int shouldBeInWrongHistory() {
			int numTimesInWH = 0;
			Vector<TimePoint> history = (currentQ2A?qHistory:aHistory);
			String st;

			if (history.size() == 0) {
				return 0;
			}

			st = ((TimePoint)history.elementAt(history.size()-1)).performance;
			for (int i=0; i<=st.length()-1; i++) {
				if (st.charAt(i) == 'X') {
					numTimesInWH++;
				} else {
					if (numTimesInWH  > 0) {
						numTimesInWH--;
					}
				}
			}

			return numTimesInWH;
		}

		private String briefTitle() {
			String str = (currentQ2A?QuestionString:AnswerString);

			if (str.length() > 50) {
				str = str.substring(0, 47)+"...";
			}

			return str;
		}

		private String lastSeen() {
			Vector<TimePoint> v = (currentQ2A?qHistory:aHistory);
			if (v.size() == 0) {
				return "never seen it";
			}
			return "Last seen "+timeAgo(today-((Word.TimePoint)v.elementAt(v.size()-1)).when);
		}

		private String timeAgo(int ago) {
			if (ago == 0) {
				return "today";
			} else if (ago == 1) {
				return "yesterday";
			} else if (ago >= 365) {
				return ((int)(ago/365) +" years ago");
			} else if (ago >= 61) {
				return ((int)(ago/30.4) +" months ago");
			} else if (ago >= 31) {
				return ("1 month ago");
			} else if (ago >= 14) {
				return ((int)(ago/7) +" weeks ago");
			} else {
				return (ago +" days ago");
			}
		}

		private String dispHist() {
			String str = "";
			TimePoint trav;
			Vector<TimePoint> history;

			if (currentQ2A) {
				history = qHistory;
			} else {
				history = aHistory;
			}

			for (int i=history.size()-1; i>=((history.size()>MAXHIST)?history.size()-MAXHIST:0); i--) {
				trav = (TimePoint)history.elementAt(i);
				str += timeAgo(today - trav.when) +"  ";
				for (int j=0; j<trav.performance.length(); j++) {
					if (trav.performance.charAt(j) == 'X') {
						str += "\u25BC";
					} else {
						str += "\u25B2";
					}
				}
				str += "\n";
			}

			return str;
		}

		public String toString() {
			String stringToReturn = AnswerString +"|"+ QuestionString +"|"+ category;
			for (int i=(qHistory.size()>MAXHIST)?qHistory.size()-MAXHIST:0; i<qHistory.size(); i++) {
				stringToReturn += "|" + (TimePoint)qHistory.elementAt(i);
			}
			stringToReturn += "|<-Q and A->";
			for (int i=(aHistory.size()>MAXHIST)?aHistory.size()-MAXHIST:0; i<aHistory.size(); i++) {
				stringToReturn += "|" + (TimePoint)aHistory.elementAt(i);
			}

			return stringToReturn;
		}

		public String iFliprExport() {
			String stringToReturn = "=\""+AnswerString +"\"\t=\""+ QuestionString +"\"\t"+ category;
			
			stringToReturn += "\t"+qHistory.size()+"\t"+scoreHistory(qHistory);
			stringToReturn += "\t"+aHistory.size()+"\t"+scoreHistory(aHistory);
			stringToReturn += "\t"+((MAXHIST-qHistory.size())+scoreHistory(qHistory)+(5-Math.min(aHistory.size(), 5))+scoreHistory(aHistory));

			return stringToReturn;
		}
		
		public int scoreHistory(Vector<TimePoint> v) {
			int count;
			
			count = 0;
			for (int i = 0; i<v.size(); i++) {
				count += v.elementAt(i).countWrong();
            }
			
			return count;
		}
		
		

	}

	public static final String QUIT = "Quit", QUITWOSAVING = "Quit without updating data", CHANGE_LIST = "Change World List", PROGRESS = "Save Progress", ENTERWORDS = "Enter new words";
	public static final String CUT = "Cut", COPY = "Copy", PASTE = "Paste", HELP = "Help", ABOUT = "About";
	public static final String STATS = "Show statistics", LEFT = "Words left to master", HIDETAGGED = "Hide tagged words", BLENDER = "Blender - shuffle definitions";
	public static final String HIGHERR = "Weight words with high error rates", REFRESH = "Refresh words not seen recently", BRANDNEW = "Learn brand new words", LATEST = "Weight latest words";
	public static final String YES = "Yes", NO = "No";
	public static final String SHOW_ANSWER = "Show", SKIP_WORD = "Skip", RIGHT = "Right", WRONG = "Wrong", UNDO = "Undo last action";
	public static final String DBFILE = "C:\\Documents and Settings\\npankrat\\My Documents\\jProjects\\Flash\\turk.database";

	private int MAXSCORE = 25;
	private int RESET_Q_SCORE = 9;
	private int RESET_A_SCORE = 7;
	private int MAXHIST = 10;
	private int CORCT_INC = 2;
	private int WRONG_INC = 4;
	private int RFRSH_INC = 7;
	private int LABEL_SIZE = 510;

	private JFrame parentFrame;
	private JRootPane pane;
	private Hashtable<String, Word> hash;
	private JLabel instrBar;
	private JLabel catBar;
	private JLabel thresholdLabel;
	private JLabel questions;
	private JLabel answers;
	private JButton showRight;
	private JButton skipWrong;
	private ButtonGroup whichButton;
	private JRadioButton questionButton;
	private JRadioButton answerButton;
	private boolean answerHasBeenShown = false;
	private Word currentWord;
	private boolean currentQ2A;
	private String questionLabel, answerLabel;
	private JDialog aboutBox = null;
	private JWindow wordStatsBox = null;
//	private JDialog enterTable = null;
	private JMenuBar menuBar;
	private JCheckBoxMenuItem menuStatsOption, menuRefreshOption, menuLeftOption, menuHighErrOption, menuBrandNewOption, menuLatestOption, menuTaggedOption, menuBlenderOption;
	private JMenuItem menuUndoItem;
	private int today;
	private int[] qScoreArray, aScoreArray;
	private int qScoreTotal, aScoreTotal;
	private int threshold;
	private Vector<String> wrongHistory;
	private boolean wrongFlag;
	private int[] previousWordAction = {-1,-1,-1,-1,-1}; // {word's hashNum, skipped?, wrong?, in wrongHistory?}
	private int questionsLeft;
//	private int answersLeft;
	private JLabel leftLabel;
//	private String lastCommand;
	private JSlider threshSlider;
	private boolean brandNewQuestions, brandNewAnswers;

	public Incipient(JFrame frame) {
		parentFrame = frame;
		makeMenu(parentFrame);
	}

	public JMenu createMenu(JMenuBar bar, String name, char mnemonic) {
		JMenu menu = (JMenu)bar.add(new JMenu(name));
		menu.setMnemonic(mnemonic);
		menu.addMenuListener(this);
		return menu;
	}

	public JMenuItem createMenuItem(JMenu menu, String name, int mnemonic, String tooltip, int accelerator, boolean checkbox) {
		JMenuItem menuItem = checkbox?(JCheckBoxMenuItem)menu.add(new JCheckBoxMenuItem(name)):(JMenuItem)menu.add(new JMenuItem(name));
		menuItem.setMnemonic(mnemonic);
		menuItem.setToolTipText(tooltip);
		if (accelerator != -1) {
			menuItem.setAccelerator(KeyStroke.getKeyStroke(accelerator, ActionEvent.CTRL_MASK));
		}
		menuItem.addActionListener(this);
		if (checkbox) {
			((JCheckBoxMenuItem)menuItem).setState(false);
		}
		return menuItem;
	}

	public void makeMenu(JFrame frame) {
		JMenu fileMenu, editMenu, optionsMenu, helpMenu;

		menuBar = new JMenuBar();

		fileMenu = createMenu(menuBar, "File", 'F');
		createMenuItem(fileMenu, ENTERWORDS, KeyEvent.VK_N, "Enter new words into the database", KeyEvent.VK_N, false);
		createMenuItem(fileMenu, CHANGE_LIST, KeyEvent.VK_C, "Change list of words to something else... something worth changing to", KeyEvent.VK_C, false);
		fileMenu.addSeparator();
		createMenuItem(fileMenu, PROGRESS, KeyEvent.VK_P, "Update database file (helpful if you might not want to save your progress on the next segment)", KeyEvent.VK_P, false);
		fileMenu.addSeparator();
		createMenuItem(fileMenu, QUITWOSAVING, KeyEvent.VK_W, "Quit without saving all changes made since opening the database", KeyEvent.VK_W, false);
		createMenuItem(fileMenu, QUIT, KeyEvent.VK_Q, "Quit after saving the updated database", KeyEvent.VK_Q, false);

		editMenu = createMenu(menuBar, "Edit", 'E');
		menuUndoItem = (JMenuItem)createMenuItem(editMenu, UNDO, KeyEvent.VK_U, "Backup one word and any undo changes to score", KeyEvent.VK_Z, false);
		createMenuItem(editMenu, SKIP_WORD, KeyEvent.VK_K, "Skip this word without logging any changes", KeyEvent.VK_A, false);
		createMenuItem(editMenu, CUT, KeyEvent.VK_T, "Cut selected text", KeyEvent.VK_X, false);
		createMenuItem(editMenu, COPY, KeyEvent.VK_C, "Copy selected text", KeyEvent.VK_C, false);
		createMenuItem(editMenu, PASTE, KeyEvent.VK_P, "Paste text", KeyEvent.VK_V, false);

		optionsMenu = createMenu(menuBar, "Options", 'O');
		menuStatsOption = (JCheckBoxMenuItem)createMenuItem(optionsMenu, STATS, KeyEvent.VK_S, "Show statistics", KeyEvent.VK_S, true);
		menuLeftOption = (JCheckBoxMenuItem)createMenuItem(optionsMenu, LEFT, KeyEvent.VK_L, "The number of remaining words above the threshold", KeyEvent.VK_3, true);
		optionsMenu.addSeparator();
		menuRefreshOption = (JCheckBoxMenuItem)createMenuItem(optionsMenu, REFRESH, KeyEvent.VK_R, "Add 1 point per week not seen", KeyEvent.VK_R, true);
		menuHighErrOption = (JCheckBoxMenuItem)createMenuItem(optionsMenu, HIGHERR, KeyEvent.VK_E, "Add 1 point per error made beyond in the last "+MAXHIST+" tries", KeyEvent.VK_E, true);
		menuLatestOption = (JCheckBoxMenuItem)createMenuItem(optionsMenu, LATEST, KeyEvent.VK_L, "Adds "+MAXHIST+" points minus the number of times the word has been seen", KeyEvent.VK_L, true);
		menuBrandNewOption = (JCheckBoxMenuItem)createMenuItem(optionsMenu, BRANDNEW, KeyEvent.VK_B, "Automatically changes threshold; ensures all words are below the initial threshold before proceeding", KeyEvent.VK_B, true);
		menuTaggedOption = (JCheckBoxMenuItem)createMenuItem(optionsMenu, HIDETAGGED, KeyEvent.VK_T, "Removes tagged (wrong) entries from pool while "+REFRESH+", "+HIGHERR+", or "+LATEST+" is checked", KeyEvent.VK_T, true);
		menuBlenderOption = (JCheckBoxMenuItem)createMenuItem(optionsMenu, BLENDER, KeyEvent.VK_D, "Shuffles definitions", KeyEvent.VK_D, true);

		helpMenu = createMenu(menuBar, "Help", 'H');
		createMenuItem(helpMenu, HELP, KeyEvent.VK_H, "Display help screen", -1, true);
		helpMenu.addSeparator();
		createMenuItem(helpMenu, ABOUT, KeyEvent.VK_A, "Display information about this program", -1, true);

		frame.setJMenuBar(menuBar);
	}

	public void menuCanceled(MenuEvent me) {}
	public void menuSelected(MenuEvent me) { pane.setEnabled(false); }
	public void menuDeselected(MenuEvent me) { pane.setEnabled(true); }

	public void keyPressed(KeyEvent ke) {System.err.println("down");}
	public void keyReleased(KeyEvent ke) {System.err.println("up");}
	public void keyTyped(KeyEvent ke) {
		char keyChar = ke.getKeyChar();

		if (pane.isEnabled()) {
			if ((keyChar == 'n' || keyChar == 'N' || keyChar == 'w') && answerHasBeenShown) {
				actionPerformed(new ActionEvent(parentFrame, 1, WRONG));
			}
		}

		if (menuStatsOption.getState()) {
			if (keyChar == KeyEvent.VK_ESCAPE) {
				menuStatsOption.setState(false);
				actionPerformed(new ActionEvent(parentFrame, 1, STATS));
			}
		}

	}

	public void actionPerformed(ActionEvent ae) {
		menuBar.setSelected(null);
		String command;
		if (ae.getActionCommand() == null) {
			command = ((MenuItem)ae.getSource()).getActionCommand();
		} else {
			command = ae.getActionCommand();
		}

		if (wordStatsBox != null) {
			menuStatsOption.setState(false);
			command = STATS;
		}

		if (command.equals(SHOW_ANSWER)) {
			showAnswer();

			showRight.setText(RIGHT);
			showRight.setActionCommand(RIGHT);
			skipWrong.setText(WRONG);
			skipWrong.setActionCommand(WRONG);
		} else if (command.equals(SKIP_WORD)) {
			pickNewWord(currentWord.hashNum, 1, -1, -1, (wrongFlag?1:0));
			showWord();

			showRight.setText(SHOW_ANSWER);
			showRight.setActionCommand(SHOW_ANSWER);
			skipWrong.setText(SKIP_WORD);
			skipWrong.setActionCommand(SKIP_WORD);
		} else if (command.equals(RIGHT) || command.equals(WRONG)) {
			boolean wasRighto = checkAnswer(command.equals(RIGHT));
			pickNewWord(currentWord.hashNum, 0, (command.equals(RIGHT)?0:1), (wasRighto?1:0), (wrongFlag?1:0));
			showWord();

			showRight.setText(SHOW_ANSWER);
			showRight.setActionCommand(SHOW_ANSWER);
			skipWrong.setText(SKIP_WORD);
			skipWrong.setActionCommand(SKIP_WORD);
		} else if (command.equals(UNDO) && previousWordAction[0] != -1 && pane.isEnabled()) {
			undoLastAction();

			showRight.setText(SHOW_ANSWER);
			showRight.setActionCommand(SHOW_ANSWER);
			skipWrong.setText(SKIP_WORD);
			skipWrong.setActionCommand(SKIP_WORD);
		} else if (command.equals(HELP)) {
			System.err.println("Help! I need somebody. Help! not just anybody...");
		} else if (command.equals(LEFT)) {
			if (menuLeftOption.getState()) {
				leftLabel.setVisible(true);
			} else {
				leftLabel.setVisible(false);
			}
		} else if (command.equals(STATS)) {
			if (menuStatsOption.getState()) {
				if (wordStatsBox == null) {
					IPanel wsPanel = new IPanel(400, 230);
					wsPanel.setLayout(new XYLayout());
					wsPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createRaisedBevelBorder(), BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Statistics for \'"+currentWord.briefTitle()+"\'")));

					wordStatsBox = new JWindow(parentFrame);
					wordStatsBox.getContentPane().add(wsPanel, BorderLayout.CENTER);

					JLabel njl = new JLabel("Score = "+(currentQ2A?currentWord.q_score:currentWord.a_score) +", "+ currentWord.lastSeen(), JLabel.CENTER);
					njl.setFont(new Font("Arial", Font.PLAIN, 20));
					wsPanel.add(njl, new XYConstraints(0, 0, 400, 30));

					JTextArea njta = new JTextArea(currentWord.dispHist());
					njta.setFont(new Font("Arial", Font.PLAIN, 12));
					njta.setBackground(answers.getBackground());
					njta.setEditable(false);
					njta.setEnabled(false);
					njta.setDisabledTextColor(Color.black);
					wsPanel.add(njta, new XYConstraints(50, 35, 300, 160));
				}
				wordStatsBox.pack();
				Point p = parentFrame.getLocationOnScreen();
				wordStatsBox.setLocation(p.x + 103, p.y + 115);
				wordStatsBox.setVisible(true);
			} else {
				wordStatsBox.dispose();
				wordStatsBox = null;
			}
		} else if (command.equals(REFRESH)) {
			calcScoreSums();
			if (!answerHasBeenShown) {
				showWord();
			} else {
				showWord();
				showAnswer();
			}
		} else if (command.equals(LATEST)) {
			calcScoreSums();
			if (!answerHasBeenShown) {
				showWord();
			} else {
				showWord();
				showAnswer();
			}
		} else if (command.equals(HIGHERR)) {
			calcScoreSums();
			if (!answerHasBeenShown) {
				showWord();
			} else {
				showWord();
				showAnswer();
			}
		} else if (command.equals(BRANDNEW)) {
			calcScoreSums();
			if (menuBrandNewOption.getState() && (brandNewQuestions || brandNewAnswers)) {
				menuRefreshOption.setState(false);
				menuHighErrOption.setState(false);
				menuLatestOption.setState(false);
				if (brandNewQuestions) {
					threshSlider.setValue(RESET_Q_SCORE);
					currentQ2A = true;
				} else if (brandNewAnswers) {
					threshSlider.setValue(RESET_A_SCORE);
				}
			}
			pickNewWord(-1, -1, -1, -1, -1);
			showWord();
		} else if (command.equals(HIDETAGGED)) {
			calcScoreSums();
			if (!answerHasBeenShown) {
				showWord();
			} else {
				showWord();
				showAnswer();
			}
		} else if (command.equals(ABOUT)) {
			if(aboutBox == null) {
				AboutPanel panel = new AboutPanel(aboutBox);
				panel.setLayout(new BorderLayout());

				aboutBox = new JDialog(parentFrame, "About Incipient Flash Cards!", false);
				aboutBox.getContentPane().add(panel, BorderLayout.CENTER);

				JPanel buttonpanel = new JPanel();
				buttonpanel.setOpaque(false);
				JButton button = (JButton) buttonpanel.add(new JButton("OK"));
				panel.add(buttonpanel, BorderLayout.SOUTH);
				button.addActionListener(new OkAction(aboutBox));
			}
			aboutBox.pack();
			Point p = parentFrame.getLocationOnScreen();
			aboutBox.setLocation(p.x + 125, p.y + 50);
			aboutBox.setVisible(true);
		} else if (command.equals(ENTERWORDS)) {
			//			if(enterTable == null) {
			//				EnterWords wordsPanel = new EnterWords(enterTable);
			//				wordsPanel.setLayout(new BorderLayout());
			//
			//				enterTable = new JDialog(parentFrame, "Add some new words", false);
//				enterTable.getContentPane().add(wordsPanel, BorderLayout.CENTER);
//
//				JPanel buttonpanel = new JPanel();
//				buttonpanel.setOpaque(false);
//				JButton button = (JButton) buttonpanel.add(new JButton("OK"));
//				wordsPanel.add(buttonpanel, BorderLayout.SOUTH);
//				button.addActionListener(new OkAction(enterTable));
//			}
//			enterTable.pack();
//			Point p = parentFrame.getLocationOnScreen();
//			enterTable.setLocation(p.x + 125, p.y + 50);
//			enterTable.show();
            } else if (command.equals(PROGRESS)) {
                saveData();
            } else if (command.equals(QUITWOSAVING)) {
                System.exit(0);
            } else if (command.equals(QUIT)) {
                saveData();
                System.exit(0);
            }
        }

        private void calcScoreSums() {
            questionsLeft = 0;
//            answersLeft = 0;
            qScoreTotal = 0;
            aScoreTotal = 0;
            qScoreArray = new int[hash.size()];
            aScoreArray = new int[hash.size()];
            Word trav;

            for (int i=0; i<hash.size(); i++) {
                trav = (Word)hash.get(i+"");
                trav.calcScores();
                if (trav.q_score >= threshold) {
                    qScoreArray[i] = trav.q_score;
                    qScoreTotal += trav.q_score;
                    questionsLeft++;
                } else {
                    qScoreArray[i] = 0;
                }
                if (trav.a_score >= threshold) {
                    aScoreArray[i] = trav.a_score;
                    aScoreTotal += trav.a_score;
//                    answersLeft++;
                } else {
                    aScoreArray[i] = 0;
                }
            }
        }

        private void pickNewWord(int prevWord, int skipped, int wasWrong, int wasInWrongHistory, int wasCyclingWrongHistory) {
            int target, count, attempt=0;
            int scoreTotal;
            int[] scoreArray;
            boolean brandNew;

            previousWordAction[0] = prevWord;
            previousWordAction[1] = skipped;
            previousWordAction[2] = wasWrong;
            previousWordAction[3] = wasInWrongHistory;
            previousWordAction[4] = wasCyclingWrongHistory;

            if (currentQ2A) {
                scoreTotal = qScoreTotal;
                scoreArray = qScoreArray;
                brandNew = brandNewQuestions;
            } else {
                scoreTotal = aScoreTotal;
                scoreArray = aScoreArray;
                brandNew = brandNewAnswers;
            }

            if  (wrongHistory.size() >= 5 || (scoreTotal == 0 && wrongHistory.size() > 0)) {
                wrongFlag = true;
            }
            if  (wrongHistory.size() == 0) {
                wrongFlag = false;
            }

            if (wrongFlag) {
                do {
                    target = (int)(Math.random()*wrongHistory.size());
                    attempt++;
                } while (wrongHistory.elementAt(target).equals(previousWordAction[0]+"") && attempt < 10);
                if (attempt < 10) {  // avoid picking the same word over and over
                    currentWord = (Word)hash.get(wrongHistory.elementAt(target));
                    return;
                }
            }

            if (scoreTotal == 0 && menuBrandNewOption.getState() && brandNew && threshSlider.getValue() == (currentQ2A?RESET_Q_SCORE:RESET_A_SCORE)) {
                threshSlider.setValue(2);
                if (currentQ2A) {
                    brandNewQuestions = false;
                    scoreTotal = qScoreTotal;
                    scoreArray = qScoreArray;
                } else {
                    brandNewAnswers = false;
                    scoreTotal = aScoreTotal;
                    scoreArray = aScoreArray;
                }
            }
            if (scoreTotal == 0 && menuBrandNewOption.getState() && brandNewAnswers && !brandNewQuestions && threshSlider.getValue() == 2) {
                answerButton.setSelected(true);
                scoreTotal = aScoreTotal;
                scoreArray = aScoreArray;
                threshSlider.setValue(RESET_A_SCORE);
            }

            while (scoreTotal == 0 && threshold > 2) {
                threshSlider.setValue(threshSlider.getValue()-1);
                if (currentQ2A) {
                    scoreTotal = qScoreTotal;
                    scoreArray = qScoreArray;
                } else {
                    scoreTotal = aScoreTotal;
                    scoreArray = aScoreArray;
                }
            }
            if (scoreTotal == 0) {
                menuTaggedOption.setState(false);
                actionPerformed(new ActionEvent(parentFrame, 1, HIDETAGGED));
                if (currentQ2A) {
                    scoreTotal = qScoreTotal;
                    scoreArray = qScoreArray;
                } else {
                    scoreTotal = aScoreTotal;
                    scoreArray = aScoreArray;
                }
            }
            if (scoreTotal == 0) {
                menuRefreshOption.setState(true);
                actionPerformed(new ActionEvent(parentFrame, 1, REFRESH));
                if (currentQ2A) {
                    scoreTotal = qScoreTotal;
                    scoreArray = qScoreArray;
                } else {
                    scoreTotal = aScoreTotal;
                    scoreArray = aScoreArray;
                }
            }
            if (scoreTotal == 0) {
                menuHighErrOption.setState(true);
                actionPerformed(new ActionEvent(parentFrame, 1, HIGHERR));
                if (currentQ2A) {
                    scoreTotal = qScoreTotal;
                    scoreArray = qScoreArray;
                } else {
                    scoreTotal = aScoreTotal;
                    scoreArray = aScoreArray;
                }
            }
            if (scoreTotal == 0) {
                menuLatestOption.setState(true);
                actionPerformed(new ActionEvent(parentFrame, 1, LATEST));
                if (currentQ2A) {
                    scoreTotal = qScoreTotal;
                    scoreArray = qScoreArray;
                } else {
                    scoreTotal = aScoreTotal;
                    scoreArray = aScoreArray;
                }
            }
            if (scoreTotal == 0) {
                threshSlider.setValue(1);
                if (currentQ2A) {
                    scoreTotal = qScoreTotal;
                    scoreArray = qScoreArray;
                } else {
                    scoreTotal = aScoreTotal;
                    scoreArray = aScoreArray;
                }
            }

            do {
                target = (int)(Math.random()*scoreTotal)+1;
                count = -1;
                do {
                    count++;
                    target -= scoreArray[count];
                } while (target > 0);
                attempt++;
            } while (count == previousWordAction[0] && attempt < 10);

            currentWord = (Word)hash.get(count+"");
        }

        public void showWord() {
            int fontSize = 40;
            String temp;

            if (currentWord == null) {
                return;
            }

            if (previousWordAction[0] == -1) {
                menuUndoItem.setEnabled(false);
            } else {
                menuUndoItem.setEnabled(true);
            }

            showRight.requestFocus();
            questions.setVisible(false);
            temp = "Say the "+(currentQ2A?questionLabel:answerLabel)+" word  ";
            for (int i=0; i<wrongHistory.size(); i++) {
                temp += ".";
            }
            instrBar.setText(temp);
            questions.setText(currentQ2A?blender(currentWord.QuestionString):currentWord.AnswerString);
            answers.setText("");
            leftLabel.setText(questionsLeft+"");

            questions.setFont(new Font("Arial", Font.PLAIN, fontSize));
            while (questions.getMaximumSize().getWidth() > LABEL_SIZE) {
                fontSize -= 2;
                questions.setFont(new Font("Arial", Font.PLAIN, fontSize));
            }
            questions.setVisible(true);

            answerHasBeenShown = false;
        }

        public void showAnswer() {
            int fontSize = 40;

            if (currentWord == null) {
                return;
            }

            showRight.requestFocus();
            answers.setVisible(false);
            answers.setText(currentQ2A?currentWord.AnswerString:blender(currentWord.QuestionString));
            answers.setFont(new Font("Arial", Font.PLAIN, fontSize));
            while (answers.getMaximumSize().getWidth() > LABEL_SIZE) {
                fontSize -= 2;
                answers.setFont(new Font("Arial", Font.PLAIN, fontSize));
            }
            answers.setVisible(true);
            answerHasBeenShown = true;
        }

        public String blender(String text) {
            String[] line, subline;
            int[] lineKey, sublineKey;
            String temp, newtext = "";
//            int within;

            if (!menuBlenderOption.getState() || text.indexOf("[")>=0 || text.indexOf(" or ")>=0 || text.indexOf("...")>=0) {
                return text;
            }
            line = text.split(";");
            lineKey = random(line.length);
            for (int i = 0; i<line.length; i++) {
                newtext += i==0?"":"; ";
                temp = line[lineKey[i]];
                while (temp.startsWith(" ")) {
                    temp = temp.substring(1);
                }
                if (temp.startsWith("to ")) {
                    newtext += "to ";
                    temp = temp.substring(3);
                }
                temp = temp.replaceAll(", ", "|");
                if (text.indexOf("(")>=0) {
//                    within = 0;
                    for (int j = 0; j<temp.length(); j++) {
                        if (temp.charAt(j) == '(') {
//                            within++;
                        } else if (temp.charAt(j) == ')') {
//                            within--;
                        } else if (temp.charAt(j) == '|') {
                            temp = temp.substring(0, j) +";"+ temp.substring(j+1);
                        }
                    }
                    temp = temp.replaceAll(";", ", ");
                }

                subline = temp.split("\\|");
                sublineKey = random(subline.length);
                for (int k = 0; k<subline.length; k++) {
                    newtext += (k==0?"":", ")+subline[sublineKey[k]];
                }

            }
            return newtext;
        }

        public int[] random(int size) {
            Vector<String> source = new Vector<String>();
            int[] keys = new int[size];
            String temp;

            for (int i=0; i<size; i++) {
                source.addElement(i+"");
            }

            for (int i=size; i>0; i=i-1) {
                temp = source.elementAt((int)(Math.random()*i));
                keys[i-1] = Integer.valueOf(temp).intValue();
                source.removeElement(temp+"");
            }

            return keys;
        }

        public boolean checkAnswer(boolean correct) {
            int scoreBefore;
            boolean wasWrong = false;

            currentWord.addTimePoint(today, correct, currentQ2A);
            if (correct) {
                wasWrong = wrongHistory.remove(currentWord.hashNum+"");
            } else {
                wrongHistory.addElement(currentWord.hashNum+"");
            }


            if (currentQ2A) {
                scoreBefore = currentWord.q_score;
                currentWord.calcScores();
                qScoreTotal -= scoreBefore - currentWord.q_score;
                qScoreArray[currentWord.hashNum] = currentWord.q_score;
                if (currentWord.q_score < threshold) {
                    calcScoreSums();
                }
            } else {
                scoreBefore = currentWord.a_score;
                currentWord.calcScores();
                aScoreTotal -= scoreBefore - currentWord.a_score;
                aScoreArray[currentWord.hashNum] = currentWord.a_score;
                if (currentWord.a_score < threshold) {
                    calcScoreSums();
                }
            }

            return wasWrong;
        }


        public void undoLastAction() {
            Vector<TimePoint> history;
            Word.TimePoint t;

            currentWord = (Word)hash.get(previousWordAction[0]+"");

            if (previousWordAction[1] != 1) {
                if (currentQ2A) {
                    history = currentWord.qHistory;
                } else {
                    history = currentWord.aHistory;
                }
                t = (Word.TimePoint)history.elementAt(history.size()-1);
                if (t.performance.length() == 1) {
                    history.removeElementAt(history.size()-1);
                } else {
                    t.performance = t.performance.substring(0, t.performance.length()-1);
                }

                if (previousWordAction[2] == 1) {
                    wrongHistory.remove(currentWord.hashNum+"");
                }
                if (previousWordAction[3] == 1) {
                    wrongHistory.addElement(currentWord.hashNum+"");
                }
            }
            calcScoreSums();
            previousWordAction[0] = -1;
            if (previousWordAction[4] == 1) {
                wrongFlag = true;
            } else {
                wrongFlag = false;
            }
            showWord();
        }

        public Component createComponents() {
            pane = new JRootPane();
            pane.setBorder(BorderFactory.createEmptyBorder(0, 0, 20, 20)); //top, left, bottom, right
            pane.setLayout(new XYLayout());

            catBar = new JLabel(" 00 Everything", JLabel.LEFT);
            catBar.setFont(new Font("Arial", Font.BOLD, 20));
            pane.add(catBar, new XYConstraints(10, 10, 250, 32));
            catBar.setBorder(BorderFactory.createEtchedBorder());

            instrBar = new JLabel("", JLabel.LEFT);
            instrBar.setFont(new Font("Arial", Font.BOLD, 16));
            pane.add(instrBar, new XYConstraints(45, 60, 300, 32));

            questions = new JLabel("", JLabel.CENTER);
            questions.setFont(new Font("Arial", Font.PLAIN, 40));
            pane.add(questions, new XYConstraints(45, 90, LABEL_SIZE, 50));

            answers = new JLabel("", JLabel.CENTER);
            answers.setFont(new Font("Arial", Font.PLAIN, 40));
            pane.add(answers, new XYConstraints(45, 160, LABEL_SIZE, 50));

            questions.setBorder(BorderFactory.createEtchedBorder());
            answers.setBorder(BorderFactory.createEtchedBorder());

            showRight = new JButton(SHOW_ANSWER);
            showRight.setFont(new Font("Arial", Font.BOLD, 16));
            showRight.addActionListener(this);
            showRight.setActionCommand(SHOW_ANSWER);

            skipWrong = new JButton(SKIP_WORD);
            skipWrong.setFont(new Font("Arial", Font.BOLD, 16));
            skipWrong.addActionListener(this);
            skipWrong.setActionCommand(SKIP_WORD);

            pane.add(showRight, new XYConstraints(170, 300, 100, 30));
            pane.add(skipWrong, new XYConstraints(325, 300, 100, 30));

            leftLabel = new JLabel("", JLabel.CENTER);
            leftLabel.setFont(new Font("Arial", Font.PLAIN, 25));
            pane.add(leftLabel, new XYConstraints(500, 10, 80, 40));
            leftLabel.setVisible(false);
            //		leftLabel.setBorder(BorderFactory.createEtchedBorder());

            whichButton = new ButtonGroup();
            answerButton = new JRadioButton("Question", false);
            whichButton.add(answerButton);
            answerButton.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent ie) {
                    if (((JRadioButton)ie.getItem()).isSelected()) {
                        currentQ2A = false;
                    } else {
                        currentQ2A = true;
                    }
                    answers.setText("");
                    wrongHistory.removeAllElements(); // switch to a check for appropriate elements??
                    if (!answerHasBeenShown) {
                        showWord();
                    } else {
                        showWord();
                        showAnswer();
                    }
                }
            });
            questionButton = new JRadioButton("Answer", true);
            currentQ2A = true;
            whichButton.add(questionButton);
            pane.add(answerButton, new XYConstraints(10, 310, 80, 15));
            pane.add(questionButton, new XYConstraints(10, 330, 80, 15));

            pane.add(new JPanel(), new XYConstraints(0, 200, 570, 89)); // dummy the right size

            threshSlider = new JSlider(JSlider.VERTICAL, 1, MAXSCORE, 2);
            threshSlider.setInverted(true);
            threshSlider.addChangeListener(new ChangeListener() {
                public void stateChanged(ChangeEvent ce) {
                    JSlider slider = (JSlider)ce.getSource();
                    threshold = slider.getValue();
                    calcScoreSums();
                    if (!answerHasBeenShown) {
                        showWord();
                    } else {
                        showWord();
                        showAnswer();
                    }
                    slider.requestFocus();
                }
            });

            thresholdLabel = new JLabel("", JLabel.LEFT);
            thresholdLabel.setFont(new Font("Arial", Font.BOLD, 20));
            pane.add(thresholdLabel, new XYConstraints(45, 140, 310, 20));
            threshSlider.addFocusListener(new FocusListener() {
                public void focusGained(FocusEvent ce) {
                    JSlider slider = (JSlider)ce.getSource();
                    if (currentWord != null) {
                        thresholdLabel.setText(" Score threshold = "+slider.getValue());
                    }
                }
                public void focusLost(FocusEvent ce) {
                    thresholdLabel.setText("");
                }
            });
            pane.add(threshSlider, new XYConstraints(10, 90, 25, 120));

            return pane;
        }

        public void saveData() {
            try {
                (new File(DBFILE)).renameTo(new File(DBFILE+".bak"));
                System.err.println("Saving database...");
                PrintWriter writer = new PrintWriter(new OutputStreamWriter(new FileOutputStream(DBFILE), "UnicodeLittle"));
                Word trav;

                writer.println(questionLabel);
                writer.println(answerLabel);

                writer.println("threshold="+threshSlider.getValue()+" "+(menuRefreshOption.getState()?"Refresh ":"")+(menuHighErrOption.getState()?"HighErrors ":"")
							 +(menuLatestOption.getState()?"Latest ":"")+(menuBrandNewOption.getState()?"BrandNew ":"")+(menuTaggedOption.getState()?"HideTagged ":"")
							 +(menuBlenderOption.getState()?"Blender ":"")+(currentQ2A?"":"Answers ")+(menuLeftOption.getState()?"Left ":""));

                for (int i=0; i<hash.size(); i++) {
                        trav = (Word)hash.get(i+"");
                        writer.println(trav.toString());
                }
                writer.close();
        } catch (IOException ioe) {
                ioe.printStackTrace();
        }
    }

	public void loadData() {
		try {
			BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(DBFILE), "UnicodeLittle"));
			StringTokenizer st, st2;
			int count = 0;
			String temp;
			Word trav;
			boolean nowQ;

			brandNewQuestions = brandNewAnswers = false;

			questionLabel = reader.readLine();
			answerLabel = reader.readLine();

			questionButton.setText(questionLabel);
			answerButton.setText(answerLabel);

			st = new StringTokenizer(reader.readLine());
			while (st.hasMoreTokens()) {
				temp = st.nextToken();
				if (temp.equals("Refresh")) {
					menuRefreshOption.setState(true);
				}
				if (temp.equals("HighErrors")) {
					menuHighErrOption.setState(true);
				}
				if (temp.equals("Latest")) {
					menuLatestOption.setState(true);
				}
				if (temp.equals("BrandNew")) {
					menuBrandNewOption.setState(true);
				}
				if (temp.equals("HideTagged")) {
					menuTaggedOption.setState(true);
				}
				if (temp.equals("Blender")) {
					menuBlenderOption.setState(true);
				}
				if (temp.startsWith("threshold=")) {
					threshold = Integer.valueOf(temp.substring(temp.indexOf("=")+1)).intValue();
					threshSlider.setValue(threshold);
				}
				if (temp.equals("Answers")) {
					currentQ2A = false;
					answerButton.setSelected(true);
				}

				if (temp.equals("Left")) {
					menuLeftOption.setState(true);
					leftLabel.setVisible(true);
				}
			}

			while (reader.ready()) {
				st = new StringTokenizer(reader.readLine(),"|");
				trav = new Word(st.nextToken(), st.nextToken(), st.nextToken(), hash.size());
				nowQ = true;
				while (st.hasMoreTokens()) {
					temp = st.nextToken();
					if (temp.equals("<-Q and A->")) {
						nowQ = false;
					} else {
						st2 = new StringTokenizer(temp, ":");
						trav.addTimeSeries(Integer.valueOf(st2.nextToken()).intValue(), st2.nextToken(), nowQ);
					}
				}
				if (trav.qHistory.size() == 0 || (trav.qHistory.size() == 1 && !((Word.TimePoint)trav.qHistory.elementAt(0)).caughtUp(RESET_Q_SCORE))) {
					brandNewQuestions = true;
				}
				if (trav.aHistory.size() == 0 || (trav.aHistory.size() == 1 && !((Word.TimePoint)trav.aHistory.elementAt(0)).caughtUp(RESET_A_SCORE))) {
					brandNewAnswers = true;
				}

				for (int i=trav.shouldBeInWrongHistory(); i>0; i--) {
					wrongHistory.addElement(count+"");
				}
				hash.put(count+"", trav);
				count++;
			}
			reader.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}

		calcScoreSums();

		if (menuBrandNewOption.getState() && (brandNewQuestions || brandNewAnswers)) {
			menuRefreshOption.setState(false);
			menuHighErrOption.setState(false);
			menuLatestOption.setState(false);
			if (brandNewQuestions) {
				threshSlider.setValue(RESET_Q_SCORE);
				currentQ2A = true;
				questionButton.setSelected(true);
			} else if (brandNewAnswers) {
				threshSlider.setValue(RESET_A_SCORE);
				currentQ2A = false;
				answerButton.setSelected(true);
			}
			calcScoreSums();
		}
	}

	public void initialize() {
		today = (int)((new Date()).getTime()/86400000);
		hash = new Hashtable<String,Word>();
		wrongHistory = new Vector<String>();
		loadData();
		System.err.println(hash.size());

		pickNewWord(-1, -1, -1, -1, -1);
		showWord();

	}

	public void exportTo_iFlipr() {
        try {
            (new File(DBFILE+".out")).renameTo(new File(DBFILE+".out"+".bak"));
            System.err.println("Saving database...");
            PrintWriter writer = new PrintWriter(new OutputStreamWriter(new FileOutputStream(DBFILE+".out"), "UnicodeLittle"));
            Word trav;

            for (int i=0; i<hash.size(); i++) {
            	trav = hash.get(i+"");
            	writer.println(trav.iFliprExport());
            }
            writer.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
	}
	
	public static void main(String[] args) {
		try {
			UIManager.setLookAndFeel("com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
		} catch (Exception e) {
			System.err.println("Failed loading LookandFeel: com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
			System.err.println(e);
			try {
				UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
			} catch (Exception e2) {
				System.err.println("Failed loading CrossPlatformLookAndFeel");
				System.err.println(e2);
			}
		}

		JFrame frame = new JFrame("Incipient Flash");
		frame.setLocation(300,200);
		Incipient inc=new Incipient(frame);
		frame.getContentPane().add(inc.createComponents(), BorderLayout.CENTER);
		//Finish setting up the frame, and show it.

		frame.addKeyListener(inc);
		frame.addWindowListener(new WindowIssues(inc));

		frame.pack();
		frame.setVisible(true);

		inc.initialize();
		inc.exportTo_iFlipr();
	}
}

class WindowIssues extends WindowAdapter {
	private Incipient inc;
	
	protected WindowIssues(Incipient inc) {
		this.inc = inc;
	}

	public void windowClosing(WindowEvent e) {
		inc.saveData();
		System.exit(0);
	}
}

class AboutPanel extends JPanel {   // can we get rid of this please?
	public static final long serialVersionUID = 1L;
	private ImageIcon aboutimage = null;
	public AboutPanel(JDialog aboutbox) {
		aboutimage = new ImageIcon("kem_goz.jpg");
		setOpaque(false);
	}

	public void paint(Graphics g) {
		aboutimage.paintIcon(this, g, 0, 0);
		super.paint(g);
	}

	public Dimension getPreferredSize() {
		return new Dimension(aboutimage.getIconWidth(), aboutimage.getIconHeight());
	}
}

class IPanel extends JPanel {
	public static final long serialVersionUID = 1L;
	private int prefX;
	private int prefY;

	public IPanel(int x, int y) {
		prefX = x;
		prefY = y;
//	    setOpaque(false);
	}

	public Dimension getPreferredSize() {
		return new Dimension(prefX, prefY);
	}

}

//class EnterWords extends JPanel {
// //	Incipient inc = null;
//	JDialog enterTable;
//    JTable      tableView;
// //    Dimension   origin = new Dimension(0, 0);
//    JScrollPane mainTable;
//    JScrollPane scrollpane;
//
//
//
//	public EnterWords(JDialog enterTable) {
// //	    this.inc = inc;
//		this.enterTable = enterTable;
// //	    setOpaque(false);
//
//		mainTable = createTable();
//	}
//
//	public void doit() {
//		enterTable.add(mainTable, BorderLayout.CENTER);
//	}
//
//
//    public JScrollPane createTable() {
//        final String[] names = {"Turkish", "Engligh", "Word Group"};
//        final String[][] data = {
//		  {"biber", "pepper", "food"},
//		  {"domates", "tomato", "food"},
//		  {"tavla", "backgammon", "game"},
//        };
//
//        // Create a model of the data.
//        TableModel dataModel = new AbstractTableModel() {
//            public int getColumnCount() { return names.length; }
//            public int getRowCount() { return data.length;}
//            public Object getValueAt(int row, int col) {return data[row][col];}
//            public String getColumnName(int column) {return names[column];}
//            public Class getColumnClass(int c) {return getValueAt(0, c).getClass();}
//	        public boolean isCellEditable(int row, int col) {return col != 5;}
//            public void setValueAt(String aValue, int row, int column) { data[row][column] = aValue; }
//         };
//
//        tableView = new JTable(dataModel);
//        tableView.setRowHeight(33);
//        scrollpane = new JScrollPane(tableView);
//        return scrollpane;
//    }
//
//}

class OkAction extends AbstractAction {
	public static final long serialVersionUID = 1L;
	private JDialog aboutBox;

    protected OkAction(JDialog aboutBox) {
        super("OkAction");
        this.aboutBox = aboutBox;
    }

    public void actionPerformed(ActionEvent e) {
        aboutBox.setVisible(false);
    }
}
