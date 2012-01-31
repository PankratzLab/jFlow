package dead;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

import stats.Anova;
import stats.LeastSquares;
import stats.LogisticRegression;

import common.*;

public class compGroups {
	public static String DB_FILE = "crf_db.dat";
	public static String[] USE = {"BirthDate", "AgeAtOnset", "AgeAtExam", "Duration", "AOO", "DurationFromAOO", "EarlyOnset45", "EarlyOnset50", "EarlyOnset60", "EarlyDisease", "LateDisease", "MidDisease", "Affected", "VPD", "CONF_PD", "Male", "Caucasian", "Hispanic", "AffFather", "AffMother", "AffParent", "parkin", "polymorph", "G2019S", "Depression", "Depressed", "DSMIV_Dx_Depression", "MajorDepression", "MinorDepression", "MMSE", "Demented", "BlessedFunctionality", "Education", "MilitaryYears", "smoked", "smoker", "alcohol", "pesticides", "OnsetWithTremor", "DominantSideFirst", "LeftSideFirst", "RightSideFirst", "BothSidesFirst", "HallucinationWithDrugs", "HallucinationsWithoutDrugs", "depressionBefore", "depressionSince", "HeadInjury", "Infection", "SchwabExaminer", "SchwabSubject", "SchwabDiff", "Hoehn&Yahr", "Bradykinesia", "Rigidity", "Instability", "PersistentAsymmetry", "RestTremor", "ProgressiveDosirder", "levodopaChorea", "levodopa5PlusYears", "Course10+", "UnexplainedMotor", "Strokes", "Encephalitis", "OculogyricCrisis", "Alzheimers", "SensoryDeficits_Apraxia", "PDfromDopaDrugs", "Remission", "unilateral3PlusYears", "SupranuclearGaze", "CerebellarSigns", "Hypotension", "NoResponseLDopa", "lesionMRI", "ldopaResponse", "PDopinion", "PD>90%", "logisticE4", "logisticE2", "APOE4count", "UPDRSliving", "UPDRSmotor", "BradykinesiaScore", "BradykinesiaSubScore", "RigiditySubScore", "PIGD_score", "Tremor_score", "PIGD_scale", "PIGD_dominant", "PIGD_intermediate", "Tremor_dominant", "SpeechScore", "RestTremorScore", "ActionTremorScore", "CombinedTremorScore"};
	public static String[] IGNORE = {"UniqueID", "FamID", "IndID", "Site", "OnsetDate", "DiagnosisDate", "ExamDate", "Dx", "Gender", "BirthMonth", "BirthDay", "BirthYear", "Race", "OtherLRRK2", "InitialSymptom", "Levodopa(mg)", "Levodopa(yrs)", "Speech", "Salivation", "Swallowing", "Handwriting", "Utensils", "Dressing", "Hygiene", "BedTurning", "Falling", "FreezeWalk", "RightHanded", "Walking", "Tremor", "Sensory", "Speech_motor", "FacialExpression", "T_face", "T_r_hand", "T_l_hand", "T_r_foot", "T_l_foot", "PT_r_hand", "PT_l_hand", "R_neck", "R_r_ue", "R_l_ue", "R_r_le", "R_l_le", "FingerTabsR", "FingerTapsL", "HandR", "HandL", "RapidAltR", "RapidAltL", "LegR", "LegL", "Arising", "Posture", "Gait", "PosturalStability", "BradykinesiaUPDRS", "Ethnicity", "Onset>20", "S_Examiner_dur", "S_Subject_dur", "HY_dur", "S_Examiner_durAOO", "S_Subject_durAOO", "HY_durAOO", "S_Examiner_regress", "S_Subject_regress", "HY_regress", "S_Examiner_regressAOO", "S_Subject_regressAOO", "HY_regressAOO", "MMSEcutoff", "MMSE_edu_adj", "ExpressionScore", "ChairScore", "PostureScore", "GaitScore", "InstabilityScore", "APOE"};

	public compGroups(String filename) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		Pattern delimiters = Pattern.compile("[\\-\\s]+");
		String trav;
		Hashtable<String,Vector<String>> hash = new Hashtable<String,Vector<String>>();
		Vector<String> v = null;
		int[] phenos;
		int count, group, correction, count2;
		String[] groups, phenoNames, dataline;
		Vector<Vector<String[]>> data;
		double[] dep;
		double[][] indep, data_doubled;
		double p;
		boolean binaryData, newGroup;
		int lowestVal = 0;

		if (!new File(filename).exists()) {
			System.err.println("Error - could not find "+filename+" in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader(filename));

		hash.put("groups", new Vector<String>());
		newGroup = true;
		while (reader.ready()) {
			trav = reader.readLine();
			line = delimiters.split(trav);
			if (line.length==1&&line[0].equals("")) {
				newGroup = true;
			} else if (newGroup) {
				if (trav.equals("groups")) {
					System.err.println("Error - 'groups' is a reserved word and cannot be the complete name of a group");
					System.exit(4);
				}
				hash.get("groups").add(trav);
				hash.put(trav, v = new Vector<String>());
				newGroup = false;
			} else {
				if (line.length==1) {
					v.add(line[0]);
				} else if (line.length==2) {
					v.add((Integer.valueOf(line[0]).intValue()*1000+Integer.valueOf(line[1]).intValue())+"");
				} else {
					System.err.println("Error - compGroups can only compare either families or individuals (so group member can only be 1-2 words):");
					System.err.println("        "+trav);
					System.exit(5);
				}
			}
		}
		v = hash.get("groups");
		groups = new String[v.size()];
		for (int i = 0; i<v.size(); i++) {
			groups[i] = v.elementAt(i);
		}
		hash.remove("groups");

		if (!new File(DB_FILE).exists()) {
			System.err.println("Error - could not find database '"+DB_FILE+"' in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader(DB_FILE));

		v = new Vector<String>();
		phenoNames = reader.readLine().split("[\\s]+");
		for (int i = 0; i<phenoNames.length; i++) {
			if (containsStr(phenoNames[i], USE)) {
				v.add(i+"");
			} else if (!containsStr(phenoNames[i], IGNORE)) {
				System.err.println("Error - Unexpected database variable name '"+phenoNames[i]+"'");
				System.err.println("        recompile with variable name added to either the USE list or the IGNORE list");
				System.err.println("        (ignoring this variable for the time being");
			}
		}
		phenos = new int[v.size()];
		for (int i = 0; i<v.size(); i++) {
			phenos[i] = Integer.valueOf(v.elementAt(i)).intValue();
		}

		data = new Vector<Vector<String[]>>();
		for (int i = 0; i<groups.length; i++) {
			data.add(new Vector<String[]>());
		}

		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			if (phenoNames.length!=line.length) {
				System.err.println("Error - different number of values ("+line.length+" versus "+phenoNames.length+" phenotypes) for "+line[0]);
				System.exit(1);
			}
			trav = line[0];
			group = -1;
			for (int i = 0; i<groups.length; i++) {
				v = hash.get(groups[i]);
				if (v.contains(trav)) {
					if (group==-1) {
						group = i;
					} else {
						System.err.println("Error - \""+trav+"\" belongs to both group "+groups[group]+" and group "+groups[i]);
						System.exit(1);
					}
				}
			}
			if (group>=0) {
				dataline = new String[phenos.length];
				for (int i = 0; i<phenos.length; i++) {
					dataline[i] = line[phenos[i]];
				}
				data.elementAt(group).add(dataline);
			}
		}
		// printData(data, groups);

		writer = new PrintWriter(new FileWriter(filename+"-comparison.xls"));
		writer.print("Trait");
		for (int i = 0; i<groups.length; i++) {
			writer.print("\t"+groups[i]);
		}
		writer.println("\tANOVA\tlinear\tt-test\tlogistic\tchi-square");

		for (int pheno = 0; pheno<phenos.length; pheno++) {
			try {
				data_doubled = new double[groups.length][];
				binaryData = true;

				count = 0;
				v = new Vector<String>();
				for (int i = 0; i<groups.length; i++) {
					correction = 0;
					for (int j = 0; j<data.elementAt(i).size(); j++) {
						dataline = data.elementAt(i).elementAt(j);
						if (!dataline[pheno].equals(".")) {
							count++;
							if (!v.contains(dataline[pheno])) {
								v.add(dataline[pheno]);
							}
						} else {
							correction++;
						}
					}
					data_doubled[i] = new double[data.elementAt(i).size()-correction];
				}
				if (v.size()>2) {
					binaryData = false;
				} else if (v.size()==2) {
					if (Math.abs(Integer.valueOf(v.elementAt(0)).intValue()-Integer.valueOf(v.elementAt(1)).intValue())!=1) {
						System.err.println("Error: Phenotype "+pheno+" ("+phenoNames[phenos[pheno]]+") has 2 values that are not sequential ("+v.elementAt(0)+" and "+v.elementAt(1)+")");
					}
					lowestVal = Integer.valueOf(v.elementAt(0)).intValue()<Integer.valueOf(v.elementAt(1)).intValue()?Integer.valueOf(v.elementAt(0)).intValue():Integer.valueOf(v.elementAt(1)).intValue();
				} else if (v.size()==1) {
					lowestVal = Integer.valueOf(v.elementAt(0)).intValue();
				} else {
					System.err.println("Warning - did you know there was no data for any group for "+phenoNames[phenos[pheno]]+"?");
				}

				dep = new double[count];
				indep = new double[count][groups.length];

				count = 0;
				for (int i = 0; i<groups.length; i++) {
					count2 = 0;
					for (int j = 0; j<data.elementAt(i).size(); j++) {
						dataline = data.elementAt(i).elementAt(j);
						if (!dataline[pheno].equals(".")) {
							dep[count] = Double.valueOf(dataline[pheno]).doubleValue();
							indep[count][i] = 1;

							data_doubled[i][count2] = Double.valueOf(dataline[pheno]).doubleValue();

							count++;
							count2++;
						}
					}
				}
				if (!binaryData) {
					Anova a = new Anova(data_doubled);
					a.oneway();
					writer.print(pheno+a.getMeans());
					writer.print("\t"+ext.formDeci(a.getSig(), 5, true)+(a.getSig()<0.05?"*":""));

					LeastSquares ls = new LeastSquares(dep, indep);
					p = ls.getFsig();
					writer.print("\t"+ext.formDeci(p, 5, true)+(p<0.05?"*":""));
					writer.print("\t        ");

					writer.print("\t        ");
					writer.print("\t        ");
					writer.print("\t"+phenoNames[phenos[pheno]]);

					writer.print(a.getStdDevs());
				}

				if (binaryData) {
					Anova a = new Anova(data_doubled);
					a.oneway();
					writer.print(pheno+a.getMeans());

					int[][] chiCounts = new int[groups.length][2];
					int[] depInt = new int[dep.length];

					for (int i = 0; i<dep.length; i++) {
						depInt[i] = (int)dep[i];
						for (int j = 0; j<indep[i].length; j++) {
							chiCounts[j][depInt[i]-lowestVal] += indep[i][j];
						}
					}

					// System.out.println(phenoNames[phenos[pheno]]);
					// for (int i=0; i<chiCounts.length; i++ ) {
					// System.out.println(chiCounts[i][0] +"\t"+
					// chiCounts[i][1]);
					// }
					// System.out.println();
					// System.out.println();

					LogisticRegression lr = new LogisticRegression(depInt, indep);
					writer.print("\t        ");
					writer.print("\t        ");
					writer.print("\t        ");

					writer.print(lr.analysisFailed()?"\tfailed":"\t"+ext.formDeci(lr.getSigs()[0], 5, true)+(lr.getSigs()[0]<0.05?"*":""));
					writer.print("\t        ");
					writer.print("\t"+phenoNames[phenos[pheno]]);

					writer.print(a.getCounts());
				}
				writer.println();

			} catch (Exception e) {
				writer.println();
			}

		}

		writer.close();

	}

	public boolean containsStr(String target, String[] list) {
		for (int i = 0; i<list.length; i++) {
			if (list[i].equals(target)) {
				return true;
			}
		}
		return false;
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		// String filename = "in.dat";
		// String filename = "Tremor-BradyRigid.dat";
		// String filename = "R98Q comparison.txt";
		// String filename = "novel_lrrks_in.dat";
		String filename = "PIGD.dat";

		String usage = "\n"+"park.compGroups requires 1 argument\n"+"   (1) filename (i.e. 'file=in.dat')\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (args.length==0) {
			System.err.println("Using default filename ('"+filename+"').");
		}
		try {
			new compGroups(filename);
			// new compGroups("in-All3.dat");
			// new compGroups("in-G2019S.dat");
			// new compGroups("in-R1514Q.dat");
			// new compGroups("in-R1514Q_G2019S.dat");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
