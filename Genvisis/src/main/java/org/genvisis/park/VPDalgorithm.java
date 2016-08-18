package org.genvisis.park;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;
import org.genvisis.db.crfDB;

public class VPDalgorithm {
	public static boolean CARES_NOT_PROGENI = false;
	public static final String CRF_A_FILE = tools.CRF_DIR+(CARES_NOT_PROGENI?"CARES_":"")+"crf_a1.csv";
	public static final String CRF_B_FILE = tools.CRF_DIR+(CARES_NOT_PROGENI?"CARES_":"")+"crf_b.csv";
	public static final String CRF_F_FILE = tools.CRF_DIR+(CARES_NOT_PROGENI?"CARES_":"")+"crf_f.csv";
	// public static final String BIRTH_DATES_FILE = tools.CRF_DIR+"ninfo1_BirthDates.csv";
	public static final String[] CRF_A_HEADER = {"SITE_NO", "FAM_NO", "INFO_DT", "SUBJ_NO", "NDOB", "QA", "QB", "Q1", "Q2A", "Q2B", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10", "Q11|1|3", "Q12", "Q13", "Q14", "Q15", "Q16", "Q17", "Q18", "Q19", "Q40", "Q41_TXT", "STAFF_NO", "ERROR_FLAG", "EXCEPT_FLAG", "DE_DT"};
	public static final String[] CRF_B_HEADER = {"SITE_NO", "FAM_NO", "INFO_DT", "SUBJ_NO", "NDOB", "Q1|0|4", "Q2|0|4", "Q3|0|4", "Q4|0|4", "Q5|0|4", "Q6|0|4", "Q7|0|4", "Q8|0|4", "Q9|0|4", "Q10|0|4", "Q11|0|4", "Q12|0|4", "Q13|0|4", "STAFF1", "Q14|0|4", "Q15|0|4", "Q16A|0|4", "Q16B|0|4", "Q16C|0|4", "Q16D|0|4", "Q16E|0|4", "Q17A|0|4", "Q17B|0|4", "Q18A|0|4", "Q18B|0|4", "Q18C|0|4", "Q18D|0|4", "Q18E|0|4", "Q19A|0|4", "Q19B|0|4", "Q20A|0|4", "Q20B|0|4", "Q21A|0|4", "Q21B|0|4", "Q22A|0|4", "Q22B|0|4", "Q23|0|4", "Q24|0|4", "Q25|0|4", "Q26|0|4", "Q27|0|4", "STAFF2", "Q28|0|100", "Q29|0|100", "Q30|0|5", "Q31", "Q31_TXT", "STAFF3", "ERROR_FLAG", "EXCEPT_FLAG", "DE_DT"};
	public static final String[] CRF_F_HEADER = {"SITE_NO", "FAM_NO", "INFO_DT", "SUBJ_NO", "NDOB", "Q1", "Q2|1|2", "Q3|1|7", "Q3_TXT", "Q4|0|1", "Q5|0|1", "Q6|0|1", "Q7|0|1", "Q8|0|1", "Q9|0|1", "Q10|0|1", "Q11|0|1", "Q12|0|1", "Q13|0|1", "Q14|0|1", "Q15|0|1", "Q16|0|1", "Q17|0|1", "Q18|0|1", "Q19|0|1", "Q20|0|1", "Q21|0|1", "Q22|0|1", "Q23|0|1", "Q24|0|1", "Q25|0|1", "Q26|0|1", "Q27|0|1", "Q28|1|4", "Q29|1|4", "Q30|0|42", "Q30_R37T", "Q30_R41T", "Q31", "Q31_TXT", "STAFF_NO", "FIELD__", "SITE_PAYMENT", "EXCEPT_PAY", "PAY_COMMENT", "ERROR_FLAG", "EXCEPT_FLAG", "DE_DT"};
	public static final String[] BIRTH_DATES_HEADER = {"FamID", "IndID", "DOB"};
	public static final String[] EXCLUSION_CRITERIA = {"UnexplainedMotor", "Strokes", "Encephalitis", "OculogyricCrisis", "Alzheimers", "SensoryDeficits_Apraxia", "PDfromDopaDrugs", "Remission", "UnilateralFor3PlusYears", "SupranuclearGaze", "CerebellarSigns", "Hypotension", "NoResponseLDopa", "LesionOnMRI"};
	public static final String[] SUPPORTING_CRITERIA = {"PersistentAsymmetry", "ProgressiveDisorder", "levodopaChorea", "levodopa5PlusYears", "Course10+"};

	public static void implementAlgorithm() {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, header, tempOfficialHeader;
		String temp, trav;
		Hashtable<String,DxInfo> hash = new Hashtable<String,DxInfo>();
		DxInfo info;
		int missing, countToVPD, countToUKBBC;

		try {
			reader = new BufferedReader(new FileReader(CRF_B_FILE));
			header = reader.readLine().split(",");
			tempOfficialHeader = new String[CRF_B_HEADER.length];
			for (int i = 0; i<tempOfficialHeader.length; i++) {
				tempOfficialHeader[i] = CRF_B_HEADER[i].split("\\|")[0];
			}
			ext.checkHeader(header, tempOfficialHeader, true);
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.indexOf("\"")>=0) {
					temp = crfDB.fixMarks(temp);
				}
				while (temp.indexOf(",,")!=-1) {
					temp = temp.replaceAll(",,", ",.,");
				}
				line = temp.split(",");

				trav = line[1]+ext.formNum(Integer.parseInt(line[3]), 3);
				if (hash.containsKey(trav)) {
					System.err.println("Error - multiple instances of "+trav+" in CRF_B");
				}
				info = new DxInfo();

				temp = line[ext.indexOfStr("Q27", header)];
				info.brady_B = temp;
				if (procDouble(temp)>0) {
					info.bradykinesia = "1";
					info.inclusion.add("Bradykinesia");
					info.brainBank_inclusion.add("Bradykinesia");
				} else if (crfDB.checkData(temp).equals(".")) {
					info.bradykinesia = ".";
				}

				for (int i = 0; i<5; i++) {
					temp = line[ext.indexOfStr("Q16"+(char)(i+65), header)];
					info.restTremor += procDouble(temp);
					info.restTremorMissing += crfDB.checkData(temp).equals(".")?1:0;
				}
				info.restT_B = info.restTremor+"";

				double actionTremor = 0;
				for (int i = 0; i<2; i++) {
					temp = line[ext.indexOfStr("Q17"+(char)(i+65), header)];
					actionTremor += procDouble(temp);
				}
				info.actionT_B = actionTremor+"";
				
				for (int i = 0; i<5; i++) {
					temp = line[ext.indexOfStr("Q18"+(char)(i+65), header)];
					info.rigidity += procDouble(temp);
					info.rigidityMissing += crfDB.checkData(temp).equals(".")?1:0;
				}
				info.rigid_B = info.rigidity+"";


				temp = line[ext.indexOfStr("Q26", header)];
				if (procDouble(temp)>0) {
					info.atLeast1++;
					info.brainBank_atLeast1++;
					info.inclusion.add("Postural Instability");
					info.brainBank_inclusion.add("Postural Instability");
				}
				if (crfDB.checkData(temp).equals(".")) {
					info.atLeast1missing++;
				}

				hash.put(trav, info);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+CRF_B_FILE+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+CRF_B_FILE+"\"");
			System.exit(2);
		}

		try {
			reader = new BufferedReader(new FileReader(CRF_A_FILE));
			header = reader.readLine().split(",");
			tempOfficialHeader = new String[CRF_A_HEADER.length];
			for (int i = 0; i<tempOfficialHeader.length; i++) {
				tempOfficialHeader[i] = CRF_A_HEADER[i].split("\\|")[0];
			}
			ext.checkHeader(header, tempOfficialHeader, true);
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.indexOf("\"")>=0) {
					temp = crfDB.fixMarks(temp);
				}

				while (temp.indexOf(",,")!=-1) {
					temp = temp.replaceAll(",,", ",.,");
				}
				line = temp.split(",");

				trav = line[1]+ext.formNum(Integer.parseInt(line[3]), 3);
				if (hash.containsKey(trav)) {
					info = hash.get(trav);
				} else {
					hash.put(trav, info = new DxInfo());
				}

				temp = line[ext.indexOfStr("Q11", header)];
				if (procDouble(temp)==1 || procDouble(temp)==2) {
					info.brainBank_atLeast3++;
					info.brainBank_supporting.add("Asymetric onset");
				}
				
				hash.put(trav, info);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+CRF_F_FILE+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+CRF_F_FILE+"\"");
			System.exit(2);
		}
				
		try {
			reader = new BufferedReader(new FileReader(CRF_F_FILE));
			header = reader.readLine().split(",");
			tempOfficialHeader = new String[CRF_F_HEADER.length];
			for (int i = 0; i<tempOfficialHeader.length; i++) {
				tempOfficialHeader[i] = CRF_F_HEADER[i].split("\\|")[0];
			}
			ext.checkHeader(header, tempOfficialHeader, true);
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.indexOf("\"")>=0) {
					temp = crfDB.fixMarks(temp);
				}
				while (temp.indexOf(",,")!=-1) {
					temp = temp.replaceAll(",,", ",.,");
				}
				line = temp.split(",");

				trav = line[1]+ext.formNum(Integer.parseInt(line[3]), 3);
				if (hash.containsKey(trav)) {
					info = hash.get(trav);
				} else {
					hash.put(trav, info = new DxInfo());
					info.atLeast1missing = 10;
				}

				temp = line[ext.indexOfStr("Q4", header)];
				if (procDouble(temp)>0) {
					info.onset_GT20 = "1";
					info.inclusion.add("Onset after 20");
				} else if (crfDB.checkData(temp).equals(".")) {
					info.onset_GT20 = ".";
					info.reasoning.add("Missing 'onset before 20' flag");
				} else {
					info.reasoning.add("Onset before 20");
				}

				temp = line[ext.indexOfStr("Q5", header)];
				info.brady_F = temp;
				if (info.bradykinesia.equals("1")||procDouble(temp)>0) {
					info.bradykinesia = "1";
					HashVec.addIfAbsent("Bradykinesia", info.inclusion);
					HashVec.addIfAbsent("Bradykinesia", info.brainBank_inclusion);
				} else if (info.bradykinesia.equals("0")||!crfDB.checkData(temp).equals(".")) {
					info.bradykinesia = "0";
					info.reasoning.add("No bradykinesia");
				}
				if (info.bradykinesia.equals(".")) {
					if (crfDB.checkData(temp).equals(".")) {
						info.reasoning.add("Bradykinesia info missing");
					}
					// info.bradykinesia = ".";
					info.bradykinesia = "0";
				}

				temp = line[ext.indexOfStr("Q29", header)];
				if (procInt(temp)==1||procInt(temp)==2) {
					info.prob_GTE50 = "1";
					info.inclusion.add("Neurologist impression >50%");
				} else if (crfDB.checkData(temp).equals(".")) {
					info.reasoning.add("Neurologist impression not noted");
					// info.prob_GTE50 = ".";
				} else {
					info.reasoning.add("Neurologist impression <50%");
				}

				temp = line[ext.indexOfStr("Q9", header)];
				info.restTremor += procDouble(temp);
				info.restTremorMissing += crfDB.checkData(temp).equals(".")?1:0;
				info.restT_F = temp;

				temp = line[ext.indexOfStr("Q6", header)];
				info.rigidity += procDouble(temp);
				info.rigidityMissing += crfDB.checkData(temp).equals(".")?1:0;
				info.rigid_F = temp;

				temp = line[ext.indexOfStr("Q7", header)];
				if (procDouble(temp)>0) {
					info.atLeast1++;
					info.brainBank_atLeast1++;
					HashVec.addIfAbsent("Postural Instability", info.inclusion);
					HashVec.addIfAbsent("Postural Instability", info.brainBank_inclusion);
				} else if (crfDB.checkData(temp).equals(".")) {
					info.atLeast1missing++;
				}

				temp = line[ext.indexOfStr("Q28", header)];
				if (procInt(temp)==1) {
					info.ldopa_GT50 = "1";
					info.atLeast2++;
					info.brainBank_atLeast3++;
					info.supporting.add("Positive response to levodopa");
					info.brainBank_supporting.add("Positive response to levodopa");
				} else if (crfDB.checkData(temp).equals(".")) {
					// info.ldopa_GT50 = ".";
					info.atLeast2missing++;
				}

				for (int i = 9; i<14; i++) {
					temp = line[ext.indexOfStr("Q"+(i==9?8:i), header)];
					if (temp.equals("1")) {
						info.atLeast2++;
						info.brainBank_atLeast3++;
						info.supporting.add(SUPPORTING_CRITERIA[i-9]);
						info.brainBank_supporting.add(SUPPORTING_CRITERIA[i-9]);
					}
					info.atLeast2missing += crfDB.checkData(temp).equals(".")?1:0;
				}

				for (int i = 0; i<14; i++) {
					temp = line[ext.indexOfStr("Q"+(i+14), header)];
					if (temp.equals("1")) {
						info.exclusion.add(EXCLUSION_CRITERIA[i]);
						if (i > 0 && i != 5) {
							info.brainBank_exclusion.add(EXCLUSION_CRITERIA[i]);
						}
					}
					info.exclusionMissing += crfDB.checkData(temp).equals(".")?1:0;
				}

				hash.put(trav, info);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+CRF_F_FILE+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+CRF_F_FILE+"\"");
			System.exit(2);
		}

		try {
			reader = tools.getNinfoReader(CARES_NOT_PROGENI?5:1);
			writer = new PrintWriter(new FileWriter(tools.CRF_DIR+(CARES_NOT_PROGENI?"CARES_":"PROGENI_")+"verifiedData.xln"));
			writer.println("UniqueID\tFamID\tIndID\tDx\tVPD\tCONF_PD\tAffected\tOnset_GT20\tBradykinesia\tBrady_B\tBrady_F\tRestT_B\tRestT_F\tActionT_B\tRigidity_B\tRigidity_F\tProbability_GTE50\tAtLeast1\tAtLeast2\tHasNoExclusionCriteria\tVerified\tNum VPD criteria\tInclusion\tSupporting\tReason NVPD\tUKBBC\tNum UKBBC criteria\tUKBBC Inclusion\tUKBBC Supporting\tReason not UKBBC");
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				trav = line[0]+ext.formNum(Integer.parseInt(line[1]), 3);
				writer.print(trav+"\t"+line[0]+"\t"+line[1]+"\t"+line[3]+"\t"+tools.isVPD(line[3])+"\t"+tools.isConfPD(line[3])+"\t"+tools.isAffected(line[3]));
				if (hash.containsKey(trav)) {
					info = hash.get(trav);
					if (info.restTremor>0) {
						info.atLeast1++;
						info.brainBank_atLeast1++;
						info.atLeast2++;
						info.inclusion.add("Rest tremor");
						info.brainBank_inclusion.add("Rest tremor");
						info.supporting.add("Rest tremor");
					}
					if (info.rigidity>0) {
						info.atLeast1++;
						info.brainBank_atLeast1++;
						info.supporting.add("Rigidity");
						info.brainBank_inclusion.add("Rigidity");
					}
					writer.print("\t"+info.onset_GT20+"\t"+info.bradykinesia+"\t"+info.brady_B+"\t"+info.brady_F+"\t"+info.restT_B+"\t"+info.restT_F+"\t"+info.actionT_B+"\t"+info.rigid_B+"\t"+info.rigid_F+"\t"+info.prob_GTE50);
					// writer.print("\t"+(info.atLeast1>=1?"1":(info.atLeast1missing
					// == 4?".":"0")));
					if (info.atLeast1>=1) {
						writer.print("\t1");
					} else {
						info.reasoning.add("No rigidity, tremor or instability");
						writer.print("\t0");
					}

					// writer.print("\t"+(info.atLeast2>=2?"1":(info.atLeast2missing
					// == 7?".":"0")));

					if (info.atLeast2>=2) {
						writer.print("\t1");
					} else {
						info.reasoning.add("No supporting criteria");
						writer.print("\t0");
					}
					// writer.print("\t"+(info.exclusion.size()>0?"0":(info.exclusionMissing
					// == EXCLUSION_CRITERIA.length?".":"1")));
					// writer.print("\t"+(info.exclusion.size()>0?"0":"1"));
					writer.print("\t"+(info.exclusion.size()>0?"0":(info.exclusionMissing==EXCLUSION_CRITERIA.length?"0":"1")));
					// writer.print("\t"+(info.onset_GT20.equals("1")&&info.bradykinesia.equals("1")&&info.prob_GTE50.equals("1")&&info.atLeast1>=1&&info.atLeast2>=2&&!((double)info.exclusionMissing
					// ==
					// EXCLUSION_CRITERIA.length)&&info.exclusion.size()==0?"1":"0"));
					countToVPD = 0;
					countToUKBBC = 0;
					if (info.onset_GT20.equals("1")) {
						countToVPD++;
					}
					if (info.bradykinesia.equals("1")) {
						countToVPD++;
						countToUKBBC++;
					}
					if (info.prob_GTE50.equals("1")) {
						countToVPD++;
					}
					if (info.atLeast1>=1) {
						countToVPD++;
					}
					if (info.brainBank_atLeast1>=1) {
						countToUKBBC++;
					}
					if (info.atLeast2>=2) {
						countToVPD++;
					}
					if (info.brainBank_atLeast3>=3) {
						countToUKBBC++;
					}
					if (info.exclusion.size()==0) {
						countToVPD++;
					}
					if (info.brainBank_exclusion.size()==0) {
						countToUKBBC++;
					}
					// writer.print("\t"+(info.onset_GT20.equals("1")&&info.bradykinesia.equals("1")&&info.prob_GTE50.equals("1")&&info.atLeast1>=1&&info.atLeast2>=2&&info.exclusion.size()==0?"1":"0"));
					writer.print("\t"+(countToVPD==6?"1":"0")+"\t"+countToVPD);
					writer.print("\t"+(info.inclusion.size()==0?".":ext.listWithCommas(Array.toStringArray(info.inclusion), true)));
					writer.print("\t"+(info.supporting.size()==0?".":ext.listWithCommas(Array.toStringArray(info.supporting), true)));
					missing = info.atLeast1missing+info.atLeast2missing+info.exclusionMissing+info.restTremorMissing+info.rigidityMissing;
					temp = (missing>20?"Missing all data"+(info.exclusion.size()+info.reasoning.size()>0?"; also ":""):(missing>10?"Missing quite a bit of data"+(info.exclusion.size()+info.reasoning.size()>0?"; also ":""):""))+(info.exclusion.size()>0?"Exclusion criteria: "+ext.listWithCommas(Array.toStringArray(info.exclusion), true)+(info.reasoning.size()>0?"; also ":""):"")+ext.listWithCommas(Array.toStringArray(info.reasoning), true);
					writer.print("\t"+(temp.equals("")?".":temp));

					writer.print("\t"+(countToUKBBC==4?"1":"0")+"\t"+countToUKBBC);
					writer.print("\t"+(info.brainBank_inclusion.size()==0?".":ext.listWithCommas(Array.toStringArray(info.brainBank_inclusion), true)));
					writer.print("\t"+(info.brainBank_supporting.size()==0?".":ext.listWithCommas(Array.toStringArray(info.brainBank_supporting), true)));
					writer.println("\t"+(info.brainBank_exclusion.size()==0?".":ext.listWithCommas(Array.toStringArray(info.brainBank_exclusion), true)));
				} else {
					writer.println(".\t.\t.\t.\t.\t.\t.\t.\t\t.\t0\t0\t0\t0\t0\t0\t.\t.\tNo CRFs\t0\t0\t.\t.\tNo CRFs");
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+"ninfo1"+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+"ninfo1"+"\"");
			System.exit(2);
		}

		/*
		 * Onset_GT20=#crf_f.csv:Q4 #BradyOR=+crf_f.csv:Q5+crf_b.csv:Q27
		 * Bradykinesia=$GTE,final:BradyOR,0.25
		 * Probability_GTE50=$LTE,crf_f.csv:Q29,2
		 * #RestTremorCriteriaList=+crf_f.csv:Q9+crf_b.csv:Q16A+crf_b.csv:Q16B+crf_b.csv:Q16C+crf_b.csv:Q16D+crf_b.csv:Q16E
		 * #RestTremorCriteria=$GTE,final:RestTremorCriteriaList,0.25
		 * #RigidityCriteriaList=+crf_f.csv:Q6+crf_b.csv:Q18A+crf_b.csv:Q18B+crf_b.csv:Q18C+crf_b.csv:Q18D+crf_b.csv:Q18E
		 * #RigidityCriteria=$GTE,final:RigidityCriteriaList,0.25
		 * #AtLeast1sum=+final:RestTremorCriteria+final:RigidityCriteria+crf_f.csv:Q7+crf_b.csv:Q26
		 * AtLeast1=$GTE,final:AtLeast1sum,1 #LDopa_GT50=$LTE,crf_f.csv:Q28,1
		 * #AtLeast2sum=+final:RestTremorCriteria+final:LDopa_GT50+crf_f.csv:Q8+crf_f.csv:Q10+crf_f.csv:Q11+crf_f.csv:Q12+crf_f.csv:Q13
		 * AtLeast2=$GTE,final:AtLeast2sum,2
		 * #ExclusionCriteriaList=+crf_f.csv:Q14+crf_f.csv:Q15+crf_f.csv:Q16+crf_f.csv:Q17+crf_f.csv:Q18+crf_f.csv:Q19+crf_f.csv:Q20+crf_f.csv:Q21+crf_f.csv:Q22+crf_f.csv:Q23+crf_f.csv:Q24+crf_f.csv:Q25+crf_f.csv:Q26+crf_f.csv:Q27
		 * HasNoExclusionCriteria=$LTE,final:ExclusionCriteriaList,0
		 * 
		 * Dx=#ninfo1.csv:IllnessStatusCode
		 * VPD=$RECODE,final:Dx,CONF_PD,1,VPD,1,NVPD,0,.
		 * CONF_PD=$EQUALS,final:Dx,CONF_PD
		 * Affected=$RECODE,final:Dx,CONF_PD,1,VPD,1,NVPD,1,RPD,1,NOEV,0,NRPD,0,NXAM,.,CONF,.,.
		 * 
		 * #VerifiedCount=+final:Onset_GT20+final:Bradykinesia+final:Probability_GTE50+final:AtLeast1+final:AtLeast2+final:HasNoExclusionCriteria
		 * Verified=$GTE,final:VerifiedCount,6
		 */

	}

	public static int procInt(String str) {
		double value = 0;

		if (crfDB.checkData(str).equals(".")||crfDB.checkData(str).equals("NA")) {
			value = 0;
		} else if (crfDB.checkData(str).equals("Y")) {
			value = 1;
		} else {
			try {
				value = Double.valueOf(str).doubleValue();
			} catch (Exception e) {
				System.err.println("Error - '"+str+"' is not a valid integer");
			}
		}

		return (int)value;
	}

	public static double procDouble(String str) {
		double value = 0;

		if (crfDB.checkData(str).equals(".")||crfDB.checkData(str).equals("NA")) {
			value = 0;
		} else if (crfDB.checkData(str).equals("Y")) {
			value = 1;
		} else {
			try {
				value = Double.parseDouble(str);
			} catch (Exception e) {
				System.err.println("Error - '"+str+"' is not a valid 'double' number");
			}
		}

		return value;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "VPDalgorithm.dat";

		String usage = "\\n"+"park.VPDalgorithm requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default))\n"+"";

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
		try {
			implementAlgorithm();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}

class DxInfo {
	public String ninfo1_dx, onset_GT20, bradykinesia, brady_B, brady_F, restT_B, actionT_B, restT_F, rigid_B, rigid_F, prob_GTE50, ldopa_GT50;

	public int atLeast1, atLeast2;

	public int brainBank_atLeast1, brainBank_atLeast3;

	public double restTremor, rigidity;

	public int atLeast1missing, atLeast2missing, restTremorMissing, rigidityMissing, exclusionMissing;

	public Vector<String> inclusion, supporting, exclusion, reasoning;

	public Vector<String> brainBank_inclusion, brainBank_supporting, brainBank_exclusion, brainBank_reasoning;

	public DxInfo() {
		// ninfo1_dx = onset_GT20 = bradykinesia = prob_GTE50 = ldopa_GT50 = "-1";
		ninfo1_dx = onset_GT20 = bradykinesia = prob_GTE50 = ldopa_GT50 = "0";
		// ninfo1_dx = onset_GT20 = bradykinesia = prob_GTE50 = ldopa_GT50 = ".";
		atLeast1 = atLeast2 = 0;
		brainBank_atLeast1 = brainBank_atLeast3 = 0;
		restTremor = rigidity = 0;
		atLeast1missing = atLeast2missing = restTremorMissing = rigidityMissing = exclusionMissing = 0;
		inclusion = new Vector<String>();
		brainBank_inclusion = new Vector<String>();
		supporting = new Vector<String>();
		brainBank_supporting = new Vector<String>();
		exclusion = new Vector<String>();
		brainBank_exclusion = new Vector<String>();
		reasoning = new Vector<String>();
	}
}
