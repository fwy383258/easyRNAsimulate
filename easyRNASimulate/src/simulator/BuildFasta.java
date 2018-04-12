package simulator;

import java.util.ArrayList;
import java.util.HashMap;

public class BuildFasta {
	
	boolean getScriptReads(ArrayList<String> out, ArrayList<String> out_ip,ArrayList<String> report, Transcript script, int count, int count_ip, int length, boolean circ_flag, boolean peak_flag) {
		StringBuffer id = new StringBuffer();
		StringBuffer base_seq = new StringBuffer();
		
		int base_sum = script.getBase_sum();
		if (base_sum < length) {
			return false;
		}
		id.append(">");
		int start_index = 0;
		int end_index = base_sum - length;
		if (script.isCirc_flag()) {
			int bj_index = script.getExon(0).getEnd() - script.getExon(0).getStart();
			if (circ_flag) {
				start_index = bj_index - length + 1;
				end_index = bj_index;
				if (start_index < 0) {
					end_index -= start_index;
					start_index = 0;
				}
				if (end_index > base_sum - length) {
					start_index -= end_index - base_sum + length;
					end_index = base_sum - length;
				}
				if (start_index < 0) {
					return false;
				}
				if (peak_flag) {
					id.append("peak:");
				}
				id.append("circ:");
			}
			else {
				int front_l = bj_index - 2 * length + 2;
				if (front_l < 0) {
					front_l = 0;
				}
				int back_l = end_index - length - bj_index;
				if (back_l < 0) {
					back_l = 0;
				}
				if (front_l + back_l == 0) {
					return false;
				}
				int rand = Method.randIntReach(0, front_l + back_l - 1);
				if (rand >= front_l) {
					rand += bj_index + 1 - front_l;
				}
				start_index = rand;
				end_index = rand + length - 1;
				if (peak_flag) {
					if (start_index - length / 4 >= 0 && end_index + length / 4 < base_sum - length) {
						rand = Method.randIntReach(- length/4, length/4);
						start_index += rand;
						end_index += rand;
					}
					id.append("peak:");
				}
				id.append("circn:");
			}
		}
		else if (peak_flag) {
			if (end_index - length > start_index) {
				int rand = Method.randIntReach(start_index, end_index - length);
				start_index = rand;
				end_index = rand + length;
			}
			id.append("peak:");
		}
		id.append(script.getId());
		id.append(':');
		int id_fix = id.length();
		int len = 0;
		if (report != null && peak_flag) {
			Exon last_exon = null;
			StringBuffer report_line = new StringBuffer();
			report_line.append(script.getId());
			report_line.append("\t");
			report_line.append(script.getExon(0).getChr_symbol());
			report_line.append("\t");
			for (int i = 0; i < script.getExons().size(); i++) {
				Exon exon = script.getExon(i);
				if (exon != last_exon && start_index < exon.getEnd() - exon.getStart() + 1 + len) {
					report_line.append(exon.getStart());
					report_line.append("\t");
					report_line.append(exon.getEnd());
					report_line.append("\t");
				}
				last_exon = exon;
				len += exon.getEnd() - exon.getStart() + 1;
				if (end_index + length < len) {
					break;
				}
			}
			report.add(report_line.toString());
		}
		this.getScriptReads(out, script, count, id, id_fix, base_seq, start_index, end_index, length);
		if (peak_flag) {
			this.getScriptReads(out_ip, script, count_ip, id, id_fix, base_seq, start_index, end_index, length);
		}
//		HashMap<Integer, Integer> count_id = new HashMap<Integer, Integer>();
//		while (count > 0) {
//			id.setLength(id_fix);
//			base_seq.setLength(0);
//			int rand = Method.randIntReach(start_index, end_index);
//			if (count_id.containsKey(rand)) {
//				count_id.put(rand, count_id.get(rand) + 1);
//			}
//			else {
//				count_id.put(rand, 1);
//			}
//			id.append(count_id.get(rand));
//			id.append(':');
//			for (int i=0; i < script.getExons().size(); i++) {
//				Exon exon = script.getExon(i);
//				int l = exon.getEnd() - exon.getStart() + 1;
//				if (rand < l) {
//					int left_length = this.buildFirstSplice(exon, rand, length, base_seq, id);
//					int j = i;
//					while (left_length > 0) {
//						j++;
//						left_length = this.addBoundSplice(base_seq, id, left_length, script.getExon(j));
//					}
//					break;
//				}
//				else {
//					rand -= l;
//				}
//			}			
//			out.add(id.toString());
//			out.add(base_seq.toString());
//			count--;
//		}
//		while (count_ip > 0) {
//			id.setLength(id_fix);
//			base_seq.setLength(0);
//			int rand = Method.randIntReach(start_index, end_index);
//			if (count_id.containsKey(rand)) {
//				count_id.put(rand, count_id.get(rand) + 1);
//			}
//			else {
//				count_id.put(rand, 1);
//			}
//			id.append(count_id.get(rand));
//			id.append(':');
//			for (int i=0; i < script.getExons().size(); i++) {
//				Exon exon = script.getExon(i);
//				int l = exon.getEnd() - exon.getStart() + 1;
//				if (rand < l) {
//					int left_length = this.buildFirstSplice(exon, rand, length, base_seq, id);
//					int j = i;
//					while (left_length > 0) {
//						j++;
//						left_length = this.addBoundSplice(base_seq, id, left_length, script.getExon(j));
//					}
//					break;
//				}
//				else {
//					rand -= l;
//				}
//			}			
//			out_ip.add(id.toString());
//			out_ip.add(base_seq.toString());
//			count_ip--;
//		}
		return true;
	}
	
	private void getScriptReads(ArrayList<String> out, Transcript script, int count, StringBuffer id, int id_fix, StringBuffer base_seq, int start_index, int end_index, int length) {
		HashMap<Integer, Integer> count_id = new HashMap<Integer, Integer>();
		while (count > 0) {
			id.setLength(id_fix);
			base_seq.setLength(0);
			int rand = Method.randIntReach(start_index, end_index);
			if (count_id.containsKey(rand)) {
				count_id.put(rand, count_id.get(rand) + 1);
			}
			else {
				count_id.put(rand, 1);
			}
			id.append(count_id.get(rand));
			id.append(':');
			for (int i=0; i < script.getExons().size(); i++) {
				Exon exon = script.getExon(i);
				int l = exon.getEnd() - exon.getStart() + 1;
				if (rand < l) {
					int left_length = this.buildFirstSplice(exon, rand, length, base_seq, id);
					int j = i;
					while (left_length > 0) {
						j++;
						left_length = this.addBoundSplice(base_seq, id, left_length, script.getExon(j));
					}
					break;
				}
				else {
					rand -= l;
				}
			}			
			out.add(id.toString());
			out.add(base_seq.toString());
			count--;
		}
	}
	
//	int buildFirstSplice(Transcript script, int offset, int length, StringBuffer base_seq, StringBuffer id) {
//		int left_length = 0;
//		if (length > 0) {
//			String exon_base_seq = exon.getBase_seq();
//			left_length = length + offset - exon_base_seq.length() ;
//			if (left_length > 0) {
//				base_seq.append(exon_base_seq.substring(offset));
//				id.append(exon.getChr_symbol());
//				id.append(':');
//				id.append(exon.getStart() + offset);
//				id.append('-');
//				id.append(exon.getEnd());
//			}
//			else {
//				base_seq.append(exon_base_seq.substring(offset, offset + length));
//				id.append(exon.getChr_symbol());
//				id.append(':');
//				id.append(exon.getStart() + offset);
//				id.append('-');
//				id.append(exon.getStart() + offset + length - 1);
//			}
//		}
//		return left_length;
//	}
	
	int buildFirstSplice(Exon exon, int offset, int length, StringBuffer base_seq, StringBuffer id) {
		int left_length = 0;
		if (length > 0) {
			String exon_base_seq = exon.getBase_seq();
			left_length = length + offset - exon_base_seq.length() ;
			if (left_length > 0) {
				base_seq.append(exon_base_seq.substring(offset));
				id.append(exon.getChr_symbol());
				id.append(':');
				id.append(exon.getStart() + offset);
				id.append('-');
				id.append(exon.getEnd());
			}
			else {
				base_seq.append(exon_base_seq.substring(offset, offset + length));
				id.append(exon.getChr_symbol());
				id.append(':');
				id.append(exon.getStart() + offset);
				id.append('-');
				id.append(exon.getStart() + offset + length - 1);
			}
		}
		return left_length;
	}
	
	int addBoundSplice(StringBuffer base_seq, StringBuffer id, int length, Exon next_exon) {
		int left_length = 0;
		if (length > 0) {
			left_length = length - next_exon.getBase_seq().length();
			if (left_length > 0) {
				base_seq.append(next_exon.getBase_seq());
				appendID(id, next_exon.getStart(), next_exon.getEnd());
			}
			else {
				base_seq.append(next_exon.getBase_seq().substring(0, length));
				appendID(id, next_exon.getStart(), next_exon.getStart() + length - 1);
			}
		}
		return left_length;
	}
	
	void addAllBoundForwardSplice(StringBuffer base_seq, StringBuffer id, int length, int exon_index, ArrayList<Exon> exon_list){
		if (length > 0) {
			int next_index = exon_list.size() - 1;
			if (exon_index < next_index) {
				next_index = Method.randIntReach(exon_index + 1, next_index);
			}
			else {
				System.out.println("Warning: Last exon forward splice");
				id.append('L');
			}
			Exon next_exon = exon_list.get(next_index);
			int left_length = length - next_exon.getBase_seq().length();
			if (left_length > 0) {
				base_seq.append(next_exon.getBase_seq());
				appendID(id, next_exon.getStart(), next_exon.getEnd());
				addAllBoundForwardSplice(base_seq, id, left_length, next_index, exon_list);
			}
			else {
				base_seq.append(next_exon.getBase_seq().substring(0, length));
				appendID(id, next_exon.getStart(), next_exon.getStart() + length - 1);
			}
		}
	}
	
	private void appendID(StringBuffer id, int start, int end) {
		id.append('_');
		id.append(start);
		id.append('-');
		id.append(end);
	}
}
